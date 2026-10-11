"""The HTML page and the .md twin must not drift apart.

Three defects motivated this file, each found by hand rather than by a test:
  · the .md emitted "| Protein change (HGVS p.) |  |" where the HTML omitted the row;
  · the frequency was stated TWICE in the .md, once in a format superseded weeks earlier
    and kept alive because the line is pickled into the record at build time;
  · the .md had no ancestry breakdown after the HTML gained one.

The surfaces are deliberately NOT identical — the .md is the long-form surface and the
page the scannable one — so these tests pin specific shared FACTS and specific rules,
not whole-document equality.
"""
import glob
import gzip
import os
import pickle
import random

import pytest

CACHE = os.environ.get("RECORDS_DIR") or "/data2/sugi-variant-cache/records"


def _records(limit_genes=25, want=None):
    files = sorted(glob.glob(os.path.join(CACHE, "*.pkl.gz")))
    if not files:
        pytest.skip(f"no cached records at {CACHE}")
    random.Random(11).shuffle(files)
    for f in files[:limit_genes]:
        try:
            with gzip.open(f, "rb") as fh:
                recs = pickle.load(fh)
        except Exception:
            continue
        for r in recs:
            if want is None or want(r):
                yield r


def _both(rec):
    os.environ.setdefault("BASE_PATH", "/variant")
    import app
    from sugivariant.enrich import disagreement_flag
    from sugivariant.render import render_body
    g, pops = app._freq_ctx(rec)
    html = app._render("variant.html", v=rec, canonical=rec.get("canonical_slug"),
                       nav="variant", jsonld="", disagreement=disagreement_flag(rec),
                       freq_g=g, freq_pops=pops)
    return html, render_body(rec), g, pops


# ── rule: frequency lives in the frequency block, not in the concordance readout ──
def test_concordance_lines_never_carry_frequency_prose():
    """Frequency is shown but never counted as concordant evidence, and stating it in
    both places is what produced the duplicate (in two different formats)."""
    from sugivariant.enrich import concordance
    for gnomad in ({"absent": True, "band": "absent from gnomAD v4.1"},
                   {"absent": False, "is_common": True, "band": "common (popmax 20%)"},
                   {"absent": False, "is_common": False, "band": "ultra-rare",
                    "faf95": 3.69e-04, "popmax": 5.36e-04}):
        c = concordance("Pathogenic", {"class": "likely_pathogenic", "score": "0.99"}, gnomad)
        assert not any("gnomAD" in ln or "filtering AF" in ln.lower() for ln in c["lines"]), \
            f"frequency prose leaked back into concordance lines: {c['lines']}"


def test_common_in_gnomad_still_raises_the_disagreement_flag():
    """Removing the prose must not remove the signal it used to accompany."""
    from sugivariant.enrich import concordance
    c = concordance("Pathogenic", None, {"absent": False, "is_common": True,
                                         "band": "common (popmax 20%)"})
    assert c["flags"] and "disagree" in (c["verdict"] or "").lower()


# ── rule: key/value rows with no value are never emitted ──
def test_kv_table_drops_rows_with_an_empty_value():
    from sugivariant.util import kv_table
    out = kv_table([("Gene", "HFE"), ("Protein change (HGVS p.)", None),
                    ("Coding change", ""), ("Type", "SNV")])
    assert "Protein change" not in out and "Coding change" not in out
    assert "HFE" in out and "SNV" in out


def test_md_emits_no_labelled_row_without_a_value():
    seen = 0
    for rec in _records():
        _, md, _, _ = _both(rec)
        for section in ("## Identity", "## For patients"):
            if section not in md:
                continue
            seg = md.split(section, 1)[1].split("\n## ", 1)[0]
            for line in seg.split("\n"):
                if not line.startswith("| ") or "---" in line:
                    continue
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                if len(cells) == 2 and cells[0] and not cells[1]:
                    pytest.fail(f"{rec.get('canonical_slug')}: empty value row -> {line!r}")
        seen += 1
        if seen >= 300:
            break
    assert seen, "no records examined"


# ── parity: the frequency facts must read the same on both surfaces ──
def test_frequency_facts_agree_across_surfaces():
    checked = 0
    for rec in _records(want=lambda r: (r.get("gnomad") or {}).get("popmax") is not None):
        html, md, g, pops = _both(rec)
        pct = f"{g['popmax'] * 100:.3g}%"
        assert pct in html, f"{rec.get('canonical_slug')}: grpmax {pct} missing from HTML"
        assert pct in md, f"{rec.get('canonical_slug')}: grpmax {pct} missing from .md"
        if g.get("faf95") is not None:
            faf = f"{g['faf95'] * 100:.3g}%"
            assert faf in html and faf in md, \
                f"{rec.get('canonical_slug')}: FAF95 {faf} not on both surfaces"
        if len(pops) > 1:
            from sugivariant.render import anc_name
            for code, _v in pops[:3]:
                assert anc_name(code) in html and anc_name(code) in md, \
                    f"{rec.get('canonical_slug')}: ancestry {code} not on both surfaces"
        checked += 1
        if checked >= 60:
            break
    if not checked:
        pytest.skip("no records with a grpmax frequency in the sample")


def test_absence_is_framed_as_pm2_on_both_surfaces():
    checked = 0
    for rec in _records(want=lambda r: (r.get("gnomad") or {}).get("absent")):
        html, md, _, _ = _both(rec)
        assert "PM2" in html, f"{rec.get('canonical_slug')}: HTML drops the PM2 framing"
        assert "PM2" in md, f"{rec.get('canonical_slug')}: .md drops the PM2 framing"
        checked += 1
        if checked >= 20:
            break
    if not checked:
        pytest.skip("no gnomAD-absent records in the sample")


def test_md_never_states_the_frequency_twice():
    """Both the live format and the legacy pickled one, which is why it went unnoticed."""
    checked = 0
    for rec in _records(want=lambda r: (r.get("gnomad") or {}).get("popmax") is not None):
        _, md, _, _ = _both(rec)
        assert "filtering AF (faf95)" not in md, \
            f"{rec.get('canonical_slug')}: legacy raw-decimal frequency line resurfaced"
        assert md.count("- Population frequency") <= 1, \
            f"{rec.get('canonical_slug')}: frequency stated more than once"
        checked += 1
        if checked >= 60:
            break
    if not checked:
        pytest.skip("no suitable records in the sample")


# ── grpmax vs raw maximum ─────────────────────────────────────────────────────
# gnomAD excludes bottlenecked groups (fin/asj/ami/remaining) from grpmax, so a group
# with a HIGHER raw frequency can sit above the grpmax row in a frequency-sorted table.
# Measured: that happens for 12.6% of variants. The badge used to read "highest", which
# made correct data look like a sorting bug (reported on kdm6a-p-arg571gln).
def test_excluded_groups_are_never_reported_as_grpmax():
    """The exclusion set is empirical — if gnomAD ever starts reporting one of these as
    grpmax, this fails and the labelling needs revisiting."""
    from sugivariant.render import GRPMAX_EXCLUDED
    seen = 0
    for rec in _records(limit_genes=40):
        g = rec.get("gnomad") or {}
        if not g.get("ancestry"):
            continue
        assert g["ancestry"] not in GRPMAX_EXCLUDED, \
            f"{rec.get('canonical_slug')}: {g['ancestry']} reported as grpmax but is in GRPMAX_EXCLUDED"
        seen += 1
    assert seen, "no records with a grpmax ancestry examined"


def test_grpmax_row_is_labelled_grpmax_not_highest():
    """A bottlenecked group outranking grpmax must read as data, not as a bug."""
    checked = 0
    for rec in _records(limit_genes=60):
        g = rec.get("gnomad") or {}
        pops = sorted(((k, float(v)) for k, v in (g.get("populations") or {}).items() if v),
                      key=lambda kv: -kv[1])
        if len(pops) < 2 or not g.get("ancestry") or pops[0][0] == g["ancestry"]:
            continue
        from sugivariant.render import GRPMAX_EXCLUDED
        if pops[0][0] not in GRPMAX_EXCLUDED:
            continue
        html, md, _, _ = _both(rec)
        for surface, name in ((html, "HTML"), (md, ".md")):
            assert "grpmax" in surface, f"{rec.get('canonical_slug')}: {name} never names grpmax"
            assert "not grpmax-eligible" in surface, \
                f"{rec.get('canonical_slug')}: {name} does not mark the excluded group"
            assert "bottlenecked" in surface, \
                f"{rec.get('canonical_slug')}: {name} does not explain why it outranks grpmax"
        assert ">highest<" not in html, "the misleading 'highest' badge is back"
        checked += 1
        if checked >= 8:
            break
    if not checked:
        pytest.skip("no variant where an excluded group outranks grpmax in the sample")


# ── SpliceAI status is stated, never implied by silence ───────────────────────
# Rendering nothing when SpliceAI has no score read as "no splice signal" — a positive
# claim made from an absence. Our release also holds only ~1 of 3 alternate alleles per
# scored position, so a canonical splice variant could sit next to a Δ0.99 donor-loss
# score for a different substitution and show nothing at all.
def test_spliceai_absence_is_stated_on_both_surfaces():
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1-plus-1g-c", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 1, "name": "X:c.1+1G>C", "hgvs_c": "c.1+1G>C",
           "review_status": "criteria provided, single submitter", "conditions": [],
           "consequence": {"label": "canonical splice-site (±1/±2)", "type": "splice"},
           "spliceai": {"status": "other_allele", "effect": "donor_loss",
                        "score": "0.99", "ref": "C", "alt": "A"}}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        # A scored position with this allele missing is NOT a gap in our data: the
        # release keeps every allele reaching delta 0.2, so this allele was evaluated
        # and fell under it. That is a prediction of no effect, and saying only
        # "not scored" threw the information away.
        # Assert the CLAIM, not the prose — the page copy is deliberately terse and
        # gets tightened; what must never drift is what it asserts.
        assert "0.2" in surface, f"{name} omits the threshold that makes absence meaningful"
        assert "0.99" in surface, f"{name} omits the same-position score"
        assert "evaluated" in surface, f"{name} does not say the allele was evaluated"
        assert "not scored for this allele" not in surface, (
            f"{name} still frames a below-threshold allele as missing data")
        assert "not a gap" in surface or "not missing" in surface, (
            f"{name} does not distinguish a below-threshold prediction from missing data")

    rec["spliceai"] = {"status": "absent"}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "not assessed" in surface, f"{name} renders silence instead of 'not assessed'"
        # Position absent has TWO readings (out of MANE scope, or every allele below
        # threshold) and we cannot tell them apart, so neither may be asserted.
        assert "not evidence against" in surface, f"{name} omits the absence caveat"


def test_legacy_spliceai_records_keep_their_readout():
    """Records pickled before the status field existed carry {effect, score} only. They
    must still render and still contribute their concordance line."""
    os.environ.setdefault("BASE_PATH", "/variant")
    from sugivariant.enrich import concordance, spliceai_scored
    legacy = {"effect": "donor_loss", "score": "0.92"}
    assert spliceai_scored(legacy) is True
    c = concordance("Pathogenic", None, None, legacy, {"phylop": 7.4})
    assert any("SpliceAI" in ln for ln in c["lines"])
    assert spliceai_scored({"status": "absent"}) is False
    assert spliceai_scored(None) is False


def test_canonical_splice_gain_label_is_flagged_incomplete():
    """Our release stores one of SpliceAI's four deltas, so ~60% of canonical ±1/±2
    variants that get a row are labelled a GAIN — which cannot be the mechanism at a
    ±1/±2 position. Say so rather than publishing it bare."""
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1-plus-2t-g", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 2, "name": "X:c.1+2T>G", "hgvs_c": "c.1+2T>G",
           "review_status": "criteria provided, single submitter", "conditions": [],
           "consequence": {"label": "canonical splice-site (±1/±2)", "type": "splice"},
           "spliceai": {"status": "scored", "effect": "donor_gain", "score": "0.87"}}
    html, _, _, _ = _both(rec)
    assert "effect label incomplete" in html


# ── CADD: a described score, never a verdict ──────────────────────────────────
# CADD is the first score we carry that covers non-missense variants allele-specifically.
# It has NO single ClinGen-calibrated threshold, and the panels that use it disagree
# (GN167/GN208 >=20.0 for UTR, GN164 >=28.1 for synonymous/indel, GN125 imports Pejaver's
# missense numbers, GN141 declines it outright). We show the score and the disagreement,
# and pick no threshold ourselves.
def test_cadd_is_shown_with_its_disagreement_and_no_threshold():
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1-plus-1g-c", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 3, "name": "X:c.1+1G>C", "hgvs_c": "c.1+1G>C",
           "review_status": "criteria provided, single submitter", "conditions": [],
           "cadd": {"phred": 29.6}}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "29.6" in surface, f"{name} omits the CADD score"
        assert "disagree" in surface or "declines" in surface, \
            f"{name} presents CADD without the panel disagreement"
    # we must not assert a verdict of our own
    for bad in ("CADD: pathogenic", "CADD supports pathogenicity", "meets PP3", "meets BP4"):
        assert bad not in html and bad not in md, f"CADD presented as a verdict: {bad!r}"


def test_cadd_absent_renders_nothing_rather_than_a_zero():
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-p-y", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 4, "name": "X:c.1A>G", "hgvs_c": "c.1A>G",
           "review_status": "criteria provided, single submitter", "conditions": [],
           "cadd": None}
    html, md, _, _ = _both(rec)
    assert "<th>CADD</th>" not in html
    assert "- CADD:" not in md


# ── regional constraint (gnomAD RMC + Gnocchi) ────────────────────────────────
#
# Both are REGIONAL: they describe the interval the variant falls in, and every variant
# in that interval shares the figure. Sitting in a table of allele-specific scores, a
# bare o/e or Z would read as a property of this variant, so both surfaces must say
# otherwise — that is the claim these tests pin.

_RMC_REC = {"canonical_slug": "x-p-pro9ala", "gene_symbol": "X", "classification": "VUS",
            "variation_id": 5, "name": "X:c.25C>G", "hgvs_c": "c.25C>G",
            "hgvs_p": "p.Pro9Ala", "review_status": "criteria provided, single submitter",
            "conditions": [],
            "rmc": {"transcript": "ENST00000357654", "oe": 0.23, "constrained": True,
                    "start": 43045681, "end": 43124096, "p": 4.8083e-05, "on_mane": True}}

_GNO_REC = {"canonical_slug": "x-c-100-plus-55a-g", "gene_symbol": "X", "classification": "VUS",
            "variation_id": 6, "name": "X:c.100+55A>G", "hgvs_c": "c.100+55A>G",
            "review_status": "criteria provided, single submitter", "conditions": [],
            "gnocchi": {"z": 4.39, "oe": 0.41, "start": 786000, "end": 787000,
                        "obs": 40.0, "exp": 97.0}}


def test_regional_constraint_is_framed_as_regional_on_both_surfaces():
    os.environ.setdefault("BASE_PATH", "/variant")
    html, md, _, _ = _both(_RMC_REC)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "0.23" in surface, f"{name} omits the RMC observed/expected"
        assert "region" in surface.lower(), f"{name} fails to say the figure is regional"
        assert "MANE" in surface, f"{name} omits which transcript the region came from"
    # never an applied code, and never presented as a property of this allele
    for bad in ("meets PM1", "meets PP3", "this variant is constrained"):
        assert bad not in html and bad not in md, f"RMC presented as a verdict: {bad!r}"


def test_gnocchi_is_shown_without_any_calibrated_threshold():
    os.environ.setdefault("BASE_PATH", "/variant")
    html, md, _, _ = _both(_GNO_REC)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "4.39" in surface, f"{name} omits the Gnocchi Z score"
        assert "calibrat" in surface.lower(), \
            f"{name} shows Gnocchi without stating that it has no ACMG calibration"
        assert "1 kb" in surface, f"{name} does not say the score describes a 1 kb window"


def test_regional_constraint_absent_renders_nothing_rather_than_a_zero():
    """The recurring Jinja trap: a missing key is Undefined, and `Undefined is not none`
    is True. Legacy records predate both datasets, so both must render nothing."""
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-p-z", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 7, "name": "X:c.2A>G", "hgvs_c": "c.2A>G",
           "review_status": "criteria provided, single submitter", "conditions": []}
    html, md, _, _ = _both(rec)                      # keys entirely ABSENT, not None
    assert "<th>Regional constraint</th>" not in html
    assert "<th>Non-coding constraint</th>" not in html
    assert "Regional missense constraint:" not in md
    assert "Non-coding constraint" not in md


def test_each_regional_dataset_is_gated_to_its_own_variant_class():
    """RMC models missense; Gnocchi models non-coding sequence. Showing either outside
    its class attaches a figure that does not describe the variant."""
    from sugivariant.enrich import is_missense, is_noncoding
    cases = [({"hgvs_p": "p.Pro309Ala"}, True, False),
             ({"hgvs_p": "p.Arg1443Ter"}, False, False),
             ({"hgvs_p": "p.Lys23fs"}, False, False),
             ({"hgvs_p": "p.Val12del"}, False, False),
             ({"hgvs_p": "p.Met1?"}, False, False),
             ({"hgvs_c": "c.100+55A>G"}, False, True),
             ({"hgvs_c": "c.-23G>A"}, False, True)]
    for rec, mis, non in cases:
        assert is_missense(rec) is mis, f"is_missense wrong for {rec}"
        assert is_noncoding(rec) is non, f"is_noncoding wrong for {rec}"


def test_rmc_declines_to_guess_when_transcripts_are_ambiguous(monkeypatch):
    """RMC regions are per-transcript and a gene has dozens. With no MANE match and more
    than one candidate, showing one would attach a region from a transcript the variant
    was never reported on."""
    from sugivariant import enrich as EN
    rows = [{"id": "ENST1:100", "transcript": "ENST1", "oe": "0.2", "constrained": "true"},
            {"id": "ENST2:100", "transcript": "ENST2", "oe": "0.9", "constrained": "false"}]
    monkeypatch.setattr(EN, "_interval_map", lambda pos, ds, c=None: rows)
    monkeypatch.setattr(EN, "_interval_attrs", lambda i, ds, c: {})
    assert EN.rmc_for("1:100", None, {}) is None            # ambiguous → nothing
    got = EN.rmc_for("1:100", "ENST2", {})                  # MANE picks the right one
    assert got["transcript"] == "ENST2" and got["on_mane"] is True
    assert got["oe"] == 0.9


def test_rmc_constrained_false_string_is_not_read_as_true(monkeypatch):
    """biobtree projects the boolean through the map chain as the STRING "false", and
    bool("false") is True. Truth-testing it would mark every unconstrained region as
    depleted of missense variation — a clinical claim, inverted."""
    from sugivariant import enrich as EN
    monkeypatch.setattr(EN, "_interval_attrs", lambda i, ds, c: {})
    for raw, expected in (("false", False), ("true", True), (False, False), (True, True)):
        monkeypatch.setattr(EN, "_interval_map", lambda pos, ds, c=None, _r=raw: [
            {"id": "E:1", "transcript": "E", "oe": "0.5", "constrained": _r}])
        got = EN.rmc_for("1:100", "E", {})
        assert got["constrained"] is expected, f"constrained={raw!r} read as {got['constrained']}"


def test_gnocchi_autosome_only_gap_is_disclosed_not_silent():
    """Gnocchi is computed on the autosomes only, so every X-linked gene's non-coding
    variants have no window. Silence there is indistinguishable from an unconstrained
    window — the same misreading we fixed for SpliceAI."""
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-9-plus-5a-g", "gene_symbol": "MECP2", "classification": "VUS",
           "variation_id": 8, "name": "MECP2:c.9+5A>G", "hgvs_c": "c.9+5A>G",
           "review_status": "criteria provided, single submitter", "conditions": [],
           "gnocchi": {"z": None, "not_covered": "X"}}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "not assessed" in surface, f"{name} falls silent instead of disclosing the gap"
        assert "autosome" in surface, f"{name} does not explain why there is no window"
        assert "not evidence" in surface, f"{name} omits that absence is not evidence"
    # and it must not print a Z of None/0
    assert "Z None" not in html and "Z 0.00" not in html


def test_interval_fetch_failure_degrades_the_build_and_is_not_cached(monkeypatch):
    """A network blip must not become a permanent "no constraint here" claim. On failure
    the accessor flags the build (so the record is never persisted) and caches nothing,
    so the next attempt retries. A genuine empty result IS cached — it is real data."""
    from sugivariant import enrich as EN

    calls = {"n": 0}

    def boom(ident, chain, cap=None):
        calls["n"] += 1
        raise RuntimeError("connection reset by peer")

    monkeypatch.setattr(EN, "map_all", boom)
    EN.begin_build()
    cache = {}
    assert EN.gnocchi_for("1:786500", cache) is None
    assert EN.rmc_for("1:786500", "ENST1", cache) is None
    assert EN.build_degraded(), "a transport failure did not degrade the build"
    assert cache == {}, "a failure was cached — the retry would be skipped"
    assert calls["n"] == 2, "the second call reused a cached failure"

    # a genuine empty (position in no interval) is real data and SHOULD be cached
    monkeypatch.setattr(EN, "map_all", lambda i, c, cap=None: [])
    EN.begin_build()
    cache = {}
    assert EN.gnocchi_for("1:786500", cache) is None
    assert cache.get("map:1:786500") == [], "a genuine empty was not cached"
    assert not EN.build_degraded(), "an empty result wrongly degraded the build"


# ── the ClinVar ref/alt unlock (BioBTree release, Oct 2026) ──────────────────
#
# Giving indels a chr:pos:ref:alt key is the change with the widest blast radius in
# this file's history: it routes ~12% of variants into every coordinate-keyed lookup
# for the first time, including the one that asserts absence from gnomAD.

def test_indel_coordinate_is_built_from_position_vcf_not_start():
    """On an indel position_vcf and start DIFFER — the VCF form carries the anchor
    base. Keying off `start` would look up a neighbouring position and read as a
    false absence, so this pins which field is used."""
    from sugivariant.enrich import variant_coordinate
    rec = {"chromosome": "17", "start": 43093803, "position_vcf": 43093802,
           "ref_allele": "CT", "alt_allele": "C", "assembly": "GRCh38",
           "variant_type": "Deletion"}
    assert variant_coordinate(rec) == "17:43093802:CT:C"


def test_non_grch38_record_gets_no_coordinate():
    """ClinVar files some variants on GRCh37 only. Those coordinates would key every
    GRCh38 dataset at the wrong locus, which is worse than having no key at all."""
    from sugivariant.enrich import variant_coordinate
    rec = {"chromosome": "17", "position_vcf": 41223222, "ref_allele": "A",
           "alt_allele": "G", "assembly": "GRCh37"}
    assert variant_coordinate(rec) is None


def test_legacy_records_without_ref_alt_still_resolve_unchanged():
    """The ~4.1M records already on disk predate ClinVar's ref/alt and have no
    position_vcf. They must keep resolving through the genomic-HGVS route."""
    from sugivariant.enrich import variant_coordinate
    legacy = {"chromosome": "6", "start": 26093233,
              "hgvs_expressions": ["NC_000006.12:g.26093233G>T"]}
    assert variant_coordinate(legacy) == "6:26093233:G:T"


def test_a_missing_indel_is_never_claimed_absent_from_gnomad(monkeypatch):
    """THE safety case for the ref/alt unlock.

    "absent from gnomAD v4.1" renders as an ACMG PM2-supporting rarity claim. For an
    indel a coordinate miss is not reliable evidence of absence: ClinVar and gnomAD
    can left-align the same change differently inside a tandem repeat. Measured on
    benign BRCA1/BRCA2/CFTR/ATM/NF1 indels the VCF triplet resolves 67% of the time —
    enough to REPORT a hit, not enough to assert a miss. An SNV miss is unchanged."""
    import sugibiobtree
    from sugivariant import enrich

    def not_found(coord, dataset, *a, **k):
        raise sugibiobtree.client.BiobtreeError(
            f"entry → HTTP 400: Entry not found identifier {coord} dataset {dataset}")

    monkeypatch.setattr(sugibiobtree, "entry", not_found)
    enrich.begin_build()
    indel = {"chromosome": "17", "start": 43093803, "position_vcf": 43093802,
             "ref_allele": "CT", "alt_allele": "C", "assembly": "GRCh38", "rsid": None}
    g = enrich.gnomad_frequency(indel)
    assert not (g or {}).get("absent"), (
        "an indel that merely failed to match gnomAD's representation was reported "
        f"as absent, which reads as an ACMG rarity claim: {g}")

    snv = {"chromosome": "17", "start": 43071205, "position_vcf": 43071205,
           "ref_allele": "A", "alt_allele": "G", "assembly": "GRCh38", "rsid": None}
    assert enrich.gnomad_frequency(snv).get("absent") is True, (
        "an SNV miss IS a reliable absence and must still be stated")


def test_coord_is_snv_reads_the_alleles_not_the_type():
    from sugivariant.enrich import coord_is_snv
    assert coord_is_snv("17:43071205:A:G") is True
    assert coord_is_snv("17:43093802:CT:C") is False
    assert coord_is_snv("17:43090936:C:CT") is False
    assert coord_is_snv("17:43071205") is None      # positional key, not a verdict
    assert coord_is_snv(None) is None


def test_spliceai_declares_its_snv_only_scope_for_an_indel():
    """Indels reach spliceai_for for the first time now that they have a coordinate.
    The release is built from a substitutions-only source, so the honest answer is
    out-of-scope — NOT the 'absent' state, which a reader may take as no signal."""
    from sugivariant.enrich import spliceai_for
    sa = spliceai_for("17:43093802:CT:C", {"by_allele": {}, "by_pos": {}})
    assert sa["status"] == "snv_only", sa
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1del", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 2, "name": "X:c.1del", "hgvs_c": "c.1del",
           "variant_type": "Deletion", "conditions": [],
           "review_status": "criteria provided, single submitter", "spliceai": sa}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "substitutions only" in surface or "substitutions-only" in surface, (
            f"{name} hides the SNV-only scope")
        assert "not evidence" in surface, f"{name} omits the absence caveat"


def test_repeat_context_is_framed_as_confidence_not_pathogenicity():
    """A repeat overlap says the call and its frequency are less reliable. It says
    NOTHING about whether the variant is damaging, and must never read that way."""
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1del", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 3, "name": "X:c.1del", "hgvs_c": "c.1del",
           "variant_type": "Deletion", "conditions": [],
           "review_status": "criteria provided, single submitter",
           "repeat": {"n": 2, "types": ["tandem_repeat"], "classes": ["Simple_repeat"],
                      "names": ["(TAACCC)n"], "tandem": True, "segdup": False}}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert ("simple repeat" in surface.lower() or "simple_repeat" in surface.lower()), (
            f"{name} omits the repeat class")
        assert "ambiguous" in surface, f"{name} omits the alignment ambiguity"
        assert "not a pathogenicity signal" in surface, (
            f"{name} fails to disclaim repeat context as a pathogenicity signal")


def test_repeat_flags_tandem_and_segdup_separately():
    """They have different consequences — ambiguous alignment vs unreliable mapping —
    so collapsing them into one 'repetitive' flag would lose the distinction. An
    interspersed SINE/LINE is neither and must set no flag."""
    from sugivariant import enrich
    rows = [{"repeat_type": "tandem_repeat", "repeat_class": "Simple_repeat",
             "repeat_name": "(TAACCC)n"}]
    enrich._interval_map = lambda pos, ds, c=None: rows
    r = enrich.repeat_for("1:10200", {})
    assert r["tandem"] is True and r["segdup"] is False

    enrich._interval_map = lambda pos, ds, c=None: [
        {"repeat_type": "segmental_dup", "repeat_class": "Segmental_duplication",
         "repeat_name": "chr17:43153705"}]
    r = enrich.repeat_for("17:43125364", {})
    assert r["segdup"] is True and r["tandem"] is False

    enrich._interval_map = lambda pos, ds, c=None: [
        {"repeat_type": "repeatmasker", "repeat_class": "SINE", "repeat_name": "AluSx"}]
    r = enrich.repeat_for("1:1000", {})
    assert r["tandem"] is False and r["segdup"] is False, (
        "an interspersed element is not a tandem repeat and implies no alignment "
        f"ambiguity: {r}")


def test_pangolin_is_shown_as_independent_and_uncalibrated():
    """Pangolin's value is being a SECOND opinion to SpliceAI. ClinGen's SVI
    calibration covers SpliceAI only, so Pangolin must carry no band and no code."""
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1-plus-1g-c", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 4, "name": "X:c.1+1G>C", "hgvs_c": "c.1+1G>C",
           "variant_type": "single nucleotide variant", "conditions": [],
           "review_status": "criteria provided, single submitter",
           "consequence": {"label": "canonical splice-site (±1/±2)", "type": "splice"},
           "pangolin": {"status": "scored", "score": 0.87, "effect": "loss",
                        "gain": 0.02, "loss": -0.87, "gene": "ENSG00000012048"}}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "0.87" in surface, f"{name} omits the Pangolin score"
        assert "independent" in surface.lower(), f"{name} omits the independence framing"
        assert "Pangolin" in surface, f"{name} does not name the predictor"


def test_new_keys_absent_entirely_render_clean():
    """The Jinja Undefined trap, which has caused three outages here: a MISSING dict
    key is Undefined, and `Undefined is not none` evaluates True. Every record on disk
    today lacks both new keys, so this is the live path, not an edge case."""
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1a-g", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 5, "name": "X:c.1A>G", "hgvs_c": "c.1A>G", "conditions": [],
           "review_status": "criteria provided, single submitter"}
    assert "pangolin" not in rec and "repeat" not in rec
    html, md, _, _ = _both(rec)
    # Both names appear in every page's provenance footer, so target the ROW itself.
    assert "<th>Pangolin</th>" not in html, "HTML invented a Pangolin row from Undefined"
    assert "<th>Repeat context</th>" not in html, "HTML invented a repeat row"
    assert "- Pangolin:" not in md, ".md invented a Pangolin row from Undefined"
    assert "- Repeat context:" not in md, ".md invented a repeat row"
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "None" not in surface, f"{name} leaked a None into the page"


def test_splice_agreement_compares_the_right_thing():
    """The first version of this compared SpliceAI's score to 0.2 — but a SpliceAI
    'scored' state already means >= 0.2, because the release is filtered there, so
    that branch always said "agrees" and the disagreement case was unreachable.

    The informative comparison is whether the two predictors agree a signal EXISTS:
    Pangolin scoring a variant SpliceAI leaves under its own threshold is a real
    disagreement between independent models, and it must be visible."""
    from sugivariant.enrich import splice_agreement
    pg = {"status": "scored", "score": 0.87, "effect": "loss"}

    both = splice_agreement({"status": "scored", "effect": "donor_loss",
                             "score": "0.99"}, pg)
    assert both["state"] == "agree", both

    disagree = splice_agreement({"status": "other_allele", "score": "0.9"}, pg)
    assert disagree["state"] == "spliceai_below", disagree
    assert "0.2" in disagree["text"]
    assert splice_agreement({"status": "absent"}, pg)["state"] == "spliceai_below"

    # Agree that there is a signal, disagree on its direction — worth surfacing.
    opposed = splice_agreement({"status": "scored", "effect": "acceptor_gain",
                                "score": "0.6"}, pg)
    assert opposed["state"] == "direction_differs", opposed
    # The text names SpliceAI's direction only: it is rendered inside the Pangolin
    # row, which already states Pangolin's, so naming both was duplication.
    assert "gain" in opposed["text"], opposed

    # Nothing to compare: no Pangolin score, or an indel neither release covers.
    assert splice_agreement({"status": "scored", "score": "0.9"},
                            {"status": "absent"}) is None
    assert splice_agreement({"status": "snv_only"}, pg) is None
    assert splice_agreement(None, pg) is None
    # A legacy record carries a score with no status field at all.
    assert splice_agreement({"effect": "donor_loss", "score": "0.92"},
                            pg)["state"] == "agree"


# ── reading the reporting threshold from biobtree's metadata (their Oct 2026 ask) ──

def test_threshold_is_parsed_from_the_published_prose():
    """biobtree publishes the filter as PROSE, not a number, so it is parsed out.
    The exact strings below are what the live /ws/meta returns."""
    from sugivariant import enrich
    enrich._DATASET_META = {
        "spliceai": {"threshold": "reported only if max SpliceAI delta (DS_AG/AL/DG/DL)"
                                  " >= 0.2; a variant absent from this dataset is below"
                                  " threshold (no predicted splice effect), not 'not assessed'"},
        "pangolin": {"threshold": "reported only if the largest Pangolin delta-score >= 0.2;"
                                  " absence = below threshold, not 'not assessed'"},
    }
    try:
        assert enrich.dataset_threshold("spliceai", 0.99) == 0.2
        assert enrich.dataset_threshold("pangolin", 0.99) == 0.2
    finally:
        enrich._DATASET_META = None


def test_threshold_falls_back_rather_than_guessing():
    """A wrong threshold would mis-state what an ABSENT score means, which is a claim
    about the variant. Prose is not a machine-readable field, so every doubtful parse
    must yield the pinned default instead of a guess."""
    from sugivariant import enrich
    cases = [
        {},                                             # dataset absent from meta
        {"threshold": ""},                              # present but empty
        {"threshold": "filtered to a meaningful signal"},   # no number at all
        {"threshold": ">= 7"},                          # not a delta in (0,1)
        {"threshold": ">= 0"},                          # degenerate
        {"threshold": ">= abc"},                        # unparseable
    ]
    for rec in cases:
        enrich._DATASET_META = {"spliceai": rec}
        try:
            assert enrich.dataset_threshold("spliceai", 0.2) == 0.2, rec
        finally:
            enrich._DATASET_META = None


def test_a_failed_meta_fetch_is_bounded_but_still_retryable(monkeypatch):
    """dataset_threshold() runs per variant and the fetch carries a 5s timeout, so an
    uncached failure would cost 5s on every variant of a build. It must be cached —
    but only briefly, or one transient outage blanks every version stamp for the life
    of the process."""
    from sugivariant import enrich
    calls = {"n": 0}

    def boom(*a, **k):
        calls["n"] += 1
        raise OSError("connection refused")

    monkeypatch.setattr("urllib.request.urlopen", boom)
    enrich._DATASET_META = None
    enrich._DATASET_META_FAILED_AT = 0.0
    try:
        assert enrich._dataset_meta() == {}
        assert calls["n"] == 1
        for _ in range(25):
            assert enrich._dataset_meta() == {}
        assert calls["n"] == 1, (
            f"the 5s timeout was paid {calls['n']} times; a build would stall")
        # ...but the window must lapse, so the stamps come back after an outage.
        enrich._DATASET_META_FAILED_AT -= (enrich._DATASET_META_RETRY_AFTER + 1)
        enrich._dataset_meta()
        assert calls["n"] == 2, "a failure was cached permanently"
    finally:
        enrich._DATASET_META = None
        enrich._DATASET_META_FAILED_AT = 0.0


# How each dataset is keyed is an assumption baked into every lookup in enrich.py.
# It used to be invisible — which is how two datasets got reported as empty when they
# were fully populated, having simply been queried the wrong way. biobtree now
# publishes it, so the assumption can be checked instead of believed.
_KEY_ASSUMPTIONS = {
    # looked up with entry() on a full chr:pos:ref:alt key
    "cadd": "chr:pos:ref:alt", "spliceai": "chr:pos:ref:alt",
    "pangolin": "chr:pos:ref:alt", "gnomad_variant": "chr:pos:ref:alt",
    "alphamissense": "chr:pos:ref:alt", "revel": "chr:pos:ref:alt",
    # ref/alt-agnostic: enrich truncates the key to chr:pos for these
    "conservation": "chr:pos",
    # interval-indexed: unreachable by entry(), so reached via the map chain
    "gnomad_gnocchi": "interval", "gnomad_rmc": "interval", "repeat": "interval",
}


def test_our_keying_assumptions_match_what_biobtree_publishes():
    """If upstream re-keys a dataset, every lookup against it starts returning
    nothing — silently, and indistinguishably from "this variant has no data". That
    is the most expensive failure mode in this codebase, so it gets a tripwire."""
    from sugivariant import enrich
    enrich._DATASET_META = None
    enrich._DATASET_META_FAILED_AT = 0.0
    meta = enrich._dataset_meta()
    if not meta:
        pytest.skip("biobtree /ws/meta unreachable")
    missing, wrong = [], []
    for group, assumed in _KEY_ASSUMPTIONS.items():
        rec = meta.get(group)
        if rec is None:
            continue                      # dataset not in this deployment; not our concern
        scheme = rec.get("keyScheme")
        if not scheme:
            missing.append(group)         # can't verify — reported, not failed
            continue
        if assumed == "interval":
            ok = "interval" in scheme.lower()
        elif assumed == "chr:pos":
            ok = scheme.startswith("chr:pos") and "ref:alt" not in scheme
        else:
            ok = scheme.strip().startswith("chr:pos:ref:alt")
        if not ok:
            wrong.append(f"{group}: we key it as {assumed!r}, biobtree says {scheme!r}")
    assert not wrong, "keying assumptions now disagree with upstream:\n  " + "\n  ".join(wrong)
    if missing:
        print(f"note: no keyScheme published for {missing} — assumption unverified")
