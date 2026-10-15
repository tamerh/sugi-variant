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
    from sugivariant.enrich import disagreement_flag, refresh_concordance
    from sugivariant.render import render_body
    refresh_concordance(rec)          # what app.py does before rendering
    g, pops = app._freq_ctx(rec)
    html = app._render("variant.html", v=rec, canonical=rec.get("canonical_slug"),
                       nav="variant", jsonld="", freq_g=g, freq_pops=pops)
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


def test_rmc_only_ever_reports_the_genes_own_mane_transcript(monkeypatch):
    """RMC regions are per-transcript and a gene has dozens, so a region from the wrong
    transcript is a figure the variant was never reported under.

    The positional route could not always tell which transcript a region belonged to and
    accepted a lone region from ANY transcript. Asking the gene's MANE transcript makes
    that unrepresentable, which is what this pins: a non-MANE region is dropped, not
    ranked."""
    from sugivariant import enrich as EN
    rows = [{"id": "ENST1:100", "transcript": "ENST1", "oe": "0.2", "constrained": "true"},
            {"id": "ENST2:100", "transcript": "ENST2", "oe": "0.9", "constrained": "false"}]
    monkeypatch.setattr(EN, "map_all", lambda *a, **k: rows)
    monkeypatch.setattr(EN, "mane_select", lambda h: {"enst": "ENST2"})
    monkeypatch.setattr(EN, "_interval_attrs",
                        lambda i, ds, c: {"start": 50, "end": 150, "oe_chisq_p": "0.01"})
    regions = EN.gene_rmc("HGNC:1", {})
    assert [r["transcript"] for r in regions] == ["ENST2"], (
        f"a region from a non-MANE transcript survived: {regions}")
    got = EN.rmc_for("1:100", regions)
    assert got["oe"] == 0.9 and got["on_mane"] is True, got
    # and a position outside every region yields nothing rather than the nearest one
    assert EN.rmc_for("1:9999", regions) is None

    # with no MANE transcript there is nothing to key on, so nothing is claimed
    monkeypatch.setattr(EN, "mane_select", lambda h: None)
    assert EN.gene_rmc("HGNC:1", {}) == []


def test_rmc_constrained_false_string_is_not_read_as_true(monkeypatch):
    """biobtree projects the boolean through the map chain as the STRING "false", and
    bool("false") is True. Truth-testing it would mark every unconstrained region as
    depleted of missense variation — a clinical claim, inverted."""
    from sugivariant import enrich as EN
    monkeypatch.setattr(EN, "_interval_attrs",
                        lambda i, ds, c: {"start": 50, "end": 150, "oe_chisq_p": "0.01"})
    monkeypatch.setattr(EN, "mane_select", lambda h: {"enst": "E"})
    for raw, expected in (("false", False), ("true", True), (False, False), (True, True)):
        monkeypatch.setattr(EN, "map_all", lambda *a, _r=raw, **k: [
            {"id": "E:1", "transcript": "E", "oe": "0.5", "constrained": _r}])
        got = EN.rmc_for("1:100", EN.gene_rmc("HGNC:1", {}))
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
    monkeypatch.setattr(EN, "mane_select", lambda h: {"enst": "ENST1"})
    assert EN.gene_rmc("HGNC:1", cache) is EN.FETCH_ERROR
    assert EN.rmc_for("1:786500", EN.FETCH_ERROR) is EN.FETCH_ERROR
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
    # `absent` means there is no SpliceAI record for the position, so the delta could be
    # anything. Reading it as "below 0.2" made the Pangolin row contradict the SpliceAI
    # row two above it ("not assessed") on 957 of 1,249 measured cases.
    absent = splice_agreement({"status": "absent"}, pg)
    assert absent["state"] == "spliceai_absent", absent
    assert "0.2" not in absent["text"], absent
    assert absent["text"] == "not in the SpliceAI release"

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


# ── mobile: nothing may force the page wider than a phone screen ──────────────
#
# Found 2026-10-12 on a real phone. The citation line prints the page's full
# canonical URL, and for an indel slug that is ~59 unbreakable characters
# (https://sugi.bio/variant/pten-c-802-minus-9-802-minus-8insc) — about 384px at
# .82rem, past a 360px viewport, so the whole page scrolled sideways. It was not a
# new bug; it only became visible because indels are what this release made worth
# looking at, and indels get the long HGVS-derived slugs. A missense slug makes a
# 40-char URL and fits, which is why it had never shown up.

def _css():
    import pathlib
    return pathlib.Path("static/css/app.css").read_text(encoding="utf-8")


def _rule(css, selector):
    import re
    m = re.search(re.escape(selector) + r"\s*\{[^}]*\}", css)
    return m.group(0) if m else ""


def test_the_citation_url_can_wrap():
    """The one container that prints a full URL must be allowed to break it."""
    css = _css()
    rule = _rule(css, ".cite-text")
    assert rule, ".cite-text rule disappeared"
    assert "overflow-wrap:anywhere" in rule or "word-break:break-word" in rule, (
        "the citation line prints the page's full canonical URL as one unbreakable "
        f"token; without wrapping a long indel slug scrolls the page: {rule}")


def test_long_slug_pages_do_not_emit_unwrappable_tokens():
    """Render a page whose slug is as long as the real ones and check the citation
    URL is present AND that its container wraps — the two have to move together."""
    os.environ.setdefault("BASE_PATH", "/variant")
    slug = "pten-c-802-minus-9-802-minus-8insc"
    rec = {"canonical_slug": slug, "gene_symbol": "PTEN", "classification": "Likely benign",
           "variation_id": 9, "name": "NM_000314.8(PTEN):c.802-9_802-8insC",
           "hgvs_c": "c.802-9_802-8insC", "variant_type": "Insertion", "conditions": [],
           "review_status": "criteria provided, single submitter"}
    html_out, _, _, _ = _both(rec)
    if slug not in html_out:
        pytest.skip("citation block not rendered in this harness")
    longest = max((w for w in html_out.split() if "sugi.bio/variant/" in w), key=len, default="")
    assert len(longest) >= 40, "expected a long canonical URL to be present"
    assert "overflow-wrap:anywhere" in _rule(_css(), ".cite-text"), (
        f"page emits a {len(longest)}-char unbreakable URL with no wrapping")


def test_nowrap_tables_release_their_labels_on_narrow_viewports():
    """A nowrap cell sets its table's minimum width, so a long label pushes the whole
    page sideways on a phone. The labels are the only thing naming each metric, so
    they wrap on small screens rather than being shortened. Numeric cells keep
    nowrap deliberately — a split number reads as two numbers."""
    css = _css()
    assert "white-space:nowrap" in _rule(css, ".vd-preds th"), (
        "expected the desktop nowrap rule to still exist")
    import re
    narrow = re.findall(r"@media\s*\(max-width:\s*720px\)\s*\{(.*?)\n\}", css, re.S)
    blob = " ".join(narrow)
    for sel in (".vd-preds th", ".sources-table .src-name-cell"):
        assert sel in blob and "white-space:normal" in blob, (
            f"{sel} is not released from nowrap on narrow viewports")
    assert "overflow:hidden" not in _rule(css, ".sources-table-wrap "), (
        "the sources wrapper clips instead of scrolling, so the longest source names "
        "are unreachable on a phone rather than merely cramped")


# ── SpliceAI masked deltas (biobtree v2.16.0) ─────────────────────────────────
#
# The masked release is what Ensembl VEP and the Broad lookup serve. On 2026-10-13
# prod returned all four masked keys for every one of PTEN's 1,851 scored variants
# with the value 0.00, while dev returned the real figures — same version, same
# commit, same last_built. Trusting the keys' presence would have rendered "masked
# predicts no effect" on every splice readout on the site, so presence is never
# taken as evidence; a non-zero value somewhere in the gene is.

def test_masked_claims_need_the_gene_to_actually_have_masked_data():
    """THE safety case for masked. Unpopulated zeros are indistinguishable from a
    genuine "no effect" on any single variant, so the gate is per gene."""
    from sugivariant.enrich import spliceai_for
    raw_only = {"effect": "acceptor_loss", "score": "0.21", "ds_ag": 0.01, "ds_al": 0.21,
                "ds_dg": 0.0, "ds_dl": 0.0,
                "ds_ag_masked": 0.0, "ds_al_masked": 0.0,
                "ds_dg_masked": 0.0, "ds_dl_masked": 0.0}
    coord = "10:87894028:T:G"

    # prod's shape: every masked delta zero across the gene -> claim nothing
    unloaded = {"by_allele": {coord: raw_only}, "by_pos": {}, "masked_loaded": False}
    sa = spliceai_for(coord, unloaded)
    assert sa.get("masked_verdict") is None, (
        "masked data is not populated for this gene, so no masked claim may be made: "
        f"{sa.get('masked_verdict')}")

    # dev's shape: the gene has real masked values somewhere, so this allele's
    # all-zero masked IS the real negative evidence
    loaded = {"by_allele": {coord: raw_only}, "by_pos": {}, "masked_loaded": True}
    sa = spliceai_for(coord, loaded)
    assert sa["masked_verdict"] == "no_effect", sa
    assert sa["masked_max"] == 0.0


def test_masked_agreement_suppresses_the_disagreement_text():
    from sugivariant.enrich import spliceai_for
    agree = {"effect": "donor_gain", "score": "0.78", "ds_ag": 0.0, "ds_al": 0.0,
             "ds_dg": 0.78, "ds_dl": 0.52,
             "ds_ag_masked": 0.0, "ds_al_masked": 0.0,
             "ds_dg_masked": 0.78, "ds_dl_masked": 0.52}
    sa = spliceai_for("X:154032202:C:T",
                      {"by_allele": {"X:154032202:C:T": agree}, "by_pos": {},
                       "masked_loaded": True})
    assert sa["masked_verdict"] == "agrees", sa
    assert sa["dominant"] == "donor gain", sa


def test_gene_spliceai_detects_unpopulated_masked_deltas(monkeypatch):
    """masked_loaded is read off rows already fetched, so it costs nothing and starts
    working by itself once upstream populates the data — no deploy needed."""
    from sugivariant import enrich
    zeros = [{"id": "1:1:A:G", "effect": "acceptor_loss", "score": "0.3",
              "ds_ag": "0.00", "ds_al": "0.30", "ds_dg": "0.00", "ds_dl": "0.00",
              "ds_ag_masked": "0.00", "ds_al_masked": "0.00",
              "ds_dg_masked": "0.00", "ds_dl_masked": "0.00"}]
    real = [dict(zeros[0], ds_al_masked="0.30")]
    monkeypatch.setattr(enrich, "map_all", lambda *a, **k: zeros)
    assert enrich.gene_spliceai("HGNC:1")["masked_loaded"] is False
    monkeypatch.setattr(enrich, "map_all", lambda *a, **k: real)
    assert enrich.gene_spliceai("HGNC:1")["masked_loaded"] is True


def test_masked_disagreement_is_visible_on_both_surfaces():
    os.environ.setdefault("BASE_PATH", "/variant")
    sa = {"status": "scored", "effect": "acceptor_loss", "score": "0.21",
          "ds_ag": 0.01, "ds_al": 0.21, "ds_dg": 0.0, "ds_dl": 0.0,
          "ds_ag_masked": 0.0, "ds_al_masked": 0.0, "ds_dg_masked": 0.0,
          "ds_dl_masked": 0.0, "masked_verdict": "no_effect", "masked_max": 0.0,
          "dominant": "acceptor loss", "dominant_score": 0.21, "dominant_dp": -10}
    rec = {"canonical_slug": "x-c-83t-g", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 11, "name": "X:c.83T>G", "hgvs_c": "c.83T>G", "conditions": [],
           "variant_type": "single nucleotide variant",
           "review_status": "criteria provided, single submitter", "spliceai": sa}
    html_out, md, _, _ = _both(rec)
    assert "masked 0" in html_out, "HTML hides that masked disagrees"
    assert "10 bp away" in html_out, "HTML omits the distance that explains it"
    assert "no splice signal" in md, ".md does not say what VEP would show"
    assert "cryptic site" in md, ".md omits why raw and masked differ"
    for surface, name in ((html_out, "HTML"), (md, ".md")):
        assert "0.21" in surface, f"{name} drops the raw score"


# ── which upstream data a gene's records were built from ──────────────────────
#
# biobtree v2.16.1 added data_db_version, which tracks the INDEXED data. last_built
# does not: on 2026-10-13 prod reported last_built 2026-10-13T02:43:47 while the
# masked SpliceAI deltas were absent, and it did not change when they landed. The
# records cache had no way to know what it was built from, which is why three
# full-corpus rebuilds happened in three days.

def test_the_stamp_reads_data_db_version_and_never_last_built(monkeypatch):
    from sugivariant import enrich
    enrich._DATASET_META = {
        "clinvar": {"data_db_version": "db_v52", "last_built": "2026-10-13T02:43:47"},
        "spliceai": {"data_db_version": "db_v52", "last_built": "2026-10-13T02:43:47"},
        "cadd": {"last_built": "2026-10-13T02:43:47"},   # no data version published
    }
    try:
        v = enrich.record_data_versions()
        assert v == {"clinvar": "db_v52", "spliceai": "db_v52"}, v
        assert not any("2026" in str(x) for x in v.values()), (
            "last_built leaked into the stamp; it does not track the data")
    finally:
        enrich._DATASET_META = None


def test_stale_datasets_separates_out_of_date_from_cannot_tell(tmp_path):
    """"Cannot tell" and "up to date" must not collapse into each other: an unstamped
    legacy file would otherwise report itself current and never be rebuilt."""
    from sugivariant import enrich
    recs = tmp_path / "BRCA1.pkl.gz"
    recs.write_bytes(b"x")
    cur = {"clinvar": "db_v52", "spliceai": "db_v52"}

    # never stamped -> cannot tell, which is NOT the same as ()
    assert enrich.stale_datasets(recs, cur) is None

    assert enrich.write_data_stamp(recs, cur) is True
    assert enrich.read_data_stamp(recs) == cur
    assert enrich.stale_datasets(recs, cur) == ()          # up to date

    moved = {"clinvar": "db_v52", "spliceai": "db_v53"}
    assert enrich.stale_datasets(recs, moved) == ("spliceai",)

    added = dict(cur, cadd="db_v2")
    assert enrich.stale_datasets(recs, added) == ("cadd",)

    # meta unreachable -> cannot tell, and refuse to write a stamp we cannot vouch for
    assert enrich.stale_datasets(recs, {}) is None
    assert enrich.write_data_stamp(recs, {}) is False


def test_the_stamp_lives_beside_the_records_not_inside_them(tmp_path):
    """A new key on every record is the Jinja Undefined trap that has caused three
    outages here, and would carry the same string 4.1M times."""
    from sugivariant import enrich
    recs = tmp_path / "TP53.pkl.gz"
    recs.write_bytes(b"x")
    enrich.write_data_stamp(recs, {"clinvar": "db_v52"})
    assert (tmp_path / "TP53.dbver.json").exists(), "stamp not written beside the records"
    assert recs.read_bytes() == b"x", "the records file was modified"


def test_rmc_uses_the_variants_true_gene_not_the_files_gene():
    """A record file can hold a variant from an overlapping gene. ClinVar files TP53's
    antisense neighbour WRAP53 under TP53's hgnc_id, so keying RMC on the id reaches
    the wrong gene and finds nothing. collect() already derives the authoritative gene
    from the HGVS name into gene_symbol, and that is what must be keyed on.

    Verified live: gene_rmc("WRAP53") recovers o/e 0.97007 on WRAP53's own MANE
    transcript ENST00000396463 for the variant at 17:7689619 — the same figure the old
    positional route found, but reported as the gene's own region rather than a guess
    flagged on_mane=False."""
    import inspect
    from sugivariant import collect
    src = inspect.getsource(collect.attach_enrichment)
    assert "EN.gene_rmc(rec.get(\"gene_symbol\")" in src, (
        "RMC is keyed on something other than the variant's true gene symbol")


def test_gene_rmc_caches_per_gene_so_a_file_pays_once(monkeypatch):
    """The whole point of the gene-level route is one fetch per gene instead of one
    positional lookup per missense variant (measured 0.36 upstream calls each)."""
    from sugivariant import enrich as EN
    calls = {"map": 0, "mane": 0}

    def fake_map(ident, chain, cap=None):
        calls["map"] += 1
        return [{"id": "E:1", "transcript": "E", "oe": "0.5", "constrained": "true"}]

    monkeypatch.setattr(EN, "map_all", fake_map)
    monkeypatch.setattr(EN, "mane_select",
                        lambda g: (calls.__setitem__("mane", calls["mane"] + 1),
                                   {"enst": "E"})[1])
    monkeypatch.setattr(EN, "_interval_attrs",
                        lambda i, ds, c: {"start": 1, "end": 1000, "oe_chisq_p": "0.01"})
    cache = {}
    for _ in range(50):
        regions = EN.gene_rmc("TP53", cache)
        assert EN.rmc_for("17:500", regions)["oe"] == 0.5
    assert calls["map"] == 1, f"the gene was fetched {calls['map']} times, not once"
    assert calls["mane"] == 1, f"MANE resolved {calls['mane']} times, not once"
    # a second gene in the same file costs exactly one more
    for _ in range(20):
        EN.gene_rmc("WRAP53", cache)
    assert calls["map"] == 2, calls


# ── AlphaMissense abstention: 'ambiguous' is not a vote ──────────────────────
# AlphaMissense makes a three-way call. Collapsing it to `== "likely_pathogenic"`
# turned every 'ambiguous' into a tolerated vote, so pages claimed "all 3 predictors
# call this tolerated" when the ACMG-weighted tool had explicitly declined to call —
# on ~2,700 pathogenic/LP variants that also carried a severity-3 disagreement flag.
_AMB = {"class": "ambiguous", "score": "0.5367", "short": "R520L", "uniprot": "Q96G74"}
_REVEL_BENIGN = {"score": 0.249, "band": "benign-leaning", "direction": "benign"}
_SAPROT_TOL = {"llr": -5.52, "damaging": False}


def _conc_ambiguous(**kw):
    from sugivariant.enrich import concordance
    args = dict(classification="Likely pathogenic", am=_AMB, gnomad=None,
                revel=_REVEL_BENIGN, saprot=_SAPROT_TOL)
    args.update(kw)
    return concordance(args.pop("classification"), args.pop("am"), args.pop("gnomad"), **args)


def test_an_ambiguous_alphamissense_is_never_counted_as_a_tolerated_vote():
    c = _conc_ambiguous()
    cons = c["consensus"]
    assert cons["total"] == 2, f"AlphaMissense abstained but was counted: {cons}"
    assert "AlphaMissense" not in cons["names"]
    assert cons["abstained"] == ["AlphaMissense"]
    assert "all 3" not in cons["summary"]


def test_nothing_agrees_with_a_call_alphamissense_did_not_make():
    joined = " | ".join(_conc_ambiguous()["lines"])
    assert "agrees with AlphaMissense" not in joined, joined
    assert "differs from AlphaMissense" not in joined, joined
    assert "no call either way" in joined


def test_the_consensus_says_which_predictors_it_is_speaking_for():
    # "both predictors" would be a claim about all three; "both scoring predictors"
    # is a claim only about the two that spoke, and the abstention is named.
    c = _conc_ambiguous()
    summary = c["consensus"]["summary"]
    # it must name the predictors it speaks for, and say which one sat out — "both
    # scoring predictors call this tolerated" was true but left the reader to work out
    # which two, from a card of up to thirteen rows.
    assert "REVEL and SaProt" in summary, summary
    # It SCORED the variant ('ambiguous'); it is indeterminate, not silent.
    assert "AlphaMissense is indeterminate" in summary, summary
    assert "AlphaMissense makes no call" not in summary, summary
    assert "AlphaMissense" not in summary.split(";")[0], "an abstainer is not a speaker"


def test_a_single_remaining_predictor_is_not_a_consensus():
    # AM ambiguous + REVEL indeterminate leaves only SaProt. One tool is not agreement.
    c = _conc_ambiguous(revel={"score": 0.547, "band": "indeterminate",
                               "direction": "indeterminate"})
    assert c["consensus"] is None, c["consensus"]


def test_an_abstention_downgrades_the_clinvar_disagreement_and_names_the_tool():
    from sugivariant.enrich import disagreement_flag
    from sugivariant.enrich import refresh_concordance
    rec = {"classification": "Likely pathogenic", "alphamissense": _AMB,
           "revel": _REVEL_BENIGN, "saprot": _SAPROT_TOL}
    refresh_concordance(rec)          # what every render path does first
    f = disagreement_flag(rec)
    assert f["category"] == "predictor_vs_clinvar"
    assert f["severity"] == 2, "the ACMG-weighted tool abstained — this is weaker evidence"
    assert "AlphaMissense ambiguous" in f["reason"], f["reason"]


def test_a_real_likely_benign_call_still_raises_the_full_severity_flag():
    # The downgrade must apply ONLY to abstentions, not to genuine discordance.
    from sugivariant.enrich import disagreement_flag
    from sugivariant.enrich import refresh_concordance
    rec = {"classification": "Likely pathogenic",
           "alphamissense": {"class": "likely_benign", "score": "0.07"},
           "revel": _REVEL_BENIGN, "saprot": _SAPROT_TOL}
    refresh_concordance(rec)
    f = disagreement_flag(rec)
    assert f["severity"] == 3, f
    assert "ambiguous" not in f["reason"]


def test_conservation_cannot_agree_with_an_absent_or_opposite_call():
    cons_hi = {"phylop": 7.585, "phastcons": 1.0, "gerp": 0.541}
    line = [l for l in _conc_ambiguous(conservation=cons_hi)["lines"] if "Conservation" in l][0]
    assert "agrees" not in line, line
    # against a genuine likely-benign call, a conserved position DIFFERS
    from sugivariant.enrich import concordance
    line2 = [l for l in concordance("Likely benign", {"class": "likely_benign", "score": "0.05"},
                                    None, conservation=cons_hi)["lines"]
             if "Conservation" in l][0]
    assert "differs from AlphaMissense" in line2, line2


def test_the_readout_is_recomputed_from_the_record_so_fixes_need_no_rebuild():
    # A record pickled before the fix carries the wrong readout. Rendering must not
    # trust it: refresh_concordance rebuilds it from the stored predictor fields.
    from sugivariant.enrich import refresh_concordance
    stale = {"classification": "Likely pathogenic", "alphamissense": _AMB,
             "revel": _REVEL_BENIGN, "saprot": _SAPROT_TOL,
             "concordance": {"lines": ["all 3 predictors call this tolerated"],
                             "verdict": "stale", "flags": [],
                             "consensus": {"n_damaging": 0, "total": 3, "unanimous": True,
                                           "names": "AlphaMissense, REVEL, SaProt",
                                           "summary": "all 3 predictors call this tolerated"}}}
    c = refresh_concordance(stale)
    assert c["consensus"]["total"] == 2
    assert stale["concordance"] is c, "must replace the record's readout in place"


# ── the abstention invariants, as a checker ──────────────────────────────────
# Written as a function rather than inline asserts for two reasons: the truth-table
# test below sweeps every predictor combination through it, and the mutation test
# feeds it the pre-fix output to prove the checker actually catches the regression.
def abstention_violations(c, am=None, revel=None):
    """Claims a concordance readout must never make. Returns a list of violations."""
    bad = []
    amb = (am or {}).get("class") == "ambiguous"
    rev_ind = (revel or {}).get("direction") == "indeterminate"
    joined = " | ".join(c.get("lines") or [])
    cons = c.get("consensus") or {}
    names, summary = cons.get("names") or "", cons.get("summary") or ""
    if amb:
        if "with AlphaMissense" in joined:
            bad.append("claims agreement/disagreement with a call AlphaMissense never made")
        if "AlphaMissense" in names:
            bad.append("consensus counts AlphaMissense as a voter after it abstained")
        if "AlphaMissense" not in (cons.get("abstained") or []) and cons:
            bad.append("consensus does not record the AlphaMissense abstention")
        if "conserved (agrees)" in joined:
            bad.append("conservation 'agrees' with an abstention")
    if rev_ind:
        if "REVEL" in names:
            bad.append("consensus counts REVEL after it returned indeterminate")
        if "indeterminate AlphaMissense" in joined:
            bad.append("nonsense REVEL/AlphaMissense agreement wording")
    # the counted total must equal the number of tools that actually made a call
    calls = sum([bool(am) and not amb, bool(revel) and not rev_ind, True])  # SaProt always calls
    if cons and cons.get("total", 0) > calls:
        bad.append(f"consensus total {cons['total']} exceeds the {calls} tools that called")
    if cons and cons.get("total", 0) < 2:
        bad.append("a consensus was formed from fewer than two predictors")
    return bad


_AM_CLASSES = [None,
               {"class": "likely_pathogenic", "score": "0.91"},
               {"class": "ambiguous", "score": "0.45"},
               {"class": "likely_benign", "score": "0.07"}]
_REVEL_STATES = [None,
                 {"score": 0.9, "band": "strongly pathogenic-leaning", "direction": "pathogenic"},
                 {"score": 0.55, "band": "indeterminate", "direction": "indeterminate"},
                 {"score": 0.2, "band": "benign-leaning", "direction": "benign"}]
_SAPROT_STATES = [None, {"llr": -12.4, "damaging": True}, {"llr": -5.5, "damaging": False}]


def test_no_predictor_combination_can_produce_a_false_agreement_claim():
    """Sweep every AlphaMissense x REVEL x SaProt state against every ClinVar call."""
    from sugivariant.enrich import concordance
    checked = 0
    for cls in ("Pathogenic", "Likely pathogenic", "Uncertain significance",
                "Conflicting classifications of pathogenicity", "Benign", None):
        for am in _AM_CLASSES:
            for revel in _REVEL_STATES:
                for saprot in _SAPROT_STATES:
                    c = concordance(cls, am, None, None,
                                    {"phylop": 7.5, "phastcons": 1.0}, revel, saprot)
                    v = abstention_violations(c, am, revel)
                    assert not v, (f"cls={cls} am={am} revel={revel} saprot={saprot}\n"
                                   + "\n".join(v) + "\n" + "\n".join(c["lines"]))
                    checked += 1
    assert checked == 6 * 4 * 4 * 3


def test_the_checker_catches_the_pre_fix_output():
    """A mutation test on the guard itself. This is verbatim what the page used to
    render for /otud5-p-arg515leu; if the checker ever stops flagging it, the test
    above has become decorative."""
    pre_fix = {"lines": ["AlphaMissense **ambiguous** (0.5367)",
                         "REVEL 0.249 (benign-leaning) — agrees with AlphaMissense",
                         "SaProt LLR -5.52 (tolerated) — agrees with AlphaMissense",
                         "→ **In-silico consensus:** all 3 predictors call this tolerated "
                         "(AlphaMissense, REVEL, SaProt)",
                         "Conservation: phyloP 7.585, GERP 0.541, phastCons 1 — conserved (agrees)"],
               "consensus": {"n_damaging": 0, "total": 3, "unanimous": True,
                             "names": "AlphaMissense, REVEL, SaProt",
                             "summary": "all 3 predictors call this tolerated"}}
    v = abstention_violations(pre_fix, _AMB, _REVEL_BENIGN)
    assert len(v) >= 4, v
    assert any("never made" in x for x in v)
    assert any("as a voter" in x for x in v)
    assert any("conservation 'agrees'" in x for x in v)


def test_when_both_calibrated_tools_abstain_the_page_says_which():
    from sugivariant.enrich import concordance
    c = concordance("Likely pathogenic", {"class": "ambiguous", "score": "0.5521"}, None, None,
                    {"phylop": 5.529, "phastcons": 1.0},
                    {"score": 0.547, "band": "indeterminate", "direction": "indeterminate"},
                    {"llr": -5.197, "damaging": False})
    assert c["consensus"] is None, "one remaining tool is not a consensus"
    # the verdict must not go silent: a missing verdict reads as missing DATA, which is
    # a different claim from a tool declining to call.
    assert c["verdict"] == ("AlphaMissense and REVEL are indeterminate here — "
                            "the available evidence is described below"), c["verdict"]


def test_a_non_missense_record_still_gets_a_verdict():
    # `abstained` is built inside the AlphaMissense branch; the verdict chain reads it
    # unconditionally, so a record with no AlphaMissense at all must not trip over it.
    from sugivariant.enrich import concordance
    c = concordance("Pathogenic", None, None, conservation={"phylop": 5.0})
    assert c["verdict"] and "no call" not in c["verdict"]


def test_the_alphamissense_class_bands_are_the_published_ones():
    """We take `am_class` from upstream rather than deriving it, so this pins the
    boundary we rely on: ambiguous spans 0.34-0.564 (AlphaMissense, Cheng et al. 2023).
    Measured over 600 genes / 93,231 scored missense: likely_benign <= 0.3400,
    ambiguous 0.3400-0.5639, likely_pathogenic >= 0.5640 — disjoint, no overlap.
    If upstream ever re-bands, the score/class pairs stop agreeing and this fails."""
    from sugivariant.enrich import concordance
    for score, cls, expect_vote in ((0.20, "likely_benign", True),
                                    (0.45, "ambiguous", False),
                                    (0.80, "likely_pathogenic", True)):
        c = concordance("Uncertain significance", {"class": cls, "score": str(score)}, None,
                        None, None, _REVEL_BENIGN, _SAPROT_TOL)
        counted = "AlphaMissense" in (c["consensus"] or {}).get("names", "")
        assert counted is expect_vote, f"{cls} ({score}) counted={counted}"


def test_every_disagreement_flag_call_site_refreshes_the_readout_first():
    """disagreement_flag() classifies whatever readout the record carries, so each
    caller must refresh it first. index.py was missed on the first pass, and it is the
    one site that PERSISTS its answer (variant.flag → the /set worklist and the
    similar-variants ⚠), where a stale category survives on disk uncorrectable."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parent.parent
    offenders = []
    # refresh_record() is the entry point; refresh_concordance() alone also satisfies the
    # requirement, since it is what disagreement_flag actually reads.
    for f in list(root.glob("*.py")) + list((root / "sugivariant").glob("*.py")):
        src = f.read_text()
        calls = [ln for ln in src.splitlines()
                 if "disagreement_flag(" in ln and "def disagreement_flag" not in ln
                 and not ln.strip().startswith("#")]
        if calls and not ("refresh_record" in src or "refresh_concordance" in src):
            offenders.append(f"{f.name}: calls disagreement_flag but never refreshes")
    assert not offenders, "\n".join(offenders)


def test_the_page_and_the_md_both_name_an_abstaining_predictor():
    """Surface parity for the abstention: whatever each surface shows, neither may
    report agreement with a call AlphaMissense did not make."""
    rec = next(_records(limit_genes=400,
                        want=lambda r: (r.get("alphamissense") or {}).get("class") == "ambiguous"
                        and (r.get("revel") or {}).get("direction") not in (None, "indeterminate")
                        and r.get("saprot")), None)
    if rec is None:
        pytest.skip("no ambiguous-AlphaMissense record in the sampled genes")
    html, md, _, _ = _both(rec)
    for name, surface in (("HTML", html), (".md", md)):
        assert "agrees with AlphaMissense" not in surface, f"{name}: false agreement claim"
        assert "differs from AlphaMissense" not in surface, f"{name}: false disagreement claim"
        assert "all 3 predictors" not in surface, f"{name}: counts an abstaining tool"
    assert "ambiguous" in html
    # the HTML tile must disclose which tool went quiet, not just the surviving tally
    assert "AlphaMissense indeterminate" in html, "the consensus tile hides the abstention"
    assert "no call either way" in md


# ── rule: SaProt goes dark for whole PROTEINS, and that must be visible ──
# SaProt is structure-aware, and AlphaFold DB's human proteome release excluded
# sequences longer than 2,700 residues. So 182 genes with scored missense -- TTN,
# BRCA2, ATM, APC, RYR1, DMD -- have NO SaProt score for any variant (289,227
# variants, 11.9% of the AlphaMissense-covered set). Before this, the row simply
# vanished, the consensus never named it, and the note under the card still told the
# reader SaProt was one of the agreement signals. Found on brca2-p-arg2336his.
def test_a_dark_saprot_is_named_as_abstaining():
    from sugivariant.enrich import concordance
    c = concordance("Pathogenic", {"class": "likely_pathogenic", "score": "0.91"}, None,
                    revel={"score": 0.8, "band": "damaging", "direction": "pathogenic"},
                    saprot=None)
    cons = c["consensus"]
    assert "SaProt" in cons["abstained"], f"a silent SaProt passed unmentioned: {cons}"
    assert "SaProt makes no call" in cons["summary"], cons["summary"]
    # and it must not be counted as a consenting voice
    assert cons["total"] == 2 and "SaProt" not in cons["names"]


def test_a_scoring_saprot_is_not_reported_as_abstaining():
    from sugivariant.enrich import concordance
    c = concordance("Pathogenic", {"class": "likely_pathogenic", "score": "0.91"}, None,
                    revel={"score": 0.8, "band": "damaging", "direction": "pathogenic"},
                    saprot={"llr": -12.0, "damaging": True})
    cons = c["consensus"]
    assert cons["abstained"] == [] and cons["total"] == 3
    assert "no call" not in cons["summary"], cons["summary"]


def test_three_quiet_predictors_still_read_as_a_sentence():
    """With SaProt named, all three can go quiet at once; the verdict used
    " and ".join, which produced "A and B and C"."""
    from sugivariant.enrich import concordance
    c = concordance("Likely pathogenic", _AMB, None,
                    revel={"score": 0.5, "band": "indeterminate",
                           "direction": "indeterminate"}, saprot=None)
    # Two reasons, kept apart: AM and REVEL returned scores that point neither way,
    # SaProt returned nothing at all.
    assert c["verdict"].startswith(
        "AlphaMissense and REVEL are indeterminate and SaProt makes no call here"), c["verdict"]


def test_the_page_shows_a_saprot_row_even_when_saprot_is_dark():
    """An absent predictor must be stated, like SpliceAI's snv_only/absent rows —
    not omitted, which is indistinguishable from a page that never looked."""
    rec = next(_records(limit_genes=600,
                        want=lambda r: (r.get("alphamissense") or {}).get("class")
                        and not r.get("saprot")), None)
    if rec is None:
        pytest.skip("every sampled AlphaMissense-scored record has a SaProt score")
    html, _md, _g, _p = _both(rec)
    assert "SaProt" in html, "a structureless protein drops SaProt off the page entirely"
    assert "not scored" in html, "SaProt's absence is not stated"
    assert "2,700 residues" in html, "the page does not say WHY SaProt is silent"
    # and the note must not claim a tool that never spoke is an agreement signal
    assert "REVEL and SaProt are independent agreement signals" not in html, \
        "the note credits a SaProt that made no call"


def test_the_note_still_credits_saprot_when_it_did_score():
    rec = next(_records(limit_genes=600,
                        want=lambda r: (r.get("alphamissense") or {}).get("class")
                        and r.get("saprot") and r.get("revel")), None)
    if rec is None:
        pytest.skip("no record with AlphaMissense, REVEL and SaProt in the sample")
    html, _md, _g, _p = _both(rec)
    assert "REVEL and SaProt are independent agreement signals" in html
    assert "not scored" not in html.split("Splicing")[0], "a scored SaProt rendered as absent"


@pytest.mark.skipif(not os.environ.get("SUGI_CORPUS_INVARIANTS"),
                    reason="corpus sweep; set SUGI_CORPUS_INVARIANTS=1")
def test_no_cached_record_renders_a_false_agreement_claim():
    """The unit tests above cover constructed inputs. This one sweeps real records,
    which is how the bug reached production in the first place: every input involved
    was individually plausible and no constructed case had been written for it."""
    from sugivariant.enrich import refresh_concordance
    checked = bad = 0
    for r in _records(limit_genes=1200):
        am = r.get("alphamissense") or {}
        if am.get("class") != "ambiguous":
            continue
        checked += 1
        v = abstention_violations(refresh_concordance(r), am, r.get("revel"))
        if v:
            bad += 1
            if bad <= 5:
                print(f"\n/{r.get('canonical_slug')}: " + "; ".join(v))
    assert checked, "no ambiguous-AlphaMissense records found — sweep proved nothing"
    assert not bad, f"{bad}/{checked} records violate the abstention invariants"


# ── "absent from gnomAD" is only sayable when we looked and found nothing ─────
# Absence is an ACMG PM2-supporting rarity claim that propagates into the meta
# description and the JSON-LD FAQ, so it must never be produced by missing data.
# Two shapes did exactly that, on 22.6% and 3.5% of pathogenic/LP records.
def test_a_gnomad_row_with_coverage_but_no_allele_count_is_not_an_absence():
    from sugivariant.enrich import refresh_frequency
    # uniform shape across all 1,602 measured cases: AN set, af/ac/popmax all absent
    rec = {"gnomad": {"source": "gnomAD v4.1", "absent": False, "af": None, "ac": None,
                      "popmax": None, "an": 1605142, "band": "absent from gnomAD v4.1"}}
    g = refresh_frequency(rec)
    assert g["unknown"] is True
    assert g["absent"] is False
    assert "absent" not in g["band"], g["band"]
    assert "1,605,142" in g["band"], "the coverage we do know should still be reported"


def test_an_empty_dbsnp_frequency_string_is_not_a_zero():
    from sugivariant.enrich import refresh_frequency
    g = refresh_frequency({"gnomad": {"source": "dbSNP/gnomAD", "frequency": "",
                                      "absent": True, "band": "absent from gnomAD"}})
    assert g["unknown"] is True and g["absent"] is False
    assert "absent" not in g["band"], g["band"]


def test_an_explicit_dbsnp_zero_is_still_a_real_absence():
    """The fix must not throw away a genuine zero: "0" is dbSNP reporting a frequency,
    "" is dbSNP reporting nothing. Only the second is missing data."""
    from sugivariant.enrich import refresh_frequency
    g = refresh_frequency({"gnomad": {"source": "dbSNP/gnomAD", "frequency": "0",
                                      "absent": True, "band": "absent from gnomAD"}})
    assert g["absent"] is True and g.get("unknown") is False


def test_a_real_coordinate_miss_on_an_snv_still_claims_absence():
    """A reliable key that found nothing IS an absence — the PM2 signal must survive."""
    from sugivariant.enrich import refresh_frequency
    g = refresh_frequency({"gnomad": {"source": "gnomAD v4.1", "absent": True,
                                      "af": None, "ac": None, "popmax": None, "an": None,
                                      "band": "absent from gnomAD v4.1"}})
    assert g["absent"] is True and "absent from gnomAD" in g["band"]


def test_a_real_frequency_is_left_untouched():
    from sugivariant.enrich import refresh_frequency
    band = "rare (popmax 0.12%, nfe)"
    g = refresh_frequency({"gnomad": {"source": "gnomAD v4.1", "absent": False, "af": 0.0012,
                                      "ac": 18, "an": 1500000, "popmax": 0.0012,
                                      "ancestry": "nfe", "band": band}})
    assert g["band"] == band and g["unknown"] is False and g["absent"] is False


def test_no_surface_calls_a_missing_frequency_an_absence():
    """Both surfaces, on a real record carrying the no-allele-count shape."""
    rec = next(_records(limit_genes=500, want=lambda r: (
        (r.get("gnomad") or {}).get("source") == "gnomAD v4.1"
        and (r.get("gnomad") or {}).get("absent") is False
        and not (r.get("gnomad") or {}).get("af")
        and (r.get("gnomad") or {}).get("popmax") is None
        and (r.get("gnomad") or {}).get("an"))), None)
    if rec is None:
        pytest.skip("no no-allele-count record in the sampled genes")
    html, md, g, _ = _both(rec)
    assert g["unknown"] is True
    for name, surface in (("HTML", html), (".md", md)):
        assert "Absent from gnomAD" not in surface, f"{name} claims absence"
        assert "absent from gnomAD" not in surface, f"{name} claims absence"
        assert "PM2-supporting" not in surface or "cannot be assessed" in surface, \
            f"{name} offers the PM2 rarity signal off a missing number"
    # and the block must still RENDER — hiding it would lose the disclosure entirely
    assert "no allele count reported" in html
    assert "cannot be assessed" in md


# ── phyloP bands: a negative score is not "moderately conserved" ─────────────
def test_a_negative_phylop_is_never_called_conserved():
    """One threshold (p >= 2) drove three different claims: everything in (-2, 2) read
    "moderately conserved position" on the non-missense branch — where conservation is
    the single ACMG-weighted signal — including phyloP -0.8 on 2,193 pathogenic/LP
    records whose GERP and phastCons both said the opposite."""
    from sugivariant.enrich import concordance, conservation_note
    for p in (-0.096, -0.808, -1.9, -2.5, 0.0):
        for weighted in (True, False):
            note = conservation_note(p, weighted)
            assert "conserved position" not in note or "not a conserved" in note, (p, note)
            assert "moderately conserved" not in note, (p, note)
    # and on a real non-missense readout, which is where it carried ACMG weight
    line = [l for l in concordance("Pathogenic", None, None, None,
                                   {"phylop": -0.808, "gerp": -2.24, "phastcons": 0.003})["lines"]
            if "Conservation" in l][0]
    assert "not a conserved position" in line, line
    assert "PP3" not in line, "a non-conserved position must not offer PP3-type support"


def test_the_two_conservation_branches_agree_on_the_same_value():
    """phyloP 1.9 read "low conservation" on a missense page and "moderately conserved"
    on a non-missense one. Both now come from one band definition."""
    from sugivariant.enrich import conservation_note
    for p in (1.9, 0.5, -0.5, 3.0, -3.0):
        a, b = conservation_note(p, True), conservation_note(p, False)
        # the only licensed difference is the PP3-type note on the weighted branch
        assert a.replace(" (PP3-type support)", "").replace("position", "").strip() == \
               b.replace("position", "").strip() or p >= 2, (p, a, b)


def test_splice_agreement_is_refreshed_from_the_stored_blocks():
    from sugivariant.enrich import refresh_record
    rec = {"spliceai": {"status": "absent"},
           "pangolin": {"status": "scored", "effect": "loss", "score": 0.35},
           "splice_agreement": {"state": "spliceai_below", "text": "SpliceAI below 0.2"}}
    refresh_record(rec)
    assert rec["splice_agreement"]["state"] == "spliceai_absent"
    assert "0.2" not in rec["splice_agreement"]["text"]


# ── "we did not look" is not "the release has nothing" ───────────────────────
def test_an_unloaded_spliceai_cache_is_not_reported_as_absence():
    """spliceai_for returned `absent` when no per-gene index had been built, turning
    "we did not look" into a claim about the release. BRCA1 c.4357+2T>A — a pathogenic
    canonical donor whose real readout is donor_loss Δ0.99 — rendered as "not assessed"."""
    from sugivariant.enrich import spliceai_for
    assert spliceai_for("17:43082402:A:T", {})["status"] == "not_loaded"
    assert spliceai_for("17:43082402:A:T", None)["status"] == "not_loaded"
    # a gene genuinely carrying no rows still yields the 3-key dict, so it is NOT not_loaded
    loaded_but_empty = {"by_allele": {}, "by_pos": None, "masked_loaded": False}
    assert spliceai_for("17:43082402:A:T", loaded_but_empty)["status"] != "not_loaded"


def test_an_unloaded_cache_makes_no_splice_agreement_claim():
    from sugivariant.enrich import splice_agreement
    pg = {"status": "scored", "effect": "loss", "score": 0.35}
    assert splice_agreement({"status": "not_loaded"}, pg) is None


def test_both_surfaces_say_we_did_not_look_rather_than_not_assessed():
    import os as _os
    _os.environ.setdefault("BASE_PATH", "/variant")
    import app
    from sugivariant.render import render_body
    from sugivariant.enrich import refresh_record, disagreement_flag
    rec = next(_records(limit_genes=60, want=lambda r: r.get("spliceai")), None)
    if rec is None:
        pytest.skip("no record with a spliceai block")
    rec = dict(rec)
    rec["spliceai"] = {"status": "not_loaded"}
    rec.pop("_sag_fresh", None); rec.pop("_conc_fresh", None)
    refresh_record(rec)
    g, pops = app._freq_ctx(rec)
    html = app._render("variant.html", v=rec, canonical=rec.get("canonical_slug"), nav="variant",
                       jsonld="", freq_g=g, freq_pops=pops)
    md = render_body(rec)
    assert "not looked up on this page" in html
    assert "not looked up on this page" in md
    for surface in (html, md):
        assert "no allele at this position is in our SpliceAI release" not in surface


def test_no_surface_lowercases_a_dataset_or_field_name():
    """Jinja's |capitalize and str.capitalize() lowercase the REST of the string, which
    turned "gnomAD v4.1 (AN 1,208,295)" into "gnomad v4.1 (an ...)" and, earlier,
    "SpliceAI" into "Spliceai". Both read as sloppy on a reference page."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parent.parent
    bad = [f"{f.name}: |capitalize" for f in (root / "templates").glob("*.html")
           if "|capitalize" in f.read_text()]
    assert not bad, bad + ["use the ucfirst filter instead"]
    rec = next(_records(limit_genes=400, want=lambda r: (r.get("gnomad") or {}).get("an")
                        and not (r.get("gnomad") or {}).get("af")
                        and (r.get("gnomad") or {}).get("popmax") is None), None)
    if rec is None:
        pytest.skip("no no-allele-count record in the sampled genes")
    html, md, _, _ = _both(rec)
    for name, surface in (("HTML", html), (".md", md)):
        for wrong in ("gnomad v4.1", "(an ", "Spliceai", "Alphamissense", "Clinvar"):
            assert wrong not in surface, f"{name} contains {wrong!r}"


# ── colour must encode the stance toward ClinVar, not internal unanimity ─────
def test_a_unanimous_tolerated_call_against_a_pathogenic_one_is_not_styled_as_agreement():
    """`unanimous` means the predictors agree WITH EACH OTHER. Styling on it alone put
    the agreement colour on "both scoring predictors call this tolerated" beneath a
    Likely pathogenic classification — green on exactly the page where the evidence
    conflicts, and the one a reader most needs to read as a conflict."""
    from sugivariant.enrich import concordance
    c = concordance("Likely pathogenic", _AMB, None, None, None, _REVEL_BENIGN, _SAPROT_TOL)
    cons = c["consensus"]
    assert cons["unanimous"] is True, "the predictors do agree with each other"
    assert cons["stance"] == "discordant", "...but they contradict the clinical call"


def test_stance_covers_every_clinvar_direction():
    from sugivariant.enrich import concordance
    dmg = ({"class": "likely_pathogenic", "score": "0.97"},
           {"score": 0.9, "band": "strongly pathogenic-leaning", "direction": "pathogenic"},
           {"llr": -12.0, "damaging": True})
    tol = (_AMB, _REVEL_BENIGN, _SAPROT_TOL)
    cases = [("Pathogenic", dmg, "concordant"), ("Pathogenic", tol, "discordant"),
             ("Benign", tol, "concordant"), ("Benign", dmg, "discordant"),
             ("Uncertain significance", dmg, "neutral"),
             ("Conflicting classifications of pathogenicity", tol, "neutral")]
    for cls, (am, rv, sp), want in cases:
        got = concordance(cls, am, None, None, None, rv, sp)["consensus"]["stance"]
        assert got == want, f"{cls} + {'damaging' if am is not _AMB else 'tolerated'}: {got} != {want}"


def test_the_rendered_page_does_not_paint_a_discordant_consensus_green():
    rec = next(_records(limit_genes=500, want=lambda r: (
        (r.get("alphamissense") or {}).get("class") == "ambiguous"
        and "athogenic" in (r.get("classification") or "")
        and "onflict" not in (r.get("classification") or "")
        # REVEL must actually SPEAK: with AlphaMissense abstaining, an indeterminate
        # REVEL leaves SaProt alone and a consensus of one is correctly not rendered.
        # Without this the test picked whichever record the sample happened to yield and
        # failed on correct behaviour (pigo-p-leu957phe) after the record set changed.
        and (r.get("revel") or {}).get("direction") not in (None, "indeterminate")
        and r.get("saprot"))), None)
    if rec is None:
        pytest.skip("no ambiguous-AM pathogenic record in the sampled genes")
    html, _, _, _ = _both(rec)
    import re
    # the consensus now renders inside the predictor group it summarises, not as a
    # paragraph below the whole card
    m = re.search(r'class="vd-consensus([^"]*)">(.*?)</td>', html, re.S)
    assert m, "consensus row not rendered"
    klass, text = m.group(1), m.group(2)
    if "tolerated" in text:
        assert "is-agree" not in klass, f"discordant consensus painted as agreement: {klass!r}"
        assert "is-flag" in klass, klass


def test_every_record_persist_site_guards_against_collapse():
    """Any code that writes a gene's .pkl.gz must consult records_collapsed() first.

    `if recs:` only rejects an EMPTY rebuild. Staging is missing the hgnc>>clinvar link
    for ~0.5% of genes, which returns ONE record for a gene with hundreds — SOX10 (537
    stored) and PALS1 (91) were written that way on 2026-10-10. The write is
    unrecoverable: the good records are gone, and because the resolver indexes the
    pickles, the next index build drops those pages to 404. The guard has to sit at the
    write, not the caller, so a new build script inherits it."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parent.parent
    offenders = []
    for f in list(root.glob("*.py")) + list((root / "sugivariant").glob("*.py")):
        src = f.read_text()
        writes = [ln for ln in src.splitlines()
                  if "pickle.dump(" in ln and not ln.strip().startswith("#")]
        if writes and "records_collapsed" not in src:
            offenders.append(f"{f.name}: pickle.dump of records without a collapse guard")
    assert not offenders, "\n".join(offenders)


# ── rule: a WITHHELD SaProt score must not read as a score that never existed ──
# SaProt is keyed on AlphaMissense's residue number. For ~170 accessions the two were
# computed on different isoforms, so past the divergence point our key lands on the
# wrong residue -- and still RESOLVES when the shifted position carries the same
# reference letter. MEN1/O00255: aligned to 149, +5 from 150 on. 8,055 such scores were
# live, and 523 triage flags were built on them. Found via brca2-p-arg2336his.
def _men1_like():
    """A gene whose accession is scored on a 5-shorter isoform: aligned to 149,
    coincidental hits after it."""
    recs = []
    def rec(pos, scored, slug):
        r = {"canonical_slug": slug, "gene_symbol": "MEN1",
             "alphamissense": {"class": "likely_benign", "score": "0.1",
                               "uniprot": "O00255", "short": f"G{pos}E"}}
        if scored:
            r["saprot"] = {"llr": -7.0, "damaging": False}
        return r
    for pos in (2, 100, 149):
        recs.append(rec(pos, True, f"men1-aligned-{pos}"))      # real scores
    for pos in (150, 151, 300):
        recs.append(rec(pos, False, f"men1-miss-{pos}"))        # key missed
    for pos in (162, 391):
        recs.append(rec(pos, True, f"men1-coincidence-{pos}"))  # WRONG residue
    return recs


def test_a_coincidental_saprot_hit_past_the_divergence_is_withheld():
    from sugivariant.enrich import drop_misaligned_saprot
    recs = _men1_like()
    n = drop_misaligned_saprot(recs)
    assert n == 2, f"expected the two coincidental hits withheld, got {n}"
    by = {r["canonical_slug"]: r for r in recs}
    for slug in ("men1-coincidence-162", "men1-coincidence-391"):
        assert by[slug].get("saprot") is None, f"{slug} kept a wrong-residue LLR"
        assert by[slug]["saprot_misaligned"] == 150, "the divergence residue is not recorded"
    # the aligned prefix is CORRECT and must survive -- this is not a blanket drop
    for slug in ("men1-aligned-2", "men1-aligned-100", "men1-aligned-149"):
        assert by[slug]["saprot"] is not None, f"{slug} lost a valid score"


def test_a_fully_covered_accession_is_untouched():
    from sugivariant.enrich import drop_misaligned_saprot
    recs = [r for r in _men1_like() if r.get("saprot")]      # hits only, no misses
    assert drop_misaligned_saprot(recs) == 0
    assert all(r.get("saprot") for r in recs)


def test_a_fully_dark_accession_is_untouched():
    """No hits at all (no AlphaFold model) -- there is no divergence to infer and
    nothing to withhold; the page already says 'not scored'."""
    from sugivariant.enrich import drop_misaligned_saprot
    recs = [r for r in _men1_like() if not r.get("saprot")]
    assert drop_misaligned_saprot(recs) == 0
    assert not any(r.get("saprot_misaligned") for r in recs)


def test_misses_above_the_top_hit_do_not_invent_a_divergence():
    """A miss past the end of the scored sequence is ordinary coverage, not a
    mismatch, and must not withhold the scores below it."""
    from sugivariant.enrich import drop_misaligned_saprot
    recs = _men1_like()
    recs = [r for r in recs
            if int(r["alphamissense"]["short"][1:-1]) <= 149
            or not r.get("saprot")]
    recs = [r for r in recs if int(r["alphamissense"]["short"][1:-1]) not in (150, 151)]
    # now: hits at 2/100/149, a single miss at 300 which is ABOVE the top hit
    assert drop_misaligned_saprot(recs) == 0
    assert all(r.get("saprot") for r in recs
               if int(r["alphamissense"]["short"][1:-1]) <= 149)


def test_the_page_says_withheld_not_absent():
    """The two reasons are different claims and must not render the same."""
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = next(_records(limit_genes=600,
                        want=lambda r: (r.get("alphamissense") or {}).get("class")
                        and not r.get("saprot")), None)
    if rec is None:
        pytest.skip("no AlphaMissense-scored record without SaProt in the sample")
    rec = dict(rec)
    rec["saprot_misaligned"] = 150
    html, _md, _g, _p = _both(rec)
    assert "withheld" in html, "a withheld score renders as if it never existed"
    assert "numbering disagrees past residue 150" in html
    assert "not scored" not in html.split("Splicing")[0], \
        "withheld and never-scored must not render identically"


# ── rule: a gene-level figure must be keyed on the VARIANT's gene, not ClinVar's id ──
# ClinVar files an overlapping gene's variant under the HOST gene's hgnc_id. Keying MANE
# Select on hgnc_id therefore returned the host gene's transcript AND then labelled the
# variant's real MANE as the deviation: a pathogenic ATM frameshift read "MANE Select
# NM_152587 · reported on NM_000051", where NM_152587 is C11orf65's and NM_000051 IS
# ATM's. 59,805 pages stated it backwards, 3,669 of them P/LP. gene_rmc was fixed for
# this exact reason ("ClinVar's hgnc_id can be the overlapping host gene's");
# mane_select was missed.
def test_mane_select_is_keyed_on_the_variants_own_gene():
    import pathlib
    src = (pathlib.Path(__file__).resolve().parent.parent
           / "sugivariant" / "collect.py").read_text(encoding="utf-8")
    assert 'EN.mane_select(rec.get("hgnc_id"))' not in src, \
        "MANE keyed on ClinVar's hgnc_id returns the HOST gene's transcript"
    assert 'rec.get("gene_symbol") or rec.get("hgnc_id")' in src, \
        "MANE must resolve from the variant's own gene symbol first"


def test_mane_is_cached_per_gene_not_per_record_file():
    """One record file holds the host gene's variants AND the overlapping genes', so a
    single file-level MANE is wrong for every gene in it but the first."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parent.parent
           / "sugivariant" / "collect.py").read_text(encoding="utf-8")
    assert 'if "mane" not in ctx' not in src, \
        "a file-level MANE cache gives the first gene's transcript to all the others"
    assert "mane_cache" in src


def test_every_gene_level_lookup_keys_on_the_gene_not_the_hgnc_id():
    """The general rule behind both: any per-gene fetch reached from a record must
    prefer the variant's own symbol. Catches the next one added."""
    import pathlib, re
    src = (pathlib.Path(__file__).resolve().parent.parent
           / "sugivariant" / "collect.py").read_text(encoding="utf-8")
    offenders = [ln.strip() for ln in src.splitlines()
                 if re.search(r'EN\.\w+\(\s*rec\.get\("hgnc_id"\)\s*\)', ln)
                 and not ln.strip().startswith("#")]
    assert not offenders, (
        "gene-level lookups keyed on ClinVar's hgnc_id, which can be the overlapping "
        "host gene's: " + "; ".join(offenders))


# ── rule: a bare SYMBOL must never be handed to a `>>hgnc` chain ──
# The chain resolves a symbol by SEARCH, and HGNC ranks PREVIOUS symbols first, so the
# first hit is routinely another gene: MET returns HGNC:20709 before MET's own
# HGNC:7029, BAP1 returns HGNC:950 only third, OPA1 returns HGNC:11957 before
# HGNC:8140. Fixing the MANE keying by passing gene_symbol straight into mane_select
# reintroduced the FLG/FGFR1 bug on 88,939 records before this guard was added.
def test_a_symbol_is_resolved_to_its_own_hgnc_id_before_any_chain():
    import pathlib, re
    src = (pathlib.Path(__file__).resolve().parent.parent
           / "sugivariant" / "enrich.py").read_text(encoding="utf-8")
    bad = [ln.strip() for ln in src.splitlines()
           if re.search(r'map_all\(\s*(gene|symbol)\s*,\s*">>hgnc', ln)
           and not ln.strip().startswith("#")]
    assert not bad, ("a bare symbol reaches a >>hgnc chain, which resolves previous "
                     "symbols first: " + "; ".join(bad))


def test_hgnc_for_symbol_prefers_the_approved_symbol():
    """Pins the behaviour rather than the implementation: an id passes through, and a
    symbol whose search ranks another gene first still resolves to its own."""
    from sugivariant.enrich import hgnc_for_symbol
    assert hgnc_for_symbol("HGNC:7029") == "HGNC:7029", "an id must pass through"
    assert hgnc_for_symbol(None) is None
    assert hgnc_for_symbol("") is None


@pytest.mark.skipif(not os.environ.get("CACHE_DIR"), reason="needs upstream")
def test_the_colliding_symbols_resolve_to_their_own_gene():
    from sugivariant.enrich import hgnc_for_symbol
    for sym, want in (("MET", "HGNC:7029"), ("BAP1", "HGNC:950"),
                      ("OPA1", "HGNC:8140"), ("FLG", "HGNC:3748")):
        got = hgnc_for_symbol(sym)
        if got is None:
            pytest.skip("upstream unavailable")
        assert got == want, f"{sym} resolved to {got}, not its own {want}"


# ── rule: the machine-readable surface must not outrun the corrections ──
# refresh_frequency() strips a stored gnomAD band that asserts an absence the data cannot
# support. The JSON-LD was built ABOVE that call, so 187,493 pages (84,456 P/LP) served
# HTML saying "no allele count reported in gnomAD v4.1" beside structured data saying "is
# absent from gnomAD v4.1" -- manufacturing PM2-supporting rarity on the surface search
# engines and language models quote.
def test_the_jsonld_is_built_after_the_record_is_refreshed():
    import pathlib, re
    src = (pathlib.Path(__file__).resolve().parent.parent / "app.py").read_text(encoding="utf-8")
    body = src[src.index("def variant_page"):]
    body = body[:body.index("\n@app.") if "\n@app." in body else len(body)]
    i_refresh = body.find("_refresh(rec)")
    i_jsonld = body.find("as_script_tag(")
    assert i_refresh != -1 and i_jsonld != -1, "variant_page no longer has both calls"
    assert i_refresh < i_jsonld, \
        "JSON-LD is built before _refresh(rec): it will carry an absence claim the page retracted"


def test_benign_likely_benign_is_not_coloured_as_uncertain():
    """The only ClinVar call that merges benign and likely-benign evidence fell through
    cls_class's default and rendered in the VUS grey on 73,495 variants."""
    os.environ.setdefault("BASE_PATH", "/variant")
    import app
    assert app.cls_class("Benign/Likely benign") == "ben"
    assert app.cls_class("Pathogenic/Likely pathogenic") == "path"
    assert app.cls_class("Uncertain significance") == "vus"


def test_no_surface_still_claims_a_pathogenic_only_corpus():
    """The build gate now produces VUS and benign too. Telling a reader their variant has
    'no pathogenic-gated page' invites a clinical inference from a stale build note."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parent.parent / "templates"
    bad = []
    for f in root.glob("*.html"):
        if "disagreements" in f.name:
            continue                      # dead templates; the routes 308-redirect
        t = f.read_text(encoding="utf-8")
        for phrase in ("pathogenic-gated",
                       "Pathogenic / Likely-pathogenic / Conflicting only"):
            if phrase in t:
                bad.append(f"{f.name}: {phrase}")
    assert not bad, "stale corpus-scope claims: " + "; ".join(bad)


def test_the_md_discloses_saprots_dark_states_like_the_html():
    """96317f5 fixed this for the HTML and left the .md silent: 97,584 AM-scored pages
    said nothing about SaProt while the footer called it an agreement signal, and 7,794
    withheld scores were flattened into 'makes no call'."""
    rec = next(_records(limit_genes=600,
                        want=lambda r: (r.get("alphamissense") or {}).get("class")
                        and not r.get("saprot")), None)
    if rec is None:
        pytest.skip("no AlphaMissense-scored record without SaProt in the sample")
    html, md, _g, _p = _both(rec)
    assert "SaProt" in html
    assert "SaProt" in md, "the .md does not mention SaProt at all"
    if rec.get("saprot_misaligned"):
        assert "withheld" in md, "a withheld score must not read as an abstention"
    else:
        assert "not scored" in md


# ── rule: a page built without its gene's context must say so ──
# _collect_one() sets rec["_partial"] for a variant ClinVar files under an overlapping
# gene, so the page is built with ctx={} and loses gene constraint, dosage, structure,
# pathways, panels, hotspot, similar-variants and the AlphaMissense percentile. The flag
# was set and read NOWHERE: 67,186 pages served the thinner set silently. SpliceAI
# already discloses its own absence this way.
def test_a_partial_page_discloses_what_was_not_looked_up():
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = next(_records(limit_genes=40), None)
    if rec is None:
        pytest.skip("no cached records")
    rec = dict(rec)
    rec["_partial"] = True
    html, md, _g, _p = _both(rec)
    for name, surface in (("HTML", html), (".md", md)):
        assert "without this gene's context" in surface, f"{name}: no disclosure"
        assert "not looked up" in surface, f"{name}: does not say what is missing"


def test_a_normal_page_carries_no_such_note():
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = next(_records(limit_genes=40), None)
    if rec is None:
        pytest.skip("no cached records")
    rec = dict(rec)
    rec.pop("_partial", None)
    html, md, _g, _p = _both(rec)
    assert "without this gene's context" not in html
    assert "without this gene's context" not in md


# ── rule: an OUTAGE must not be read as the evidence that triggers withholding ──
# drop_misaligned_saprot infers the divergence point from which positions missed, and D
# is the LOWEST miss. So one failed lookup inside the aligned prefix withholds every
# correct score at or above it, for the whole accession. biobtree staging restarted
# during the 2026-10-14 alias recovery (18,982 lookups) and surfaced that saprot_for was
# flattening a transport failure to the same None a genuine miss produces.
def test_a_failed_saprot_lookup_does_not_withhold_the_accession():
    from sugivariant.enrich import drop_misaligned_saprot, FETCH_ERROR
    recs = _men1_like()
    by = {r["canonical_slug"]: r for r in recs}
    # the miss at 150 -- the one that sets D -- was really an outage, not a mismatch
    by["men1-miss-150"]["saprot"] = FETCH_ERROR
    assert drop_misaligned_saprot(recs) == 0, \
        "a transient fetch failure was mistaken for the divergence point"
    for slug in ("men1-aligned-2", "men1-aligned-100", "men1-aligned-149",
                 "men1-coincidence-162", "men1-coincidence-391"):
        assert by[slug].get("saprot") is not None, f"{slug} withheld on an outage"
        assert not by[slug].get("saprot_misaligned")


def test_the_fetch_error_sentinel_is_never_persisted():
    """It informs the decision; it must not reach a pickle, where every reader
    expects a dict or None."""
    from sugivariant.enrich import drop_misaligned_saprot, FETCH_ERROR
    recs = _men1_like()
    recs[0]["saprot"] = FETCH_ERROR
    drop_misaligned_saprot(recs)
    for r in recs:
        assert r.get("saprot") is not FETCH_ERROR, "the sentinel leaked into a record"
        assert r.get("saprot") is None or isinstance(r["saprot"], dict)


def test_saprot_for_distinguishes_a_failed_fetch_from_a_genuine_miss():
    from sugivariant import enrich as EN
    real = EN._coord_entry
    try:
        EN._coord_entry = lambda *a, **k: EN.FETCH_ERROR
        assert EN.saprot_for("O00255", "G150E") is EN.FETCH_ERROR, \
            "a failed fetch must not return the None a genuine miss returns"
        EN._coord_entry = lambda *a, **k: None
        assert EN.saprot_for("O00255", "G150E") is None, \
            "a genuine negative must stay None"
    finally:
        EN._coord_entry = real


def test_scoping_the_divergence_pass_leaves_other_accessions_exactly_as_they_were():
    """An INCREMENTAL repair re-fetches only the accessions that gained scores. A group
    whose scores are already withheld presents them as misses, so an unscoped re-run
    would read its own withholding as fresh coverage evidence and ratchet D downward on
    groups nothing changed for."""
    from sugivariant.enrich import drop_misaligned_saprot
    recs = _men1_like()
    drop_misaligned_saprot(recs)                       # first, full pass
    settled = {r["canonical_slug"]: (r.get("saprot"), r.get("saprot_misaligned"))
               for r in recs}
    # a second accession, untouched by the repair, is added to the same file
    other = {"canonical_slug": "other-1", "gene_symbol": "MEN1",
             "alphamissense": {"class": "likely_benign", "score": "0.1",
                               "uniprot": "P99999", "short": "G7E"},
             "saprot": {"llr": -3.0, "damaging": False}}
    recs.append(other)
    assert drop_misaligned_saprot(recs, only={"P99999"}) == 0
    for r in recs:
        if r["canonical_slug"] == "other-1":
            continue
        assert (r.get("saprot"), r.get("saprot_misaligned")) == settled[r["canonical_slug"]], \
            f"{r['canonical_slug']} moved although it was out of scope"


# ── rule: a predictor the page SHOWS a score for is not a predictor that stayed silent ──
# The consensus tile sits directly above the evidence card, so a reader sees both at
# once. "REVEL no call" over a card row reading "REVEL 0.636 (indeterminate)" is a
# contradiction on one screen, and only one of the two can be true. Three states, not
# two: scored-and-directional (votes), scored-but-indeterminate (named, does not vote),
# unscored (named, does not vote). Found by the user on ush2a-p-pro1490his, where
# AlphaMissense + SaProt read 2/2 damaging and the tile called REVEL silent.
def test_an_indeterminate_predictor_is_not_described_as_silent():
    rec = next(_records(limit_genes=400,
                        want=lambda r: (r.get("revel") or {}).get("direction") == "indeterminate"
                        and (r.get("alphamissense") or {}).get("class")
                        in ("likely_pathogenic", "likely_benign")
                        and r.get("saprot")), None)
    if rec is None:
        pytest.skip("no indeterminate-REVEL record in the sampled genes")
    html, md, _, _ = _both(rec)
    score = str((rec["revel"])["score"])
    for name, surface in (("HTML", html), (".md", md)):
        if score not in surface:
            continue                      # that surface does not show the value at all
        assert "REVEL no call" not in surface, \
            f"{name}: tile calls REVEL silent while showing REVEL {score}"
        assert "REVEL makes no call" not in surface, \
            f"{name}: readout calls REVEL silent while showing REVEL {score}"


def test_the_two_reasons_a_predictor_does_not_vote_are_phrased_differently():
    """One tool returned a score pointing neither way, the other returned nothing.
    Both are out of the tally; they are not the same claim, and the readout that
    mentions them must not merge them into one."""
    from sugivariant.enrich import concordance
    c = concordance("Pathogenic", {"class": "ambiguous", "score": "0.44"}, None,
                    revel={"score": 0.636, "band": "indeterminate",
                           "direction": "indeterminate"},
                    saprot=None)
    assert "AlphaMissense and REVEL are indeterminate" in c["verdict"], c["verdict"]
    assert "SaProt makes no call" in c["verdict"], c["verdict"]
    assert "REVEL makes no call" not in c["verdict"], c["verdict"]
    # and REVEL's actual value is still on the card, which is why calling it silent
    # would contradict what the reader can see
    assert any("REVEL 0.636" in ln for ln in c["lines"]), c["lines"]


def test_a_lone_voting_predictor_does_not_claim_the_others_were_silent():
    """AlphaMissense alone agreeing is reported as one predictor, and the indeterminate
    REVEL beside it keeps its value on the card rather than being called silent."""
    from sugivariant.enrich import concordance
    c = concordance("Pathogenic", {"class": "likely_pathogenic", "score": "0.91"}, None,
                    revel={"score": 0.636, "band": "indeterminate",
                           "direction": "indeterminate"},
                    saprot=None)
    assert c["consensus"] is None, "one voting tool is not a consensus"
    assert any("REVEL 0.636" in ln for ln in c["lines"]), c["lines"]
    joined = " ".join(c["lines"] + [c.get("verdict") or ""])
    assert "REVEL makes no call" not in joined, joined


def test_an_indeterminate_and_an_unscored_tool_are_both_still_named():
    """Splitting the phrasing must not drop either from the readout."""
    from sugivariant.enrich import concordance
    c = concordance("Pathogenic", {"class": "likely_pathogenic", "score": "0.91"}, None,
                    revel={"score": 0.5, "band": "indeterminate",
                           "direction": "indeterminate"},
                    saprot={"llr": -12.0, "damaging": True})
    cons = c["consensus"]
    assert cons["total"] == 2 and cons["n_damaging"] == 2
    assert cons["indeterminate"] == ["REVEL"] and cons["unscored"] == []
    assert set(cons["abstained"]) == {"REVEL"}, "the union must still list it"
    assert "REVEL is indeterminate" in cons["summary"], cons["summary"]


# ── rule: a dataset gates rebuilds only if a record actually reads it ──
# RECORD_DATASETS is stamped beside each gene so an upstream release rebuilds only the
# genes it affects. A name in it that nothing fetches inverts that: it makes every gene
# look stale for free. `gtex` sat there unfetched, and when upstream moved db_v52 ->
# db_v53 it marked all 18,810 genes for a ~90h rebuild over data no record holds.
def test_record_datasets_are_all_actually_fetched():
    import pathlib
    import re
    from sugivariant.enrich import RECORD_DATASETS
    root = pathlib.Path(__file__).resolve().parent.parent
    src = "\n".join(p.read_text() for p in (root / "sugivariant").glob("*.py"))
    src += (root / "app.py").read_text()
    # Drop the declaration itself, or every name "finds" itself -- and drop COMMENTS,
    # or the note explaining why gtex was removed counts as a use of gtex. A dataset
    # mentioned only in prose is exactly the case this test exists to catch.
    src = re.sub(r"RECORD_DATASETS = \([^)]*\)", "", src, flags=re.S)
    src = "\n".join(ln.split("#", 1)[0] for ln in src.splitlines())
    unused = [g for g in RECORD_DATASETS if not re.search(re.escape(g), src)]
    assert not unused, (
        "RECORD_DATASETS names a dataset the package never fetches, so every gene "
        f"reads as stale on it for nothing: {unused}")


def test_no_evidence_table_tooltip_sits_in_a_cell_the_containment_rule_misses():
    """A tooltip in the evidence table must be anchored to a cell, not to the word.

    .vd-preds-wrap is a scroll container (overflow-x:auto, so the table can scroll
    instead of the page on a narrow phone). An absolutely positioned bubble whose
    containing block is INSIDE that wrapper is both clipped by it and counted in its
    scroll width -- so a 22rem tip opened from a word near the right-hand end of the
    wide third column stuck out past the card and put a horizontal scrollbar under a
    card that otherwise fits. Found on jak3-c-2351-minus-15a-g.

    app.css pins the bubble to `td:not(.vd-num)` and caps it at the cell's width. A
    tooltip added to a vd-num cell would fall outside that rule and bring the scrollbar
    back, silently -- so this checks the MARKUP, which is what can drift, as well as
    that the rule is still there. The header column is deliberately exempt: it sits at
    the card's left edge, where a left-anchored bubble has room.
    """
    import pathlib
    import re
    css = (pathlib.Path(__file__).resolve().parent.parent
           / "static" / "css" / "app.css").read_text()
    for sel in (".vd-preds td:not(.vd-num){position:relative}",
                ".vd-preds td:not(.vd-num) .ti{position:static}",
                ".vd-preds td:not(.vd-num) .ti::after{max-width:min(22rem,100%)}"):
        assert sel in css, f"the evidence-table tooltip containment rule is gone: {sel}"

    checked = 0
    for rec in _records(limit_genes=12, want=lambda r: r.get("spliceai") or r.get("cadd")):
        html, _, _, _ = _both(rec)
        if 'class="vd-preds-wrap"' not in html:
            continue
        table = html.split('class="vd-preds-wrap"', 1)[1].split("</table>", 1)[0]
        for attrs, body in re.findall(r"<td([^>]*)>(.*?)</td>", table, re.S):
            if 'class="ti"' in body:
                checked += 1
                assert "vd-num" not in attrs, (
                    f"{rec.get('canonical_slug')}: a tooltip sits in a vd-num cell, which "
                    "the containment rule does not cover -- it will widen the wrapper's "
                    "scroll area and put a scrollbar under the card")
    if not checked:
        pytest.skip("no evidence-table tooltip in the sampled genes")


def test_a_bracketed_spliceai_absence_is_stated_as_evidence_not_as_a_gap():
    """A position inside the gene's scored span must not be called "not assessed".

    SpliceAI's release scores every possible substitution inside a gene's body and keeps
    only the alleles reaching Δ0.2, so what it holds for a gene is one contiguous interval
    sampled sparsely. A position between the lowest and highest record we hold for the
    gene was therefore evaluated, and the absence of every alternate allele puts all four
    deltas of each under 0.2 — a prediction of no splicing effect, which is evidence.

    Found on JAK3 c.2351-15A>G: the base one away carries acceptor_loss 0.77, so the
    intron was plainly scored, and the page still said "not assessed". The remaining
    `absent` state keeps its caveat and is checked alongside, because the whole point is
    that the two are now told apart rather than merged.
    """
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-2351-minus-15a-g", "gene_symbol": "X",
           "classification": "Likely benign", "variation_id": 2, "name": "X:c.2351-15A>G",
           "hgvs_c": "c.2351-15A>G",
           "review_status": "criteria provided, single submitter", "conditions": [],
           "spliceai": {"status": "below_threshold", "nearest_dp": 1,
                        "scored_positions": 1138, "span_lo": 17825373,
                        "span_hi": 17847949, "threshold": 0.2}}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "0.2" in surface, f"{name} omits the threshold that makes absence meaningful"
        assert "not assessed" not in surface, (
            f"{name} still discloses a gap for a position the release evaluated")
        assert "evaluated" in surface, f"{name} does not say the position was evaluated"
        assert "not a gap" in surface or "not missing" in surface, (
            f"{name} does not distinguish a below-threshold prediction from missing data")
        # the evidence the claim rests on, so a reader can check it rather than trust it
        assert "1 bp" in surface or "1,138" in surface, (
            f"{name} asserts the position was evaluated without showing what says so")

    rec["spliceai"] = {"status": "absent"}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "not assessed" in surface, (
            f"{name}: an unbracketed absence must still disclose the gap")
        assert "not evidence against" in surface, f"{name} omits the absence caveat"


def test_an_unknown_spliceai_status_degrades_to_the_gap_wording_not_to_silence():
    """Records live on disk; the code ships separately. Neither order may blank the row.

    Record files are rewritten in place by transform passes, and the prod image is built
    and deployed on its own schedule, so a record can carry a status the running code has
    no branch for — records upgraded before the code ships, or the code rolled back after.
    An unmatched status used to fall off the end of the chain and render NOTHING, which is
    the one thing the SpliceAI row exists to prevent: silence reads as "no splice signal".
    Fall back to the honest gap wording instead.
    """
    os.environ.setdefault("BASE_PATH", "/variant")
    rec = {"canonical_slug": "x-c-1-plus-1g-c", "gene_symbol": "X", "classification": "VUS",
           "variation_id": 3, "name": "X:c.1+1G>C", "hgvs_c": "c.1+1G>C",
           "review_status": "criteria provided, single submitter", "conditions": [],
           "spliceai": {"status": "a_status_from_a_later_build"}}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "SpliceAI" in surface, f"{name} drops the SpliceAI row for an unknown status"
        assert "not assessed" in surface, f"{name} renders silence for an unknown status"

    # a legacy record — a score and no status at all — must NOT be swept into that
    # fallback: it HAS a score, so "not assessed" would be the lie the chain prevents
    rec["spliceai"] = {"effect": "donor_loss", "score": "0.99"}
    html, md, _, _ = _both(rec)
    assert "not assessed" not in md, ".md calls a legacy scored record unassessed"
    assert "not assessed" not in html, "HTML calls a legacy scored record unassessed"
