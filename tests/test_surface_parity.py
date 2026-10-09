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
        assert "not scored for this allele" in surface, f"{name} hides the missing allele"
        assert "0.99" in surface, f"{name} omits the same-position score"
        assert "not evidence against" in surface, f"{name} omits the absence caveat"

    rec["spliceai"] = {"status": "absent"}
    html, md, _, _ = _both(rec)
    for surface, name in ((html, "HTML"), (md, ".md")):
        assert "not assessed" in surface, f"{name} renders silence instead of 'not assessed'"


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
