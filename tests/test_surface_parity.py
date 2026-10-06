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
