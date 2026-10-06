"""Render every cached record through BOTH output surfaces and assert nothing raises.

This is the regression gate for the Jinja-Undefined bug class, which has caused two
production 500 outages: a missing dict key is `Undefined`, and `Undefined is not none`
evaluates TRUE, so a guard like `{% if g.popmax is not none %}` walks straight into
`format(Undefined * 100)` and raises. The records cache is never auto-invalidated, so
records on disk predate current fields and the shapes genuinely drift (commit b598202).

Unit tests over pure functions cannot catch this — only rendering a real record can.
The cache is optional (CI may not have it), so the test skips when it is absent.
"""
import glob
import gzip
import os
import pickle
import random

import pytest

CACHE = os.environ.get("RECORDS_DIR") or "/data2/sugi-variant-cache/records"
MAX_GENES = int(os.environ.get("RENDER_TEST_GENES", "40"))     # keep it quick by default


def _records():
    files = sorted(glob.glob(os.path.join(CACHE, "*.pkl.gz")))
    if not files:
        pytest.skip(f"no cached records at {CACHE}")
    random.Random(7).shuffle(files)
    for f in files[:MAX_GENES]:
        try:
            with gzip.open(f, "rb") as fh:
                recs = pickle.load(fh)
        except Exception:
            continue
        for r in recs:
            yield os.path.basename(f), r


def test_every_cached_record_renders_both_surfaces():
    os.environ.setdefault("BASE_PATH", "/variant")
    import app
    from sugivariant.enrich import disagreement_flag
    from sugivariant.render import render_body

    n = 0
    failures = []
    for gene, rec in _records():
        n += 1
        slug = rec.get("canonical_slug")
        try:
            app._render("variant.html", v=rec, canonical=slug, nav="variant",
                        disagreement=disagreement_flag(rec))
        except Exception as e:
            failures.append(f"HTML {gene} {slug}: {type(e).__name__}: {e}")
        try:
            render_body(rec)
        except Exception as e:
            failures.append(f"MD   {gene} {slug}: {type(e).__name__}: {e}")
        if len(failures) >= 5:
            break

    assert n, "no records found to render"
    assert not failures, f"{len(failures)} render failure(s) over {n} records:\n" + "\n".join(failures)


def test_gnomad_fallback_shape_renders_a_frequency():
    """The dbSNP fallback has no popmax/af. It must still show its frequency, and must
    NOT be labelled grpmax (a global frequency is not a group max) — see a320452."""
    os.environ.setdefault("BASE_PATH", "/variant")
    import re

    import app
    from sugivariant.enrich import disagreement_flag

    row = re.compile(r"gnomAD v4\.1</th><td[^>]*>([^<]*)</td><td[^>]*>([^<]*)")
    checked = 0
    for _gene, rec in _records():
        g = rec.get("gnomad")
        if not (isinstance(g, dict) and g.get("source") == "dbSNP/gnomAD" and not g.get("absent")):
            continue
        html = app._render("variant.html", v=rec, canonical=rec.get("canonical_slug"),
                           nav="variant", disagreement=disagreement_flag(rec))
        m = row.search(html)
        if not m:
            continue
        value, label = m.group(1).strip(), m.group(2).strip()
        assert value != "—", f"{rec.get('canonical_slug')}: frequency present but rendered as a dash"
        assert "grpmax" not in label, f"{rec.get('canonical_slug')}: global frequency mislabelled {label!r}"
        checked += 1
        if checked >= 25:
            break
    if not checked:
        pytest.skip("no fallback-shape records in the sampled genes")


def test_faf95_is_surfaced_and_is_the_named_ba1_bs1_input():
    """FAF95 (grpmax filtering AF) is ClinGen's designated BA1/BS1 input per the gnomAD
    v4 guidance v3.0 (June 2025). We stored it as `faf` and displayed only faf99; both
    surfaces must now show FAF95, including for records cached before the `faf95` key
    existed (the dict shape on disk is never auto-migrated)."""
    os.environ.setdefault("BASE_PATH", "/variant")
    import re

    import app
    from sugivariant.enrich import disagreement_flag
    from sugivariant.render import render_body

    row = re.compile(r"gnomAD v4\.1</th><td[^>]*>[^<]*</td><td[^>]*>(.*?)</td>", re.S)
    checked = 0
    for _gene, rec in _records():
        g = rec.get("gnomad") or {}
        faf95 = g.get("faf95") if g.get("faf95") is not None else g.get("faf")
        if not isinstance(g, dict) or g.get("absent") or faf95 is None:
            continue
        html = app._render("variant.html", v=rec, canonical=rec.get("canonical_slug"),
                           nav="variant", disagreement=disagreement_flag(rec))
        m = row.search(html)
        if m:
            assert "FAF95" in m.group(1), f"{rec.get('canonical_slug')}: FAF95 missing from the row"
        body = render_body(rec)
        if "Population frequency" in body:
            assert "FAF95" in body, f"{rec.get('canonical_slug')}: FAF95 missing from the .md"
        checked += 1
        if checked >= 20:
            break
    if not checked:
        pytest.skip("no records with a published FAF in the sampled genes")


# ── population-frequency block ────────────────────────────────────────────────
def _freq_render(rec):
    import app
    from sugivariant.enrich import disagreement_flag
    g, pops = app._freq_ctx(rec)
    return app._render("variant.html", v=rec, canonical=rec.get("canonical_slug"),
                       nav="variant", jsonld="", disagreement=disagreement_flag(rec),
                       freq_g=g, freq_pops=pops), g, pops


def _mk(gnomad):
    return {"canonical_slug": "x-p-y", "gene_symbol": "X", "classification": "VUS",
            "review_status": "criteria provided, single submitter", "conditions": [],
            "gnomad": gnomad}


def test_freq_block_hides_the_table_for_a_single_ancestry_group():
    """A one-row ancestry table only repeats the "Highest group" figure above it, and a
    third of variants report exactly one group."""
    os.environ.setdefault("BASE_PATH", "/variant")
    html, _, pops = _freq_render(_mk({
        "popmax": 5.93e-06, "ancestry": "nfe", "af": 4.34e-06, "ac": 7, "an": 1614062,
        "faf95": 2.47e-06, "faf99": 1.59e-06, "band": "very rare", "absent": False,
        "source": "gnomAD v4.1", "populations": {"nfe": 5.93e-06}}))
    assert len(pops) == 1
    assert 'class="pf"' in html and "Highest group" in html
    assert "pf-tbl" not in html, "single-ancestry variants must not render a one-row table"


def test_freq_block_shows_the_table_for_multiple_groups():
    os.environ.setdefault("BASE_PATH", "/variant")
    html, _, pops = _freq_render(_mk({
        "popmax": 5.36e-04, "ancestry": "eas", "af": 1.69e-05, "ac": 27, "an": 1601778,
        "ac_grpmax": 24, "an_grpmax": 44812, "faf95": 3.69e-04, "faf99": 3.14e-04,
        "band": "ultra-rare", "absent": False, "source": "gnomAD v4.1",
        "populations": {"eas": 5.36e-04, "amr": 1.67e-05, "nfe": 8.56e-07}}))
    assert len(pops) == 3
    assert "pf-tbl" in html and "East Asian" in html and "European (non-Finnish)" in html
    assert "24/44,812" in html, "grpmax must carry its AC/AN"
    assert "ClinGen BA1/BS1 input" in html


def test_freq_block_states_absence_rather_than_omitting_it():
    """Absence from gnomAD is PM2-supporting evidence, so it must be stated, not blank."""
    os.environ.setdefault("BASE_PATH", "/variant")
    html, _, _ = _freq_render(_mk({"absent": True, "is_common": False, "popmax": None,
                                   "band": "absent from gnomAD v4.1", "source": "gnomAD v4.1"}))
    assert "pf-absent" in html and "Absent from gnomAD v4.1" in html
    assert "pf-tbl" not in html


def test_freq_block_handles_the_legacy_dbsnp_fallback_shape():
    """Records cached before gnomad_for() emitted `af` carry only the raw `frequency`
    string. Without deriving af from it the whole block silently disappears."""
    os.environ.setdefault("BASE_PATH", "/variant")
    html, g, _ = _freq_render(_mk({"frequency": "0.00537", "absent": False,
                                   "is_common": False, "band": "low-frequency (gnomAD MAF 0.00537)",
                                   "source": "dbSNP/gnomAD"}))
    assert g["af"] == 0.00537, "af must be derived from the legacy frequency key"
    assert "Global frequency" in html and "via dbSNP" in html
    assert "ClinGen BA1/BS1 input" not in html, "no FAF exists on the fallback path"
    assert "pf-tbl" not in html, "no ancestry breakdown exists on the fallback path"
