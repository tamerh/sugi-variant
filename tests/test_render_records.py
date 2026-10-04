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
