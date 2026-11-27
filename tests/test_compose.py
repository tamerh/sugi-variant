"""Phase 3: composing /gene/X from the two corpora.

The mechanism under test is narrow on purpose: remove the Atlas page furniture the
composed page supplies itself, append the Atlas markdown as a section, and prove the
Atlas prose itself is untouched (the browser gate,
`migration/verify_gene_composed.py`, compares the post-JS `.entity-body` DOM).
"""
from __future__ import annotations

import gzip
import pickle

import pytest

import app as A
from sugivariant import atlas_corpus as AC
from sugivariant import atlas_web as ATW
from sugivariant import compose
from sugivariant import render as VR

FRAGMENT = (
    '<article class="entity-page" data-entity-type="gene">'
    '<a class="atlas-feedback" href="mailto:x">Send feedback</a>'
    '<div class="container container-wide">'
    '<p class=breadcrumb><a href=/atlas/>Atlas</a> / <a href=/atlas/gene/>Gene</a> / BRCA1</p>'
    '<header class=entity-header><h1 class=entity-title>BRCA1</h1>'
    '<span class=entity-kind-badge>gene</span></header>'
    '<div class=entity-layout>'
    '<aside class=entity-toc-aside><p>On this page</p><nav id=TableOfContents>'
    '<ul><li><a href=#summary>Summary</a></li></ul></nav></aside>'
    '<button type=button id=entity-toc-reopen>≡</button>'
    '<div class=entity-main>'
    '<details class=entity-toc-inline><summary>On this page</summary></details>'
    '<div class=entity-body>'
    '<h2 id=summary>Summary</h2><p>BRCA1 is a <b>tumour</b> suppressor.</p>'
    '<table><tr><th>A</th></tr><tr><td>1</td></tr></table>'
    '</div></div></div>'
    '<footer class=entity-footer><div class=entity-footer-about>'
    '<h3>About this page</h3><p>Built by Sugi Atlas.</p></div></footer>'
    '</div></article>'
)

ATLAS_MD = (
    "# BRCA1\n\n"
    "> Gene · canonical: https://sugi.bio/atlas/gene/BRCA1/\n"
    "> Atlas v1.16.0 · generated 2026-11-26T18:30:44+00:00\n\n"
    "## Summary\n\nBRCA1 is a tumour suppressor.\n\n"
    "### Function\n\nDNA repair.\n"
)


def test_strip_removes_the_furniture_and_keeps_the_prose():
    out = compose.strip_page_furniture(FRAGMENT)
    for gone in ("breadcrumb", "entity-header", "entity-kind-badge",
                 "entity-toc-aside", "entity-toc-inline", "entity-toc-reopen",
                 "container container-wide"):
        assert gone not in out, gone
    # The content and its provenance stay exactly as built.
    assert "class=entity-body" in out
    assert "BRCA1 is a <b>tumour</b> suppressor." in out
    assert out.count("<table>") == 1
    assert "entity-footer" in out and "About this page" in out


def test_demote_markdown_drops_the_h1_and_nests_the_sections():
    out = compose.demote_markdown(ATLAS_MD)
    assert not out.startswith("# BRCA1")
    assert out.startswith("> Gene · canonical:")
    assert "\n## Summary" not in out
    assert "\n### Summary" in out
    assert "\n#### Function" in out
    # Provenance lines are kept: they name the build the prose came from.
    assert "Atlas v1.16.0 · generated" in out


@pytest.fixture()
def corpus(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    (cache / "gene").mkdir(parents=True)
    rec = {"schema": "sugi-atlas-entity/1", "type": "gene", "slug": "BRCA1",
           "title": "BRCA1 gene", "canonical": "https://sugi.bio/atlas/gene/BRCA1/",
           "body_html": FRAGMENT, "body_md": ATLAS_MD, "sections": [],
           "frag_sha": "a", "md_sha": "b", "bytes": len(FRAGMENT)}
    with gzip.open(cache / "gene" / "BRCA1.pkl.gz", "wb") as f:
        pickle.dump(rec, f)
    import sqlite3
    index = tmp_path / "atlas_index.db"
    conn = sqlite3.connect(index)
    conn.executescript(
        "CREATE TABLE entity(type TEXT, slug TEXT, path TEXT, title TEXT, canonical TEXT,"
        " frag_sha TEXT, md_sha TEXT, bytes INTEGER, sections TEXT, PRIMARY KEY(type,slug));"
        "CREATE TABLE alias(key TEXT, type TEXT, slug TEXT, name TEXT);"
        "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);")
    conn.execute("INSERT INTO entity VALUES(?,?,?,?,?,?,?,?,?)",
                 ("gene", "BRCA1", "gene/BRCA1.pkl.gz", "BRCA1 gene",
                  rec["canonical"], "a", "b", len(FRAGMENT), "[]"))
    conn.commit()
    conn.close()
    keys = ("ROOT", "REDIRECTS", "CACHE", "INDEX", "ENABLED", "RELEASE",
            "RELEASE_MODE", "_lru", "_lru_bytes", "_tls", "_redirects", "_sha_cache")
    saved = {k: getattr(ATW, k) for k in keys}
    ATW.configure(AC.describe(cache, index))
    yield
    for k, v in saved.items():
        setattr(ATW, k, v)


def test_gene_block_returns_the_stripped_atlas_half(corpus):
    blk = compose.gene_block("BRCA1")
    assert blk is not None
    assert "entity-header" not in blk["html"]
    assert "class=entity-body" in blk["html"]
    assert blk["md"].startswith("> Gene · canonical:")
    assert compose.gene_block("NOPE") is None


def test_composed_template_has_one_identity_and_the_atlas_assets(corpus):
    ctx = {
        "gene": "BRCA1", "total": 2, "groups": [], "flagged": 0, "counts": {},
        "dgroups": {}, "cats": [], "meta": {}, "dcounts": {}, "sv_rows": [],
        "sv_total": 0, "cap": 200, "base": "", "nav": "variant",
        "atlas": compose.gene_block("BRCA1"),
    }
    html = A.env.get_template("gene_composed.html").render(**ctx)
    assert html.count("<h1") == 1
    assert "Gene context" in html
    assert "atlas-site.css" in html          # the fragment's own stylesheet
    assert "atlas-entity.js" in html         # the runtime enhancement
    assert "entity-header" not in html       # no second identity
    assert "BRCA1 is a <b>tumour</b> suppressor." in html


def test_gene_md_gains_the_context_section():
    ctx = {"gene": "BRCA1", "total": 2, "groups": [], "flagged": 0, "counts": {},
           "dgroups": {}, "cats": [], "meta": {}, "dcounts": {}, "sv_rows": [],
           "sv_total": 0, "cap": 200}
    md = VR.render_gene_md(ctx, atlas_md=compose.demote_markdown(ATLAS_MD))
    assert md.count("\n# ") == 0 and md.startswith("# BRCA1")
    assert "## Gene context {#gene-context}" in md
    assert "### Summary" in md
    assert A._md_section(md, "Gene context") is not None
    # Without the context the twin is the hub's, unchanged.
    plain = VR.render_gene_md(ctx)
    assert "Gene context" not in plain
