"""Phase 1a: the app must answer Atlas URLs exactly as the static site does.

These are URL-semantics tests, not content tests: trailing slash, 301 table, 404, the
.md twin, ETag/304 and path traversal. They are the gate that has to pass before any
traffic or any renderer is involved.
"""
import json

import pytest
from starlette.testclient import TestClient

import app as A
from sugivariant import atlas_web as ATW


@pytest.fixture()
def site(tmp_path, monkeypatch):
    root = tmp_path / "atlas"
    (root / "gene" / "BRCA1").mkdir(parents=True)
    (root / "gene" / "BRCA1" / "index.html").write_text(
        "<html><body><article class=\"entity-page\">BRCA1 PAGE</article></body></html>")
    (root / "gene" / "BRCA1" / "index.md").write_text("# BRCA1\n\nmarkdown body\n")
    (root / "disease" / "x-syndrome").mkdir(parents=True)
    (root / "disease" / "x-syndrome" / "index.html").write_text("<html>X</html>")
    (root / "index.html").write_text("<html>atlas home</html>")
    (root / "about").mkdir()
    (root / "about" / "index.html").write_text("<html>about</html>")

    rules = tmp_path / "redirects.json"
    rules.write_text(json.dumps({"rules": [
        {"from": "/atlas/disease/old-name/", "to": "/atlas/disease/x-syndrome/"},
        {"from": "/atlas/disease/old-name", "to": "/atlas/disease/x-syndrome/"},
    ]}))
    monkeypatch.setattr(ATW, "ROOT", root)
    monkeypatch.setattr(ATW, "REDIRECTS", rules)
    monkeypatch.setattr(ATW, "ENABLED", True)
    monkeypatch.setattr(ATW, "_redirects", None)
    return TestClient(A.app)


def test_entity_page_is_served_with_an_etag(site):
    r = site.get("/atlas/gene/BRCA1/")
    assert r.status_code == 200
    assert "BRCA1 PAGE" in r.text
    assert r.headers["content-type"].startswith("text/html")
    assert r.headers["etag"]
    assert "s-maxage" in r.headers["cache-control"]


def test_missing_entity_is_a_real_404(site):
    assert site.get("/atlas/gene/NOPE/").status_code == 404
    assert site.get("/atlas/drug/nope/").status_code == 404


def test_trailing_slash_is_canonical(site):
    r = site.get("/atlas/gene/BRCA1", follow_redirects=False)
    assert r.status_code == 301
    assert r.headers["location"].endswith("/atlas/gene/BRCA1/")
    assert site.get("/atlas", follow_redirects=False).headers["location"].endswith("/atlas/")


def test_markdown_twin_is_text_markdown(site):
    r = site.get("/atlas/gene/BRCA1/index.md")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/markdown")
    assert r.text.startswith("# BRCA1")


def test_rename_table_301s_to_its_successor(site):
    for path in ("/atlas/disease/old-name/", "/atlas/disease/old-name"):
        r = site.get(path, follow_redirects=False)
        assert r.status_code == 301
        assert r.headers["location"].endswith("/atlas/disease/x-syndrome/")


def test_etag_revalidation_returns_304(site):
    etag = site.get("/atlas/gene/BRCA1/").headers["etag"]
    r = site.get("/atlas/gene/BRCA1/", headers={"If-None-Match": etag})
    assert r.status_code == 304


def test_path_traversal_is_refused(site):
    for path in ("/atlas/../../etc/passwd", "/atlas/gene/../../app.py"):
        assert site.get(path).status_code in (404, 400)


def test_disabled_deployment_does_not_serve_atlas(monkeypatch):
    monkeypatch.setattr(ATW, "ENABLED", False)
    assert TestClient(A.app).get("/atlas/gene/BRCA1/").status_code == 404


def test_site_level_pages_work(site):
    assert site.get("/atlas/about/").status_code == 200
    assert "about" in site.get("/atlas/about/").text
    assert site.get("/atlas/").status_code == 200


# ── Phase 1b: the pickle corpus ────────────────────────────────────────────────
def _corpus(tmp_path, monkeypatch):
    import gzip
    import pickle
    import sqlite3
    cache = tmp_path / "cache"
    (cache / "gene").mkdir(parents=True)
    rec = {
        "schema": "sugi-atlas-entity/1", "type": "gene", "slug": "BRCA1",
        "title": "BRCA1 gene", "canonical": "https://sugi.bio/atlas/gene/BRCA1/",
        "body_html": '<article class=entity-page>BRCA1 <a href=/atlas/gene/TP53/>TP53</a>'
                     ' <a href=https://example.org>x</a></article>',
        "body_md": "# BRCA1\n", "sections": [{"id": "summary", "title": "Summary"}],
        "frag_sha": "x", "md_sha": "y", "bytes": 120,
    }
    with gzip.open(cache / "gene" / "BRCA1.pkl.gz", "wb") as f:
        pickle.dump(rec, f)
    db = tmp_path / "index.db"
    conn = sqlite3.connect(db)
    conn.executescript(
        "CREATE TABLE entity(type TEXT, slug TEXT, path TEXT, title TEXT, canonical TEXT,"
        " frag_sha TEXT, md_sha TEXT, bytes INTEGER, sections TEXT, PRIMARY KEY(type,slug));")
    conn.execute("INSERT INTO entity VALUES(?,?,?,?,?,?,?,?,?)",
                 ("gene", "BRCA1", "gene/BRCA1.pkl.gz", rec["title"], rec["canonical"],
                  "x", "y", 120, "[]"))
    conn.commit()
    conn.close()
    monkeypatch.setattr(ATW, "CACHE", cache)
    monkeypatch.setattr(ATW, "INDEX", db)
    monkeypatch.setattr(ATW, "_lru", {})
    monkeypatch.setattr(ATW, "_lru_bytes", 0)
    monkeypatch.setattr(ATW, "_tls", __import__("threading").local())
    return rec


def test_entity_record_loads_from_the_pickle(tmp_path, monkeypatch):
    rec = _corpus(tmp_path, monkeypatch)
    assert ATW.cache_enabled() is True
    assert ATW.has_entity("gene", "BRCA1") is True
    assert ATW.has_entity("gene", "NOPE") is False
    assert ATW.entity("gene", "BRCA1")["body_html"] == rec["body_html"]


def test_rewrite_prefixes_only_internal_atlas_links():
    html = '<a href=/atlas/gene/TP53/>t</a> <a href="https://x/atlas/">e</a> <a href=#s>s</a>'
    out = ATW.rewrite_links(html, "/variant")
    assert "href=/variant/atlas/gene/TP53/" in out
    assert "https://x/atlas/" in out          # external untouched
    assert "href=#s" in out                   # anchor untouched
    assert ATW.rewrite_links(html, "") == html


def test_excerpt_strips_markup():
    rec = {"body_html": "<h2>Summary</h2><p>BRCA1 is a <b>tumour</b> suppressor.</p>",
           "title": "t"}
    assert ATW.excerpt(rec) == "Summary BRCA1 is a tumour suppressor."
