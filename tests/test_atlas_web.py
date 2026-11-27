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
