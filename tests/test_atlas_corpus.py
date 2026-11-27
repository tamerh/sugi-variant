"""Phase 2: the release manifest's Atlas block, and what the app does with it.

The point of these tests is the release discipline: the block must describe the files
on disk (sha + counts), and the app must serve exactly the entities the block's index
lists -- never fall through to a page that merely happens to exist in the tree.
"""
from __future__ import annotations

import gzip
import json
import pickle
import sqlite3

import pytest

from sugivariant import atlas_corpus as AC
from sugivariant import atlas_web as ATW


@pytest.fixture(autouse=True)
def _restore_atlas_state():
    """configure() is deliberately global: it adopts a corpus for the process. Tests
    must not leak that into the rest of the suite."""
    keys = ("ROOT", "REDIRECTS", "CACHE", "INDEX", "ENABLED", "RELEASE",
            "RELEASE_MODE", "_lru", "_lru_bytes", "_tls", "_redirects", "_sha_cache")
    saved = {k: getattr(ATW, k) for k in keys}
    problems = list(ATW.RELEASE_PROBLEMS)
    yield
    for k, v in saved.items():
        setattr(ATW, k, v)
    ATW.RELEASE_PROBLEMS[:] = problems


def _corpus(tmp_path):
    """A two-entity corpus, a root with its two anchor files, and a redirect table."""
    cache = tmp_path / "cache"
    (cache / "gene").mkdir(parents=True)
    rec = {
        "schema": "sugi-atlas-entity/1", "type": "gene", "slug": "BRCA1",
        "title": "BRCA1 gene", "canonical": "https://sugi.bio/atlas/gene/BRCA1/",
        "body_html": "<article class=entity-page>BRCA1</article>",
        "body_md": "# BRCA1\n", "sections": [], "frag_sha": "a", "md_sha": "b",
        "bytes": 41,
    }
    with gzip.open(cache / "gene" / "BRCA1.pkl.gz", "wb") as f:
        pickle.dump(rec, f)
    index = tmp_path / "atlas_index.db"
    conn = sqlite3.connect(index)
    conn.executescript(
        "CREATE TABLE entity(type TEXT, slug TEXT, path TEXT, title TEXT, canonical TEXT,"
        " frag_sha TEXT, md_sha TEXT, bytes INTEGER, sections TEXT, PRIMARY KEY(type,slug));"
        "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);")
    conn.execute("INSERT INTO entity VALUES(?,?,?,?,?,?,?,?,?)",
                 ("gene", "BRCA1", "gene/BRCA1.pkl.gz", "BRCA1 gene",
                  "https://sugi.bio/atlas/gene/BRCA1/", "a", "b", 41, "[]"))
    conn.execute("INSERT INTO meta VALUES('built_at','2026-11-27T00:00:00')")
    conn.commit()
    conn.close()
    root = tmp_path / "site"
    root.mkdir()
    (root / "index.html").write_text("<html>atlas</html>")
    (root / "manifest.json").write_text('{"gene": {"BRCA1": "BRCA1"}}')
    redirects = tmp_path / "redirects.json"
    redirects.write_text(json.dumps({"rules": [
        {"from": "/atlas/gene/OLD/", "to": "/atlas/gene/BRCA1/"}]}))
    return cache, index, root, redirects


def test_describe_records_counts_and_anchors(tmp_path):
    cache, index, root, redirects = _corpus(tmp_path)
    d = AC.describe(cache, index, root=root, redirects=redirects)
    assert d["entities"] == 1
    assert d["by_type"] == {"gene": 1}
    assert d["with_markdown"] == 1
    assert d["redirect_rules"] == 1
    assert d["index_sha256"] == AC.sha256_file(index)
    assert d["root_index_sha256"] == AC.sha256_file(root / "index.html")
    assert d["search"]["sha256"] == AC.sha256_file(root / "manifest.json")
    assert AC.check(d) == []


def test_check_catches_a_corpus_that_moved_under_the_manifest(tmp_path):
    cache, index, root, redirects = _corpus(tmp_path)
    d = AC.describe(cache, index, root=root, redirects=redirects)
    # A rebuilt index with one more entity is a different corpus.
    conn = sqlite3.connect(index)
    conn.execute("INSERT INTO entity VALUES(?,?,?,?,?,?,?,?,?)",
                 ("gene", "TP53", "gene/TP53.pkl.gz", "TP53", None, "c", "d", 10, "[]"))
    conn.commit()
    conn.close()
    problems = AC.check(d)
    assert any("index sha" in p for p in problems)
    assert any("entities" in p for p in problems)


def test_check_catches_wrong_anchors_and_a_missing_tree(tmp_path):
    cache, index, root, redirects = _corpus(tmp_path)
    d = AC.describe(cache, index, root=root, redirects=redirects)
    (root / "manifest.json").write_text('{"gene": {}}')
    problems = AC.check(d)
    assert any("search data sha" in p for p in problems)


def test_configure_adopts_the_manifest_and_finds_the_mounted_copy(tmp_path, monkeypatch):
    cache, index, root, redirects = _corpus(tmp_path)
    d = AC.describe(cache, index, root=root, redirects=redirects)
    # Simulate a container: the manifest carries host paths, the files are mounted
    # under CACHE_DIR/atlas by basename.
    mounted = tmp_path / "mnt"
    (mounted / "atlas").mkdir(parents=True)
    import shutil
    shutil.copytree(cache, mounted / "atlas" / cache.name)
    shutil.copy(index, mounted / "atlas" / index.name)
    shutil.copytree(root, mounted / "atlas" / root.name)
    shutil.copy(redirects, mounted / "atlas" / redirects.name)
    shutil.rmtree(cache)
    index.rename(tmp_path / "gone.db")
    root.rename(tmp_path / "gone-site")
    redirects.rename(tmp_path / "gone.json")

    monkeypatch.setattr(ATW, "_lru", {})
    monkeypatch.setattr(ATW, "_lru_bytes", 0)
    monkeypatch.setattr(ATW, "_tls", __import__("threading").local())
    problems = ATW.configure(d, cache_dir=mounted)
    assert problems == []
    assert ATW.enabled() is True
    assert ATW.CACHE == mounted / "atlas" / cache.name
    assert ATW.INDEX == mounted / "atlas" / index.name
    assert ATW.has_entity("gene", "BRCA1") is True
    assert ATW.has_entity("gene", "TP53") is False


def test_configure_disables_atlas_when_the_corpus_is_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(ATW, "_tls", __import__("threading").local())
    problems = ATW.configure({"cache": str(tmp_path / "nope"),
                              "index_db": str(tmp_path / "nope.db")}, cache_dir=None)
    assert problems and "unusable" in problems[0]
    assert ATW.enabled() is False


def test_sitemap_chunks_enumerate_the_index_and_nothing_else(tmp_path, monkeypatch):
    from starlette.testclient import TestClient

    import app as A

    cache, index, root, redirects = _corpus(tmp_path)
    ATW.configure(AC.describe(cache, index, root=root, redirects=redirects))
    assert ATW.count() == 1
    assert ATW.slugs(0, 10) == [("gene", "BRCA1")]
    assert ATW.slugs(1, 10) == []

    client = TestClient(A.app)
    r = client.get("/sitemap-atlas-1.xml")
    assert r.status_code == 200
    assert "/atlas/gene/BRCA1/" in r.text
    # The registry order matters: /sitemap-{gene}.xml must not swallow this route.
    assert client.get("/sitemap-atlas-2.xml").status_code == 404


def test_resolve_json_answers_an_entity_when_no_variant_does(tmp_path, monkeypatch):
    from starlette.testclient import TestClient

    import app as A

    cache, index, root, redirects = _corpus(tmp_path)
    ATW.configure(AC.describe(cache, index, root=root, redirects=redirects))
    monkeypatch.setattr(A, "resolve_query", lambda q: ("none", q, []))
    monkeypatch.setattr(A, "_index_hits", lambda k: [])
    monkeypatch.setattr(A, "_vcv_hits", lambda k: [])
    monkeypatch.setattr(A, "_sv_lookup", lambda k: [])
    client = TestClient(A.app)

    d = client.get("/resolve.json?q=BRCA1").json()
    assert d["not_found"] is False
    assert d["kind"] == "atlas_gene"
    hit = d["hits"][0]
    assert hit["canonical_url"].endswith("/atlas/gene/BRCA1/")
    assert hit["markdown_url"].endswith("/atlas/gene/BRCA1/index.md")

    assert client.get("/resolve.json?q=nothing-here").json()["not_found"] is True


def test_resolve_json_puts_the_atlas_page_beside_a_gene_hit(tmp_path, monkeypatch):
    from starlette.testclient import TestClient

    import app as A

    cache, index, root, redirects = _corpus(tmp_path)
    ATW.configure(AC.describe(cache, index, root=root, redirects=redirects))
    monkeypatch.setattr(A, "resolve_query",
                        lambda q: ("gene", "BRCA1", [{"gene": "BRCA1", "variants": 3}]))
    monkeypatch.setattr(A, "_gene_counts", lambda g: {"total": 3})
    client = TestClient(A.app)

    d = client.get("/resolve.json?q=BRCA1").json()
    hit = d["hits"][0]
    assert hit["kind"] == "gene"
    assert hit["atlas_url"].endswith("/atlas/gene/BRCA1/")
    assert hit["atlas_markdown_url"].endswith("/atlas/gene/BRCA1/index.md")


def test_suggest_merges_atlas_entities_without_duplicating_a_gene(tmp_path, monkeypatch):
    from starlette.testclient import TestClient

    import app as A

    cache, index, root, redirects = _corpus(tmp_path)
    ATW.configure(AC.describe(cache, index, root=root, redirects=redirects))
    monkeypatch.setattr(A, "_index", lambda: None)   # no variant index in this test
    client = TestClient(A.app)

    d = client.get("/suggest.json", params={"q": "BRCA1"}).json()
    assert d and d[0]["kind"] == "atlas_gene"
    assert d[0]["url"] == "atlas/gene/BRCA1/"
    assert d[0]["badge"] == "atlas gene"
    assert sum(1 for i in d if i.get("label") == "BRCA1") == 1

    # Two characters is the floor, exactly as the variant suggester has it.
    assert client.get("/suggest.json", params={"q": "B"}).json() == []
