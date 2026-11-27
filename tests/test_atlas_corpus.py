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
