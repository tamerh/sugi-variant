"""Snapshot manifest: read-only audit of records/ + index.db.

No biotree, no app, no serving path. These tests pin the two properties that
matter before the manifest is allowed near the live cache:

  * it describes exactly what is on disk and in the index (counts, missing,
    extra, stale), and
  * it never mutates either one.
"""
import gzip
import json
import os
import pickle
import sqlite3

from sugivariant.snapshot import (
    build_manifest, read_index, read_stamp, rebuild_queue, scan_records,
    serve_manifest, verify_manifest,
)


def _records(root, genes):
    """genes = {GENE: (n_records, stamp_or_None)}"""
    os.makedirs(root, exist_ok=True)
    for gene, (n, stamp) in genes.items():
        with gzip.open(os.path.join(root, f"{gene}.pkl.gz"), "wb") as f:
            pickle.dump([{"variation_id": i, "gene_symbol": gene,
                          "canonical_slug": f"{gene.lower()}-{i}"} for i in range(n)], f)
        if stamp is not None:
            with open(os.path.join(root, f"{gene}.dbver.json"), "w") as f:
                json.dump(stamp, f)


def _index(path, gene_counts):
    conn = sqlite3.connect(path)
    conn.executescript(
        "CREATE TABLE gene_meta(gene TEXT PRIMARY KEY, n INTEGER);"
        "CREATE TABLE variant(vcv INTEGER PRIMARY KEY, gene TEXT, slug TEXT, flag TEXT);"
        "CREATE TABLE alias(key TEXT, vcv INTEGER);")
    vcv = 1
    for gene, n in gene_counts.items():
        conn.execute("INSERT INTO gene_meta VALUES(?,?)", (gene, n))
        for i in range(n):
            conn.execute("INSERT INTO variant VALUES(?,?,?,?)",
                         (vcv, gene, f"{gene.lower()}-{i}", None))
            vcv += 1
    conn.commit()
    conn.close()


def _fixture(tmp_path):
    rec = tmp_path / "records"
    _records(rec, {
        "AAA": (2, {"clinvar": "v1"}),
        "BBB": (1, {"clinvar": "v0"}),       # stale
        "CCC": (0, None),                    # empty, unstamped, orphan
        "FFF": (1, {"clinvar": "v1"}),       # pickle behind index (index 3)
        "GGG": (4, {"clinvar": "v1"}),       # surplus (index 2)
    })
    idx = tmp_path / "index.db"
    _index(idx, {"AAA": 2, "BBB": 1, "DDD": 5, "EEE": 0, "FFF": 3, "GGG": 2})
    return str(rec), str(idx)


def test_build_manifest_describes_disk_and_index(tmp_path):
    rec, idx = _fixture(tmp_path)
    manifest, report = build_manifest(
        rec, idx, count_records=True, current_versions={"clinvar": "v1"},
        release="v1.0.0", released_at="2026-11-06")

    assert manifest["schema"] == "sugi-variant-release/1"
    assert manifest["release"] == "v1.0.0"
    assert manifest["released_at"] == "2026-11-06"
    state = manifest["corpora"]["small_variant"]["gene_state"]
    assert set(state) == {"AAA", "BBB", "CCC", "DDD", "EEE", "FFF", "GGG"}

    # the three differences that matter
    assert report["index_only"] == ["DDD", "EEE"]
    assert report["orphan"] == ["CCC"]                # records, no current index rows
    assert report["records_only"] == []
    assert report["missing"] == ["DDD"]
    assert report["known_empty"] == ["EEE"]          # n=0: no page expected
    assert report["empty"] == ["CCC"]
    assert report["stale"] == [{"gene": "BBB", "datasets": ["clinvar"]}]
    assert report["unstamped"] == ["CCC"]
    assert {m["gene"] for m in report["mismatch"]} == {"FFF", "GGG"}

    assert state["AAA"]["status"] == "ok"
    assert state["BBB"]["status"] == "stale"
    assert state["CCC"]["status"] == "orphan"
    assert state["DDD"]["status"] == "missing"
    assert state["EEE"]["status"] == "known_empty"
    assert state["FFF"]["status"] == "mismatch"
    assert state["GGG"]["status"] == "mismatch"
    assert state["AAA"]["n"] == 2
    assert state["DDD"]["n"] is None

    # served corpus excludes orphan and missing files
    corpus = manifest["corpora"]["small_variant"]
    assert corpus["genes"] == 4
    assert corpus["record_files"] == 5
    assert corpus["orphan_files"] == 1
    assert corpus["known_empty"] == 1
    assert corpus["missing"] == 1
    assert corpus["records"] == 8               # AAA 2 + BBB 1 + FFF 1 + GGG 4


def test_serve_manifest_is_compact_and_excludes_unservable(tmp_path):
    rec, idx = _fixture(tmp_path)
    manifest, _ = build_manifest(rec, idx, count_records=True, release="v1.0.0",
                                 released_at="2026-11-06", biobtree_version="v2.20.0")
    serve = serve_manifest(manifest)
    assert serve["schema"] == "sugi-variant-release-serve/1"
    assert serve["release"] == "v1.0.0"
    assert serve["released_at"] == "2026-11-06"
    assert serve["biobtree_version"] == "v2.20.0"
    genes = serve["corpora"]["small_variant"]["genes"]
    assert set(genes) == {"AAA", "BBB", "FFF", "GGG"}
    assert "CCC" not in genes      # orphan
    assert "DDD" not in genes      # missing record
    assert "EEE" not in genes      # known-empty
    assert genes["AAA"]["path"] == "AAA.pkl.gz"


def test_rebuild_queue_targets_only_missing_variants(tmp_path):
    rec, idx = _fixture(tmp_path)
    manifest, _ = build_manifest(rec, idx, count_records=True)
    q = rebuild_queue(manifest)
    assert [e["gene"] for e in q["queue"]] == ["DDD", "FFF"]
    assert [e["pages"] for e in q["queue"]] == [5, 2]
    assert q["totals"]["pages_to_add"] == 7
    assert q["totals"]["rows_to_add"] == 7
    assert q["totals"]["by_reason"] == {"no_record": 1, "pages_behind": 1}
    assert q["uncounted"] == []


def test_rebuild_queue_separates_alias_only_and_slug_drift(tmp_path):
    """Three failure modes a row count cannot see:

    HHH: three VCVs share one slug -> 0 missing pages, 2 alias-only rows.
    III: same row count on both sides but different slugs -> 1 missing page.
    JJJ: the index carries a coordinate slug whose VCV already has a pickle page
         under an alias -> 0 missing pages, 1 stale slug to fix by reindex.
    """
    rec = tmp_path / "records"
    _records(rec, {"HHH": (1, {"clinvar": "v1"}), "III": (1, {"clinvar": "v1"}),
                   "JJJ": (1, {"clinvar": "v1"})})
    idx = tmp_path / "index.db"
    conn = sqlite3.connect(idx)
    conn.executescript(
        "CREATE TABLE gene_meta(gene TEXT PRIMARY KEY, n INTEGER);"
        "CREATE TABLE variant(vcv INTEGER PRIMARY KEY, gene TEXT, slug TEXT, flag TEXT);"
        "CREATE TABLE alias(key TEXT, vcv INTEGER);")
    for gene, n in {"HHH": 3, "III": 1, "JJJ": 1}.items():
        conn.execute("INSERT INTO gene_meta VALUES(?,?)", (gene, n))
    conn.executemany("INSERT INTO variant VALUES(?,?,?,?)", [
        (0, "HHH", "hhh-0", None),       # the one the pickle holds
        (1, "HHH", "hhh-0", None),       # alias-only duplicates
        (2, "HHH", "hhh-0", None),
        (3, "III", "iii-9", None),       # pickle has iii-0
        (4, "JJJ", "5:100:g:a", None),   # stale slug; alias names jjj-0
    ])
    conn.executemany("INSERT INTO alias VALUES(?,?)",
                     [("jjj-0", 4), ("5:100:g:a", 4)])
    conn.commit()
    conn.close()

    manifest, _ = build_manifest(str(rec), str(idx), count_records=True)
    q = rebuild_queue(manifest)
    assert [e["gene"] for e in q["queue"]] == ["III"]
    assert q["queue"][0]["pages"] == 1
    assert q["queue"][0]["rows"] == 0
    assert [e["gene"] for e in q["reindex"]] == ["JJJ"]
    assert q["reindex"][0]["pages_stale"] == 1
    assert [e["gene"] for e in q["alias_only"]] == ["HHH"]
    assert q["alias_only"][0]["alias_only"] == 2
    assert q["totals"]["pages_to_add"] == 1
    assert q["totals"]["reindex_genes"] == 1
    assert q["totals"]["stale_pages_to_reindex"] == 1
    assert q["totals"]["alias_only_genes"] == 1
    assert q["totals"]["alias_only_rows"] == 2

    state = manifest["corpora"]["small_variant"]["gene_state"]
    assert state["HHH"]["pages_missing"] == 0
    assert state["HHH"]["alias_only"] == 2
    assert state["III"]["pages_missing"] == 1
    assert state["III"]["missing_slugs"] == ["iii-9"]
    assert state["III"]["extra_slugs"] == ["iii-0"]
    assert state["JJJ"]["pages_missing"] == 0
    assert state["JJJ"]["pages_stale"] == 1
    assert state["JJJ"]["missing_slugs"] == []
    assert state["JJJ"]["stale_slugs"] == ["5:100:g:a"]


def test_rebuild_queue_does_not_guess_when_counts_are_missing(tmp_path):
    rec, idx = _fixture(tmp_path)
    manifest, _ = build_manifest(rec, idx, count_records=False)
    q = rebuild_queue(manifest)
    # DDD is still a definite no_record gap; the rest need a count
    assert [e["gene"] for e in q["queue"]] == ["DDD"]
    assert set(q["uncounted"]) == {"AAA", "BBB", "FFF", "GGG"}


def test_build_manifest_does_not_write_to_records_or_index(tmp_path):
    rec, idx = _fixture(tmp_path)
    before = {p: os.stat(os.path.join(rec, p)).st_mtime_ns for p in os.listdir(rec)}
    manifest, _ = build_manifest(rec, idx, count_records=True)
    after = {p: os.stat(os.path.join(rec, p)).st_mtime_ns for p in os.listdir(rec)}
    assert before == after
    # read_index leaves the db untouched (no journal/WAL side files)
    assert not any(p.startswith("index.db-") for p in os.listdir(os.path.dirname(idx)))
    assert manifest["corpora"]["small_variant"]["index_records"] == 13


def test_verify_detects_a_deleted_record(tmp_path):
    rec, idx = _fixture(tmp_path)
    manifest, _ = build_manifest(rec, idx, count_records=True)
    os.remove(os.path.join(rec, "BBB.pkl.gz"))
    res = verify_manifest(manifest, rec, idx, count_records=True)
    assert res["ok"] is False
    assert any(e["gene"] == "BBB" and e["error"] == "file-missing" for e in res["errors"])


def test_verify_detects_a_changed_count(tmp_path):
    rec, idx = _fixture(tmp_path)
    manifest, _ = build_manifest(rec, idx, count_records=True)
    # replace AAA's pickle with a shorter one
    with gzip.open(os.path.join(rec, "AAA.pkl.gz"), "wb") as f:
        pickle.dump([{"variation_id": 1}], f)
    res = verify_manifest(manifest, rec, idx, count_records=True)
    assert res["ok"] is False
    assert any(e["gene"] == "AAA" and e["error"] == "size-changed" for e in res["errors"])


def test_scan_reads_stamps_and_can_sample_counts(tmp_path):
    rec, _ = _fixture(tmp_path)
    found = scan_records(rec, count_sample=1)
    assert set(found) == {"AAA", "BBB", "CCC", "FFF", "GGG"}
    assert found["AAA"]["stamp"] == {"clinvar": "v1"}
    assert found["CCC"]["stamp"] is None
    # exactly one gene was counted by the deterministic sample
    assert sum(1 for r in found.values() if r["n"] is not None) == 1


def test_read_stamp_is_none_when_absent(tmp_path):
    p = tmp_path / "X.pkl.gz"
    p.write_bytes(b"not a pickle")
    assert read_stamp(str(p)) is None


def test_read_index_is_read_only(tmp_path):
    idx = tmp_path / "index.db"
    _index(idx, {"AAA": 2})
    out = read_index(str(idx))
    assert out["gene_meta"] == {"AAA": 2}
    assert out["variant_counts"] == {"AAA": 2}
    assert out["total"] == 2
