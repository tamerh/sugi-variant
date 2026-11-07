"""A duplicate VCV's dead URL must resolve to the page that serves the variant.

When the build collapses two ClinVar VCVs into one page, only the survivor is
written; the dropped VCV keeps its index row and its aliases, so its URL 404s
even though the same variant is served under the survivor's slug. The merge pass
re-points those keys and deletes the dead row. A row whose slug nothing serves is
left alone -- that is a real gap and deleting it would hide it.
"""
import importlib.util
from pathlib import Path

from sugivariant.index import index_gene, open_db

_SPEC = importlib.util.spec_from_file_location(
    "merge_stale_slugs",
    Path(__file__).resolve().parents[1] / "scripts/passes/merge_stale_slugs.py")
MS = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(MS)


def _rec(vcv, slug, hgvs_p="p.Leu144=", hgvs_c="c.432G>A",
         coordinate="21:34880633:C:T"):
    return {"variation_id": str(vcv), "gene_symbol": "GENE1", "canonical_slug": slug,
            "slugs": [slug, "gene1-c-432g-a"], "hgvs_p": hgvs_p, "hgvs_c": hgvs_c,
            "coordinate": coordinate,
            "classification": "Likely benign", "review_status": "reviewed by expert panel"}


def _dead_row(conn, vcv, slug, gene="GENE1", hgvs_p=None, hgvs_c=None,
              coordinate="21:34880635:G:A"):
    conn.execute(
        "INSERT OR REPLACE INTO variant VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
        (vcv, gene, slug, hgvs_p, hgvs_c, "Likely benign", None, 2, None,
         coordinate, None, None))
    conn.execute("INSERT OR IGNORE INTO alias VALUES(?,?)", (slug, vcv))


def test_dead_duplicate_row_resolves_to_the_survivor(tmp_path):
    db = str(tmp_path / "index.db")
    conn = open_db(db)
    index_gene(conn, "GENE1", recs=[_rec(2181528, "gene1-p-leu144")])
    _dead_row(conn, 2879097, "21:34880635:g:a", hgvs_p="p.Leu144=", hgvs_c="c.430C>T")
    conn.commit()

    merged, kept = MS.merge_gene(conn, "GENE1", [_rec(2181528, "gene1-p-leu144")],
                                 apply=True)
    assert [m["vcv"] for m in merged] == [2879097]
    assert kept == []
    assert [r["vcv"] for r in conn.execute(
        "SELECT vcv FROM alias WHERE key='21:34880635:g:a'")] == [2181528]
    assert conn.execute("SELECT COUNT(*) FROM variant WHERE vcv=2879097").fetchone()[0] == 0
    conn.close()


def test_a_row_nothing_serves_is_left_alone(tmp_path):
    db = str(tmp_path / "index.db")
    conn = open_db(db)
    index_gene(conn, "GENE1", recs=[_rec(1, "gene1-p-leu144")])
    _dead_row(conn, 2, "gene1-p-arg145gly", hgvs_p="p.Arg145Gly", hgvs_c="c.433C>G",
              coordinate="21:34880634:G:T")
    conn.commit()

    merged, kept = MS.merge_gene(conn, "GENE1", [_rec(1, "gene1-p-leu144")], apply=True)
    assert merged == []
    assert [k["vcv"] for k in kept] == [2]
    assert conn.execute("SELECT COUNT(*) FROM variant WHERE vcv=2").fetchone()[0] == 1
    conn.close()


def test_gene_meta_follows_the_rows_actually_filed(tmp_path):
    db = str(tmp_path / "index.db")
    conn = open_db(db)
    index_gene(conn, "GENE1", recs=[_rec(1, "gene1-p-leu144")])
    _dead_row(conn, 2, "21:34880635:g:a", hgvs_p="p.Leu144=")
    conn.execute("UPDATE gene_meta SET n=2 WHERE gene='GENE1'")   # as an older index counted it
    conn.commit()
    assert conn.execute("SELECT n FROM gene_meta WHERE gene='GENE1'").fetchone()[0] == 2

    MS.merge_gene(conn, "GENE1", [_rec(1, "gene1-p-leu144")], apply=True)
    assert conn.execute("SELECT n FROM gene_meta WHERE gene='GENE1'").fetchone()[0] == 1
    conn.close()


def test_dry_run_changes_nothing(tmp_path):
    db = str(tmp_path / "index.db")
    conn = open_db(db)
    index_gene(conn, "GENE1", recs=[_rec(1, "gene1-p-leu144")])
    _dead_row(conn, 2, "21:34880635:g:a", hgvs_p="p.Leu144=")
    conn.commit()

    merged, _kept = MS.merge_gene(conn, "GENE1", [_rec(1, "gene1-p-leu144")], apply=False)
    assert [m["vcv"] for m in merged] == [2]
    assert conn.execute("SELECT COUNT(*) FROM variant WHERE vcv=2").fetchone()[0] == 1
    conn.close()
