"""Resolution index: coordinate normalization + key→VCV lookup across all three
arrival keys. Builds an in-memory index from synthetic records (no biobtree)."""
from sugivariant import index as IX


def _rec(vcv, gene, slugs, hgvs_p, hgvs_c, rsid, coord, cls="Pathogenic"):
    return {"variation_id": vcv, "gene_symbol": gene, "canonical_slug": slugs[0],
            "slugs": slugs, "hgvs_p": hgvs_p, "hgvs_c": hgvs_c, "classification": cls,
            "review_status": "criteria provided, single submitter", "rsid": rsid,
            "coordinate": coord, "conditions": [{"name": "Cowden syndrome 1"}],
            "concordance": {"consensus": None, "flags": []}}


def _db():
    conn = IX.open_db(":memory:")
    IX.index_gene(conn, "PTEN", recs=[
        _rec("7829", "PTEN", ["pten-p-arg130gln", "pten-c-389g-a"],
             "p.Arg130Gln", "c.389G>A", "rs121909229", "10:87933148:G:A")])
    return conn


def test_norm_coordinate():
    assert IX.norm_coordinate("10:87933148:G:A") == "10:87933148:G:A"
    assert IX.norm_coordinate("chr10-87933148-g-a") == "10:87933148:G:A"
    assert IX.norm_coordinate("PTEN R130Q") is None
    assert IX.norm_coordinate("rs1") is None


def test_lookup_all_three_key_schemes():
    c = _db()
    for key in ["rs121909229",            # rsID
                "10:87933148:g:a",         # coordinate (lookup lowercases)
                "pten-r130q",              # 1-letter short
                "pten-p-arg130gln",        # canonical p. slug
                "pten-c-389g-a",           # c. slug
                "pten-arg130gln"]:         # bare p. form
        hits = IX.lookup(c, key)
        assert len(hits) == 1 and hits[0]["slug"] == "pten-p-arg130gln", key
    assert IX.lookup(c, "pten-r130w") == []          # different residue → no hit


def test_reindex_is_idempotent():
    c = _db()
    IX.index_gene(c, "PTEN", recs=[
        _rec("7829", "PTEN", ["pten-p-arg130gln", "pten-c-389g-a"],
             "p.Arg130Gln", "c.389G>A", "rs121909229", "10:87933148:G:A")])
    assert IX.stats(c)["variants"] == 1                # no duplicate rows
    assert len(IX.lookup(c, "rs121909229")) == 1


# ── gene_meta ↔ variant consistency ───────────────────────────────────────────
# gene_meta used to be keyed on the ENUMERATED gene while rows are filed under
# r["gene_symbol"] (taken from the HGVS transcript). The two diverged, producing
# 1,067 "phantom" genes that claimed 46,944 variants they had no rows for — each a
# public URL that fell through to an unbounded live build and was advertised by
# type-ahead. These tests pin the invariant.
def _mkdb():
    import sqlite3

    from sugivariant.index import SCHEMA
    c = sqlite3.connect(":memory:")
    c.row_factory = sqlite3.Row
    c.executescript(SCHEMA)
    return c


def _add(c, vcv, gene):
    c.execute("INSERT INTO variant VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
              (vcv, gene, f"slug{vcv}", None, None, "VUS", "criteria provided", 1,
               None, None, None, None))


def test_reconcile_gene_meta_fixes_phantoms_counts_and_gaps():
    from sugivariant.index import gene_known_empty, reconcile_gene_meta
    c = _mkdb()
    for vcv, gene in [(1, "TTN"), (2, "TTN"), (3, "TTN"), (4, "NEWGENE")]:
        _add(c, vcv, gene)
    c.executemany("INSERT INTO gene_meta VALUES(?,?)",
                  [("TTN-AS1", 6410),      # phantom: claims rows it does not have
                   ("TTN", 2),             # undercount: rows were filed here from elsewhere
                   ("GONEGENE", 0)])       # legitimately enumerated-empty
    c.commit()

    phantom, added, fixed = reconcile_gene_meta(c)
    assert (phantom, added, fixed) == (1, 1, 1)

    meta = dict(c.execute("SELECT gene, n FROM gene_meta").fetchall())
    assert meta["TTN-AS1"] == 0          # phantom zeroed
    assert meta["TTN"] == 3              # count corrected
    assert meta["NEWGENE"] == 1          # missing row added
    assert meta["GONEGENE"] == 0         # negative marker preserved

    # the phantom is now gated, a real gene is not, and an unseen gene may still build
    assert gene_known_empty(c, "TTN-AS1") is True
    assert gene_known_empty(c, "GONEGENE") is True
    assert gene_known_empty(c, "TTN") is False
    assert gene_known_empty(c, "NEVERSEEN") is False


def test_index_gene_writes_meta_for_the_symbol_rows_are_filed_under():
    """A record whose gene_symbol differs from the enumerated gene must not leave the
    enumerated name claiming variants it has no rows for."""
    from sugivariant.index import index_gene
    c = _mkdb()
    recs = [{"variation_id": 11, "gene_symbol": "TTN", "canonical_slug": "ttn-p-x1y",
             "classification": "VUS", "review_status": "criteria provided", "conditions": []}]
    index_gene(c, "TTN-AS1", recs)          # enumerated as TTN-AS1, filed under TTN
    meta = dict(c.execute("SELECT gene, n FROM gene_meta").fetchall())
    assert meta.get("TTN") == 1, meta
    assert meta.get("TTN-AS1") == 0, f"enumerated name must not claim rows: {meta}"


def test_live_index_gene_meta_matches_variant_rows():
    """Invariant on the real index, when present: no gene_meta row may claim a count
    that the variant table does not actually hold."""
    import os
    import pathlib
    import sqlite3

    import pytest
    db = os.environ.get("INDEX_DB") or "/data/sugi-variant/cache/index.db"
    if not pathlib.Path(db).exists():
        pytest.skip(f"no index at {db}")
    c = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    bad = c.execute("SELECT COUNT(*) FROM gene_meta m WHERE m.n != "
                    "(SELECT COUNT(*) FROM variant v WHERE v.gene = m.gene)").fetchone()[0]
    missing = c.execute("SELECT COUNT(*) FROM (SELECT DISTINCT gene FROM variant v WHERE NOT "
                        "EXISTS (SELECT 1 FROM gene_meta m WHERE m.gene = v.gene))").fetchone()[0]
    c.close()
    assert bad == 0, f"{bad} gene_meta rows disagree with the variant table"
    assert missing == 0, f"{missing} genes have variants but no gene_meta row"
