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
