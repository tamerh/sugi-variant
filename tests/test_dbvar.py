"""dbVar structural-variant prototype: record parsing, aliases, interval lookup.

Pure SQLite; no biotree and no network. The build/collect path is exercised in
the live prototype, not here.
"""
from sugivariant.dbvar import (
    _insert, aliases, open_db, lookup, overlapping, record_from_attrs, render_md,
    should_build,
)


def _rec(nssv, chrom, start, end, svtype="DEL", sig="Pathogenic", vcv="VCV000980940.1"):
    return record_from_attrs(nssv, {"Dbvar": {
        "chromosome": chrom, "start": start, "end": end, "svtype": svtype,
        "svlen": str(end - start), "clinical_significance": sig,
        "origin": "germline", "phenotype_ids": ["MONDO:0011929"],
        "clinvar_accessions": [vcv], "region_id": "nsv4728187",
        "description": "1p36 deletion",
    }})


def _db(tmp_path):
    conn = open_db(str(tmp_path / "sv.db"))
    _insert(conn, _rec("NSSV1", "1", 10001, 5532775))
    _insert(conn, _rec("NSSV2", "1", 10001, 2649537, svtype="DUP", sig="Benign",
                       vcv="VCV002501007.2"))
    _insert(conn, _rec("NSSV3", "9", 72521608, 72838297, svtype="DEL",
                       vcv="VCV003000001.1"))
    conn.commit()
    return conn


def test_record_from_attrs_parses_the_dbvar_shape():
    rec = _rec("NSSV1", "1", 10001, 5532775)
    assert rec["schema"] == "sugi-variant-sv/1"
    assert rec["chrom"] == "1" and rec["start"] == 10001 and rec["end"] == 5532775
    assert rec["svtype"] == "DEL"
    assert rec["phenotypes"] == ["MONDO:0011929"]
    assert rec["clinvar"] == ["VCV000980940.1"]
    assert rec["nsv"] == "nsv4728187"


def test_aliases_cover_identity_region_and_clinvar():
    keys = aliases(_rec("NSSV1", "1", 10001, 5532775))
    assert "nssv1" in keys
    assert "nsv4728187" in keys
    assert "1:10001-5532775" in keys
    assert "vcv000980940.1" in keys
    assert "980940" in keys


def test_lookup_and_interval_overlap(tmp_path):
    conn = _db(tmp_path)
    assert [r["nssv"] for r in lookup(conn, "NSSV1")] == ["NSSV1"]
    assert [r["nssv"] for r in lookup(conn, "1:10001-5532775")] == ["NSSV1"]
    assert [r["nssv"] for r in lookup(conn, "980940")] == ["NSSV1"]
    # NSSV2 overlaps; NSSV3 is on another chromosome
    hits = overlapping(conn, "1", 2_000_000, 2_100_000)
    assert {r["nssv"] for r in hits} == {"NSSV1", "NSSV2"}
    assert overlapping(conn, "1", 6_000_000, 6_100_000) == []
    assert [r["nssv"] for r in overlapping(conn, "9", 72_600_000, 72_700_000)] == ["NSSV3"]


def test_insert_is_idempotent(tmp_path):
    conn = _db(tmp_path)
    _insert(conn, _rec("NSSV1", "1", 10001, 5532775))
    conn.commit()
    assert len(lookup(conn, "nssv1")) == 1
    assert len(overlapping(conn, "1", 10001, 10002)) == 2


def test_gate_keeps_germline_acmg_classes_only():
    assert should_build({"significance": "Pathogenic"})
    assert should_build({"significance": "Uncertain significance"})
    assert should_build({"significance": "conflicting data from submitters"})
    for skip in ("not provided", "drug response", "risk factor", "evidence_only",
                 "Tier I - Strong", "protective", ""):
        assert not should_build({"significance": skip}), skip


def test_render_md_states_the_snv_predictor_boundary():
    md = render_md(_rec("NSSV1", "1", 10001, 5532775), ["TMC1"])
    assert "# dbVar NSSV1" in md
    assert "Pathogenic" in md
    assert "TMC1" in md
    assert "do not apply" in md
    assert "5,522,774 bp" in md
