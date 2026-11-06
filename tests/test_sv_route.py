"""The /sv/<key> routes over a prototype dbVar store.

One module-scoped temporary sv.db so the app's thread-local connection is reused
across tests; the records are synthetic and need no biotree.
"""
import sqlite3

import pytest
from starlette.testclient import TestClient

import app
from sugivariant.dbvar import _insert, open_db, record_from_attrs


def _rec(nssv, chrom, start, end, vcv, svtype="DEL"):
    return record_from_attrs(nssv, {"Dbvar": {
        "chromosome": chrom, "start": start, "end": end, "svtype": svtype,
        "svlen": str(end - start), "clinical_significance": "Pathogenic",
        "origin": "germline", "phenotype_ids": ["MONDO:0011929"],
        "clinvar_accessions": [vcv], "region_id": "nsv1",
        "description": "synthetic deletion",
    }})


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    db = tmp_path_factory.mktemp("svdb") / "sv.db"
    conn = open_db(str(db))
    _insert(conn, _rec("NSSV1", "1", 10001, 5532775, "VCV000980940.1"))
    _insert(conn, _rec("NSSV2", "1", 10001, 2649537, "VCV002501007.2", "DUP"))
    conn.execute("INSERT INTO sv_gene VALUES('TMC1','NSSV1')")
    conn.commit()
    conn.close()
    app.SV_DB = db
    return TestClient(app.app)


def test_html_page_has_the_sv_fields(client):
    r = client.get("/sv/NSSV1")
    assert r.status_code == 200
    assert "dbVar NSSV1" in r.text
    assert "Pathogenic" in r.text
    assert "TMC1" in r.text
    assert "do not apply" in r.text


def test_lookup_by_vcv_and_exact_coordinate(client):
    for key in ("VCV000980940.1", "NSSV1", "nsv1"):
        r = client.get(f"/sv/{key}")
        assert r.status_code == 200, key
    assert "dbVar NSSV1" in client.get("/sv/VCV000980940.1").text


def test_interval_query_returns_a_list(client):
    r = client.get("/sv/1:2000000-2100000")
    assert r.status_code == 200
    assert "NSSV1" in r.text and "NSSV2" in r.text


def test_markdown_twin(client):
    r = client.get("/sv/NSSV1.md")
    assert r.status_code == 200
    assert r.text.startswith("# dbVar NSSV1")
    assert "TMC1" in r.text


def test_unknown_is_404(client):
    assert client.get("/sv/NSSV999").status_code == 404
    assert client.get("/sv/NSSV999.md").status_code == 404
