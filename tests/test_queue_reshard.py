"""The queue runner must tell a legitimate carrier reshard from a real collapse.

A gene pickle can carry another gene's ClinVar records (the "carrier-gene defect"),
so dropping them shrinks the file without losing anything. The collapse guard
exists to catch *real* loss, so it must stay armed unless every dropped variant is
served by the gene it actually belongs to.
"""
import gzip
import importlib.util
import pickle
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "queue_rebuild", Path(__file__).resolve().parents[1] / "scripts/passes/queue_rebuild.py")
QR = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(QR)


def _write(records_dir, gene, rows):
    path = Path(records_dir) / f"{gene}.pkl.gz"
    with gzip.open(path, "wb") as f:
        pickle.dump(rows, f)


def _rec(vcv, gene, hgvs_p=None, hgvs_c=None, coord=None):
    return {"variation_id": vcv, "gene_symbol": gene, "hgvs_p": hgvs_p,
            "hgvs_c": hgvs_c, "coordinate": coord}


def test_no_stored_file_is_not_a_reshard(tmp_path):
    assert QR.carrier_reshard(str(tmp_path), "GENE1", [_rec("1", "GENE1")])[0] is False


def test_carrier_rows_covered_by_their_own_pickle(tmp_path):
    _write(tmp_path, "HOST", [_rec("1", "HOST")] + [_rec(str(v), "OTHER") for v in range(10, 15)])
    _write(tmp_path, "OTHER", [_rec(str(v), "OTHER") for v in range(10, 15)])
    safe, dropped, uncovered = QR.carrier_reshard(
        str(tmp_path), "HOST", [_rec("1", "HOST")])
    assert safe is True
    assert dropped == 5
    assert uncovered == {}


def test_uncovered_carrier_rows_refuse(tmp_path):
    _write(tmp_path, "HOST", [_rec("1", "HOST")] + [_rec("10", "OTHER"), _rec("11", "OTHER")])
    _write(tmp_path, "OTHER", [_rec("10", "OTHER")])   # 11 is not there
    safe, dropped, uncovered = QR.carrier_reshard(
        str(tmp_path), "HOST", [_rec("1", "HOST")])
    assert safe is False
    assert dropped == 2
    assert uncovered == {"OTHER": 1}


def test_missing_carrier_pickle_refuses(tmp_path):
    _write(tmp_path, "HOST", [_rec("1", "HOST"), _rec("10", "GONE")])
    safe, _dropped, uncovered = QR.carrier_reshard(
        str(tmp_path), "HOST", [_rec("1", "HOST")])
    assert safe is False
    assert uncovered == {"GONE": 1}


def test_same_protein_change_counts_as_served(tmp_path):
    # Two VCVs, different nucleotides, one protein change: the build makes one page.
    _write(tmp_path, "HOST", [_rec("1", "HOST"), _rec("10", "OTHER", hgvs_p="p.Pro45=")])
    _write(tmp_path, "OTHER", [_rec("20", "OTHER", hgvs_p="p.Pro45=")])
    safe, _dropped, uncovered = QR.carrier_reshard(
        str(tmp_path), "HOST", [_rec("1", "HOST")])
    assert safe is True
    assert uncovered == {}


def test_same_coding_change_under_a_neighbour_gene_is_served(tmp_path):
    # LRTOMT c.529C>T is the same variant the neighbouring annotation calls
    # LRRC51 c.529C>T; one page exists, so the other VCV is not a lost page.
    _write(tmp_path, "HOST", [_rec("1", "HOST"), _rec("10", "OTHER", hgvs_c="c.529C>T")])
    _write(tmp_path, "OTHER", [_rec("20", "ELSEWHERE", hgvs_c="c.529C>T")])
    safe, _dropped, uncovered = QR.carrier_reshard(
        str(tmp_path), "HOST", [_rec("1", "HOST")])
    assert safe is True
    assert uncovered == {}


def test_a_different_variant_at_the_same_residue_is_not_served(tmp_path):
    _write(tmp_path, "HOST", [_rec("1", "HOST"),
                              _rec("10", "OTHER", hgvs_c="c.5833C>T", hgvs_p="p.Gln1945Ter")])
    _write(tmp_path, "OTHER", [_rec("20", "OTHER", hgvs_c="c.5834C>T", hgvs_p="p.Gln1945Arg")])
    safe, _dropped, uncovered = QR.carrier_reshard(
        str(tmp_path), "HOST", [_rec("1", "HOST")])
    assert safe is False
    assert uncovered == {"OTHER": 1}


def test_a_real_collapse_still_refuses(tmp_path):
    _write(tmp_path, "HOST", [_rec(str(v), "HOST") for v in range(1, 6)])
    safe, dropped, uncovered = QR.carrier_reshard(
        str(tmp_path), "HOST", [_rec("1", "HOST")])
    assert safe is False
    assert dropped == 4
    assert uncovered == {"HOST": 4}
