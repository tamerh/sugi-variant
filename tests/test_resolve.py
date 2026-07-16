"""Resolver match-keys: a GENE + HGVS/AA query in any shape must normalize to the
same key the target record exposes. Pure (no biobtree / no server)."""
from sugivariant.slug import _norm_hgvs
from app import _hgvs_keys, _RSID_RE


PTEN_R130Q = {"hgvs_p": "p.Arg130Gln", "hgvs_c": "c.389G>A"}
ACTA1_P309A = {"hgvs_p": "p.Pro309Ala", "hgvs_c": "c.925C>G"}


def _matches(rec, query_hgvs):
    return _norm_hgvs(query_hgvs) in _hgvs_keys(rec)


def test_all_query_shapes_hit_the_same_record():
    # protein: 3-letter (with/without p.), 1-letter short
    assert _matches(PTEN_R130Q, "p.Arg130Gln")
    assert _matches(PTEN_R130Q, "Arg130Gln")
    assert _matches(PTEN_R130Q, "R130Q")
    # coding: with/without c.
    assert _matches(PTEN_R130Q, "c.389G>A")
    assert _matches(PTEN_R130Q, "389G>A")
    # the bare 3-letter demand shape
    assert _matches(ACTA1_P309A, "pro309ala")
    assert _matches(ACTA1_P309A, "P309A")


def test_no_false_match():
    assert not _matches(PTEN_R130Q, "R130W")     # different residue
    assert not _matches(PTEN_R130Q, "c.388G>A")  # different position


def test_rsid_pattern():
    assert _RSID_RE.match("rs121909229") and _RSID_RE.match("RS80357906")
    assert not _RSID_RE.match("rs") and not _RSID_RE.match("pten-p-arg130gln")


def test_non_missense_has_no_short_key():
    rec = {"hgvs_p": "p.Lys328del", "hgvs_c": "c.983_985del"}
    keys = _hgvs_keys(rec)
    assert _norm_hgvs("c.983_985del") in keys
    assert _norm_hgvs("p.Lys328del") in keys  # still matchable by the p. form
