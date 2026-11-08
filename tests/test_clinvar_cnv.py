"""ClinVar CNVs: the clinical record for copy-number variants dbVar does not carry.

The parser is the whole risk here — an interval placed on the wrong assembly, or a
CNV silently treated as a small variant, would put a wrong coordinate in front of a
clinical reader. These tests pin what is accepted and what is refused.
"""
from sugivariant.clinvar_cnv import cnv_from_attrs, vcv_id


def _attrs(name, kind="copy number loss", cls="Pathogenic", **extra):
    d = {"name": name, "type": kind, "germline_classification": cls,
         "phenotype_ids": ["MONDO:0000001"]}
    d.update(extra)
    return {"Clinvar": d}


def test_grch38_name_is_placed():
    rec = cnv_from_attrs(5020085, _attrs("GRCh38/hg38 7q11.23(chr7:73293482-74732632)x1"))
    assert rec["chrom"] == "7"
    assert (rec["start"], rec["end"]) == (73293482, 74732632)
    assert rec["svtype"] == "DEL"
    assert rec["nssv"] == "VCV005020085"
    assert rec["source"] == "clinvar"
    assert rec["clinvar"] == ["VCV005020085"]


def test_copy_number_gain_is_a_duplication():
    rec = cnv_from_attrs(1, _attrs("GRCh38/hg38 7q11.23(chr7:100-200)x3",
                                   kind="copy number gain"))
    assert rec["svtype"] == "DUP"


def test_hg19_only_is_refused():
    """A GRCh37 interval is 1-2 Mb away from its hg38 equivalent; placing it would be
    a wrong coordinate, so it is dropped rather than guessed."""
    assert cnv_from_attrs(395455, _attrs("GRCh37/hg19 7q11.23(chr7:73470470-73498327)x1")) is None


def test_structured_coordinates_win_when_the_name_has_none():
    rec = cnv_from_attrs(5038037, _attrs("GRCh38/hg38 2q22.1-22.3", start=140732866,
                                         stop=144690027, chromosome="2"))
    assert (rec["chrom"], rec["start"], rec["end"]) == ("2", 140732866, 144690027)


def test_non_cnv_types_are_not_claimed():
    assert cnv_from_attrs(9, _attrs("NM_000001.1(G):c.1A>G", kind="single nucleotide variant")) is None


def test_a_classification_outside_the_gate_is_refused():
    assert cnv_from_attrs(10, _attrs("GRCh38/hg38 1p36(chr1:100-200)x1",
                                     cls="drug response")) is None


def test_vcv_ids_are_zero_padded_to_clinvar_width():
    assert vcv_id(5020085) == "VCV005020085"
    assert vcv_id(17662) == "VCV000017662"
