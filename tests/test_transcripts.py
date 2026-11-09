"""The transcript block: one row per accession, roles named, nothing inferred.

The data is already on every record (`hgvs_expressions`); what matters is that we
group it by transcript, keep the highest version, label MANE Select and
ClinVar-reported, and never invent a protein change per transcript (a c. string
cannot be translated without that transcript's CDS).
"""
from sugivariant.present import ROLE_MANE, ROLE_PLUS, ROLE_REPORTED, transcript_view

RUNX1 = {
    "hgvs_expressions": [
        "LRG_482:g.1190396T>C",
        "NC_000021.8:g.36171613A>G",
        "NC_000021.9:g.34799316A>G",
        "NG_011402.2:g.1190396T>C",
        "NM_001001890.3:c.871T>C",
        "NM_001754.4:c.952T>C",
        "NM_001754.5:c.952T>C",
    ],
    "mane": {"mrna": "NM_001754", "clinvar_transcript": "NM_001754",
             "clinvar_on_mane": True},
}


def test_one_row_per_accession_and_highest_version_wins():
    tv = transcript_view(RUNX1)
    accessions = [r["accession"] for r in tv["rows"]]
    assert accessions == ["NM_001754", "NM_001001890"]
    mane = tv["rows"][0]
    assert (mane["version"], mane["versions"], mane["also"]) == ("5", ["4"], ".4")
    assert mane["change"] == "c.952T>C"


def test_roles_name_the_transcript_not_a_rank():
    tv = transcript_view(RUNX1)
    assert tv["rows"][0]["roles"] == [ROLE_MANE, ROLE_REPORTED]
    assert tv["rows"][1]["roles"] == []          # alternative: no role claimed


def test_mane_is_first_even_when_another_accession_sorts_lower():
    rec = dict(RUNX1, mane={"mrna": "NM_001001890", "clinvar_transcript": "NM_001754"})
    tv = transcript_view(rec)
    assert tv["rows"][0]["accession"] == "NM_001001890"
    assert ROLE_MANE in tv["rows"][0]["roles"]
    assert ROLE_REPORTED in tv["rows"][1]["roles"]


def test_genomic_expressions_are_not_transcript_rows():
    tv = transcript_view(RUNX1)
    assert all(r["accession"].startswith("NM_") for r in tv["rows"])
    genomic = {(g["accession"], g["change"]) for g in tv["genomic"]}
    assert ("NC_000021.9", "g.34799316A>G") in genomic
    assert ("LRG_482", "g.1190396T>C") in genomic
    assert ("NG_011402.2", "g.1190396T>C") in genomic


def test_lrg_transcript_and_lrg_locus_are_told_apart():
    rec = {"hgvs_expressions": ["LRG_482t1:c.430C>T", "LRG_482:g.1190396T>C"]}
    tv = transcript_view(rec)
    assert [r["accession"] for r in tv["rows"]] == ["LRG_482t1"]
    assert [g["accession"] for g in tv["genomic"]] == ["LRG_482"]


def test_the_expressions_fallback_derives_no_protein_change():
    """Old records carry only coding changes; nothing may be invented from them."""
    tv = transcript_view(RUNX1)
    for row in tv["rows"]:
        assert row["protein"] is None
        assert row["consequence"] is None
        assert "p." not in row["change"]


STRUCTURED = {
    "hgvs_expressions": ["NM_001754.5:c.430C>T"],
    "hgvs_transcripts": [
        {"transcript": "NM_001754.5", "hgvs_c": "NM_001754.5:c.430C>T",
         "hgvs_p": "NP_001745.2:p.Leu144=", "consequence": "synonymous variant",
         "so_id": "SO:0001819", "mane_select": True, "mane_plus_clinical": False},
        {"transcript": "NM_004992.4", "hgvs_c": "NM_004992.4:c.349C>T",
         "hgvs_p": "NP_004983.1:p.Leu117=", "consequence": "synonymous variant",
         "mane_select": False, "mane_plus_clinical": True},
    ],
    "transcript_count": 44,
    "mane": {"mrna": "NM_001754", "clinvar_transcript": "NM_001754"},
}


def test_structured_rows_carry_the_upstream_protein_and_consequence():
    tv = transcript_view(STRUCTURED)
    mane = tv["rows"][0]
    assert mane["protein"] == "p.Leu144="          # verbatim, prefix stripped only
    assert mane["consequence"] == "synonymous variant"
    assert mane["roles"] == [ROLE_MANE, ROLE_REPORTED]
    plus = tv["rows"][1]
    assert plus["roles"] == [ROLE_PLUS]


def test_omitted_transcripts_are_counted_not_implied():
    tv = transcript_view(STRUCTURED)
    assert tv["others"] == 42                     # 44 annotated, 2 anchored
    assert "44" in tv["note"] and "anchored" in tv["note"]



def test_genomic_only_record_has_no_transcript_rows():
    rec = {"hgvs_expressions": ["NC_000006.12:g.26093233G>T"]}
    tv = transcript_view(rec)
    assert tv["rows"] == []
    assert tv["genomic"] == [{"accession": "NC_000006.12", "change": "g.26093233G>T"}]
    assert tv["note"] is None


def test_no_expressions_means_no_block():
    assert transcript_view({"hgvs_expressions": []}) is None
    assert transcript_view({}) is None


def test_gene_without_mane_says_so():
    rec = {"hgvs_expressions": ["NR_123456.1:n.10A>G"]}
    tv = transcript_view(rec)
    assert tv["mane_missing"] is True
    assert tv["rows"][0]["roles"] == []


def test_note_names_the_reported_transcript():
    tv = transcript_view(RUNX1)
    assert "NM_001754" in tv["note"]
    assert "MANE" in tv["note"]
    assert "do not choose" in tv["note"]


def test_the_page_block_and_md_line_name_the_versions(monkeypatch):
    """Versions live in one page-level block (Atlas shape), not the site footer, and
    the .md twin carries the same facts."""
    import app as A
    from pathlib import Path

    monkeypatch.setattr(A, "RELEASE_MODE", True)
    monkeypatch.setattr(A, "_RELEASE_VERSION", "v1.5.0")
    monkeypatch.setattr(A, "_RELEASE_DATE", "2026-11-07")
    monkeypatch.setattr(A, "_RELEASE_META", {"biobtree_version": "v2.20.0"})
    assert A.release_line_md() == ("_Release v1.5.0 · BioBTree v2.20.0 · "
                                   "last updated 2026-11-07_")
    root = Path(A.__file__).resolve().parent
    footer = (root / "templates/base.html").read_text()
    block = (root / "templates/_sources_line.html").read_text()
    assert "Release <strong>" not in footer          # nothing in the main footer
    assert "About this page" in block and "Sources" in block
    for field in ("Generated", "Sugi Variant", "BioBTree", "ClinVar"):
        assert field in block
