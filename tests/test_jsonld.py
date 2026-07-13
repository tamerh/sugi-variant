"""Variant JSON-LD + FAQPage — the AI-citation payload. Pure (no biobtree):
given an enriched record, emit a valid schema.org @graph whose FAQPage answers
the literal 'Is GENE HGVS pathogenic?' question verbatim."""
import json
from sugivariant import jsonld as VJ

_REC = {
    "canonical_slug": "acta1-p-pro309ala", "variation_id": "4685499",
    "gene_symbol": "ACTA1", "hgnc_id": "HGNC:129", "hgvs_p": "p.Pro309Ala",
    "hgvs_c": "c.925C>G", "rsid": "rs1559", "hgvs_expressions": ["NM_001100.4:c.925C>G"],
    "classification": "Pathogenic", "review_status": "criteria provided, single submitter",
    "submitter_count": 2,
    "conditions": [{"mondo_id": "MONDO:1", "name": "nemaline myopathy"}],
    "alphamissense": {"class": "likely_pathogenic", "score": "0.91"},
    "gnomad": {"band": "absent from gnomAD", "absent": True},
    "digest": {"name": "nemaline myopathy", "inheritance": ["Autosomal dominant"]},
}


def test_graph_is_valid_and_has_both_nodes():
    g = VJ.build_jsonld(_REC, {"generated_at": "2026-07-07T00:00:00+00:00"})
    json.dumps(g)                                        # serializable
    types = [n["@type"] for n in g["@graph"]]
    assert "MedicalEntity" in types and "FAQPage" in types
    ent = next(n for n in g["@graph"] if n["@type"] == "MedicalEntity")
    ids = {i["propertyID"]: i["value"] for i in ent["identifier"]}
    assert ids["ClinVar"] == "VCV4685499" and ids["dbSNP"] == "rs1559"
    assert ids["HGVS.p"] == "p.Pro309Ala"
    assert ent["dateModified"].startswith("2026-07-07")
    assert any("clinvar" in c for c in ent["citation"])


def test_faq_answers_the_pathogenicity_question():
    faq = next(n for n in VJ.build_jsonld(_REC)["@graph"] if n["@type"] == "FAQPage")
    qs = {q["name"]: q["acceptedAnswer"]["text"] for q in faq["mainEntity"]}
    lead_q = "Is ACTA1 p.Pro309Ala (c.925C>G) pathogenic?"
    assert lead_q in qs and "Pathogenic" in qs[lead_q]
    # only-answerable questions emitted: condition + frequency + inheritance + in-silico
    assert any("condition" in q for q in qs)
    assert any("common" in q for q in qs)
    assert any("inherited" in q for q in qs)


def test_faq_omits_unanswerable_questions():
    bare = {"canonical_slug": "g-p-x", "variation_id": "1", "gene_symbol": "G",
            "hgvs_p": "p.Ala1Val", "classification": "Pathogenic", "submitter_count": 1,
            "conditions": [], "review_status": ""}
    faq = next(n for n in VJ.build_jsonld(bare)["@graph"] if n["@type"] == "FAQPage")
    qs = [q["name"] for q in faq["mainEntity"]]
    assert any("pathogenic" in q.lower() for q in qs)     # always answerable
    assert not any("common" in q for q in qs)             # no gnomad → omitted
    assert not any("condition" in q for q in qs)          # no conditions → omitted
