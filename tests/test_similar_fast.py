"""Serve-time similar-variant mesh: the fast index must match the old semantics.

The record build no longer computes `similar`; the app builds `similar_index(recs)`
once per gene and calls `similar_variants_from_index()` per page. These tests pin
the scoring, the tie-break, the relation label and the self-exclusion.
"""
from sugivariant.enrich import (
    similar_index, similar_variants, similar_variants_from_index,
)


def _rec(slug, hgvs_p="p.Arg10Gln", conds=(), vtype="single nucleotide variant",
         review="criteria provided, single submitter", cls="Pathogenic", gene="G"):
    return {"canonical_slug": slug, "gene_symbol": gene, "hgvs_p": hgvs_p,
            "hgvs_c": "c.29G>A", "variant_type": vtype,
            "conditions": [{"mondo_id": m} for m in conds],
            "review_status": review, "classification": cls}


def test_relation_order_and_self_exclusion():
    me = _rec("g-p-arg10gln", conds=["A"])
    same_res = _rec("g-p-arg10gln-x", conds=["B"],
                    review="reviewed by expert panel", cls="Likely pathogenic")
    same_cond = _rec("g-p-arg20gln", hgvs_p="p.Arg20Gln", conds=["A"])
    same_type = _rec("g-p-arg30gln", hgvs_p="p.Arg30Gln", conds=["C"])
    unrelated = _rec("g-p-arg40gln", hgvs_p="p.Arg40Gln", conds=["D"],
                     vtype="deletion")
    recs = [me, same_res, same_cond, same_type, unrelated]
    out = similar_variants_from_index(me, similar_index(recs))
    slugs = [x["slug"] for x in out]
    assert "g-p-arg10gln" not in slugs
    assert slugs == ["g-p-arg10gln-x", "g-p-arg20gln", "g-p-arg30gln"]
    assert [x["relation"] for x in out] == [
        "same residue", "same condition", "same type"]


def test_tie_break_is_stars_then_interest_then_slug():
    me = _rec("me")
    low = _rec("a-low", hgvs_p="p.Gly99Asp", cls="Benign",
               review="no assertion criteria provided")
    high = _rec("z-high", hgvs_p="p.Gly99Asp", cls="Pathogenic",
                review="reviewed by expert panel")
    mid = _rec("m-mid", hgvs_p="p.Gly99Asp", cls="Uncertain significance")
    out = similar_variants_from_index(me, similar_index([me, low, high, mid]))
    assert [x["slug"] for x in out] == ["z-high", "m-mid", "a-low"]
    assert out[0]["stars"] == 3


def test_limit_is_eight_and_wrapper_matches():
    me = _rec("me")
    recs = [me] + [_rec(f"g-p-gly{i}asp", hgvs_p="p.Gly99Asp") for i in range(12)]
    fast = similar_variants_from_index(me, similar_index(recs))
    slow = similar_variants(me, recs)
    assert len(fast) == 8
    assert [x["slug"] for x in fast] == [x["slug"] for x in slow]
    assert [x["relation"] for x in fast] == ["same type"] * 8


def test_same_residue_scores_more_than_same_condition():
    me = _rec("me", conds=["A"])
    residue = _rec("res", conds=["B"])
    condition = _rec("cond", hgvs_p="p.Arg20Gln", conds=["A"])
    out = similar_variants_from_index(me, similar_index([me, condition, residue]))
    assert [x["slug"] for x in out] == ["res", "cond"]
