"""Collector logic that shapes the page's primary condition + submitter metric.
Pure functions — no biobtree. Gate the benchmark-2026-07 fixes before any regen."""
from sugivariant.collect import _order_conditions


def _conds(*names):
    return [{"mondo_id": f"MONDO:{i}", "name": n} for i, n in enumerate(names)]


def test_order_conditions_leads_with_vcep_disease():
    # TP53 R248Q case: arbitrary MONDO order buries the VCEP disease at the end.
    conds = _conds("gastric cancer", "breast carcinoma", "Li-Fraumeni syndrome")
    vcep = [{"disease": "Li-Fraumeni syndrome", "assertion": "Pathogenic"}]
    out = _order_conditions(conds, vcep)
    assert out[0]["name"] == "Li-Fraumeni syndrome"
    # the rest keep their original relative order (stable)
    assert [c["name"] for c in out[1:]] == ["gastric cancer", "breast carcinoma"]


def test_order_conditions_case_insensitive_match():
    conds = _conds("gastric cancer", "li-fraumeni SYNDROME")
    vcep = [{"disease": "Li-Fraumeni syndrome"}]
    assert _order_conditions(conds, vcep)[0]["name"] == "li-fraumeni SYNDROME"


def test_order_conditions_noop_without_vcep():
    # Most genes have no VCEP → order is untouched (and somatic no-framing preserved).
    conds = _conds("acute myeloid leukemia", "myelodysplastic syndrome")
    assert _order_conditions(conds, []) == conds
    assert _order_conditions(conds, None) == conds


def test_order_conditions_single_or_empty():
    assert _order_conditions([], [{"disease": "x"}]) == []
    one = _conds("only one")
    assert _order_conditions(one, [{"disease": "other"}]) == one
