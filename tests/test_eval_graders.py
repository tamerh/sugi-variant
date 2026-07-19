"""Offline unit tests for the deterministic eval graders (tools/eval/graders.py)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools", "eval"))
import graders as G


def test_base_code():
    assert G.base_code("PM2_Supporting") == "PM2"
    assert G.base_code("PVS1") == "PVS1"
    assert G.base_code("PS4_Moderate") == "PS4"


def test_criteria_prf_base_vs_strict():
    truth = ["PM2_Supporting", "PS3", "PM3", "PVS1"]
    pred = ["PM2", "PS3", "PM3"]                  # missing PVS1; PM2 sans modifier
    base = G.criteria_prf(pred, truth, mode="base")
    assert base["tp"] == 3 and base["fn"] == 1 and base["fp"] == 0
    assert abs(base["recall"] - 0.75) < 1e-9 and base["precision"] == 1.0
    strict = G.criteria_prf(pred, truth, mode="strict")   # 'PM2' != 'PM2_Supporting'
    assert strict["tp"] == 2 and "PM2" not in {G.base_code(c) for c in []}
    # perfect match
    perfect = G.criteria_prf(truth, truth, mode="base")
    assert perfect["f1"] == 1.0


def test_classification_match():
    assert G.classification_match("Pathogenic", "Pathogenic")["exact"]
    m = G.classification_match("Likely pathogenic", "Pathogenic")
    assert not m["exact"] and m["within_one"]          # adjacent bins
    far = G.classification_match("Benign", "Pathogenic")
    assert not far["within_one"]
    vus = G.classification_match("Uncertain significance", "VUS")
    assert vus["exact"]


def test_af_match():
    assert G.af_match(1.0e-5, 1.1e-5, rel_tol=0.25)
    assert not G.af_match(1.0e-3, 1.0e-5)
    assert G.af_match(None, None)                       # both absent
    assert not G.af_match(None, 0.01)


def test_extract_citations():
    text = "See PMID: 9450897 and pmid 2574002; also doi:10.1038/s41586-023-06291-2."
    c = G.extract_citations(text)
    assert "9450897" in c["pmids"] and "2574002" in c["pmids"]
    assert any(d.startswith("10.1038/") for d in c["dois"])


def test_harmful_flags():
    assert G.harmful_flags("Benign", "Pathogenic")["wrong_direction"]
    assert not G.harmful_flags("Likely pathogenic", "Pathogenic")["wrong_direction"]
    assert not G.harmful_flags("Uncertain significance", "Pathogenic")["wrong_direction"]
