"""Deterministic variant-enrichment logic: missense parsing, frequency banding,
submitter consensus, and the cross-source concordance verdict. Pure functions —
no biobtree. These gate the derived claims before any regen."""
from sugivariant.enrich import (missense_short, protein_position, _freq_band,
                                  submitter_consensus, concordance, residue_hotspot)


def test_missense_short_and_position():
    assert missense_short("p.Pro309Ala") == "P309A"
    assert missense_short("p.Arg14Ter") == "R14*"
    assert missense_short("p.Lys328del") is None      # not a missense → AM n/a
    assert missense_short("p.Asn117fs") is None
    assert protein_position("p.Pro309Ala") == 309
    assert protein_position("p.Lys328del") == 328


def test_freq_band():
    assert _freq_band("", True) == "absent from gnomAD"
    assert _freq_band("0.00005", False).startswith("ultra-rare")
    assert _freq_band("0.0005", False).startswith("rare")
    assert _freq_band("0.27", False).startswith("common")


def test_submitter_consensus():
    uni = submitter_consensus([{"classification": "Pathogenic"}] * 3)
    assert uni["n"] == 3 and "unanimous" in uni["verdict"]
    conf = submitter_consensus([{"classification": "Pathogenic"}] * 2
                               + [{"classification": "Uncertain significance"}])
    assert "conflicting" in conf["verdict"] and "2 Pathogenic" in conf["verdict"]
    assert submitter_consensus([]) is None


def test_concordance_agree_and_flag():
    # AM (now a dict) likely-pathogenic → the concordant predictor; absence is
    # shown but NOT counted as concordant (audit P2)
    ok = concordance("Pathogenic", {"class": "likely_pathogenic", "score": "0.99"},
                     {"absent": True, "is_common": False, "band": "absent from gnomAD v4.1"})
    assert "concordant" in ok["verdict"].lower() and not ok["flags"]
    assert any("PM2" in ln for ln in ok["lines"])          # rarity framed as PM2-support only
    # common in gnomAD → disagreement flagged
    bad = concordance("Pathogenic", None,
                      {"absent": False, "is_common": True, "band": "common (popmax 20%)"})
    assert bad["flags"] and "disagree" in bad["verdict"].lower()
    # AlphaMissense likely-benign vs ClinVar Pathogenic → discordant flag
    disc = concordance("Pathogenic", {"class": "likely_benign", "score": "0.1"}, None)
    assert disc["flags"]
    # absence ALONE is not concordant evidence (no in-silico predictor) → no verdict
    rare_only = concordance("Pathogenic", None, {"absent": True, "band": "absent from gnomAD v4.1"})
    assert rare_only["verdict"] is None and not rare_only["flags"]


def test_concordance_conservation_line():
    c = concordance("Pathogenic", None, None, None,
                    {"phylop": 8.5, "phastcons": 1.0})
    assert any("highly conserved" in ln and "PP3" in ln for ln in c["lines"])
    assert "concordant" in (c["verdict"] or "").lower()     # conservation is the concordant predictor


def test_plain_summary_wording_and_condition():
    from sugivariant.enrich import plain_summary
    rec = {"gene_symbol": "TP53", "classification": "Pathogenic",
           "review_status": "reviewed by expert panel", "submitter_count": 46,
           "conditions": [{"name": "Li-Fraumeni syndrome"}], "concordance": {"flags": []}}
    s = plain_summary(rec)
    # benchmark 2026-07: no "clinical labs"; leads with the (now VCEP-ordered) condition
    assert "clinical lab" not in s
    assert "46 submitters" in s
    assert "Li-Fraumeni syndrome" in s
    # singular submitter
    rec["submitter_count"] = 1
    assert "1 submitter," in plain_summary(rec) and "submitters" not in plain_summary(rec)


def test_disagreement_flag():
    from sugivariant.enrich import disagreement_flag
    def rec(cls, consensus=None, flags=None):
        return {"classification": cls,
                "concordance": {"consensus": consensus, "flags": flags or []}}
    # pathogenic but predictors lean tolerated → predictor_vs_clinvar
    d = disagreement_flag(rec("Pathogenic",
        {"n_damaging": 0, "total": 2, "unanimous": True, "summary": "all 2 predictors call this tolerated"}))
    assert d["category"] == "predictor_vs_clinvar" and d["severity"] == 3
    # pathogenic + AlphaMissense-benign flag, too few predictors for a consensus → still flagged
    d2 = disagreement_flag(rec("Pathogenic", None, ["AlphaMissense predicts likely-benign"]))
    assert d2["category"] == "predictor_vs_clinvar"
    # conflicting ClinVar but predictors unanimous → resolves_conflicting
    d3 = disagreement_flag(rec("Conflicting classifications of pathogenicity",
        {"n_damaging": 3, "total": 3, "unanimous": True, "summary": "all 3 predictors call this damaging"}))
    assert d3["category"] == "resolves_conflicting" and d3["severity"] == 2
    # predictors split among themselves → predictors_split
    d4 = disagreement_flag(rec("Conflicting classifications of pathogenicity",
        {"n_damaging": 1, "total": 3, "unanimous": False, "summary": "1/3 damaging — mixed"}))
    assert d4["category"] == "predictors_split"
    # concordant pathogenic (all predictors damaging) → no flag
    assert disagreement_flag(rec("Pathogenic",
        {"n_damaging": 3, "total": 3, "unanimous": True, "summary": "all 3 predictors call this damaging"})) is None
    # nothing scored → no flag
    assert disagreement_flag(rec("Pathogenic", None)) is None


def test_disagreement_flag_non_missense_lof():
    from sugivariant.enrich import disagreement_flag
    # conflicting + predicted-LoF in a LoF-intolerant gene → non-missense QC flag
    r = {"classification": "Conflicting classifications of pathogenicity",
         "concordance": {"consensus": None, "flags": []},
         "lof_context": {"label": "nonsense (stop-gain)", "haploinsufficient": True,
                         "constrained": False, "lof_disease_gene": True}}
    d = disagreement_flag(r)
    assert d["category"] == "lof_resolves_conflicting"
    assert "haploinsufficiency" in d["reason"]
    # same LoF but neither haploinsufficient nor constrained → no flag (signal weak)
    r["lof_context"] = {"label": "nonsense (stop-gain)", "lof_disease_gene": False}
    assert disagreement_flag(r) is None
    # LoF-disease gene but ClinVar pathogenic (not conflicting) → no flag (nothing to resolve)
    r2 = dict(r); r2["classification"] = "Pathogenic"
    r2["lof_context"] = {"label": "frameshift", "lof_disease_gene": True, "haploinsufficient": True}
    assert disagreement_flag(r2) is None


def test_residue_hotspot():
    idx = {309: [{"hgvs_p": "p.Pro309Ala", "label": "ACTA1 p.Pro309Ala", "slug": "a", "classification": "Pathogenic"},
                 {"hgvs_p": "p.Pro309Leu", "label": "ACTA1 p.Pro309Leu", "slug": "b", "classification": "Pathogenic"}]}
    hs = residue_hotspot("p.Pro309Ala", idx)
    assert hs["position"] == 309 and len(hs["others"]) == 1        # excludes itself
    assert residue_hotspot("p.Pro309Ala", {309: [{"hgvs_p": "p.Pro309Ala", "slug": "a"}]}) is None


# ── Batch 3 derived logic ────────────────────────────────────────────────────
def test_am_percentile():
    from sugivariant.enrich import am_percentile
    am = {f"V{i}": ("x", str(i / 100)) for i in range(100)}   # scores 0.00..0.99
    p = am_percentile("0.98", am)
    assert p["n"] == 100 and p["top_pct"] <= 3          # 0.98 is near the top
    assert am_percentile("0.5", {"a": ("x", "0.1")}) is None   # too few to rank


def test_structural_context():
    from sugivariant.enrich import structural_context
    intervals = [{"type": "helix", "desc": "", "begin": 300, "end": 320},
                 {"type": "region of interest", "desc": "actin-binding", "begin": 100, "end": 150}]
    sc = structural_context("p.Pro309Ala", intervals)
    assert sc["position"] == 309 and any("helix" in f for f in sc["features"])
    assert structural_context("p.Pro200Ala", intervals) is None   # 200 in no interval


def test_similar_variants_ranks_same_residue_first():
    from sugivariant.enrich import similar_variants
    me = {"canonical_slug": "g-p-pro9ala", "hgvs_p": "p.Pro9Ala", "gene_symbol": "G",
          "variant_type": "single nucleotide variant", "conditions": [{"mondo_id": "MONDO:1"}]}
    recs = [me,
            {"canonical_slug": "g-p-pro9leu", "hgvs_p": "p.Pro9Leu", "gene_symbol": "G",
             "variant_type": "x", "conditions": []},                         # same residue (score 3)
            {"canonical_slug": "g-p-arg5his", "hgvs_p": "p.Arg5His", "gene_symbol": "G",
             "variant_type": "y", "conditions": [{"mondo_id": "MONDO:1"}]}]   # same condition (score 2)
    sim = similar_variants(me, recs)
    assert sim[0]["slug"] == "g-p-pro9leu"          # same-residue ranks first
    assert me["canonical_slug"] not in [s["slug"] for s in sim]


def test_submission_timeline():
    from sugivariant.enrich import submission_timeline
    tl = submission_timeline([{"classification": "Pathogenic", "date": "2002-01-01"},
                              {"classification": "Pathogenic", "date": "2018-06-01"}])
    assert tl["first"] == 2002 and tl["last"] == 2018 and tl["stable"]
    tl2 = submission_timeline([{"classification": "Pathogenic", "date": "2010-01-01"},
                               {"classification": "Uncertain significance", "date": "2020-01-01"}])
    assert not tl2["stable"]


def test_plain_summary_calibrated():
    from sugivariant.enrich import plain_summary
    # strong review (2★) pathogenic → "considered disease-causing", variant condition
    strong = plain_summary({"gene_symbol": "ACTA1", "classification": "Pathogenic",
                            "review_status": "criteria provided, multiple submitters, no conflicts",
                            "submitter_count": 3, "conditions": [{"name": "nemaline myopathy"}]})
    assert "considered disease-causing" in strong and "3 submitters" in strong
    assert "nemaline myopathy" in strong
    # 1★ single-submitter pathogenic → hedged "limited review"
    weak = plain_summary({"gene_symbol": "G", "classification": "Pathogenic",
                          "review_status": "criteria provided, single submitter",
                          "submitter_count": 1, "conditions": []})
    assert "limited review" in weak
    # in-silico-discordant pathogenic → "predictors disagree"
    disc = plain_summary({"gene_symbol": "G", "classification": "Pathogenic",
                          "review_status": "criteria provided, multiple submitters, no conflicts",
                          "concordance": {"flags": ["AlphaMissense predicts likely-benign"]},
                          "submitter_count": 4, "conditions": []})
    assert "predictors disagree" in disc
    # likely-pathogenic keeps the hedge; conflicting reads as disagreement
    assert "likely disease-causing" in plain_summary(
        {"gene_symbol": "G", "classification": "Likely pathogenic", "conditions": []})
    assert "disagree" in plain_summary(
        {"gene_symbol": "G", "classification": "Conflicting classifications of pathogenicity",
         "conditions": []})


# ── gnomAD v4.1 (by coordinate) + SpliceAI (both went live 2026-07-09) ────────
def test_variant_coordinate_picks_grch38():
    from sugivariant.enrich import variant_coordinate
    # both assemblies present; must pick the one whose pos == ClinVar start (GRCh38)
    rec = {"chromosome": "1", "start": 229432266,
           "hgvs_expressions": ["NC_000001.10:g.229568013G>C",   # GRCh37 — must NOT win
                                "NC_000001.11:g.229432266G>C"]}   # GRCh38 — matches start
    assert variant_coordinate(rec) == "1:229432266:G:C"
    # no genomic SNV HGVS (indel) → None
    assert variant_coordinate({"chromosome": "1", "start": 5,
                               "hgvs_expressions": ["NM_x:c.10del"]}) is None


def test_gnomad_band_uses_popmax():
    from sugivariant.enrich import _gnomad_band
    assert _gnomad_band(0.06, 0.08, "nfe").startswith("common")
    assert "too common" in _gnomad_band(0.06, 0.08, "nfe")
    assert _gnomad_band(1e-5, 2e-5, "afr").startswith("very rare")
    assert _gnomad_band(None, None, None).startswith("absent")


def test_spliceai_for_threshold():
    from sugivariant.enrich import spliceai_for
    cache = {"1:100:C:T": {"effect": "acceptor_gain", "score": "0.61"},
             "1:200:A:G": {"effect": "donor_loss", "score": "0.05"}}
    assert spliceai_for("1:100:C:T", cache)["effect"] == "acceptor_gain"
    assert spliceai_for("1:200:A:G", cache) is None      # below 0.2 → not surfaced
    assert spliceai_for("1:999:C:T", cache) is None      # not in cache


# ── #2 batch: PharmGKB / CIViC gating + variant-landscape aggregate ──────────
def test_pharmgkb_for_from_cache():
    from sugivariant.enrich import pharmgkb_for
    cache = {"rs4244285": {"annotations": [{"drugs": "clopidogrel"}, {"drugs": "prasugrel"}],
                           "clinical": [{"chemicals": "clopidogrel", "level": "1A", "type": "Efficacy"}]}}
    p = pharmgkb_for("rs4244285", cache)
    assert "clopidogrel" in p["drugs"] and p["clinical"][0]["level"] == "1A"
    assert pharmgkb_for("rs999", cache) is None
    assert pharmgkb_for(None, cache) is None


def test_civic_gate_blocks_non_cancer_genes():
    from sugivariant.enrich import civic_for
    # has_civic False → no lookup at all (the gate), returns None
    assert civic_for("12345", False) is None


def test_gene_variant_landscape():
    from sugivariant.enrich import gene_variant_landscape
    recs = [{"classification": "Pathogenic", "variant_type": "single nucleotide variant",
             "hgvs_p": "p.Arg48His"},
            {"classification": "Pathogenic", "variant_type": "single nucleotide variant",
             "hgvs_p": "p.Arg48Cys"},                       # same residue 48 → recurrent
            {"classification": "Likely pathogenic", "variant_type": "Deletion",
             "hgvs_p": "p.Lys300del"}]
    lc = gene_variant_landscape(recs)
    assert lc["n"] == 3
    assert lc["by_type"]["single nucleotide variant"] == 2
    assert (48, 2) in lc["recurrent"]                       # residue 48 hit twice
    assert lc["residue_span"] == (48, 300)


# ── MaveDB functional-assay join (staging: hgvs_pro now in projection) ────────
def test_mavedb_for_dedups_by_score_set():
    from sugivariant.enrich import mavedb_for
    cache = {"p.Ala302Thr": [
        {"score": "-1.2", "score_set": "urn:mavedb:81-a-1", "license": "CC0", "title": "SGE"},
        {"score": "0.1", "score_set": "urn:mavedb:81-a-1", "license": "CC0", "title": "SGE"},  # dup set
        {"score": "-0.9", "score_set": "urn:mavedb:81-a-2", "license": "CC0", "title": "HDR"}]}
    m = mavedb_for("p.Ala302Thr", cache)
    assert len(m) == 2                                   # one row per score set
    assert {r["score_set"] for r in m} == {"urn:mavedb:81-a-1", "urn:mavedb:81-a-2"}
    assert mavedb_for("p.Gly1Arg", cache) is None
    assert mavedb_for(None, cache) is None


# ── mechanism narrative (variant → gene → pathway → disease; deterministic) ───
def test_mechanism_narrative_gates():
    from sugivariant.enrich import mechanism_narrative
    rec = {"gene_symbol": "PTEN", "hgvs_p": "p.Arg173Cys", "variant_type": "single nucleotide variant",
           "structural": {"features": ["the phosphatase domain"]},
           "alphamissense": {"class": "likely_pathogenic", "score": "0.99"},
           "conditions": [{"name": "Cowden syndrome 1"}]}
    pw = {"disease_pathways": [{"id": "R-HSA-5674404", "name": "PTEN Loss of Function in Cancer"}],
          "top_function": "protein tyrosine phosphatase activity"}
    s = mechanism_narrative(rec, pw)
    assert "alters PTEN at p.Arg173Cys" in s and "phosphatase domain" in s
    # Reactome pathway is stated as a FACT (annotation), not fabricated as "the mechanism of"
    # the condition (benchmark 2026-07). Per-variant AlphaMissense/SpliceAI live in the
    # Computational card, not this gene-function narrative.
    assert "PTEN Loss of Function in Cancer" in s and "In Reactome" in s
    assert "likely pathogenic" not in s          # AM prediction is NOT re-stated here
    assert "mechanism linked to" not in s        # no fabricated causal disease claim
    assert "thought to act" not in s             # (that label is in the render header)
    # NO anchor (no disease pathway, no MF term) → no narrative (never hand-waves)
    assert mechanism_narrative(rec, {"disease_pathways": [], "top_function": None}) is None


# ── REVEL (agreement signal, ClinGen tiers) + GERP (2026-07-13) ──────────────
def test_revel_band_is_descriptive_not_acmg():
    # §8: REVEL must NOT be labelled with ACMG codes (invites additive PP3/BP4 misuse)
    from sugivariant.enrich import _revel_band
    assert _revel_band(0.955) == "strongly pathogenic-leaning"
    assert _revel_band(0.80) == "moderately pathogenic-leaning"
    assert _revel_band(0.70) == "pathogenic-leaning"
    assert _revel_band(0.40) == "indeterminate"
    assert _revel_band(0.20) == "benign-leaning"
    assert _revel_band(0.10) == "moderately benign-leaning"
    assert _revel_band(0.005) == "strongly benign-leaning"
    # no ACMG code ever leaks out of the band
    for s in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0):
        assert "PP3" not in _revel_band(s) and "BP4" not in _revel_band(s)


def test_revel_is_agreement_not_a_second_vote():
    from sugivariant.enrich import concordance
    am = {"class": "likely_pathogenic", "score": "0.99"}
    revel = {"score": 0.955, "band": "strongly pathogenic-leaning", "direction": "pathogenic"}
    c = concordance("Pathogenic", am, None, revel=revel)
    # REVEL is shown as agreeing, but does NOT add a second concordant predictor
    assert any("REVEL" in ln and "agrees with AlphaMissense" in ln for ln in c["lines"])
    assert "1 independent predictor" in c["verdict"]        # AM only, not 2
    # no ACMG code on the REVEL line (§8)
    assert not any("PP3" in ln or "BP4" in ln for ln in c["lines"])
    # a benign REVEL vs pathogenic AM → "differs from" (still not a counted flag)
    d = concordance("Pathogenic", am, None,
                    revel={"score": 0.1, "band": "moderately benign-leaning", "direction": "benign"})
    assert any("differs from AlphaMissense" in ln for ln in d["lines"])


def test_molecular_consequence():
    from sugivariant.enrich import molecular_consequence as mc
    assert mc({"hgvs_p": "p.Gln592fs"})["type"] == "frameshift"
    assert mc({"hgvs_p": "p.Ser1146Ter"})["type"] == "nonsense"
    assert mc({"hgvs_p": "p.Met1?"})["type"] == "start_loss"
    assert mc({"hgvs_c": "c.588+1G>C"})["type"] == "splice"       # canonical +1
    assert mc({"hgvs_c": "c.2312-1G>A"})["type"] == "splice"      # canonical -1
    assert mc({"hgvs_c": "c.1100+13A>G"}) is None                 # deep intronic, NOT canonical
    assert mc({"hgvs_p": "p.Lys328del"})["type"] == "inframe_indel"
    assert mc({"hgvs_p": "p.Lys328del"})["lof"] is False
    # nonsense/frameshift/splice are all predicted-LoF
    assert mc({"hgvs_p": "p.Gln592fs"})["lof"] is True
    assert mc({"hgvs_p": "p.Arg130Gln"}) is None                  # plain missense → typed by AM, not here


def test_lof_context_uses_haploinsufficiency_and_constraint():
    from sugivariant.enrich import lof_context
    GERM = {"conditions": [{"name": "Cowden syndrome 1"}]}   # germline (non-somatic) condition
    # pLI-constrained + germline → lof_disease_gene via constraint
    lc = lof_context({"consequence": {"label": "nonsense (stop-gain)", "lof": True, "type": "nonsense"},
                      "gene_context": {"constraint": {"loeuf": "0.12", "pli": "0.99"}}, **GERM})
    assert lc["constrained"] is True and lc["lof_disease_gene"] is True
    assert lc["nmd_caveat"] is True                         # nonsense → NMD caveat flagged
    # ClinGen haploinsufficiency=3 catches a LoF-disease gene pLI MISSES (LDLR-like)
    lc2 = lof_context({"consequence": {"label": "frameshift", "lof": True, "type": "frameshift"},
                       "gene_context": {"constraint": {"loeuf": "1.1", "pli": "0.0"},
                                        "dosage": {"haplo": "3"}}, **GERM})
    assert lc2["haploinsufficient"] is True and lc2["lof_disease_gene"] is True
    # §8: a SOMATIC/myeloid condition → mechanism NOT asserted even in a haploinsufficient gene
    lc_som = lof_context({"consequence": {"label": "frameshift", "lof": True, "type": "frameshift"},
                          "gene_context": {"dosage": {"haplo": "3"}},
                          "conditions": [{"name": "Acute myeloid leukemia"}]})
    assert lc_som["haploinsufficient"] is True and lc_som["lof_disease_gene"] is False
    # neither haploinsufficient nor constrained → not a LoF-disease gene
    lc3 = lof_context({"consequence": {"label": "frameshift", "lof": True, "type": "frameshift"},
                       "gene_context": {"constraint": {"loeuf": "1.2", "pli": "0.01"},
                                        "dosage": {"haplo": "0"}}, **GERM})
    assert lc3["lof_disease_gene"] is False
    # in-frame (lof False) → no context
    assert lof_context({"consequence": {"label": "in-frame indel", "lof": False}}) is None


def test_clingen_criteria_guards():
    from sugivariant.enrich import clingen_criteria
    assert clingen_criteria(None) is None       # no id → no fabrication, no biobtree call
    assert clingen_criteria("") is None


def test_am_isoform_mismatch():
    from sugivariant.enrich import am_isoform_mismatch
    # AM's transcript numbers the residue 309; ClinVar p.HGVS says 334 → mismatch flag
    mm = am_isoform_mismatch({"short": "P309A"}, "p.Pro334Ala")
    assert mm == {"am": "P309A", "clinvar": "P334A"}
    # same residue → no flag
    assert am_isoform_mismatch({"short": "R130Q"}, "p.Arg130Gln") is None
    # missing data → no flag (never fabricate)
    assert am_isoform_mismatch(None, "p.Arg130Gln") is None
    assert am_isoform_mismatch({"short": "R130Q"}, None) is None
    assert am_isoform_mismatch({"short": None}, "p.Arg130Gln") is None


def test_nonmissense_banner_is_explicit_not_reassuring():
    from sugivariant.enrich import concordance
    # non-missense (am=None), no conservation → the panel must say predictors DON'T APPLY
    c = concordance("Pathogenic", None, None)
    assert any("do not apply" in ln for ln in c["lines"])
    # and it must never read as a benign/absent-predictor signal
    assert not any("not scored" in ln for ln in c["lines"])


def test_conservation_gerp_and_nonmissense_counts():
    from sugivariant.enrich import concordance
    # non-missense (am=None): conservation is the primary computational signal → counts
    c = concordance("Pathogenic", None, None, None, {"phylop": 7.4, "gerp": 4.9, "phastcons": 1})
    assert any("GERP 4.9" in ln for ln in c["lines"])
    assert "concordant" in (c["verdict"] or "").lower()
    # missense: conservation is only an agreement check (not counted → AM alone drives)
    m = concordance("Pathogenic", {"class": "likely_pathogenic", "score": "0.9"}, None, None,
                    {"phylop": 7.4, "gerp": 4.9, "phastcons": 1})
    assert any("agrees" in ln for ln in m["lines"])
    assert "1 independent predictor" in m["verdict"]         # AM only


# ── SaProt (orthogonal protein-LM predictor) + in-silico consensus (2026-07-13) ─
def test_saprot_agreement_and_consensus():
    from sugivariant.enrich import concordance
    am = {"class": "likely_pathogenic", "score": "0.99"}
    revel = {"score": 0.955, "band": "strongly pathogenic-leaning", "direction": "pathogenic"}
    esm = {"llr": -14.5, "damaging": True}
    c = concordance("Pathogenic", am, None, revel=revel, saprot=esm)
    # SaProt shown as an agreement signal (raw LLR), not a counted vote
    assert any("SaProt LLR -14.5" in ln and "agrees with AlphaMissense" in ln for ln in c["lines"])
    assert "1 independent predictor" in c["verdict"]        # AM alone carries the ACMG weight
    # consensus summary across the 3 independent missense predictors
    assert any("consensus" in ln and "all 3" in ln and "damaging" in ln for ln in c["lines"])
    # a tolerated SaProt vs pathogenic AM → "differs from" + mixed consensus
    d = concordance("Pathogenic", am, None,
                    saprot={"llr": -0.1, "damaging": False})
    assert any("differs from AlphaMissense" in ln for ln in d["lines"])
    assert any("consensus" in ln and "1/2" in ln for ln in d["lines"])
