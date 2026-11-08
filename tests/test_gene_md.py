"""The gene hub's Markdown twin: exact counts, capped lists that say so.

A resolver hit of kind "gene" had no page to fetch, which made a generic client report
a fetch error for a valid answer. The twin closes that, and the risk it introduces is
the opposite one: a list so long it is useless, or a cap that silently hides variants.
"""

from sugivariant.render import render_gene_md


def _ctx():
    rows = [{"slug": f"brca1-p-v{i}", "hgvs_p": f"p.Val{i}Ala", "hgvs_c": None,
             "classification": "Pathogenic", "stars": 3, "primary_condition": "Breast cancer",
             "flag": None} for i in range(80)]
    return {
        "gene": "BRCA1", "total": 11591, "flagged": 4521, "cap": 200,
        "groups": [{"cls": "Pathogenic", "total": 2506, "shown": rows},
                   {"cls": "Uncertain significance", "total": 7000, "shown": rows}],
        "dgroups": {"predictor_vs_clinvar": rows}, "dcounts": {"predictor_vs_clinvar": 2221},
        "dcap": 50, "cats": ["predictor_vs_clinvar"],
        "meta": {"predictor_vs_clinvar": {"label": "Predictors vs ClinVar"}},
        "sv_rows": [{"nssv": "NSSV15120701", "svtype": "DEL", "significance": "Pathogenic"}],
        "sv_total": 628, "sv_cap": 100,
    }


def test_summary_states_the_exact_count():
    md = render_gene_md(_ctx(), url="https://sugi.bio/variant/gene/BRCA1")
    assert md.startswith("# BRCA1 — 11,591 ClinVar germline variants")
    assert "| Variants in this corpus | 11,591 |" in md
    assert "| Pathogenic | 2,506 |" in md
    assert "https://sugi.bio/variant/gene/BRCA1" in md


def test_capped_lists_say_how_many_they_omitted():
    md = render_gene_md(_ctx())
    assert "showing 50 of 2,506" in md
    assert md.count("| [BRCA1 p.Val") == 100          # 50 per group, two groups
    assert "_Showing 20 of 2,221._" in md             # disagreement cap


def test_the_sections_are_addressable():
    from app import _md_headings, _md_section
    md = render_gene_md(_ctx())
    assert _md_headings(md) == ["Summary", "Variants", "Evidence disagreements",
                                "Structural variants", "Sources"]
    assert _md_section(md, "Summary").startswith("## Summary")
    assert _md_section(md, "Variants").count("###") == 2
    assert _md_section(md, "Nope") is None


def test_a_gene_with_no_flagged_or_sv_rows_still_renders():
    ctx = _ctx()
    ctx.update({"dgroups": {}, "dcounts": {}, "sv_rows": [], "sv_total": 0})
    md = render_gene_md(ctx)
    assert "## Evidence disagreements" not in md
    assert "## Structural variants" not in md
    assert "## Summary" in md


def test_counts_state_the_basis_before_the_number():
    """"How many pathogenic BRCA1 variants" must have one defensible answer. The unit
    is the page, not the ClinVar VCV, and the block has to say so."""
    ctx = _ctx()
    ctx["counts"] = {
        "count_basis": "distinct variant pages: one page per protein change, from "
                       "ClinVar VCVs whose HGVS transcript gene is BRCA1",
        "clinvar_build": "db_v56", "release": "v1.8.4",
        "total": 11591, "by_class": {"Pathogenic": 2506, "Benign": 719},
        "unclassified": 0, "counts_are_disjoint": True,
        "structural_variants": {"total": 628, "separate_corpus": True},
        "excluded": ["structural variants and CNVs: a separate corpus"],
    }
    md = render_gene_md(ctx)
    assert "## Counts {#counts}" in md
    assert "| Count basis |" in md and "one page per protein change" in md
    assert "| ClinVar build | db_v56 |" in md
    assert "| Total variants | 11,591 |" in md
    assert "| Classes sum to total | yes |" in md
    assert "| Structural variants (separate corpus, not in total) | 628 |" in md
    assert "**Not counted here:**" in md


def test_a_non_disjoint_count_says_so():
    ctx = _ctx()
    ctx["counts"] = {"count_basis": "x", "total": 100, "by_class": {"A": 10},
                     "unclassified": 0, "counts_are_disjoint": False,
                     "structural_variants": {"total": 0}, "excluded": []}
    assert "| Classes sum to total | no |" in render_gene_md(ctx)
