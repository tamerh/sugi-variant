"""Variant collector — a ClinVar variation id → structured page record.

Anchor dataset = ClinVar (enumerated gene-first via >>hgnc>>clinvar; see the gene
§6 collector). Per-variant `entry(vid, "clinvar")` gives HGVS/rsID/coords/
classification + the per-submitter table + phenotype MONDO ids; cross-links add
the linked conditions and (when the gene has a ClinGen VCEP) the expert-panel
ACMG assertion. Positional in-silico predictions (AlphaMissense/SpliceAI) are
optional annotations, joined by substitution/position — never page generators.
"""
from sugibiobtree import entry, map_all
from sugivariant.slug import parse_hgvs, variant_slugs
from sugivariant import enrich as EN

# Classifications we build pages for (Phase 1 gate: the high-value + most-searched
# slice; the benign long tail is excluded — see VARIANT_PAGES_SPEC.md §3).
BUILD_CLASSES = {
    "pathogenic", "likely pathogenic", "pathogenic/likely pathogenic",
    "conflicting classifications of pathogenicity",
}


def _attrs(e):
    a = e.get("Attributes") or {}
    return a.get("Clinvar") or a.get("ClinVar") or {}


def should_build(classification):
    return (classification or "").strip().lower() in BUILD_CLASSES


def _order_conditions(conditions, vcep):
    """Lead with the ClinGen VCEP / expert-panel disease so the primary condition
    is the DEFINING one, not an arbitrary first >>clinvar>>mondo edge. conditions[0]
    drives the Summary headline, the mechanism narrative, the patient digest routing
    and the condition links, so the order here propagates everywhere (benchmark
    2026-07: TP53 Li-Fraumeni hotspots headlined 'gastric cancer' with Li-Fraumeni
    buried at 24/28). Stable otherwise; a no-op when there is no VCEP (most genes)."""
    vcep_diseases = {(c.get("disease") or "").strip().lower()
                     for c in (vcep or []) if c.get("disease")}
    if not vcep_diseases or len(conditions) < 2:
        return conditions
    # False (0) sorts before True (1); sort is stable → VCEP disease(s) first, rest as-is
    return sorted(conditions,
                  key=lambda c: (c.get("name") or "").strip().lower() not in vcep_diseases)


def collect(variation_id):
    """Full page record for a ClinVar variation id, or None if it's not a
    buildable classification / can't resolve a slug."""
    v = _attrs(entry(variation_id, "clinvar"))
    if not v:
        return None
    cls = v.get("germline_classification")
    if not should_build(cls):
        return None

    gene = v.get("gene_symbol")
    c_form, p_form = parse_hgvs(v.get("name"))
    canonical, slugs = variant_slugs(gene, c_form, p_form)
    if not canonical:
        return None   # no parseable HGVS → can't make a deterministic URL

    # Linked conditions (MONDO → disease pages). phenotype_ids mixes OMIM/MedGen/
    # MONDO; the >>clinvar>>mondo edge gives the clean linkable set + names.
    conditions = [{"mondo_id": m.get("id"), "name": m.get("name")}
                  for m in map_all(variation_id, ">>clinvar>>mondo") if m.get("id")]

    # ClinGen VCEP expert-panel ACMG assertion (authority tier above raw ClinVar);
    # only a few hundred genes have a VCEP → usually empty.
    vcep = [{"id": c.get("id"), "assertion": c.get("assertion") or c.get("classification"),
             "panel": c.get("vcep") or c.get("panel"), "disease": c.get("disease")}
            for c in map_all(variation_id, ">>clinvar>>clingen_variant") if c.get("id")]
    # Attach the full applied-ACMG-code record (authority tier) per VCEP assertion —
    # one entry() call, and only for the ~3% of variants that HAVE a VCEP.
    for c in vcep:
        c["criteria"] = EN.clingen_criteria(c["id"])

    conditions = _order_conditions(conditions, vcep)

    subs = v.get("submissions") or []
    # Distinct submitters, not raw submission rows — a lab with multiple SCVs counts
    # once (benchmark 2026-07: '57 submissions' was shown as '57 clinical labs').
    n_submitters = len({s.get("submitter_name") for s in subs if s.get("submitter_name")})
    return {
        "variation_id": variation_id,
        "gene_symbol": gene,
        "hgnc_id": v.get("hgnc_id"),
        "name": v.get("name"),
        "hgvs_c": c_form,
        "hgvs_p": p_form,
        "hgvs_expressions": v.get("hgvs_expressions") or [],
        "rsid": v.get("dbsnp_id"),
        "chromosome": v.get("chromosome"),
        "start": v.get("start"),
        "stop": v.get("stop"),
        "assembly": v.get("assembly"),
        "variant_type": v.get("type"),
        "classification": cls,
        "review_status": v.get("review_status"),
        "last_evaluated": v.get("last_evaluated"),
        "submissions": [{"submitter": s.get("submitter_name"),
                         "classification": s.get("classification"),
                         "review_status": s.get("review_status"),
                         "method": s.get("method_type"),
                         "date": s.get("date_last_evaluated")} for s in subs],
        "submitter_count": n_submitters or len(subs),
        "conditions": conditions,
        "vcep": vcep,
        "phenotype_list": v.get("phenotype_list") or [],
        "canonical_slug": canonical,
        "slugs": slugs,
    }


def _is_pathogenic(cls):
    c = (cls or "").lower()
    return "pathogenic" in c and "conflict" not in c


def build_position_index(recs):
    """{protein_position: [{hgvs_p, label, slug, classification}]} over the gene's
    PATHOGENIC variants — powers the same-residue hotspot context."""
    idx = {}
    for r in recs:
        pos = EN.protein_position(r.get("hgvs_p"))
        if pos is None or not _is_pathogenic(r.get("classification")):
            continue
        idx.setdefault(pos, []).append({"hgvs_p": r.get("hgvs_p"),
                                        "label": f"{r['gene_symbol']} {r['hgvs_p']}",
                                        "slug": r["canonical_slug"],
                                        "classification": r["classification"]})
    return idx


def attach_enrichment(rec, ctx=None):
    """Add the deterministic enrichment layers to a collected record, in place:
    AlphaMissense (from the per-gene `am` cache), gnomAD frequency, cross-source
    concordance, submitter consensus, same-residue hotspot. ctx carries the
    per-gene caches so this stays cheap across a gene's variants."""
    ctx = ctx or {}
    coord = EN.variant_coordinate(rec)
    rec["coordinate"] = coord
    # All three are coordinate-keyed (entry-by-coordinate, not map — audit Tier 1/2)
    am = EN.alphamissense_for(coord)              # by genomic key → fixes isoform-mismatch misses
    gnomad = EN.gnomad_frequency(rec)             # by coordinate via entry() → fixes false-Absent
    conservation = EN.conservation_for(coord)     # phyloP / GERP / phastCons
    revel = EN.revel_for(coord)                   # ensemble missense predictor (agreement signal)
    # SaProt — protein-LM predictor, keyed by the uniprot:protein_variant AM gives us
    saprot = EN.saprot_for(am.get("uniprot"), am.get("short")) if am else None
    rec["alphamissense"] = am
    rec["am_isoform_mismatch"] = EN.am_isoform_mismatch(am, rec.get("hgvs_p"))
    rec["gnomad"] = gnomad
    rec["conservation"] = conservation
    rec["revel"] = revel
    rec["saprot"] = saprot
    rec["spliceai"] = EN.spliceai_for(coord, ctx.get("spliceai") or {})
    rec["consensus"] = EN.submitter_consensus(rec.get("submissions") or [])
    rec["concordance"] = EN.concordance(rec.get("classification"), am, gnomad,
                                        rec.get("spliceai"), conservation, revel, saprot)
    rec["hotspot"] = EN.residue_hotspot(rec.get("hgvs_p"), ctx.get("positions"))
    # Batch 3 — per-variant derived (in-memory from the per-gene caches)
    am_map = ctx.get("am") or {}          # per-gene AM distribution (for the percentile)
    rec["am_percentile"] = (EN.am_percentile(am["score"], am_map) if (am and am_map) else None)
    rec["structural"] = EN.structural_context(rec.get("hgvs_p"), (ctx.get("structure") or {}).get("intervals"))
    rec["similar"] = EN.similar_variants(rec, ctx.get("recs") or [])
    rec["timeline"] = EN.submission_timeline(rec.get("submissions") or [])
    # patient digest routed via the VARIANT's own condition (audit P1) — cached
    # per condition-MONDO; None for somatic/unannotated conditions (no germline
    # framing). plain_summary is confidence-calibrated (audit P2).
    rec["digest"] = EN.condition_digest(rec.get("conditions"), ctx.setdefault("digest_cache", {}))
    rec["plain"] = EN.plain_summary(rec)
    rec["pathways"] = ctx.get("pathways")
    rec["mechanism"] = EN.mechanism_narrative(rec, ctx.get("pathways"))
    rec["pharmgkb"] = EN.pharmgkb_for(rec.get("rsid"), ctx.get("pharmgkb") or {})
    rec["civic"] = EN.civic_for(rec.get("variation_id"), ctx.get("has_civic"))
    rec["mavedb"] = EN.mavedb_for(rec.get("hgvs_p"), ctx.get("mavedb") or {})
    # per-gene context is shared (attach the references for the renderer)
    rec["gene_context"] = ctx.get("gene_context")
    rec["structure"] = ctx.get("structure")
    rec["panels"] = ctx.get("panels")
    # patient condition links (GARD + trials), routed via the variant's condition
    conds = rec.get("conditions") or []
    rec["condition_links"] = (EN.condition_links(conds[0].get("mondo_id"), ctx.setdefault("trials_cache", {}))
                              if conds else {})
    return rec


def enumerate_gene(hgnc_id, cap_pages=None):
    """[(variation_id, classification, name)] for a gene's ClinVar variants that
    pass the build gate — the enumeration a variant build fans over.

    cap_pages=None (default) paginates to completion: the variant corpus MUST NOT
    truncate a gene's buildable set — the old 60-page cap silently dropped ~half of
    BRCA2 (3,072 seen vs 11,602 real). Pass an int only to bound it deliberately
    (e.g. a fast smoke test)."""
    out = []
    for r in map_all(hgnc_id, ">>hgnc>>clinvar", cap=cap_pages):
        if should_build(r.get("germline_classification")):
            out.append((r.get("id"), r.get("germline_classification"), r.get("name")))
    return out
