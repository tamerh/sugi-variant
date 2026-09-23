"""Deterministic markdown renderer for a variant record — NO model. Mirrors the
gene/drug/disease renderers: every fact verbatim from the collected record.
"""
import re

from sugivariant import links
from sugivariant.util import table

# A delins/ins/dup can carry a huge inserted sequence (satellite DNA, etc.). For
# DISPLAY, collapse a long run to a short preview + length so it doesn't blow out
# the page title/heading. The full HGVS is still shown (wrapped) in the identity
# section, and resolution is unaffected (this is presentation only).
_LONG_INS_RE = re.compile(r"((?:delins|ins|dup))([A-Za-z]{25,})")


def short_hgvs(s, keep=10):
    if not s:
        return s
    m = _LONG_INS_RE.search(s)
    if not m:
        return s
    seq = m.group(2)
    nt = set(seq.upper()) <= set("ACGTN")
    n, unit = (len(seq), "bp") if nt else (len(seq) // 3, "aa")
    return s[:m.start(2)] + f"{seq[:keep]}…[{n} {unit}]"

# Data sources + attribution. AlphaMissense (CC BY-NC-SA 4.0) and REVEL (ODbL) legally
# REQUIRE attribution; the rest are credited as good practice. Everything here is
# usable in this FREE/non-commercial product (see docs/variant-product-strategy).
# Single source of truth — the HTML template reads the same list via env.globals.
# (name, license/credit, biobtree group). The group joins to the live per-dataset
# build date (enrich.dataset_versions) so provenance shows the HONEST version.
DATA_SOURCES = [
    ("ClinVar", "NCBI — public domain", "clinvar"),
    ("gnomAD v4.1", "Broad Institute", "gnomad_variant"),
    ("AlphaMissense", "Cheng et al. 2023, © Google DeepMind — CC BY-NC-SA 4.0", "alphamissense"),
    ("REVEL", "Ioannidis et al. 2016 — ODbL", "revel"),
    ("SaProt", "Su et al. 2023 — MIT", "saprot"),
    ("SpliceAI", "Illumina — free for non-commercial use", "spliceai"),
    ("Conservation (phyloP/phastCons/GERP)", "UCSC", "conservation"),
    ("ClinGen", "ClinGen", "clingen"),
    ("UniProt", "CC BY 4.0", "uniprot"),
    ("Reactome / GO", "CC BY 4.0", "reactome"),
    ("MaveDB", "per-record license", "mavedb"),
    ("CIViC", "CC0", "civic_variant"),
    ("Orphanet", "Orphanet", "orphanet"),
]

# name -> (citation, homepage) for the About sources table (mirrors the Sugi Atlas table).
SOURCE_REFS = {
    "ClinVar":                              ("Landrum et al., NAR 2018",            "https://www.ncbi.nlm.nih.gov/clinvar/"),
    "gnomAD v4.1":                          ("Chen et al., Nature 2024",            "https://gnomad.broadinstitute.org/"),
    "AlphaMissense":                        ("Cheng et al., Science 2023",          "https://github.com/google-deepmind/alphamissense"),
    "REVEL":                                ("Ioannidis et al., AJHG 2016",         "https://sites.google.com/site/revelgenomics/"),
    "SaProt":                               ("Su et al., ICLR 2024",                "https://github.com/westlake-repl/SaProt"),
    "SpliceAI":                             ("Jaganathan et al., Cell 2019",        "https://github.com/Illumina/SpliceAI"),
    "Conservation (phyloP/phastCons/GERP)": ("Pollard et al. 2010; Davydov et al. 2010", "https://genome.ucsc.edu/"),
    "ClinGen":                              ("Rehm et al., NEJM 2015",              "https://clinicalgenome.org/"),
    "UniProt":                              ("UniProt Consortium, NAR 2025",        "https://www.uniprot.org/"),
    "Reactome / GO":                        ("Milacic et al. 2024; Ashburner et al. 2000", "https://reactome.org/"),
    "MaveDB":                               ("Esposito et al., Genome Biol 2019",   "https://www.mavedb.org/"),
    "CIViC":                                ("Griffith et al., Nat Genet 2017",     "https://civicdb.org/"),
    "Orphanet":                             ("Rath et al., 2012",                   "https://www.orpha.net/"),
}


def source_refs():
    """[(name, citation, url)] for the About sources table."""
    return [(name, *SOURCE_REFS.get(name, ("", ""))) for name, _lic, _grp in DATA_SOURCES]


def data_provenance():
    """[(name, license, built_date)] — sources joined to their live biobtree build
    date (None if unknown). The auditability payload: what, under what licence, as of when."""
    from sugivariant.enrich import dataset_versions
    ver = dataset_versions()
    return [(name, lic, ver.get(group)) for name, lic, group in DATA_SOURCES]


def data_asof():
    """Freshest dataset build date (YYYY-MM) across sources, or None — the single
    'data current as of' stamp for the per-page provenance line."""
    dates = [b for _n, _lic, b in data_provenance() if b]
    return max(dates)[:7] if dates else None


def attribution_md():
    return "; ".join(f"{n} ({lic}" + (f", built {b}" if b else "") + ")"
                     for n, lic, b in data_provenance())

# ClinVar review status → gold-star tier (the standard 0-4 confidence scale).
_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
}


def _review(review_status):
    """Named review-status tier for markdown prose (no ★ rating metaphor)."""
    from sugivariant.enrich import review_tier
    return review_tier(review_status)["label"]


def _label(v):
    """Human page label: 'ACTA1 p.Gln248Lys (c.742C>A)'. Long insertions are
    collapsed for display (see short_hgvs)."""
    g = v.get("gene_symbol") or ""
    p, c = short_hgvs(v.get("hgvs_p")), short_hgvs(v.get("hgvs_c"))
    core = f"{g} {p}" if p else f"{g} {c}" if c else g
    return core + (f" ({c})" if p and c else "")


def declarative(v):
    """Answer-first lead — the self-contained factual block an assistant lifts to
    answer 'is GENE pChange pathogenic?'. Includes the in-silico + population
    differentiators when present."""
    label = _label(v)
    cls = v.get("classification") or "classified"
    gene = v.get("gene_symbol")
    rs = v.get("review_status") or ""
    lead = f"**{label}** is classified **{cls}** in {gene} (ClinVar, {_review(rs)}"
    n_sub = v.get("submitter_count") or 0   # same count as At-a-glance (audit P2b)
    if n_sub:
        lead += f", {n_sub} submitter" + ("s" if n_sub != 1 else "")
    if v.get("last_evaluated"):
        lead += f"; last evaluated {v['last_evaluated']}"
    lead += ")."
    tail = []
    am = v.get("alphamissense")
    if am and am.get("class"):
        tail.append(f"AlphaMissense: {am['class'].replace('_', ' ')} ({am['score']})")
    g = v.get("gnomad")
    if g:
        tail.append(g["band"])
    cond = (v.get("conditions") or [{}])[0].get("name")
    if cond:
        tail.append(f"associated with {cond}")
    return lead + (" " + "; ".join(s[0].upper() + s[1:] for s in tail) + "." if tail else "")


def declarative_plain(v):
    """The lead sentence, markdown-stripped — for the frontmatter/meta
    description (SEO + AI-citation snippet)."""
    import re
    s = declarative(v)
    s = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", s)   # link → label
    s = re.sub(r"[*`_]", "", s)
    return re.sub(r"\s+", " ", s).strip()


def _patient_zone(v):
    """The 'For patients & families' section — plain-language + condition digest,
    every line reference-framed, ending in a counselor CTA."""
    d = v.get("digest") or {}   # the VARIANT's-condition digest (germline only)
    L = ["", "## For patients & families {#patients}", "",
         v.get("plain") or "", "",
         "*Reference information, not medical advice, and not a prediction for any "
         "individual. A genetic counselor or clinician can interpret what it means "
         "for you and your family.*"]
    facts = []
    # Inheritance/onset/prevalence come ONLY from the variant's condition digest
    # (audit P1) — so somatic/acquired conditions, which have no germline Orphanet
    # entry, never get germline framing.
    if d.get("inheritance"):
        facts.append(("Typical inheritance (of this condition)", ", ".join(d["inheritance"][:3])))
    if d.get("onset"):
        facts.append(("Age of onset (of this condition)", ", ".join(d["onset"])))
    if d.get("prevalence"):
        facts.append(("Prevalence", d["prevalence"]))
    cl = v.get("condition_links") or {}
    if cl.get("trial_count"):
        facts.append(("Clinical trials",
                      f"{cl['trial_count']} registered for this condition — see ClinicalTrials.gov"))
    if d and v.get("panels"):    # panel line only alongside a germline condition
        n = len(v["panels"])
        facts.append(("Diagnostic panels",
                      f"{v['gene_symbol']} is a diagnostic-grade gene on "
                      f"{n} Genomics England panel" + ("s" if n != 1 else "")))
    if facts:
        L += ["", table(["", ""], facts)]
    # top symptoms — of the variant's OWN condition (audit P1), only when present
    if d.get("phenotypes"):
        L += ["", f"**Commonly reported features of {d.get('name','this condition')}** "
              "(across patients with this condition; presentation varies — "
              "frequencies from Orphanet):", ""]
        L += [f"- {p['term']}" + (f" — {p['freq']}" if p.get("freq") else "")
              for p in d["phenotypes"][:8]]
    # registry + counselor CTA
    cta = []
    if cl.get("gard"):
        cta.append(f"[Condition overview & support (NIH GARD)](https://rarediseases.info.nih.gov/diseases/{cl['gard']}/index)")
    cta.append("[Find a genetic counselor (NSGC)](https://findageneticcounselor.nsgc.org/)")
    L += ["", "**Support & next steps:** " + " · ".join(cta) + "."]
    return L


def _gene_context_zone(v):
    gc = v.get("gene_context") or {}
    con, dos, val = gc.get("constraint"), gc.get("dosage"), gc.get("validity")
    if not (con or dos or val):
        return []
    L = ["", "## Gene constraint & dosage {#gene-context}", ""]
    bits = []
    if con and con.get("loeuf"):
        bits.append(f"LOEUF {con['loeuf']}, pLI {con.get('pli')}, missense-Z {con.get('mis_z')}")
    if dos and (dos.get("haplo") or dos.get("triplo")):
        bits.append(f"ClinGen dosage — haploinsufficiency {dos.get('haplo')}, triplosensitivity {dos.get('triplo')} (0–3 scale)")
    if bits:
        L.append(f"**{v.get('gene_symbol')}** population constraint: " + "; ".join(bits)
                 + ". *Lower LOEUF and higher pLI indicate a gene that tolerates loss-of-function "
                 "poorly; ClinGen dosage is the curated haploinsufficiency call.*")
    if val:
        L += ["", "Curated gene–disease validity (ClinGen):", ""]
        L += [f"- **{x['disease']}** — {x['classification']} ({x['moi']})" for x in val]
    return L


def _protein_zone(v):
    st = v.get("structural")
    struct = v.get("structure") or {}
    pdb, af = struct.get("pdb") or [], struct.get("alphafold")
    if not (st or pdb or af):
        return []
    L = ["", "## Protein context {#protein}", ""]
    if st and st.get("features"):
        L.append(f"Residue **{st['position']}** lies in " + "; ".join(st["features"]) + ".")
    if af and af.get("plddt"):
        L.append(f"\nAlphaFold model confidence (whole protein): pLDDT {af['plddt']}, "
                 f"{round(float(af['frac_high'])*100)}% of residues very-high."
                 if af.get("frac_high") else f"\nAlphaFold pLDDT {af['plddt']}.")
    if pdb:
        mut = [p for p in pdb if any(w in (p.get("title") or "").lower()
                                     for w in ("mutant", "variant"))]
        line = f"\n**{len(pdb)} experimental structure(s)** (PDB): " + ", ".join(
            f"[{p['id']}](https://www.rcsb.org/structure/{p['id']})" for p in pdb[:6])
        if mut:
            line += f" — incl. mutant structure {mut[0]['id']}"
        L.append(line + ".")
    return L


def _civic_zone(v):
    c = v.get("civic")
    if not c or not c.get("evidence"):
        return []
    L = ["", "## Cancer therapy associations (CIViC) {#civic}", "",
         "Curated clinical evidence for this variant from CIViC "
         f"(profile *{c.get('name')}*):", "",
         table(["Disease", "Therapies", "Type", "Level", "Significance"],
               [(e.get("disease"), e.get("therapies"), e.get("type"),
                 e.get("level"), e.get("significance")) for e in c["evidence"]]),
         "*CIViC evidence levels A–E (A strongest). Therapy associations are "
         "context-specific — not a treatment recommendation.*"]
    return L


def _num(s):
    """Round a raw assay score string for display (keep 3 decimals)."""
    try:
        return f"{float(s):.3f}"
    except (TypeError, ValueError):
        return s


def _mechanism_zone(v):
    """'What this gene does' — the variant→gene→pathway→disease story. Explicitly
    GENE-level (honest framing), Reactome disease-pathway callout as the headline."""
    pw = v.get("pathways") or {}
    if not (pw.get("pathways") or v.get("mechanism")):
        return []
    g = v.get("gene_symbol", "the gene")
    L = ["", "## What this gene does {#gene-function}", "",
         f"*Describes the biological roles of {g} (the gene this variant disrupts) — "
         "not effects measured for this specific variant. A damaging variant is expected "
         "to affect these functions and pathways; the degree depends on the variant.*", ""]
    if v.get("mechanism"):
        L += [f"**How this variant is thought to act:** {v['mechanism']}", ""]
    dp = pw.get("disease_pathways") or []
    if dp:
        L.append(f"**Disease mechanism (Reactome):** {g} loss or alteration is curated in "
                 + ", ".join(f"[{p['name']}](https://reactome.org/content/detail/{p['id']})"
                             for p in dp[:3]) + ".")
    pws = pw.get("pathways") or []
    if pws:
        L += ["", "### Pathways affected {#pathways}", "",
              f"{g} participates in **{len(pws)} Reactome pathway"
              + ("s" if len(pws) != 1 else "") + "** (⚕ = disease pathway):", "",
              table(["Pathway (Reactome)", "Evidence"],
                    [(("⚕ " if p["is_disease"] else "")
                      + f"[{p['name']}](https://reactome.org/content/detail/{p['id']})",
                      p["evidence"]) for p in pws[:10]])]
    mf, bp = pw.get("go_mf") or [], pw.get("go_bp") or []
    if mf or bp:
        tier = "experimentally supported" if pw.get("go_experimental") else "annotated"
        L += ["", f"### Molecular function & processes (GO — {tier}) {{#go}}", ""]
        if mf:
            L.append("**Molecular function:** " + ", ".join(g_["name"] for g_ in mf[:6]) + ".")
        if bp:
            L.append(("\n" if mf else "") + "**Biological process:** "
                     + ", ".join(g_["name"] for g_ in bp[:6]) + ".")
    return L


def _mavedb_zone(v):
    m = v.get("mavedb")
    if not m:
        return []
    return ["", "## Functional evidence (MaveDB) {#mavedb}", "",
            f"Measured in **{len(m)} multiplexed functional assay"
            + ("s" if len(m) != 1 else "") + "** (deep mutational scanning):", "",
            table(["Assay (MaveDB score set)", "Raw score", "License"],
                  [(f"[{r.get('title') or r['score_set']}](https://www.mavedb.org/score-sets/{r['score_set']}/)",
                    _num(r["score"]), r.get("license")) for r in m]),
            "*Raw per-assay scores — sign and scale differ between assays; interpret "
            "within each score set. Experimental functional evidence (an ACMG PS3/BS3 "
            "input), not a classification.*"]


def _pharmgkb_zone(v):
    p = v.get("pharmgkb")
    if not p or not (p.get("clinical") or p.get("drugs")):
        return []
    L = ["", "## Pharmacogenomics (PharmGKB) {#pgx}", ""]
    if p.get("clinical"):
        L += ["Clinical drug-response annotations:", "",
              table(["Drug(s)", "Evidence level", "Category", "Phenotype"],
                    [(c.get("chemicals"), c.get("level"), c.get("type"), c.get("phenotypes"))
                     for c in p["clinical"]])]
    elif p.get("drugs"):
        L.append("Drugs with reported response associations: " + ", ".join(p["drugs"]) + ".")
    L.append("*PharmGKB levels 1A/1B (highest) → 4. Discuss any medication "
             "decision with your prescriber — reference information, not advice.*")
    return L


def render_body(v, jsonld_tag=""):
    L = ["## Summary", "", declarative(v), ""]
    # At a glance
    L.append("**At a glance:** "
             + " · ".join(filter(None, [
                 v.get("classification"),
                 f"review: {_review(v.get('review_status'))}",
                 (f"rsID {v['rsid']}" if v.get("rsid") else None),
                 (f"{v['submitter_count']} submitter"
                  + ("s" if v.get("submitter_count") != 1 else "")
                  if v.get("submitter_count") else None),
                 (f"{len(v['conditions'])} linked condition"
                  + ("s" if len(v['conditions']) != 1 else "")
                  if v.get("conditions") else None),
                 (f"AlphaMissense {v['alphamissense']['class'].replace('_', ' ')}"
                  if v.get("alphamissense") and v["alphamissense"].get("class") else None),
                 (v["gnomad"]["band"] if v.get("gnomad") else None),
             ])))

    if jsonld_tag:
        L += ["", jsonld_tag]

    # Patient zone — high on the page (high human value + citable)
    L += _patient_zone(v)

    # Identity
    L += ["", "## Identity {#identity}", "",
          table(["Field", "Value"], [
              ("Gene", links.maybe_link(v.get("gene_symbol"),
                                        links.gene_url(symbol=v.get("gene_symbol"), hgnc_id=v.get("hgnc_id")))),
              ("Protein change (HGVS p.)", v.get("hgvs_p")),
              ("Coding change (HGVS c.)", v.get("hgvs_c")),
              ("dbSNP", (f"[{v['rsid']}](https://www.ncbi.nlm.nih.gov/snp/{v['rsid']}/)"
                         if v.get("rsid") else None)),
              ("Variant type", v.get("variant_type")),
              ("Location", (f"chr{v['chromosome']}:{v['start']}-{v['stop']} ({v['assembly']})"
                            if v.get("chromosome") else None)),
              ("ClinVar", f"[VCV{v['variation_id']}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{v['variation_id']}/)"),
          ])]
    exprs = v.get("hgvs_expressions") or []
    if exprs:
        L.append("\n**All HGVS expressions:** " + ", ".join(f"`{e}`" for e in exprs))

    # Computational & population evidence — the cross-source concordance readout.
    # Renders when there's a concordance verdict OR a molecular-consequence read
    # (non-missense LoF variants often have no coordinate → no verdict, but the
    # consequence/LoF-context is still the key computational signal for them).
    conc = v.get("concordance") or {}
    if conc.get("verdict") or v.get("lof_context") or v.get("consequence"):
        L += ["", "## Computational & population evidence {#evidence}", ""]
        if conc.get("verdict"):
            L += [f"**Concordance:** {conc['verdict']}.", ""]
            L += [f"- {ln}" for ln in conc.get("lines", [])]
        mm = v.get("am_isoform_mismatch")
        if mm:
            L.append(f"- ⚠ **Isoform check:** AlphaMissense is numbered on a different transcript "
                     f"(AM **{mm['am']}** vs ClinVar **{mm['clinvar']}**) — the AlphaMissense read here "
                     "may be for a different residue; verify against the ClinVar transcript.")
        lof = v.get("lof_context")
        if lof:
            con = ", ".join(filter(None, [
                f"LOEUF {lof['loeuf']}" if lof.get("loeuf") is not None else None,
                f"pLI {lof['pli']}" if lof.get("pli") is not None else None]))
            if lof.get("lof_disease_gene"):     # germline + (haploinsufficient or constrained)
                if lof.get("haploinsufficient"):
                    # cite constraint numbers ONLY when they also support (else they contradict
                    # the haploinsufficiency claim — e.g. ASXL1 pLI≈0)
                    basis = ("ClinGen curates **sufficient evidence for haploinsufficiency** in this gene"
                             + (f" (constraint {con})" if con and lof.get("constrained") else ""))
                else:
                    basis = f"the gene is **loss-of-function-intolerant** ({con})"
                tail = f" {basis}, supporting a loss-of-function disease mechanism."
            elif con:
                tail = f" Gene population constraint: {con}."
            else:
                tail = ""
            if lof.get("nmd_caveat"):
                tail += (" For a truncating variant the loss-of-function impact is position-dependent — "
                         "C-terminal / last-exon truncations may escape nonsense-mediated decay.")
            L.append(f"- Molecular consequence: **{lof['label']}** — a predicted loss-of-function "
                     f"variant.{tail} *Descriptive; not an applied PVS1 code.*")
        elif v.get("consequence"):
            L.append(f"- Molecular consequence: **{v['consequence']['label']}**.")
        pctl = v.get("am_percentile")
        if pctl and pctl["top_pct"] <= 10:      # only when it's a genuine standout
            L.append(f"- AlphaMissense ranks this among the **top {pctl['top_pct']}%** "
                     f"most-pathogenic-predicted substitutions in {v.get('gene_symbol')} "
                     f"(of {pctl['n']:,} modeled).")
        L.append("\n*Computational predictors are not independent (ClinGen SVI): "
                 "AlphaMissense (Cheng et al. 2023) carries the ACMG weight for missense "
                 "and conservation (phyloP/GERP/phastCons) for non-missense; REVEL "
                 "(Pejaver-2022 calibrated strength bands) and SaProt (Su et al. 2023, structure-"
                 "aware, ClinVar-independent protein language model, raw LLR) are shown as orthogonal "
                 "agreement signals, not additive. gnomAD v4.1 grpmax with the ClinGen-recommended "
                 "filtering allele frequency (faf) is the BA1/BS1/PM2 metric; thresholds are "
                 "disease-specific. Predictions, not a clinical determination.*")

    # Gene ACMG context + protein/structural context + mechanism/pathways
    L += _gene_context_zone(v)
    L += _protein_zone(v)
    L += _mechanism_zone(v)

    # Clinical significance + consensus + per-submitter table
    cons = v.get("consensus")
    L += ["", "## Clinical significance {#significance}", "",
          f"**{v.get('classification')}** — review status {_review(v.get('review_status'))} "
          f"*({v.get('review_status')})*"
          + (f", last evaluated {v['last_evaluated']}." if v.get("last_evaluated") else ".")]
    if cons:
        L.append(f"\n**Submitter consensus:** {cons['verdict']}.")
    tl = v.get("timeline")
    if tl and tl["n"] > 1:
        span = f"{tl['first']}" if tl["first"] == tl["last"] else f"{tl['first']}–{tl['last']}"
        L.append(f"\n**Submission history:** {tl['n']} submissions ({span}); "
                 + ("classifications have been stable." if tl["stable"]
                    else "classifications have differed over time."))
    vcep = v.get("vcep") or []
    if vcep:
        c = vcep[0]
        L.append(f"\n**ClinGen expert panel:** {c.get('assertion')} — {c.get('panel')} "
                 f"({c.get('disease')}). *Expert-panel curated (ACMG); the highest "
                 "ClinVar review tier.*")
        cr = c.get("criteria")
        if cr:
            if cr.get("codes_met"):
                L.append("**Applied ACMG criteria (met):** "
                         + ", ".join(f"`{x}`" for x in cr["codes_met"]) + ".")
            if cr.get("codes_not_met"):
                L.append("*Considered but not met:* "
                         + ", ".join(f"`{x}`" for x in cr["codes_not_met"]) + ".")
            if cr.get("moi"):
                L.append(f"Inheritance: {cr['moi']}.")
            if cr.get("summary"):
                L.append(f"\n> **Curation rationale (ClinGen, verbatim):** {cr['summary']}")
            prov = [x for x in (cr.get("guideline"),
                                f"approved {cr['approval_date']}" if cr.get("approval_date") else None,
                                f"[ClinGen Evidence Repository]({cr['erepo']})" if cr.get("erepo") else None) if x]
            if prov:
                L.append("*Source: " + " · ".join(prov) + ".*")
            L.append("*These are ClinGen's applied codes, reported verbatim — Sugi Variant does "
                     "not re-tally or combine them into its own classification.*")
    subs = v.get("submissions") or []
    if subs:
        L += ["", "### Submitter classifications {#submitters}", "",
              table(["Submitter", "Classification", "Review", "Method", "Date"],
                    [(s.get("submitter"), s.get("classification"),
                      s.get("review_status"), s.get("method"), s.get("date")) for s in subs])]

    # Conditions
    conds = v.get("conditions") or []
    if conds:
        L += ["", "## Associated conditions {#conditions}", "",
              ", ".join(links.maybe_link(c.get("name") or c.get("mondo_id"),
                                         links.disease_url(mondo_id=c.get("mondo_id"), name=c.get("name")))
                        for c in conds) + "."]
    extra = [p for p in (v.get("phenotype_list") or []) if p]
    if extra:
        L.append("\n**Reported phenotype names (ClinVar):** " + ", ".join(extra[:10]) + ".")

    # Same-residue hotspot context (deterministic, from the gene's own P/LP set)
    hs = v.get("hotspot")
    if hs and hs.get("others"):
        others = hs["others"]
        L += ["", "## Same-residue context {#hotspot}", "",
              f"Residue **{hs['position']}** carries **{len(others)} other pathogenic "
              f"ClinVar variant" + ("s" if len(others) != 1 else "")
              + "** — a recurrently-mutated position: "
              + ", ".join(f"[{o['label']}](/variant/{o['slug']})" for o in others[:12])
              + ". *Positional co-occurrence of independent ClinVar records, not "
              "functional proof.*"]

    # Functional (MaveDB) + cancer-therapy (CIViC) + PGx (PharmGKB) — conditional
    L += _mavedb_zone(v)
    L += _civic_zone(v)
    L += _pharmgkb_zone(v)

    # Similar variants (internal mesh)
    sim = v.get("similar") or []
    if sim:
        L += ["", "## Similar variants {#similar}", "",
              "Related " + v.get("gene_symbol", "") + " variant pages (same residue, "
              "condition, or type): "
              + ", ".join(f"[{s['label']}](/variant/{s['slug']})" for s in sim) + "."]

    L += ["", "### Data sources & attribution {#sources}", "", attribution_md(), "",
          f"*Primary record: NCBI ClinVar (variation {v['variation_id']}). Classifications "
          "reflect these databases as of the page's build date and may change. Free "
          "reference use — not medical advice; consult the primary submitters and a clinician.*"]
    return "\n".join(L)
