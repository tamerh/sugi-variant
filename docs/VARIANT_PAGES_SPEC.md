# Design Spec — `/atlas/variant/` entity type

Status: DESIGN (2026-07-07). Flagship item from the 2026-07-07 sweep (see
IMPROVEMENT_ROADMAP.md). Grounded in live biobtree probes (ACTA1, BRCA1, TP53,
CFTR, KCNQ1) + the GSC demand signal (memory `atlas-traffic-is-mostly-bots`).

## 1. Why (the demand)
~97% of Search-Console query volume is machine-shaped `GENE "HGVS"` lookups
(`acta1 "pro309ala"`, `"c.925c>g" acta1`), clustered on ACTA1/DACT1/ASXL1/NAA10/
RPL10. Atlas surfaces GENE pages for these and gets ~0 clicks because there is no
per-variant page. A variant page is the *exact object* that demand + the
AI-citation channel asks for, and — unlike the shelved pathway type — it is a
UNIQUE record (per-submitter calls, VCEP assertion, phenotype set, predictions)
that exists on no gene/disease page at that granularity.

## 2. Data model (verified reachable, deterministic)
Anchor dataset = **ClinVar**. Enumerated GENE-FIRST (HGVS is NOT independently
searchable in biobtree — `search("p.Pro334Arg")` → 0), exactly like gene §6:
- Enumerate: `map_all(hgnc_id, ">>hgnc>>clinvar")` → per-gene variation IDs
  (paginated; ACTA1 624, BRCA1 10k+, TP53 3.9k, CFTR 6k, KCNQ1 3k).
- Full record: `entry(variation_id, "clinvar")` →
  - `hgvs_expressions` (NM_/NC_/NG_/LRG_ c./g.), p. form in `name`
    (`c.1001C>G (p.Pro334Arg)`), dbSNP `rs…`, chrom/start/stop, GRCh38, var type.
  - `germline_classification` (Pathogenic/…), `review_status` (stars),
    `last_evaluated`.
  - `submissions[]` — per-submitter {submitter, classification, review_status,
    method, date}. ← the differentiator no aggregator shows cleanly.
  - `phenotype_list` + `phenotype_ids` incl. **MONDO** → mesh to disease pages.
  - `xrefs`: entrez, orphanet, medgen, mim, ensembl, mondo, dbsnp,
    **clingen_variant**, pubmed.
- Cross-links (all live): `>>clinvar>>mondo` (conditions), `>>clinvar>>
  clingen_variant` (VCEP ACMG assertion — authority tier above raw ClinVar),
  gene from `gene_symbol`/`hgnc_id`.
- Positional annotations (join by protein-substitution / genomic position, NOT
  page generators): **AlphaMissense** (`>>transcript>>alphamissense`), **SpliceAI**
  (`>>hgnc>>spliceai`) — these are grids; annotate the variant, never enumerate it.

## 3. Scope / scale (MUST gate — full ClinVar ≈ 2–3M is not buildable)
P/LP is ~25–40% of a gene's ClinVar set, concentrated in well-studied genes.
Phased gate on `germline_classification`:
- **Phase 1 (ship first): Pathogenic + Likely-pathogenic + Conflicting** →
  est. ~30–80k pages (6–12× is NOT triggered; ~1.5× current corpus). Highest
  value + most searched.
- **Phase 2: + expert-panel-reviewed VUS** (review_status ≥ "reviewed by expert
  panel") → adds the high-authority VUS without the benign long tail.
- **Never build**: Benign/Likely-benign, and single-submitter no-assertion VUS
  (the mass of the long tail; low value, huge cost).
Only for genes already in the corpus (29k gene pages). Build cost is the real
constraint (per-variant `entry()` call) — use the sharded/incremental build
(`batch.py` supports it); consider caching ClinVar entries.

## 4. Page contract (H2 zones — frozen, like other types)
1. **Summary** — declarative lead: "ACTA1 p.Pro334Arg (c.1001C>G) is a
   Pathogenic variant in ACTA1, reviewed by expert panel, associated with
   {condition}." + at-a-glance (classification, stars, rsID, condition count).
2. **Identity** — all HGVS expressions (c./p./g.), dbSNP rsID, genomic coords
   (both assemblies), variant type, gene link.
3. **Clinical significance** — current classification + review-status stars +
   last-evaluated; the **per-submitter table**; **ClinGen VCEP** ACMG assertion
   when present (highest tier).
4. **Conditions** — linked disease pages (via MONDO) + phenotype list.
5. **Computational predictions** — AlphaMissense (for the substitution) + SpliceAI
   (if in a splice region), clearly labelled as in-silico.
6. **Evidence / citations** — PubMed; outbound ClinVar/dbSNP/gnomAD links.
7. **Related** — sibling variants in the same gene/condition (mesh).

## 5. URL / slug design (the make-or-break — HGVS not searchable)
Variants are enumerated (not looked up by HGVS), so the slug MUST be
deterministic from `GENE` + normalized HGVS, and cover BOTH query shapes:
- Canonical: `/atlas/variant/acta1-p-pro334arg/` (protein form).
- Alias: `/atlas/variant/acta1-c-1001c-g/` (coding form) → same page (Hugo
  alias / redirect, or dual-emit).
Normalization: lowercase; `p.`/`c.`→`p-`/`c-`; strip `>`→`-`; `>`,`(`,`)`,` `→`-`;
collapse repeats. Deterministic + reversible-enough for an agent to guess from
`GENE "HGVS"`. Store both slugs on the record; the per-gene index lists both.
Handle collisions (multiple HGVS → same variation id) by preferring the MANE/
canonical-transcript c. form for the alias.

## 6. Mesh / cross-linking
- **variant → gene** (always), **variant → disease** (MONDO), **variant →
  drug** (only via CIViC predictive on cancer genes — sparse, enrich-only).
- **gene → variants**: a per-gene variant INDEX hub `/atlas/gene/<SYM>/variants/`
  (crawlable enumeration — agents can't guess 80k URLs, they follow index pages).
  Also promote gene §6 to link its top pathogenic variants to their pages.
- **disease → variants**: group `>>mondo` variants by condition on disease pages.
- Reverse-edge index (`reverse_edges.json`) already handles the mechanics.

## 7. Build architecture
Architecture tax is largely PAID DOWN — the pathway build generalized the 4th
entity type. To add `variant`:
- `src/atlas/page/links.py`: add `"variant"` to `_ENTITY_TYPES` + `_MANIFEST`/
  `_CANON`/`REVERSE_LABEL`/`_REVERSE_ORDER`.
- `src/atlas/variant/` vertical mirroring `src/atlas/pathway/`: `anchors.py`
  (resolve variation id → record), `collect.py`, `render.py`, `slug.py`.
- The gene §6 collector (`src/atlas/gene/sections/s06_variants.py`) is essentially
  the variant collector already — factor the per-variant logic into the shared
  path.
- `batch.py` / `build_corpus.py`: variant enumeration + sharded build + the
  per-gene index page generation.
- Frontend (`../sugi-atlas-web`): a `layouts/atlas/` variant single template
  (reuse baseof.html; add the answer-shaped `<title>` case:
  `"{GENE} {HGVS}: pathogenicity, conditions & evidence"`), a variant sitemap,
  and the `/atlas/gene/<SYM>/variants/` index layout.

## 8. AI-citation payload
Emit `facts.json` per variant (roadmap C5): {gene, hgvs_p, hgvs_c, rsid,
classification, review_stars, conditions[mondo], vcep_assertion, predictions}.
This is the maximally-citable object for the agent/pipeline traffic that IS the
demand. Pairs with the QAPage/JSON-LD Tier-1 work.

## 9. Rollout (phased)
1. Prototype the collector + slug on ONE gene (ACTA1, 233 P/LP) → validate record
   completeness, slug determinism, both-alias resolution, mesh to gene/disease.
2. Build the ACTA1/DACT1/ASXL1/NAA10/RPL10 demand-cluster genes → measure GSC
   impact on the `GENE "HGVS"` queries (does CTR move?).
3. Phase-1 gate across the corpus (~30–80k) if the pilot converts.
4. Add per-gene index hubs + facts.json + Hugo titles.

## 10. Risks / open questions
- **Build cost**: per-variant `entry()` at 30–80k pages — needs caching /
  sharding; measure a pilot gene's wall-clock first.
- **Slug collisions & HGVS variety**: multiple transcripts → multiple c. forms;
  pick MANE canonical for the alias, list the rest on the page.
- **Contract test**: the integration harness needs `variant` H2/H3 ids added.
- **Sitemap size**: 80k more URLs — split the variant sitemap.
- **Does CTR actually move?** The demand is agent/pipeline (they may fetch the
  `.md`/facts.json, not "click"). Measure via server logs + the AI-citation
  channel, not just SERP CTR (per the strategic reality).
- **Non-ClinVar variants**: GWAS/CIViC variants are out of scope for v1 (ClinVar
  is the anchor; CIViC enriches cancer-gene variant pages only).

---
## Enrichment plan (from the 2026-07-07 3-agent deep-mine)

BATCH 1 — SHIPPED (a0cf813): AlphaMissense (protein_variant join), gnomAD
frequency + band (rsID), cross-source concordance verdict (+ honest discordance
flags), submitter-consensus stat, same-residue hotspot, answer-first cited lead
+ provenance/disclaimer.

BATCH 2 — STAGED (verified reachable, join keys confirmed):
- UniProt residue features `>>hgnc>>uniprot>>ufeature` (cache per gene, join by
  protein position) — "residue is in X domain / active site / PTM / a UniProt
  disease variant". High interpretive value. Effort M.
- CIViC per-variant `clinvar>>civic_variant>>civic_evidence` (clean join via
  clinvar_ids xref) — cancer therapy/prognosis. Sparse (cancer genes). Effort S.
- PharmGKB variant annotations `dbsnp>>pharmgkb_var_annotation` (by rsID, piggy-
  backs the gnomAD dbsnp call) — drug-response + PMID + sentence. Sparse. Effort S.
- SpliceAI `>>hgnc>>spliceai` (join by chr:pos:ref:alt) — splice-region only. Effort M.
- QAPage/FAQPage JSON-LD keyed to the literal "Is GENE pChange pathogenic?" +
  the answer-first lead verbatim — the top AI-citation lever. Effort S.
- Evidence-strength scorecard (ACMG-flavored, NOT a clinical call) — LAST, highest
  framing risk; unit-test input→expected before regen. Effort S.
- PubMed link list `clinvar>>pubmed` (link-only; titles not stored). Effort S.

DEAD ENDS (checked): CADD/REVEL/SIFT/PolyPhen/dbNSFP/conservation not in corpus;
MANE/specific-transcript only reaches gene ENSG; PubMed titles/abstracts not
stored; pharmgkb_clinical is star-allele-keyed (not variant-joinable).

DIFFERENTIATION (competitive scan): the synthesized variant pages (VarSome,
Franklin) are JS/login-walled → invisible to AI crawlers; the crawlable ones
(ClinVar, gnomAD) are single-source + un-narrated. Nobody occupies
"crawlable + multi-source + answer-first + narrated" — that's the seam. Defer
visibly to ClinGen VCEP / high-star ClinVar when present (being the source that
correctly points to the higher authority is itself citation-worthy).

---
## Batch 3+ enrichment (2026-07-07 persona deep-mine, 4 agents)

BIG INSIGHT: most high-value data is ALREADY in biobtree, unwired. Two personas
(researcher / patient) → restructure page into a technical zone + a
"For patients & families" zone. ROUTING TRAP (verified): the ClinVar→MONDO edge
lands on a narrow, poorly-annotated subtype (0 HPO/prevalence/trials); route
patient data via the GENE's best Orphanet disease (`>>hgnc>>orphanet`, reuse
disease/anchors.py:296-325 selection) and the PARENT MONDO (`>>mondo>>mondoparent`)
for trials.

BUILDABLE NOW — PATIENT zone (all reachable, need strong per-section disclaimers):
- Condition digest: inheritance (`>>hgnc>>gencc` moi_title / `>>hgnc>>clingen_gene_validity` moi+Definitive / Orphanet inheritance), age of onset + prevalence + HPO-with-frequencies (Orphanet entry via the gene's best disease). One entry() call. TOP.
- "Am I alone?" scale line — FREE, from the already-enumerated gene variant set (233 other P/LP in ACTA1) + condition gene_count.
- "Talk to a genetic counselor" CTA + GARD registry link (`>>mondo>>gard`) — trivial, highest safety value.
- Clinical trials for the condition via PARENT MONDO (`>>mondo>>clinical_trials`, reuse s13).
- PanelApp diagnostic-panel line (`>>hgnc>>panelapp_gene`; ACTA1 green on Congenital myopathy) — actionability proxy.
- Plain-language one-line summary (deterministic template, NOT LLM).

BUILDABLE NOW — RESEARCHER / structural:
- Residue structural context — ufeature INTERVAL overlap over the p.-position ("residue 309 in an α-helix / α-actinin binding region / PTM site"). ufeature staged; go beyond the flat list. TOP structural pick.
- Paralog residue-equivalence — aligned residue in a near-identical paralog (ACTA1↔ACTC1 ~98.7%) carries a known P/LP variant. Reuses build_position_index; gate on high identity. Atlas-internal cross-link.
- PDB/AlphaFold line — structures exist (+ mutant title-match), global pLDDT confidence. NOTE: per-residue pLDDT NOT reachable (external-only).

BUILDABLE NOW — gene-context block (one call/gene each, cache in ctx):
- gnomAD constraint (pLI/LOEUF/mis_z) + ClinGen dosage HI/TS (`>>hgnc>>clingen_dosage`) + ClinGen gene-validity — ACMG-relevant "this gene doesn't tolerate damage".
- GWAS Catalog via rsID (`rs..>>dbsnp>>gwas`) — the ONLY signal for non-coding/common variants.

BUILDABLE NOW — derived analyses / features (mostly zero new calls, aggregate in-memory recs):
- Gene variant-landscape panel (textual lollipop + type breakdown + top recurrent residues) on the index + compact on pages.
- "Similar variants" internal mesh (same residue/type/condition/domain → top ~6 links) — thickens crawl mesh.
- AlphaMissense gene-percentile framing ("top 3% of 2,499 modeled") + residue tolerance — RAISE the am cache cap from 200 to full pagination.
- Variant→gene-drug actionability bridge (gene-level, honestly framed; reuse gene page drug data via manifest).
- Evidence/submission timeline (from submission dates). Per-gene index stats header.

BIOBTREE-TEAM ASKS (ingestion): (1) dbNSFP LICENSE-CLEAN subset ONLY — phyloP,
GERP, SIFT, MetaRNN, PrimateAI, ESM1b (NOT REVEL/CADD/PolyPhen — academic-only /
non-redistributable) → adds conservation (biggest gap) + orthogonal predictors to
harden concordance beyond AM-alone. (2) gnomAD v4 per-variant detail (per-ancestry
AF, AC/AN, popmax, flags, regional constraint — dbsnp has only single global freq).
(3) MaveDB (now CC0) functional-assay scores, variant-keyed (PS3/BS3-grade; sparse).

EXTERNAL-ONLY (outbound link, not computed): MaveDB/ProtVar (until ingested),
per-residue pLDDT/conservation, ACMG SF actionable-gene list (bundle a static
table), GTEx eQTL. DEAD ENDS: DMS/MaveDB, REVEL/CADD/SIFT/dbNSFP, per-residue
conservation, corum/string_interaction, clingen_dosage=0 for ACTA1, InterPro has
no coordinates.
