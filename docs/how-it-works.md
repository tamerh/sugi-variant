# How Sugi Variant works

Sugi Variant assembles a **deterministic reference page** for every clinically significant
human variant by cross-joining a dozen authoritative sources, resolving any variant to a
stable URL, and reporting the evidence under a strict set of interpretation rules. This page
explains the data, the resolver, and those rules. To use it, see
[getting-started.md](getting-started.md); for the family, see [README.md](README.md).

Every identifier is grounded through **BioBTree**, so gene, transcript, protein and disease
IDs all reconcile to the same canonical entities.

## Corpus

The small-variant corpus is the ClinVar germline slice: every buildable classification —
**Pathogenic**, **Likely pathogenic**, **Conflicting classifications of pathogenicity**,
**Uncertain significance** and **Benign/Likely benign** — across protein-coding, **ncRNA**
(`n.`) and **mitochondrial** (`m.`) loci (currently ~4.08M indexed variants across ~18,600
genes).

A separate **structural-variant corpus** comes from dbVar: clinical deletions,
duplications, insertions, inversions and CNVs, gated to the same canonical germline classes
(currently 82,672 SVs). dbVar's somatic/oncology/drug-response labels are out of scope for
now and are not mixed into the germline pages.

## Sources

Each fact on a page traces to one source, carried with its dataset build date:

| Domain | Source |
|---|---|
| Clinical classification, review status, submitters | **ClinVar** |
| Structural variants (type, size, breakpoints, origin) | **dbVar** (NCBI) |
| Expert-panel applied ACMG codes | **ClinGen VCEP** (reported verbatim, not re-tallied) |
| Gene constraint, dosage, gene–disease validity | **gnomAD** · **ClinGen** · **GenCC** |
| Population frequency | **gnomAD v4.1** |
| Missense pathogenicity | **AlphaMissense** (ACMG-weighted) · **REVEL** · **SaProt** (agreement signals) |
| Conservation | phyloP · phastCons · GERP |
| Splicing | **SpliceAI** |
| Protein structure/annotation | **UniProt** · **AlphaFold** · PDB |
| Mechanism / pathways | **Reactome** · **GO** |
| Functional assays | **MaveDB** |
| Somatic/therapeutic | **CIViC** |
| Disease context | **Orphanet** (inheritance, onset, prevalence) |

## Resolution

Variants are addressable three ways, all resolved through a precomputed **SQLite index**
(one row per variant + an alias table), so any lookup is O(1) with no per-request build:

- **Gene + HGVS** — protein or coding, in any shape (`PTEN R130Q`, `PTEN p.Arg130Gln`,
  `PTEN c.389G>A`)
- **dbSNP rsID** — `rs121909229`
- **GRCh38 coordinate** — `chr:pos:ref:alt`

Because HGVS isn't free-text-searchable upstream, each page's URL is a **deterministic slug**
derived from `GENE` + normalized HGVS (e.g. `pten-p-arg130gln`); the coding form is a
same-page alias. Over-long insertions are hashed so URLs stay bounded.

## Structural variants

Structural variants are a separate corpus, not more rows in the small-variant records:

- **Identity** is the dbVar `NSSV…` accession; the parent region is `nsv…`.
- **Location** is an interval. The corpus has its own SQLite store with an R-tree index, so
  a point or locus query returns every SV whose span covers it.
- **Genes are by overlap.** dbVar carries no gene assignment; the build intersects the
  interval with Ensembl gene coordinates. A large CNV can overlap many genes.
- **Page.** Type, size, breakpoints, clinical significance, origin, phenotypes, ClinVar
  accessions, overlapped genes and the dbVar description.
- **Not the SNV evidence.** AlphaMissense, REVEL, SaProt and SpliceAI are single-residue
  predictors; the SV page states that they do not apply rather than showing an empty score.
- **Same gate.** Only the canonical germline ACMG classes are built; dbVar's
  somatic/oncology/drug-response labels are a separate future corpus.

A gene page lists the SVs overlapping it, capped at the 100 most clinically significant.

## Interpretation rules

Sugi Variant **describes** authoritative evidence — it does not re-derive it. The rules that
keep it honest:

- **No self-computed classification, score, verdict, or PVS1.** ClinVar's classification and
  ClinGen's applied codes are shown as-is.
- **Predictors are not additive.** Exactly one calibrated tool carries the ACMG weight —
  **AlphaMissense** for missense, **conservation** for non-missense. REVEL and SaProt are
  shown as *agreement signals*, never independent votes. The in-silico consensus is
  descriptive.
- **Discordance is a QC flag, not a reclassification.** When the predictors and the clinical
  call disagree, the page flags it for review; it never claims to resolve the conflict.
- **Condition names are verbatim** from ClinVar/MONDO.
- **ClinVar's 0–4 gold-star review tier** is ClinVar's own signal, surfaced as a named badge
  plus stars.
- **Deterministic** — no model-written page text. Same inputs → same page, every build.

## Layout

Each variant page has three interchangeable layouts — **Dashboard** (default), **Datasheet**
and **Card** — over the same facts. Gene pages (`/gene/<SYMBOL>`) list a gene's variants by
classification and carry its evidence-disagreement section. A **set** page (`/set`) compares
several variants at once with a review worklist.
