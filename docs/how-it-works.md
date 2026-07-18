# How Sugi Variant works

Sugi Variant assembles a **deterministic reference page** for every clinically significant
human variant by cross-joining a dozen authoritative sources, resolving any variant to a
stable URL, and reporting the evidence under a strict set of interpretation rules. This page
explains the data, the resolver, and those rules. To use it, see
[getting-started.md](getting-started.md); for the family, see [README.md](README.md).

Every identifier is grounded through **BioBTree**, so gene, transcript, protein and disease
IDs all reconcile to the same canonical entities.

## Corpus

The corpus is the **pathogenic-gated** slice of ClinVar: every variant classified
**Pathogenic**, **Likely pathogenic**, or **Conflicting classifications of pathogenicity** —
**463,882** variants across **6,537** genes, including protein-coding, **ncRNA** (`n.`) and
**mitochondrial** (`m.`) loci. (Benign and VUS variants are out of scope by design.)

## Sources

Each fact on a page traces to one source, carried with its dataset build date:

| Domain | Source |
|---|---|
| Clinical classification, review status, submitters | **ClinVar** |
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
