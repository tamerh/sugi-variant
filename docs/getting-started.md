# Getting started

Sugi Variant is a website — there's nothing to install. This page covers how to find a
variant, the URL scheme, and the set/gene views. For the method behind the pages, see
[how-it-works.md](how-it-works.md).

Base URL: **https://sugi.bio/variant**

## Find a variant

Search from the home page with any of:

- **Gene + change** — `PTEN R130Q`, `PTEN p.Arg130Gln`, `BRCA2 c.1000C>T`
- **dbSNP rsID** — `rs121909229`
- **GRCh38 coordinate** — `17:7673803:C:T`

Press **Enter** to open a single variant. Autocomplete suggests as you type — pick a
**variant** to *add* it to the box (build a set), or a **gene** to open its page.

## URLs

Every page has a stable, shareable URL:

| Page | URL |
|---|---|
| Variant | `/pten-p-arg130gln` (slug = gene + normalized HGVS) |
| Variant by rsID | `/rs121909229` (redirects to the canonical slug) |
| Gene | `/gene/PTEN` |
| Evidence disagreements (per gene) | `/gene/PTEN#disagreements` |
| Set | `/set?v=slug1,slug2,…` |

Layout is a query param: append `?view=datasheet` or `?view=card` to any variant page
(default is the dashboard).

## Reading a variant page

- **Top tiles** — classification, ClinVar review tier (★), AlphaMissense call, gnomAD band,
  in-silico consensus, conditions — each labelled with its source.
- **Computational evidence** — AlphaMissense (the ACMG-weighted call) with REVEL, SaProt,
  conservation and SpliceAI as agreement signals.
- **ClinGen expert panel** — applied ACMG codes, verbatim, when a VCEP has curated the gene.
- **What this gene does** — a gene-level mechanism note (never a per-variant claim).
- **Disease context** — inheritance, onset, prevalence and linked conditions.

## Sets

Paste or build a list of variants to compare them on one page:

- **Build** at `/set` — paste `GENE change`, rsIDs or coordinates (one per line or
  comma-separated). Or type a comma-separated list into the main search box.
- **Get** a summary (class distribution, flagged-for-review count), a **review worklist** of
  the discordant ones, and a sortable/filterable table. The URL *is* the set — nothing is
  uploaded or stored.

## Gene pages

`/gene/<SYMBOL>` lists a gene's variants grouped by classification, and — near the top —
its **evidence-disagreement** section: the variants where the computational evidence doesn't
line up with the clinical call, a QC signal for a second look.

## Scope reminder

Reference and research use only — **not a diagnostic device and not medical advice**.
Interpretation for an individual must be done by a qualified clinician or genetic counselor.
