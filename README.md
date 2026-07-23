# Sugi Variant — a per-variant genetic reference

**Live: <https://sugi.bio/variant>**

Sugi Variant assembles one reference view per germline human variant by deterministically
mining the BioBTree knowledge graph. Anchored on ClinVar, each view brings together the
clinical classification and condition, gnomAD population frequency and gene constraint,
calibrated in-silico predictors (AlphaMissense, REVEL, SaProt, conservation, and SpliceAI for
splice regions), and — where a ClinGen expert panel has curated the variant — its applied ACMG
criteria, verbatim. Every fact is shown with its source and dataset build date, across all
variant classes (including non-coding and mitochondrial), and it surfaces the variants where the
independent signals disagree, as a quality-control signal for review.

It is a server-rendered FastAPI + Jinja app over the self-contained `sugivariant` package, with a
disk/ETag page cache. Nothing on a view is written by a language model: each is composed
deterministically from primary databases, so it is reproducible and every number traces to its
source and version. The method and evaluation are described in the preprint, served at `/method`.

## Documentation

- [Getting started](docs/getting-started.md) — using the site: search, URLs, sets, gene pages.
- [How it works](docs/how-it-works.md) — the method behind the pages.
- [Development](docs/development.md) — running locally, dependencies, and the code layout.
