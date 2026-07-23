# Sugi Variant — a per-variant genetic reference

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

## Run

    export ATLAS_BIOBTREE=http://<biobtree-host>:9291   # a running BioBTree REST API
    python -m uvicorn app:app --host 127.0.0.1 --port 8000

Then open <http://127.0.0.1:8000/>. Optional environment variables:

- `BASE_PATH` — URL prefix when served under a subpath behind a reverse proxy (e.g. `/variant`).
- `CACHE_DIR` — page/records cache and the prebuilt resolution `index.db` (default `./cache`).
- `GENE_CACHE_MB` — in-process LRU budget for the hottest genes (default `1500`).

## Dependencies

    pip install fastapi "uvicorn[standard]" httpx jinja2

Plus the shared `sugibiobtree` BioBTree client (one source of truth with Sugi Atlas), installed
from the `sugi-biobtree` repository.

## Layout

- `app.py` — the FastAPI app (routes, caching, identifier resolution).
- `sugivariant/` — the science package (collect / enrich / render / index / slug).
- `templates/`, `static/` — Jinja views and assets.
- `static/preprint/` — the built preprint, served at `/method`.
