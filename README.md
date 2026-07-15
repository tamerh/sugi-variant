# Sugi Variant — per-variant genetic reference (prototype)

Separate product spun out of Sugi Atlas. Dynamic FastAPI + Jinja2 + ISR renderer over
the self-contained `sugivariant` science package (`sugivariant.build.enriched_records`).
Server-rendered HTML (crawlable/AI-citable — the differentiator), disk/ETag-cached.
Sibling of Sugi Predict; reuses `atlas.css`.

## 👉 New here? Read `HANDOVER.md` first.
It has the full context: why the product exists, the domain rules you must not break
(ClinGen SVI predictor framing, patient-safety), biobtree quirks, corpus scale, and the
TODO roadmap. Design spec: `docs/VARIANT_PAGES_SPEC.md`. Audit history:
`docs/variant-audit-2026-07.md`.

## Run (dev)
    ./run.sh                      # uvicorn :8013, BASE_PATH=/variant
Then, behind nginx: https://sugi.bio/variant/pten-p-arg130gln
Local (no proxy): http://127.0.0.1:8013/pten-p-arg130gln

Env: `bioyoda` (`python`) — fastapi + the shared
`sugibiobtree` client. biobtree REST API must be up at localhost:9291.

## Deploy mapping (sugi.bio/variant)
nginx `location /variant/ { proxy_pass http://127.0.0.1:8013/; }` (see
the deploy repo/nginx-sugi.bio.conf). Trailing slash strips /variant/, so the
app serves root paths; BASE_PATH=/variant generates the public /variant/… links.

## Status
Prototype: verdict-card + evidence-panel view. TODO: search, protein lollipop SVG,
JSON-LD/.md twin, sitemap, slug→VCV index, Docker/systemd.

## Dependency: sugibiobtree (shared client)
The variant science uses the shared `sugibiobtree` biobtree client (one source of
truth with Sugi Atlas). Every env that runs this app must have it installed:

    pip install -e the sugi-biobtree repo

Dev: installed editable in the `bioyoda` env. When containerizing, the image must
`pip install` sugibiobtree (editable during dev, or pinned/vendored for prod).
