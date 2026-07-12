# Sugi Variant — dev web app (prototype)

Dynamic FastAPI + Jinja2 + ISR renderer over the deterministic variant enrichment
(shared with Sugi Atlas: `atlas.variant.enriched_records`). Server-rendered HTML
(crawlable/AI-citable), disk/ETag-cached. Sibling of Sugi Predict; reuses `atlas.css`.

## Run (dev)
    ./run.sh                      # uvicorn :8013, BASE_PATH=/variant
Then, behind nginx: https://sugi.bio/variant/pten-p-arg173cys/
Local (no proxy): http://127.0.0.1:8013/pten-p-arg173cys/

Env: bioyoda (has fastapi + the atlas deps). ATLAS_SRC points at the atlas src.

## Deploy mapping (sugi.bio/variant)
nginx `location /variant/ { proxy_pass http://127.0.0.1:8013/; }` (see
the deploy repo/nginx-sugi.bio.conf). Trailing slash strips /variant/, so the
app serves root paths; BASE_PATH=/variant generates the public /variant/… links.

## Status
Prototype: verdict-card + evidence-panel view. TODO: search, protein lollipop SVG,
JSON-LD/.md twin, sitemap, slug→VCV index, Docker/systemd.
