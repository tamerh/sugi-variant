# Development

## Run

    export ATLAS_BIOBTREE=http://<biobtree-host>:9291   # a running BioBTree REST API
    python -m uvicorn app:app --host 127.0.0.1 --port 8000

Then open <http://127.0.0.1:8000/>. Optional environment variables:

- `BASE_PATH` — URL prefix when served under a subpath behind a reverse proxy (e.g. `/variant`).
- `CACHE_DIR` — page/records cache and the prebuilt resolution `index.db` (default `./cache`).
- `GENE_CACHE_MB` — in-process LRU budget for the hottest genes (default `1500`).

## Tests

    python -m pytest -q

Some tests read the built record corpus, which is far too large to ship with the code and
is not in this repository. They take its location from the environment and skip when it is
absent, so a fresh clone runs the rest of the suite without any setup:

- `RECORDS_DIR` — directory of per-gene `<GENE>.pkl.gz` record files.
- `INDEX_DB` — the prebuilt resolution `index.db`.

## Operations

Running an instance (deploying the image, transforming records in place, backups, previews)
is driven by a single script kept in a private deployment repository together with the
container and reverse-proxy configuration, because it necessarily names hosts, ports,
container names and absolute paths. Nothing here depends on it: this repository holds the
application and the method, and the commands above are all that is needed to run it.

## Dependencies

    pip install fastapi "uvicorn[standard]" httpx jinja2

Plus the shared `sugibiobtree` BioBTree client (one source of truth with Sugi Atlas), installed
from the `sugi-biobtree` repository.

## Layout

- `app.py` — the FastAPI app (routes, caching, identifier resolution).
- `sugivariant/` — the science package (collect / enrich / render / index / slug).
- `templates/`, `static/` — Jinja views and assets.
- `static/preprint/` — the built preprint, served at `/method`.
