# Migration tooling (Atlas → Variant app)

Phase 0 of [docs/atlas-migration.md](../docs/atlas-migration.md): freeze what the live
Atlas site is, so every later phase proves parity instead of asserting it.

Nothing here touches production. Outputs go to `/data2/atlas-migration/`.

## 1. Oracle snapshot — `atlas_oracle.py`

```bash
python migration/atlas_oracle.py snapshot \
    --root /var/www/sugi.bio/atlas \
    --out /data2/atlas-migration/oracle-20261127.jsonl

python migration/atlas_oracle.py verify \
    --root /var/www/sugi.bio/atlas \
    --snapshot /data2/atlas-migration/oracle-20261127.jsonl
```

Measured 2026-11-27, release live at the time:

| | |
|---|---|
| entity pages (HTML) | **54,360** under a type, of which 54,351 are entities |
| — genes | 29,336 |
| — diseases | 18,596 |
| — drugs | 6,419 |
| site-level pages | 9 (`/atlas/about/`, the type indexes, …) |
| Markdown twins | 54,353 |
| largest page | `/atlas/gene/TP53/` 304,543 B |
| snapshot size | 108,713 lines, 28 MB, 6.7 s to build |

`verify` re-walks the tree and reports added / removed / changed by sha256 — the gate
for "did the corpus change under us".

## 2. Redirect table — `nginx_redirects.py`

```bash
python migration/nginx_redirects.py extract \
    --conf /data/ext/deploy/sugi-disease-rename-redirects.conf \
    --conf /data/ext/deploy/sugi-gene-atlas-redirects.conf \
    --conf /data/ext/deploy/sugi-disease-atlas-redirects.conf \
    --out /data2/atlas-migration/redirects.json

python migration/nginx_redirects.py verify \
    --table /data2/atlas-migration/redirects.json \
    --oracle /data2/atlas-migration/oracle-20261127.jsonl
```

Measured: **363 rules** (every `location = … { return 301 …; }`), **363/363 land on a
page that exists in the oracle, 0 misses**. The verification is offline on purpose — a
live probe from this host is blocked by Cloudflare bot protection, and an offline
target check catches the failure that matters (a rename pointing at a page that no
longer exists).

The table is the input for the Phase 1a test suite: every rule must 301 and the target
must 200, in the app.

## 3. Fragment method — decision

The app will store the **body fragment**, not a re-rendered page. Two ways to get it:

- **Preferred (Phase 1b):** refactor `hugo/layouts/atlas/single.html` so its article
  body is a partial (`partials/entity-page.html`), then add a `fragment` output format
  that renders only that partial. One small change in the Atlas web repo, no parsing.
- **Bootstrap (Phase 1a):** extract the fragment from the existing HTML. Verified safe:
  every page has **exactly one** `<article class="entity-page">…</article>`, it is 93%
  of the page bytes (TP53: 304 KB page, 283 KB article), there is **no inline CSS**
  inside it, and the only scripts are two `<script type="application/ld+json">` blocks
  (JSON-LD, which the app can also take from `entity.jsonld`).

The bootstrap extraction is what Phase 1a uses; the layout refactor is what makes it
permanent and removes the HTML parsing.

## 4. Passthrough verifier — `verify_atlas_passthrough.py`

Phase 1a serves the tree verbatim from the app (`sugivariant/atlas_web.py`, enabled by
`ATLAS_ROOT`) and this script proves it:

```bash
ATLAS_ROOT=/data2/atlas-migration/atlas-html \
ATLAS_REDIRECTS=/data2/atlas-migration/redirects.json \
  python -m uvicorn app:app --port 8095 &

python migration/verify_atlas_passthrough.py --base http://127.0.0.1:8095 \
    --oracle /data2/atlas-migration/oracle-20261127.jsonl \
    --redirects /data2/atlas-migration/redirects.json --sample 300
```

Result 2026-11-27:

| check | result |
|---|---|
| byte parity (sha256) | **54,360 / 54,360** identical |
| trailing-slash 301 | 25 / 25 |
| `.md` twins | 25 / 25 identical, `text/markdown` |
| `/atlas/*` redirects | **218 / 218**, targets 200 |
| `/biobtree/*` legacy targets | 145 / 145 resolve (rules stay in nginx) |

The source tree is a frozen copy at `/data2/atlas-migration/atlas-html` so the gate does
not depend on what production serves while we work.

## 5. Pickle corpus — `build_atlas_corpus.py`

```bash
python migration/build_atlas_corpus.py \
    --html /data2/atlas-migration/atlas-html \
    --out  /data2/atlas-migration/atlas-cache \
    --index /data2/atlas-migration/atlas_index.db --apply
```

| | |
|---|---|
| entities | **54,351** (29,336 genes, 18,596 diseases, 6,419 drugs) |
| without a fragment | 0 |
| corpus size | **1.1 GB** (from 5.7 GB of HTML), 95 s |
| per entity | `body_html` fragment, `body_md`, title, canonical, sections, shas |

The app reads it when `ATLAS_CACHE` + `ATLAS_INDEX` are set; site-level pages still
come from the tree. Parity gate:

```bash
python migration/verify_atlas_entities.py --base http://127.0.0.1:8114 \
    --index /data2/atlas-migration/atlas_index.db --sample 3000 --threads 8
# entity content: 3,039/3,039 identical · RESULT: PASS
```

That gate compares the article bytes, which is necessary but not sufficient: Atlas
enhances the entity body in the *page shell* (collapsible sections, TOC chevrons and
filter, Copy CSV on tables). The app carries that script verbatim in
`static/js/atlas-entity.js`; `verify_atlas_render.py` loads both pages in headless
Chromium, runs the script, and compares the post-JS DOM:

```bash
# server on :8115 with BASE_PATH="" so links match the static tree
LD_LIBRARY_PATH=/data/miniconda3/envs/bioyoda/lib:/data/miniconda3/envs/biobtree/lib:\
/data/miniconda3/pkgs/libgbm-1.0.7-h54a6638_0/lib:/data/miniconda3/pkgs/xorg-libxshmfence-1.3.3-hb9d3cd8_0/lib \
python migration/verify_atlas_render.py --base http://127.0.0.1:8115 --sample 25
# rendered article: 26/26 identical · RESULT: PASS
```

(The `LD_LIBRARY_PATH` is only because this host lacks Chromium's system libs; they are
in the conda envs.)

## 6. Release mode — both corpora from one manifest

```bash
# what `variant.sh release-manifest` now emits (Atlas parts optional)
PYTHONPATH=. python -m sugivariant.snapshot build \
    --records /data2/sugi-variant-cache-dev/records \
    --index   /data2/sugi-variant-cache-dev/index.db \
    --sv-db   /data2/sugi-variant-cache-dev/sv.db \
    --atlas-index /data2/atlas-migration/atlas_index.db \
    --atlas-cache /data2/atlas-migration/atlas-cache \
    --atlas-root  /data2/atlas-migration/atlas-html \
    --atlas-redirects /data2/atlas-migration/redirects.json \
    --out /tmp/full.json --serve-out /tmp/release-atlas.json \
    --release v0.0.0-test --released-at 2026-11-27
```

`corpora.atlas` carries the paths, 54,351 entities by type, `fragment_bytes`,
`corpus_bytes`, `with_markdown`, `built_at`, the 363 redirect rules and sha256 anchors
for the index, the site home page, the search data and the table. The entity list is not
in the manifest: the index is the membership.

Serve it (no `ATLAS_*` environment needed — the manifest drives it):

```bash
VARIANT_SERVE=release VARIANT_RELEASE=/tmp/release-atlas.json \
CACHE_DIR=/data2/sugi-variant-cache-dev INDEX_DB=/data2/sugi-variant-cache-dev/index.db \
python -m uvicorn app:app --port 8124
```

Release discipline, all exercised on dev:

| case | result |
|---|---|
| manifest matches the files | 200, `X-Sugi-Variant-Release` on every response |
| index sha in the manifest is wrong | **503** on `/atlas` with the reason, variants unaffected |
| entity in the tree but not in the released index (KCNA2 removed from a copy) | **404** — no fall-through to the tree |
| `.md` twin, 301 rename, `/atlas/about` | unchanged |

Installed layout for a prod cache: `$CACHE_DIR/atlas/{atlas_index.db, atlas-cache/,
atlas-html/, redirects.json}`. Manifest paths are host paths, so a container resolves a
missing one by basename under `CACHE_DIR/atlas` (the rule `sv.db` already follows).

## What is not done here

- No release integration: the manifest, Pagefind, sitemaps and `/resolve.json` entity
  hits are Phase 2.
- No changes to the Atlas build, to nginx, or to production.
