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

## What is not done here

- No pickles yet (Phase 1b); the app still reads the HTML tree.
- No fragment/shell rendering and no CSS merge (Phase 1b).
- No changes to the Atlas build, to nginx, or to production.
