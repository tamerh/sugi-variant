# Bringing Sugi Atlas into the Variant app

Status: plan, not started. Branch `feat/atlas-in-variant`, built and previewed on dev
only; production is untouched until Phase 4.

## 1. The decision

**One app, two corpora, one release.**

Atlas becomes a first-class corpus in the Variant service, exactly like
`small_variant` and `structural_variant`: same manifest, same process, same cache
headers, same `.md` twins, same provenance block, same rollback.

Three things stay deliberately separate:

- **Atlas's content generation** — its passes, collectors, section renderers and
  prose. Not one line changes in how a page is *built*.
- **Atlas's content model** — long sectioned entities, not Variant's evidence rows
  and `present.py` states. Forcing them together is a rewrite, not a migration.
- **Variant's record build** — unchanged.

What changes is only the *sink* and the *serving*: instead of `index.md` → Hugo →
54k files, the build writes one record per entity, and the app serves it.

## 2. Why

- **Maintenance.** Today: two template languages, two headers, two searches, two
  release mechanisms, "did the symlink swap?". After: one of each. The Hugo build
  and the 54k-file tree disappear.
- **Composition.** A `/gene/X` can show Atlas sections *and* the variant list because
  both are records in the same process — impossible with static files without JS or
  build-time coupling.
- **Release discipline.** Atlas inherits the manifest, the headers, the "not in the
  current release" semantics, and rollback-as-a-pointer.
- **Fewer artifacts.** The `dist` Markdown tree, the Hugo `public` tree and the
  corpus tarballs are three representations of the same 54k pages.

Cache invalidation is **not** an argument either way: URLs are stable in both models,
so a release changes content behind an unchanged URL in both. The only difference is
what a cache miss costs (file read vs 10–50 ms render), handled with `s-maxage` +
`stale-while-revalidate` + content-hash ETags.

## 3. Today

```
ATLAS (static)                                VARIANT (dynamic)
┌───────────────────────────────┐        ┌────────────────────────────────┐
│ Atlas passes (collectors)     │        │ corpus refresh (shards)        │
│   └─► dist/atlas/*/index.md   │        │   ├─ records/*.pkl.gz   (19k)  │
│            │                  │        │   ├─ index.db     (4.2M rows)  │
│            ▼                  │        │   └─ sv.db        (83k SVs)    │
│         HUGO build            │        │            │                   │
│            ▼                  │        │       release.json             │
│ /var/www/sugi.bio/atlas       │        │            ▼                   │
│  54,361 HTML files            │        │ variant-web :8014 (uvicorn)    │
└───────────────┬───────────────┘        └───────────────┬────────────────┘
                └──────────► nginx (files | proxy) ◄─────┘
                                    ▼
                                Cloudflare
```

Inventory measured on this host:

| | |
|---|---|
| Atlas entities | 29,336 genes + 18,596 diseases + 6,419 drugs = **54,351** |
| Atlas source | `/data/sugi-atlas/dist/atlas/<type>/<slug>/index.md` + `entity.jsonld` |
| Atlas oracle | `/var/www/sugi.bio/atlas/…` — **54,361** rendered HTML pages |
| Redirects today | ~110 rename rules in nginx snippets (+ gene legacy rules) |
| Variant corpus | 4,222,508 pages, 19,047 gene files, 83,042 SVs |

The oracle is the key asset: it gives a page-for-page test target for the renderer.

## 4. Target

```
        ┌──────────────────────────── BUILD TIME ────────────────────────┐
        │  Atlas passes (UNCHANGED)          Variant refresh (UNCHANGED) │
        │  └─ gene/disease/drug              └─ records + index + sv.db  │
        │        │                                     │                 │
        │        ▼                                     ▼                 │
        │  atlas-cache/gene/BRCA1.pkl.gz        variant cache           │
        │  atlas_index.db                              │                 │
        │        └──────────────► ONE release.json ◄───────────────────┘
        │                         corpora.atlas            (54k)         │
        │                         corpora.small_variant    (19k)         │
        │                         corpora.structural_variant (83k)       │
        └──────────────────────────────┬─────────────────────────────────┘
                                       ▼
        ┌──────────────────────────── SERVE TIME ────────────────────────┐
        │                 ONE FastAPI app  (one image)                   │
        │  /atlas/{type}/{slug}/     ← atlas pickles                     │
        │  /variant/{slug}           ← variant pickles                   │
        │  /gene/{SYM}               ← composes both (Phase 3)           │
        │  /sv/{key}                 ← sv.db                             │
        │  /resolve.json  /sitemap-*.xml  /about                         │
        │  shared shell · About/Sources block · .md twins                │
        └──────────────────────────────┬─────────────────────────────────┘
                                       ▼
                                  nginx ──► Cloudflare
                              /atlas  /variant → same app
```

Always available, off to the side:

```
┌────────────────────────────────────────────────────┐
│ NIGHTLY STATIC EXPORT (same renderer)              │
│   flat HTML + smoke diff of ~50 pages              │
│   = tested cold standby, not a folder that rots    │
└────────────────────────────────────────────────────┘
```

## 5. The Atlas record

```
atlas-cache/gene/BRCA1.pkl.gz
├── type: gene            slug: BRCA1
├── title                 canonical_url
├── body_html             # the fragment the build already renders; no re-render
├── body_md               # the .md twin, same source
├── sections[]            # {id, title, anchor} for the TOC and ?section=
├── jsonld
└── content_hash          built_at            release
```

`body_html` is a **fragment**, not a full page: the app wraps it in the shared
header/footer and provenance block. That means Atlas keeps its rendering, and the app
keeps its shell — neither side re-implements the other.

`atlas_index.db` mirrors `index.db`'s role: `(type, slug) → path, content_hash,
title`, used for hit-vs-404, sitemaps and the resolver.

### Why a pickle, not 54k HTML files

- atomic per-entity replace (no directory sync, no half-updated site)
- one manifest → rollback is a pointer
- one LRU with its own budget (`ATLAS_CACHE_MB`, separate from `GENE_CACHE_MB`)
- `.md` and HTML from one record, so they cannot drift
- the same release and rollback semantics already trusted for Variant

## 6. Non-negotiables

These are the conditions that make the merge safe; they are deliverables, not cleanup.

1. **404 semantics.** An unknown entity is `HTTP 404` with the not-found page — never
   a 200 that says "not found" (soft-404). A test asserts the status code per route.
   An entity present in the corpus but absent from the release is also 404.
2. **301s move with the routes.** The ~110 rename rules and the gene legacy rules must
   keep working: either stay in nginx ahead of the app, or move in-app with a table
   and a test that every rule lands on a 200. Forgetting them 404s renamed entities
   with no gate to catch it.
3. **Trailing slash.** Atlas serves `/atlas/gene/BRCA1/`; Variant serves
   `/variant/<slug>` without a slash. One app must normalise deliberately (301 to the
   canonical form) or URLs duplicate.
4. **Render diff per release.** For entities whose `content_hash` changed, diff the new
   HTML against the previous release's and fail if the delta isn't explained. This is
   the gate that catches a template regression (it would have caught the sitemap
   revert).
5. **Nightly export + smoke diff.** ~50 representative pages (one per entity type, plus
   the largest and the emptiest) compared against the live app. An unexercised standby
   is not a standby.
6. **Pagefind per release.** Search is a build artifact in both worlds; regenerate and
   serve it. "One app" must not mean "search got worse".
7. **Per-corpus cache budgets.** A burst of large Atlas gene pages must not evict the
   variant corpus.
8. **Versioned cache keys / manifest pointer** so a release never purges the previous
   one and rollback is a pointer plus a restart.

## 7. Phases

### Phase 0 — inventory and gates — **done 2026-11-27**

Tooling in `migration/` (`atlas_oracle.py`, `nginx_redirects.py`), outputs in
`/data2/atlas-migration/`.

- **Oracle** — `oracle-20261127.jsonl`, 108,713 lines (54,360 HTML + 54,353 Markdown),
  sha256 + size + mtime per page, 6.7 s to build. `verify` re-walks and reports
  added/removed/changed: the "did the corpus move under us?" gate.
- **Redirects** — `redirects.json`, **363 rules**, **363/363 land on a page present in
  the oracle, 0 misses**. Verified offline (a live probe is blocked by Cloudflare);
  offline catches the failure that matters, a rename pointing at a page that is gone.
- **Fragment** — every page has exactly one `<article class="entity-page">…</article>`,
  93% of the page bytes, no inline CSS, only JSON-LD scripts inside. Phase 1a extracts
  it as a bootstrap; Phase 1b replaces extraction with a `partials/entity-page.html` +
  `fragment` output format in the Atlas web repo, which is the permanent form.
- **Exit met:** snapshot + table + method committed; nothing served.

### Phase 1a — routing and semantics, no renderer — **done 2026-11-27**

- App serves `/atlas/{type}/{slug}/`; body temporarily from today's rendered HTML.
- 404/301/trailing-slash normalisation; `ETag` from content hash;
  `Cache-Control: public, s-maxage=86400, stale-while-revalidate=604800`.
- `.md` twin at `/atlas/{type}/{slug}/index.md` (`text/markdown`).
- Response headers: release, version, released-at, biobtree.
- Side-by-side diff vs nginx on a sample: status, content-type, canonical, title, H1,
  body hash.
- **Exit met, measured against the oracle on a frozen copy of the tree** (5.7 GB at
  `/data2/atlas-migration/atlas-html`, read-only):
  - **54,360 / 54,360 pages byte-identical** (sha256 vs the oracle)
  - **218 / 218 `/atlas/*` redirects** 301 to their successor, and every successor is 200
  - **145 / 145 `/biobtree/*` legacy targets** resolve (the rules themselves stay in nginx)
  - 25/25 trailing-slash 301s, 25/25 `.md` twins identical and `text/markdown`
  - 9 tests in `tests/test_atlas_web.py` (404, 301, traversal, ETag/304, disabled mode)
- Enabled only when `ATLAS_ROOT` is set, so a deployment that does not opt in is
  unchanged. HTML responses inherit the app's ETag middleware; `.md` carries its own.
- Preview: `dev.sugi.bio/variant/atlas/…` (dev basic-auth), no nginx change. Internal
  links inside the stored HTML are absolute `/atlas/…`, so navigation leaves the
  preview for the live site — expected until the fragment/shell phase.

### Phase 1b — the pickle corpus — **done 2026-11-27**

- Atlas build writes `body_html` fragment + metadata + `content_hash` per entity into
  `atlas-cache/<type>/<slug>.pkl.gz`, plus `atlas_index.db`.
- The app serves from pickles; Hugo leaves the serve path.
- Diff against the oracle: structure (headings, tables, links) then text; fix to
  parity on the sample, then a full 54k overnight diff.
- **Exit met:** built from the frozen tree with `migration/build_atlas_corpus.py`:
  **54,351 entities, 0 without a fragment, 1.1 GB** (from 5.7 GB), 95 s.
  `atlas_index.db` carries type/slug/title/canonical/fragment sha/md sha/sections.
  The app renders the stored fragment inside the shared shell (new `atlas_entity.html`),
  loads the site stylesheet (`static/css/atlas-site.css`) on Atlas pages only, and
  rewrites `/atlas/…` links through `BASE` so the dev prefix survives.
  Gate (`migration/verify_atlas_entities.py`): **3,039 of 3,039 sampled entities serve
  the exact fragment** (sha256 vs the record built from the static tree) and the `.md`
  twin matches. Site-level pages (about, type indexes) still come from the tree.
  Two bugs found and fixed by the gate: a SQLite connection shared across threads
  (`InterfaceError` under the concurrent verifier, surfacing as sporadic 404/500) and a
  plain `str` route return that FastAPI JSON-encoded.
- **The fragment is not the whole page: Atlas enhances it in the shell.** The static page
  carries an 8.7 KB inline script that rewrites `.entity-body` at runtime — H2 sections
  become `<details>` (only the first open) with an Expand/Collapse-all button, the TOC
  gets chevrons, open/closed state and a filter, hash links open collapsed sections, and
  every table is wrapped with a Copy CSV button. My first gate compared the `<article>`
  bytes and passed while the app page rendered as one flat document, because that script
  lives outside the fragment. Fixed by copying it verbatim to
  `static/js/atlas-entity.js` (loaded only on entity pages). **Any Atlas shell script that
  touches the entity body belongs with the fragment, and a fragment-only gate cannot see
  its absence.**
  Gate (`migration/verify_atlas_render.py`, headless Chromium): both pages are loaded, the
  script runs, and the post-JS `<article>` DOM is compared. **26 of 26 identical**,
  including the five largest entities, with the same enhancement stats (sections, which
  one is open, table wrappers, Copy CSV buttons, TOC chevrons, filter, Expand all) and no
  console errors. Only the TOC's sticky `top` is normalised, because that is a measurement
  of the shell the app owns.

### Phase 2 — one manifest and one search — **done 2026-11-27**

- **Manifest.** `release.json` grows `corpora.atlas`: cache/index/root/redirects paths,
  54,351 entities by type, fragment and on-disk bytes, `with_markdown`, `built_at`, the
  363 redirect rules and sha256 anchors for the index, the site home, the search data and
  the table. No entity list: the index is the membership, and 54k entries in the manifest
  would only add something that can drift from it.
- **Release mode.** The app adopts the block, resolves a missing path by basename under
  `CACHE_DIR/atlas` (a container mounts the cache elsewhere; same rule as `sv.db`), and
  re-checks the anchors at startup. A mismatch is a 503 on `/atlas` naming the reason,
  not stale pages. Membership is the released index: a page that exists in the tree but
  is not in the release 404s. Proved with a subset index (KCNA2 removed): 404 for KCNA2,
  200 for BRCA1, with the KCNA2 tree file still on disk.
- **Sitemaps.** `/sitemap.xml` advertises `/sitemap-atlas-1.xml` and `-2.xml` (50k per
  urlset, the sitemaps.org cap): 50,000 + 4,351 = 54,351, exactly the index. `lastmod` is
  the build date from the block, not a file mtime that a mount would rewrite.
- **Resolver.** `/resolve.json` answers disease and drug names, which have no variant page
  at all (`kind: atlas_disease` / `atlas_drug`), and puts `atlas_url` beside a gene hit so
  one request covers both corpora without turning a gene into "ambiguous". The entity
  index now carries Atlas's own name and synonym maps: **245,986 names → 216,642 alias
  rows**, so "carasil syndrome" and "CARASIL2" resolve without loading 13 MB of JSON.
- **Search — the plan was wrong here.** Atlas has no Pagefind: its search is a
  client-side match over a 13 MB `manifest.json`. Inventing Pagefind would have been a
  content/behaviour change, not a migration. What the migration must do is not *lose*
  search, so: `/suggest.json` now merges Atlas entities after the variant suggestions (a
  drug or disease name used to return nothing at all), and `/atlas/search/` with its
  `manifest.json` still ships per release. Full-text search is a separate decision.
- **Open, and deliberately not fixed here:** the app header carries no search box — search
  lives on the home page and the error page — so an Atlas page currently has no search
  control at all, where the Atlas header had one. Adding a site-wide search box is a
  navigation change that touches every Variant page, i.e. the umbrella/nav decision of
  Phase 3, and it must not be slipped in under "Phase 2". The endpoint it would use is
  ready and tested (`/suggest.json` returns variants and Atlas entities), and the render
  gate drives the real box end to end (type `imatinib`, keyboard-pick, land on
  `/atlas/drug/imatinib/`).
- **Exit met:** one manifest describes the small-variant, structural-variant and Atlas
  corpora; release and rollback exercised on dev — a bad index sha gives 503 with the
  reason, a subset index gives a real 404 for the missing entity, a good manifest serves
  the page, its `.md` twin, the 301 table, `/atlas/about` and `X-Sugi-Variant-Release`.

### Phase 3 — compose `/gene/X` (≈1 week)

- Atlas gene sections + Variant variant list + counts in one page, one provenance
  block, one `.md` twin.
- Build the mechanism now; decide the canonical URL with the umbrella/naming decision.
- **Exit:** the composed page passes the render-diff gate and the 404/301 suite.

### Phase 4 — standby and cutover (3–5 days)

- Nightly static export from the app renderer + smoke diff.
- Cutover (owner-run): nginx `/atlas` → app; static tree kept read-only for rollback.
- Delete Hugo and the tree **fully** once parity and the standby have been clean for a
  release or two — "we kept both paths" is how maintenance doubles instead of halving.
- Monitoring: 5xx rate, p95, cache hit ratio.

## 8. Environments

- **Preview:** app on `:8013`, `https://dev.sugi.bio/variant/atlas/…`.
  nginx strips `/variant/`, so the app receives `/atlas/…`, and because `BASE=/variant`
  the links it emits route back through the same prefix. The merged app is previewable
  with **zero nginx changes**.
- `test.sugi.bio` is not usable for this: its root is the static staging dir
  (`/data2/sugi-staging`, symlink-swapped by `stage.sh`), which suits static releases.
- **Prod** (`variant-web:8014`, `/var/www/sugi.bio`) is read-only to us until Phase 4.
  The only prod step is the nginx location for `/atlas`, owner-run at cutover.

## 9. Risks

| Risk | Mitigation |
|---|---|
| Soft-404s (SEO) | status-code tests per route; corpus-vs-release distinction |
| Lost 301s | redirect table + test; nginx rules can stay ahead of the app |
| Renderer drift | body fragment stored, not re-rendered; per-release render diff |
| Atlas pages are heavy (BRCA1 ≈ 256 KB) | gzip per record, separate LRU budget, ETag/304 |
| Release coupling (one image) | corpora stay independently versioned; code deploys are canary-first |
| Standby rots | nightly export + smoke diff gate |
| Search regression | Pagefind rebuilt per release |

## 10. Out of scope (deliberate)

- **Unifying the content models.** Atlas keeps its sectioned prose; Variant keeps its
  evidence rows. One runtime, two builders.
- **Changing Atlas URLs.** `/atlas/{type}/{slug}/` stays byte-identical through the
  migration; any change is a later, redirect-backed decision.
- **Umbrella naming.** Whether the product is "Sugi Atlas" with a Variant section, or
  two labels under one site, is decided after the migration works. "Sugi Variant" is
  cited by the preprint/DOI, so the label is not removed casually.
- **Server-side search.** Pagefind stays; replacing it is a separate decision.

## 11. Open decisions

1. Fragment source: Atlas layout emits `.Content`, or extraction from rendered HTML?
   (Prefer the layout; extraction only to bootstrap.)
2. `atlas_index.db` vs a 54k-entry JSON block in the manifest for hit-vs-404.
3. Whether the Atlas builder writes pickles directly, or writes Markdown and a
   separate packer converts — the latter keeps the Atlas repo untouched in Phase 1b.
4. Canonical `/gene/X` (Phase 3) vs keeping `/atlas/gene/X/` canonical.
