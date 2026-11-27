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

### Phase 0 — inventory and gates (1–2 days)

- Oracle snapshot: `path → sha256 + size` for all 54,361 rendered pages.
- Machine-readable table of every 301 (disease renames, gene legacy) and the
  trailing-slash rule, with a test that each lands on a 200.
- Decide and pin how the body fragment is obtained: a Hugo/Atlas layout that emits
  `.Content` (clean, one-line-ish), or extraction from the existing HTML (no Atlas
  change, but brittle). Prefer the layout output; use extraction only to bootstrap.
- **Exit:** snapshot + redirect table + fragment method committed; nothing served.

### Phase 1a — routing and semantics, no renderer (2–3 days)

- App serves `/atlas/{type}/{slug}/`; body temporarily from today's rendered HTML.
- 404/301/trailing-slash normalisation; `ETag` from content hash;
  `Cache-Control: public, s-maxage=86400, stale-while-revalidate=604800`.
- `.md` twin at `/atlas/{type}/{slug}/index.md` (`text/markdown`).
- Response headers: release, version, released-at, biobtree.
- Side-by-side diff vs nginx on a sample: status, content-type, canonical, title, H1,
  body hash.
- **Exit:** every sampled URL matches; all redirect rules pass; fully reversible.

### Phase 1b — the pickle corpus (≈1 week)

- Atlas build writes `body_html` fragment + metadata + `content_hash` per entity into
  `atlas-cache/<type>/<slug>.pkl.gz`, plus `atlas_index.db`.
- The app serves from pickles; Hugo leaves the serve path.
- Diff against the oracle: structure (headings, tables, links) then text; fix to
  parity on the sample, then a full 54k overnight diff.
- **Exit:** zero unexplained diffs; `.md` twins byte-identical to today's.

### Phase 2 — one manifest and one search (3–5 days)

- `release.json` grows `corpora.atlas` (type → entity → path + hash). 54k entries is a
  few MB: load lazily or put it behind `atlas_index.db`.
- Pagefind rebuilt per release and served by the app.
- One sitemap index: atlas entity sitemaps (50k chunks) + variant gene sitemaps.
- `/resolve.json` resolves gene / disease / drug entities.
- **Exit:** one manifest describes everything; release + rollback exercised end to end
  on dev.

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
