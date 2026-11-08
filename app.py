"""Sugi Variant — the web. A FastAPI + Jinja2 + ISR renderer over the self-contained
`sugivariant` package (decoupled from Sugi Atlas). Server-rendered
HTML (crawlable / AI-citable — the whole differentiator), disk-cached after first
hit. Sibling of the Sugi Predict app; reuses the atlas.css design system.

Run: uvicorn app:app --host 0.0.0.0 --port 8013
"""
import hashlib
import json
import os
import pathlib
import re
import sys

from fastapi import FastAPI, Request
from fastapi.responses import (HTMLResponse, Response, RedirectResponse,
                               PlainTextResponse, JSONResponse)
from fastapi.staticfiles import StaticFiles
from starlette.exceptions import HTTPException as StarletteHTTPException
from jinja2 import Environment, FileSystemLoader, select_autoescape

# The variant science lives in the local, self-contained sugivariant package.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sugibiobtree import map_all                          # noqa: E402
from sugivariant.build import enriched_records            # noqa: E402
from sugivariant.collect import collect                   # noqa: E402
from sugivariant.render import _label as variant_label    # noqa: E402
from sugivariant.render import meta_description as _meta_desc  # noqa: E402
from sugivariant.render import page_title as _page_title  # noqa: E402
from sugivariant.slug import _norm, _norm_hgvs            # noqa: E402
from sugivariant.enrich import refresh_record as _refresh  # noqa: E402
from sugivariant.enrich import records_collapsed  # noqa: E402
from sugivariant.present import (  # noqa: E402
    evidence_groups, evidence_view, transcript_view)
from sugivariant.enrich import (review_stars, review_tier, tier_label,  # noqa: E402
                                missense_short, disagreement_flag,
                                DISAGREEMENT_CATEGORIES)
from sugivariant.enrich import (similar_index as _build_similar_index,  # noqa: E402
                                similar_variants_from_index as _similar_from_index)
from sugivariant.enrich import (begin_build as _begin_build,   # noqa: E402
                               build_degraded as _build_degraded)
from sugibiobtree.client import BiobtreeError          # noqa: E402
from sugivariant import index as IX                       # noqa: E402
from sugivariant import dbvar as SV                       # noqa: E402

ROOT = pathlib.Path(__file__).parent
BASE = os.environ.get("BASE_PATH", "").rstrip("/")
CACHE_DIR = pathlib.Path(os.environ.get("CACHE_DIR") or (ROOT / "cache"))
SERVE_MODE = (os.environ.get("VARIANT_SERVE") or "auto").strip().lower()
RELEASE_MODE = SERVE_MODE in ("release", "snapshot")
RELEASE_PATH = pathlib.Path(os.environ.get("VARIANT_RELEASE")
                             or os.environ.get("VARIANT_SNAPSHOT")
                             or (CACHE_DIR / "release.json"))

env = Environment(loader=FileSystemLoader(str(ROOT / "templates")),
                  autoescape=select_autoescape(["html"]))
env.globals["base"] = BASE
env.globals["release"] = None
env.globals["released_at"] = None
env.globals["css_v"] = lambda: str(int(max(
    (ROOT / "static" / "css" / f).stat().st_mtime for f in ("app.css", "atlas.css"))))

# ── classification → colour class + display + star rendering (template helpers) ─
_CLS_CLASS = {
    "pathogenic": "path", "likely pathogenic": "lpath",
    "pathogenic/likely pathogenic": "path",
    "conflicting classifications of pathogenicity": "conflict",
    "uncertain significance": "vus", "likely benign": "lben", "benign": "ben",
    # Without this the only call that MERGES benign and likely-benign evidence fell
    # through to the default and rendered in the VUS grey -- the same colour as
    # "uncertain" -- on 73,495 variants, their gene-hub pills and their /set tiles.
    "benign/likely benign": "ben",
    # dbVar calls a submitter conflict "conflicting data from submitters"; render it
    # in the same conflict colour as ClinVar's "Conflicting classifications".
    "conflicting data from submitters": "conflict",
}


def cls_class(classification):
    return _CLS_CLASS.get((classification or "").strip().lower(), "vus")


def stars(review_status):
    # named review-status tier (not a ★ rating) — {n, label, cls}
    return review_tier(review_status)


from sugivariant.render import short_hgvs as _short_hgvs   # noqa: E402
env.globals.update(label=variant_label, cls_class=cls_class, stars=stars, hgvs_disp=_short_hgvs,
                   page_title=_page_title, meta_description=_meta_desc)
from sugivariant import links as _links   # noqa: E402
env.globals["atlas_gene"] = lambda sym: _links.gene_url(symbol=sym)
env.globals["atlas_disease"] = lambda name: _links.disease_url(name=name)
# Data-source attribution (AlphaMissense CC BY 4.0 + REVEL ODbL require it) — single
# source of truth in render.py so HTML and the markdown twin can't drift.
from sugivariant.render import (data_provenance, source_refs, data_asof,  # noqa: E402
                               data_asof_range, source_asof, anc_name, one_in,
                               GRPMAX_EXCLUDED)
env.globals["data_provenance"] = data_provenance
env.globals["source_refs"] = source_refs
env.globals["data_asof"] = data_asof
env.globals["data_asof_range"] = data_asof_range
env.globals["source_asof"] = source_asof
env.globals["tier_label"] = tier_label
# `v.spliceai` is ALWAYS a dict carrying a status now, so {% if v.spliceai %} is true for
# "absent" and "snv_only" too — which had the page claim SpliceAI "carries the signal" on
# variants where it was never assessed. Templates must ask whether it actually scored.
from sugivariant.enrich import spliceai_scored as _spliceai_scored  # noqa: E402
env.globals["spliceai_scored"] = _spliceai_scored
# Jinja's |capitalize lowercases the REST of the string, which turned "no allele count
# reported in gnomAD v4.1 (AN 1,208,295)" into "...gnomad v4.1 (an ...)". Same trap as
# str.capitalize() lowercasing "SpliceAI". Upper-case the first character only.
from sugivariant.render import _ucfirst as _uc1  # noqa: E402
env.filters["ucfirst"] = _uc1

# ClinGen dosage: 0-3 is the evidence scale; 30/40 are category codes. One definition,
# in render.py, so the HTML and the .md twin cannot drift apart (they did: the .md showed
# "haploinsufficiency 30 ... (0-3 scale)" on 459,525 pages).
from sugivariant.render import hi_display as _hi_display  # noqa: E402
env.globals["hi_display"] = _hi_display

# Inline markdown-bold in generated verdict text: **x** → <strong>x</strong>.
# `mdbold` returns safe HTML (rest escaped); `mdstrip` yields plain text (for
# headings / meta where we just want the words, no markers).
import markupsafe as _ms                                  # noqa: E402
_MD_BOLD = re.compile(r"\*\*(.+?)\*\*")
env.filters["mdbold"] = lambda s: _ms.Markup(
    _MD_BOLD.sub(r"<strong>\1</strong>", str(_ms.escape(s)))) if s else ""
env.filters["mdstrip"] = lambda s: _MD_BOLD.sub(r"\1", s or "")

app = FastAPI(title="Sugi Variant")
# Serve at ROOT: when deployed behind a reverse proxy under a subpath, the proxy
# strips the prefix, so routes + static are unprefixed. BASE_PATH is used ONLY in
# templates ({{ base }}) to generate the public prefixed links.
app.mount("/static", StaticFiles(directory=str(ROOT / "static")), name="static")


# ── preprint (Method page) ─────────────────────────────────────────────────────
# The /method page inlines the make4ht HTML export of the preprint. Build flow:
# in the preprint repo run `make html`, then copy main.html / main.pdf / main*.svg
# into static/preprint/ and main.css into static/preprint/main.raw.css, and scope the
# CSS under .preprint-doc. Here we lift the <body> inner HTML once at import and
# rewrite the relative SVG figure refs to the served static path. The make4ht global
# CSS is NOT loaded site-wide — the scoped copy is linked only on the method page.
def _preprint_body():
    html = (ROOT / "static" / "preprint" / "main.html").read_text(encoding="utf-8")
    body = html.split("<body>", 1)[1].split("</body>", 1)[0]
    # main0x/1x/….svg -> {base}/static/preprint/…  (keeps make4ht's single-quote form)
    body = re.sub(r"src='(main[0-9]+x\.svg)'", rf"src='{BASE}/static/preprint/\1'", body)
    return body


try:
    PREPRINT_BODY = _preprint_body()
except FileNotFoundError:
    PREPRINT_BODY = ""   # export not built yet; the template shows a short fallback

# ── record resolution (in-memory + persistent per-gene record cache) ───────────
# A gene's enriched records are expensive to build (a whale like BRCA2 ≈ minutes),
# so once built they're persisted to disk (gzip pickle) and loaded on later hits —
# this survives restarts and stops the cold-build Cloudflare timeout.
#
# The in-memory cache is a byte-budgeted LRU so a crawler walking all ~18.8k genes cannot grow
# it without bound (the full corpus would be ~13 GB in RAM). Evicted genes stay on disk and
# reload in ms; the disk cache is the durable, unbounded backing store. Tune with GENE_CACHE_MB.
import collections, threading                          # noqa: E402
_GENE_CACHE = collections.OrderedDict()                # gene -> records, LRU order (newest last)
_GENE_CACHE_BYTES = 0
_GENE_CACHE_BUDGET = int(os.environ.get("GENE_CACHE_MB", "1500")) * 1024 * 1024
_GENE_CACHE_LOCK = threading.Lock()
# Derived similar-variant indexes, one per loaded gene. Tied to the gene LRU:
# evicted with its records and rebuilt on the next page. Never written to disk.
_SIMILAR_INDEX = {}
_SIMILAR_LOCK = threading.Lock()
# Read-only release mode (VARIANT_SERVE=release): {GENE: {path, status, stale}}
# and the snapshot id used for the X-Sugi-Variant-Release header. Loaded once at startup.
_RELEASE = None
_RELEASE_STATE = {}
_RELEASE_ID = None
_RELEASE_DATE = None
_RELEASE_VERSION = None
_RELEASE_META = {}
# Measured deserialized-RAM cost per enriched record (~15 KB by RSS delta); rounded up so the
# byte budget is a safe upper bound on real memory (GENE_CACHE_MB ~= actual MB the cache uses).
_BYTES_PER_VARIANT = 16000
_RECORDS_DIR = CACHE_DIR / "records"
_SLUG_SEP = re.compile(r"-[pcnm]-")   # p./c./n.(ncRNA)/m.(mito)


def _cache_get(g):
    """Return cached records for gene g (marking it most-recently-used), or None."""
    with _GENE_CACHE_LOCK:
        recs = _GENE_CACHE.get(g)
        if recs is not None:
            _GENE_CACHE.move_to_end(g)
        return recs


def _cache_put(g, recs):
    """Insert gene g's records and evict least-recently-used genes until under the byte budget.
    Evicted genes remain on disk (records/*.pkl.gz) and reload from there on a later hit."""
    global _GENE_CACHE_BYTES
    sz = len(recs) * _BYTES_PER_VARIANT
    with _GENE_CACHE_LOCK:
        if g in _GENE_CACHE:
            _GENE_CACHE_BYTES -= len(_GENE_CACHE[g]) * _BYTES_PER_VARIANT
        _GENE_CACHE[g] = recs
        _GENE_CACHE.move_to_end(g)
        _GENE_CACHE_BYTES += sz
        while _GENE_CACHE_BYTES > _GENE_CACHE_BUDGET and len(_GENE_CACHE) > 1:
            ev_gene, ev = _GENE_CACHE.popitem(last=False)  # drop the least-recently-used gene
            _GENE_CACHE_BYTES -= len(ev) * _BYTES_PER_VARIANT
            with _SIMILAR_LOCK:
                _SIMILAR_INDEX.pop(ev_gene, None)


def _gene_of(slug):
    m = _SLUG_SEP.search(slug)
    return slug[:m.start()] if m else None


import threading                                        # noqa: E402
_GENE_LOCKS, _LOCKS_LOCK = {}, threading.Lock()


def _gene_lock(g):
    with _LOCKS_LOCK:
        return _GENE_LOCKS.setdefault(g, threading.Lock())


def _similar_index_for(gene, recs):
    """The per-gene prep index for similar-variant links, built once and cached.

    Tied to the gene LRU: _cache_put drops it with the gene, and _records drops it
    when a new record list replaces the old one. Never persisted.
    """
    g = gene.upper()
    with _SIMILAR_LOCK:
        idx = _SIMILAR_INDEX.get(g)
    if idx is not None:
        return idx
    with _gene_lock(g):
        with _SIMILAR_LOCK:
            idx = _SIMILAR_INDEX.get(g)
        if idx is None:
            idx = _build_similar_index(recs)
            with _SIMILAR_LOCK:
                _SIMILAR_INDEX[g] = idx
        return idx


def _with_similar(r, gene, recs):
    """Attach the serve-time similar list to a record exactly once per process."""
    if not r.get("_sim_fresh"):
        r["similar"] = _similar_from_index(r, _similar_index_for(gene, recs))
        r["_sim_fresh"] = True
    return r


_RELEASE_SERVABLE = ("ok", "stale", "mismatch")


def _release_servable(gene):
    s = _RELEASE_STATE.get((gene or "").upper())
    return bool(s and s.get("path") and s.get("status") in _RELEASE_SERVABLE)


def _records_release(g):
    """Snapshot mode: load the pickle named by the manifest, or return nothing.

    No biotree, no enriched_records(), no collect(), no persist. A gene the
    manifest does not list (or lists as missing/known-empty/orphan) is simply
    not in this snapshot.
    """
    recs = _cache_get(g)
    if recs is not None:
        return recs
    s = _RELEASE_STATE.get(g)
    if not s or s.get("status") not in _RELEASE_SERVABLE or not s.get("path"):
        return []
    import gzip
    import pickle
    with _gene_lock(g):
        recs = _cache_get(g)
        if recs is not None:
            return recs
        try:
            with gzip.open(_RECORDS_DIR / s["path"], "rb") as f:
                recs = pickle.load(f)
            with _SIMILAR_LOCK:
                _SIMILAR_INDEX.pop(g, None)
            _cache_put(g, recs)
            return recs
        except Exception:
            return []


def _records(gene):
    g = gene.upper()
    recs = _cache_get(g)
    if recs is not None:
        return recs
    if RELEASE_MODE:
        return _records_release(g)
    import gzip
    import pickle
    path = _RECORDS_DIR / f"{g}.pkl.gz"
    # per-gene lock: only ONE thread builds a given gene; concurrent hits wait and
    # then get the cached result (no duplicate multi-minute whale builds).
    with _gene_lock(g):
        recs = _cache_get(g)
        if recs is not None:
            return recs
        if path.exists():
            try:
                with gzip.open(path, "rb") as f:
                    recs = pickle.load(f)
                with _SIMILAR_LOCK:
                    _SIMILAR_INDEX.pop(g, None)   # new list: drop any old prep index
                _cache_put(g, recs)          # LRU insert (may evict cold genes; they stay on disk)
                return recs
            except Exception:
                pass
        _begin_build()
        recs = enriched_records(g) or []
        with _SIMILAR_LOCK:
            _SIMILAR_INDEX.pop(g, None)       # rebuilt list: old index is stale
        if _build_degraded():
            # An upstream fetch FAILED during this build, so the records under-report
            # evidence (a missing gnomAD read used to render as "absent from gnomAD" —
            # an ACMG PM2-supporting claim — and a missing predictor silently drops the
            # variant out of the disagreement worklist). Serve it once, but never cache
            # or persist it: a transient outage must not freeze into the corpus.
            for r in recs:
                r["_degraded"] = True
            return recs
        _cache_put(g, recs)
        # Don't persist empty/invalid-gene results, nor a set that COLLAPSES against
        # what is already stored: an upstream gap can return one record for a gene with
        # hundreds, and that write is unrecoverable (see records_collapsed).
        collapsed = records_collapsed(path, len(recs)) if recs else None
        if collapsed:
            print(f"refusing to persist {g}: {len(recs)} records would replace "
                  f"{collapsed} stored", file=sys.stderr, flush=True)
        if recs and not collapsed:
            try:
                _RECORDS_DIR.mkdir(parents=True, exist_ok=True)
                # PID-scoped tmp: _gene_lock is per-process but all workers share
                # /cache/records, so a fixed .tmp name can interleave between workers.
                tmp = path.with_suffix(f".{os.getpid()}.tmp")
                with gzip.open(tmp, "wb") as f:
                    pickle.dump(recs, f)
                tmp.replace(path)    # atomic
                # Stamp which upstream DATA versions these records were built from, so
                # the next release can rebuild only the genes it affects. Written after
                # the records, never before: a stamp without its records would claim a
                # gene is current when it is not.
                from sugivariant.enrich import write_data_stamp
                write_data_stamp(path)
            except Exception:
                pass
    return recs


def _collect_one(vcv):
    """Build a single variant record straight from its VCV (index hit) — no per-gene
    build. Used as a resolution fallback and for the whale fast path."""
    if RELEASE_MODE:
        return None                     # snapshot serving never builds
    from sugivariant.collect import attach_enrichment
    try:
        _begin_build()
        rec = collect(str(vcv))
        if rec:
            attach_enrichment(rec, {})              # single-variant enrichment (no per-gene caches)
            rec["_partial"] = True
            if _build_degraded():
                rec["_degraded"] = True
            return rec
    except Exception:
        pass
    return None


def _resolve(slug):
    gene = _gene_of(slug)
    recs = _records(gene) if gene else []
    for r in recs:
        if r["canonical_slug"] == slug or slug in (r.get("slugs") or []):
            return _with_similar(r, gene, recs)
    # Not in the gene's own records — resolve via the index alias table. This is also
    # how a coordinate-shaped URL resolves at all: `21:34880635:g:a` has no gene prefix
    # to slice, so _gene_of() returns None and the index row is the only thing that
    # knows which gene pickle serves it.
    ix = _index()
    if ix:
        hits = IX.lookup(ix, slug)
        if hits:
            vcv = str(hits[0]["vcv"])
            if not any(str(r.get("variation_id")) == vcv for r in recs):
                # The prefix _gene_of() sliced is not the gene that serves this row
                # (a coordinate slug like 21:34880635:g:a slices to "21"). Trust the
                # index row's gene and load that pickle instead.
                hit_gene = (hits[0].get("gene") or "").lower()
                if hit_gene and hit_gene != gene:
                    gene = hit_gene
                    recs = _records(gene)
            for r in recs:                          # legacy/alias-only slug (e.g. pre-cap delins URL)
                if str(r.get("variation_id")) == vcv:
                    return _with_similar(r, gene, recs)
            # The VCV isn't under this gene at all — ClinVar's sort-gene differs from the
            # HGVS-derived gene (overlapping loci, e.g. MT-ATP8 filed under MT-ATP6). The
            # per-gene build will never contain it, so build the single variant directly.
            if RELEASE_MODE:
                return None             # snapshot serving never builds
            return _collect_one(vcv)
    return None


_FAST_PATH_MIN = 1500   # genes bigger than this get a single-variant fast build on a cold hit


def _resolve_page(slug):
    """Resolve a variant for its page WITHOUT a multi-minute whale build on the request
    path. If the gene is already cached (memory/disk) or small, do the full build (rich
    page). If it's a big gene not yet built, build just THIS variant (~1-2s), serve it,
    and warm the full gene in the background so the next hit is complete."""
    gene = _gene_of(slug)
    if not gene:
        # A coordinate-shaped URL carries no gene prefix. The index row still knows
        # which pickle serves it, so let _resolve() take it from there.
        return _resolve(slug) if _index() else None
    if RELEASE_MODE:
        return _resolve(slug)           # read the pickle or return None; no build
    g = gene.upper()
    if g in _GENE_CACHE or (_RECORDS_DIR / f"{g}.pkl.gz").exists():
        return _resolve(slug)
    ix = _index()
    if not ix or IX.gene_count(ix, g) <= _FAST_PATH_MIN:
        return _resolve(slug)                       # small/medium gene → full build is fast
    hits = IX.lookup(ix, slug)
    if not hits:
        return _resolve(slug)
    rec = _collect_one(hits[0]["vcv"])
    if not rec:
        return _resolve(slug)
    threading.Thread(target=lambda: _records(g), daemon=True).start()   # warm the full gene
    return rec


def release_line_md():
    """`Release v1.5.0 · BioBTree v2.20.0 · last updated 2026-11-07` for the .md twins.

    One place so the HTML footer, the variant .md and the SV .md cannot drift. The
    biotree version is the annotation release the corpus was built against, so a
    reader can tell what the evidence dates are as well as when we published.
    """
    if not (RELEASE_MODE and _RELEASE_VERSION):
        return ""
    parts = [f"Release {_RELEASE_VERSION}"]
    bb = (_RELEASE_META or {}).get("biobtree_version")
    if bb:
        parts.append(f"BioBTree {bb}")
    if _RELEASE_DATE:
        parts.append(f"last updated {_RELEASE_DATE}")
    return "_" + " · ".join(parts) + "_"


def _render(tpl, **ctx):
    """Render a template, deriving the variant page's presentation states here rather
    than at each call site.

    `ui` holds every semantic state the page styles on (see sugivariant/present.py).
    Computing it here means a caller cannot render the page without it and get an
    unstyled or, worse, a differently-styled page — the same footgun that let
    index.py classify a flag from a stale readout. It is pure and adds no lookups.
    """
    if tpl == "variant.html" and "ui" not in ctx and ctx.get("v") is not None:
        ctx["ui"] = evidence_view(ctx["v"], ctx.get("freq_g"))
        ctx["ui"]["groups"] = evidence_groups(ctx["v"])
        ctx["ui"]["transcripts"] = transcript_view(ctx["v"])
    return env.get_template(tpl).render(**ctx)


# ── persistent resolution index (optional; falls back to the live path if absent) ─
INDEX_DB = pathlib.Path(os.environ.get("INDEX_DB") or (CACHE_DIR / "index.db"))
SV_DB = pathlib.Path(os.environ.get("SV_DB") or (CACHE_DIR / "sv.db"))
_tls = threading.local()


def _load_release():
    """Load the compact (or full) release manifest for VARIANT_SERVE=release."""
    global _RELEASE, _RELEASE_STATE, _RELEASE_ID, _RELEASE_DATE, _RELEASE_VERSION, _RELEASE_META, SV_DB
    import json
    if not RELEASE_PATH.exists():
        raise RuntimeError(f"VARIANT_SERVE=release but no release manifest at {RELEASE_PATH}")
    if not INDEX_DB.exists():
        raise RuntimeError(f"VARIANT_SERVE=release but no index at {INDEX_DB}")
    data = json.loads(RELEASE_PATH.read_text(encoding="utf-8"))
    corpus = ((data.get("corpora") or {}).get("small_variant") or {})
    raw = corpus.get("genes")
    if raw is None:                       # full manifest from `snapshot.py build`
        raw = corpus.get("gene_state") or {}
    state = {}
    for g, s in raw.items():
        if not isinstance(s, dict):
            continue
        state[g.upper()] = {
            "path": s.get("path"),
            "status": s.get("status"),
            "stale": bool(s.get("stale", "stale" in (s.get("issues") or []))),
        }
    _RELEASE = data
    _RELEASE_STATE = state
    _RELEASE_ID = data.get("release") or data.get("snapshot_id") or data.get("created_at")
    _RELEASE_DATE = data.get("released_at") or (data.get("created_at") or "")[:10] or None
    _RELEASE_VERSION = data.get("variant_version") or _RELEASE_ID
    _RELEASE_META = data
    sv = ((data.get("corpora") or {}).get("structural_variant") or {})
    if sv.get("db"):
        # The manifest records the path the build saw, which is the HOST path. Inside
        # the prod container the same cache is mounted at CACHE_DIR, so a missing
        # absolute path is resolved by name against this process's cache directory.
        # Without this, release mode in a container silently loses the SV corpus.
        db = pathlib.Path(sv["db"])
        if not db.exists():
            alt = CACHE_DIR / db.name
            if alt.exists():
                db = alt
        SV_DB = db
    env.globals["release"] = _RELEASE_VERSION
    env.globals["released_at"] = _RELEASE_DATE
    env.globals["biobtree_version"] = (_RELEASE_META or {}).get("biobtree_version")
    rd = corpus.get("records_dir")
    if rd:
        try:
            if pathlib.Path(rd).resolve() != _RECORDS_DIR.resolve():
                print(f"warning: release records_dir {rd} != {_RECORDS_DIR}",
                      file=sys.stderr, flush=True)
        except Exception:
            pass
    print(f"Sugi Variant: release {_RELEASE_VERSION} ({_RELEASE_DATE}); "
          f"{len(state):,} genes; biobtree {data.get('biobtree_version') or '?'}",
          flush=True)


if RELEASE_MODE:
    _load_release()



def _index():
    """A thread-local, read-only connection to the SQLite resolution index, or None if it
    hasn't been built. Per-thread connections let the request threads read the index in true
    parallel — SQLite allows unlimited concurrent readers across separate connections — instead
    of serialising on one shared connection. Read-only (mode=ro); serving never writes."""
    if not INDEX_DB.exists():
        return None
    conn = getattr(_tls, "index", None)
    if conn is None:
        conn = IX.open_ro(str(INDEX_DB))
        _tls.index = conn
    return conn


def _sv_index():
    """Thread-local read-only connection to the dbVar structural-variant store, or None."""
    if not SV_DB.exists():
        return None
    conn = getattr(_tls, "sv", None)
    if conn is None:
        conn = SV.open_ro(str(SV_DB))
        _tls.sv = conn
    return conn


def _sv_row(row):
    r = dict(row)
    for k in ("phenotypes", "clinvar"):
        try:
            r[k] = json.loads(r.get(k) or "[]")
        except Exception:
            r[k] = []
    return r


def _sv_lookup(key):
    conn = _sv_index()
    if not conn:
        return []
    rows = [_sv_row(r) for r in SV.lookup(conn, key)]
    if rows:
        return rows
    # A locus request (chr:start-end) that is not itself an SV coordinate: return
    # everything overlapping it. Exact SV coordinates are handled by the alias above.
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", (key or "").strip())
    if m:
        return [_sv_row(r) for r in SV.overlapping(
            conn, m.group(1), int(m.group(2)), int(m.group(3)))]
    return []


def _sv_genes(nssv):
    conn = _sv_index()
    return SV.genes_for_sv(conn, nssv) if conn else []


def _sv_size(v):
    try:
        return f"{abs(int(v.get('svlen'))):,} bp"
    except (TypeError, ValueError):
        return None


def _clinvar_url(acc):
    a = (acc or "").strip()
    if a.upper().startswith("VCV"):
        num = a[3:].split(".")[0].lstrip("0") or "0"
        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{num}/"
    if a.upper().startswith("RCV"):
        return f"https://www.ncbi.nlm.nih.gov/clinvar/{a}/"
    return None


env.globals["sv_size"] = _sv_size
env.globals["clinvar_url"] = _clinvar_url


def _ix_rec(row):
    """Adapt an index row to the record shape templates/redirects expect."""
    return {"canonical_slug": row["slug"], "gene_symbol": row["gene"],
            "hgvs_p": row["hgvs_p"], "hgvs_c": row["hgvs_c"],
            "classification": row["classification"], "review_status": row["review_status"],
            "rsid": row["rsid"], "coordinate": row["coordinate"],
            "variation_id": row["vcv"]}


def _index_hits(key):
    """[adapted records] for a normalized key via the index, or None if no index."""
    ix = _index()
    return [_ix_rec(r) for r in IX.lookup(ix, key)] if ix else None


def _vcv_hits(key):
    """Small variants by ClinVar VCV accession (VCV000007829 / VCV000007829.1)."""
    m = re.match(r"^VCV0*(\d+)(?:\.\d+)?$", (key or "").strip(), re.I)
    ix = _index()
    if not m or not ix:
        return []
    row = ix.execute("SELECT * FROM variant WHERE vcv=?", (int(m.group(1)),)).fetchone()
    return [_ix_rec(row)] if row else []


# ── query resolution (rsID reverse-map + gene-first HGVS match) ─────────────────
# Three arrival keys are supported (biobtree keys on identifiers, not free-text
# HGVS): an rsID resolves directly via dbSNP reverse-map; a GENE + HGVS/AA query
# resolves gene-first (enumerate + slug-match, all on our side). See the biobtree
# key-scheme note — raw HGVS is not searchable, so we normalize it to OUR slug.
_RSID_RE = re.compile(r"^rs\d+$", re.I)
_GENE_QUERY_RE = re.compile(r"^[A-Za-z][A-Za-z0-9-]{1,14}$")


def _hgvs_keys(rec):
    """Normalized HGVS match-keys for a record: the p./c. forms (with and without
    the leading operator, since the demand comes bare — `acta1 "pro309ala"`) and
    the 1-letter short (`R130Q`). All run through the same _norm_hgvs the slugs use."""
    keys = set()
    for form in (rec.get("hgvs_p"), rec.get("hgvs_c")):
        if form:
            keys.add(_norm_hgvs(form))
            keys.add(_norm_hgvs(re.sub(r"^[pc]\.", "", form)))
    short = missense_short(rec.get("hgvs_p"))
    if short:
        keys.add(_norm_hgvs(short))
    return keys


def resolve_rsid(rsid):
    """rsID → the pathogenic-gated ClinVar page records it maps to (0, 1 or many —
    an rsID is position-level and can cover several variations). Index first (no
    biobtree call), else the live dbSNP reverse-map."""
    hits = _index_hits(rsid.lower())
    if hits:
        return hits
    if RELEASE_MODE:
        return []           # release serving is read-only; no live dbSNP fallback
    try:
        variations = map_all(rsid, ">>dbsnp>>clinvar", cap=None)
    except Exception:
        return []
    out, seen = [], set()
    for v in variations:
        vid = v.get("id")
        if not vid or vid in seen:
            continue
        seen.add(vid)
        rec = collect(vid)          # gates should_build + attaches the slug
        if rec:
            out.append(rec)
    return out


def resolve_query(q):
    """(kind, term, [records]) for a free-text search: an rsID, or a GENE + HGVS/AA
    query resolved gene-first. Records are light (no enrichment) — enough to
    redirect or disambiguate; the page itself builds on the redirect."""
    q = (q or "").strip()
    if not q:
        return ("empty", q, [])
    if _RSID_RE.match(q):
        return ("rsid", q.lower(), resolve_rsid(q.lower()))

    # A bare gene symbol is a real query, not an unparseable one: send it to the
    # gene hub instead of "not found". It answers "which variant do you mean" with
    # the page that lists them.
    if len([t for t in re.split(r"[\s,]+", q) if t]) == 1 and _GENE_QUERY_RE.match(q):
        ix = _index()
        g = q.upper()
        n = IX.gene_count(ix, g) if ix else 0
        if n:
            return ("gene", g, [{"kind": "gene", "gene": g, "variants": n}])

    coord = IX.norm_coordinate(q)
    if coord:
        # A coordinate has NO biobtree→ClinVar edge, so only the index can resolve
        # it — no live fallback is possible.
        return ("coordinate", q, _index_hits(coord) or [])

    # GENE + HGVS/AA. Splitting bare tokens (ACTA1 vs pro309ala, and genes like
    # SLC2A1 that also carry digit-letter runs) can't be classified reliably, so
    # don't guess — try each non-HGVS token as the gene and keep the one that hits.
    toks = [t for t in re.split(r"[\s,]+", q) if t]
    gene_candidates = [t for t in toks if "." not in t and ">" not in t]
    if len(toks) < 2 or not gene_candidates:
        return ("unparsed", q, [])
    # Index-first: build the normalized key per gene candidate (no gene build).
    for g in gene_candidates:
        key = f"{_norm(g)}-{_norm_hgvs(' '.join(t for t in toks if t != g))}"
        hits = _index_hits(key)
        if hits:
            return ("query", q, hits)
    # Live fallback (gene-first build + slug-match) for un-indexed genes.
    tried_gene = False
    for g in gene_candidates:
        recs = _records(g)
        if not recs:
            continue
        tried_gene = True
        qkey = _norm_hgvs(" ".join(t for t in toks if t != g))
        hits = [r for r in recs if qkey in _hgvs_keys(r)]
        if hits:
            return ("query", q, hits)
    return ("query" if (tried_gene or _index()) else "nogene", q, [])


def _resolution_response(kind, term, hits):
    """One hit → redirect to its page; else render the results/disambiguation page."""
    if len(hits) == 1:
        if hits[0].get("kind") == "gene":
            return RedirectResponse(f"{BASE}/gene/{hits[0]['gene']}", status_code=307)
        return RedirectResponse(f"{BASE}/{hits[0]['canonical_slug']}", status_code=307)
    return HTMLResponse(_render("results.html", kind=kind, term=term, hits=hits))


# ── disagreement browse view ────────────────────────────────────────────────────
# The one signal no competitor surfaces (benchmark 2026-07): variants where the
# evidence doesn't line up. Descriptive QC, not reclassification (variant.sh §8).
FEATURED_GENES = ["BRCA1", "BRCA2", "TP53", "PTEN", "MLH1", "MSH2",
                  "LDLR", "SCN1A", "KCNQ1", "ASXL1"]
CAT_META = {
    "predictor_vs_clinvar": {
        "label": "Predictors vs ClinVar", "tone": "flag",
        "blurb": "ClinVar calls these pathogenic, but the computational predictors "
                 "lean tolerated — the highest-value review flag."},
    "resolves_conflicting": {
        "label": "Predictors agree, ClinVar conflicts", "tone": "info",
        "blurb": "ClinVar submitters conflict, but the independent predictors "
                 "unanimously agree — a QC review flag, not a reclassification."},
    "lof_resolves_conflicting": {
        "label": "Predicted LoF vs a conflicting call", "tone": "info",
        "blurb": "ClinVar submitters conflict, and it's a predicted loss-of-function "
                 "change in a gene where LoF causes disease — flagged for review, not "
                 "resolved (verify NMD-escape for C-terminal truncations). Non-missense."},
    "vus_predictors_lean": {
        "label": "Predictors lean, ClinVar uncertain", "tone": "info",
        "blurb": "ClinVar classifies these uncertain (VUS), but the independent predictors "
                 "unanimously lean one way — a triage signal for review, not a reclassification."},
    "predictors_split": {
        "label": "Predictors split", "tone": "warn",
        "blurb": "The independent predictors disagree with each other — read with care."},
}


def _grouped_flags(recs):
    """Group a gene's flagged records by disagreement category, most-relevant first
    within each group (pathogenic → higher review stars → slug)."""
    groups = {cat: [] for cat in DISAGREEMENT_CATEGORIES}
    for r in recs:
        _refresh(r)                   # pickled readouts may predate a fix to them
        f = disagreement_flag(r)
        if f:
            groups[f["category"]].append({"rec": r, "flag": f})
    for items in groups.values():
        items.sort(key=lambda it: (-review_stars(it["rec"].get("review_status")),
                                   it["rec"]["canonical_slug"]))
    return groups


# ── ETag html-cache middleware (mirrors Sugi Predict: dynamic but CDN-cacheable) ─
@app.middleware("http")
async def _cache_html(request, call_next):
    resp = await call_next(request)
    if request.method == "GET" and resp.status_code == 200 \
            and resp.headers.get("content-type", "").startswith("text/html"):
        body = b"".join([chunk async for chunk in resp.body_iterator])
        etag = '"' + hashlib.md5(body).hexdigest() + '"'
        h = dict(resp.headers)
        for k in ("content-length", "content-type"):
            h.pop(k, None)
        h["ETag"] = etag
        # CDN-cacheable: browsers revalidate via ETag (max-age=0), but Cloudflare
        # caches for a day (s-maxage) and serves stale while revalidating — so a page
        # built once is served fast from the edge instead of rebuilding the gene.
        if "no-store" in (resp.headers.get("cache-control") or ""):
            h["Cache-Control"] = "no-store"          # degraded/incomplete render — never cache it
        else:
            h["Cache-Control"] = "public, max-age=0, s-maxage=86400, stale-while-revalidate=604800, stale-if-error=86400"
        if request.headers.get("if-none-match") == etag:
            return Response(status_code=304, headers=h)
        return HTMLResponse(content=body, headers=h)
    return resp


@app.middleware("http")
async def _release_header(request, call_next):
    """Name the snapshot on every response, not just pages, so a deploy is auditable."""
    resp = await call_next(request)
    if RELEASE_MODE and _RELEASE_ID:
        resp.headers["X-Sugi-Variant-Release"] = str(_RELEASE_ID)
        if _RELEASE_VERSION:
            resp.headers["X-Sugi-Variant-Version"] = str(_RELEASE_VERSION)
        if _RELEASE_DATE:
            resp.headers["X-Sugi-Variant-Released-At"] = str(_RELEASE_DATE)
        bb = (_RELEASE_META or {}).get("biobtree_version")
        if bb:
            resp.headers["X-Sugi-BioBTree-Version"] = str(bb)
    return resp


@app.exception_handler(StarletteHTTPException)
async def _err(request, exc):
    # A 405 here can only ever be a non-GET method on a page route: every route in this
    # app is @app.get, there is no form, no <input> and no POST endpoint on any page,
    # and the JS does no fetch, no pushState and registers no service worker. So a 405
    # is never a real "you may not do that" -- it is a URL that renders perfectly well
    # over GET, reached with the wrong verb, and showing a styled error page for it is
    # the wrong answer.
    #
    # Observed 2026-10-15: mobile Firefox POSTed /variant/jak3-c-2351-minus-15a-g on a
    # browser refresh, with that same page as the referer, and got "Method Not Allowed".
    # A reload re-issues the method of the history entry it is reloading, so that entry
    # had to have been CREATED by a POST, and nothing we serve can create one. What can
    # is the Cloudflare managed challenge, which is active on /variant/ and /predict/
    # (cf-mitigated: challenge) and is served at the requested URL itself, so clearing
    # it can leave a POST entry in history for a page URL. Two such POSTs exist in the
    # whole access log, so this is rare -- but it costs four lines to never show it.
    #
    # 303 is the specific code that tells a browser to redo the request as GET (302
    # may preserve the method). The no-store is belt and braces and does NOT reach the
    # edge: the sugi.bio nginx vhost sets `proxy_hide_header Cache-Control` on
    # /variant/ and stamps its own "no-cache, stale-if-error=86400" on every response,
    # deliberately, so a deploy is visible without a hard refresh. Cloudflare therefore
    # never sees this header. It is sent anyway because it is the correct thing for the
    # app to say, and because Cloudflare caches GET and HEAD only, so a 303 answering a
    # POST is not a cacheable response in the first place.
    if exc.status_code == 405 and request.method not in ("GET", "HEAD"):
        return RedirectResponse(str(request.url), status_code=303,
                                headers={"Cache-Control": "no-store"})
    return HTMLResponse(_render("error.html", code=exc.status_code,
                                detail=exc.detail), status_code=exc.status_code)


@app.exception_handler(BiobtreeError)
async def _upstream_down(request, exc):
    """Upstream (BioBTree) is unreachable or erroring after its retries. That is a
    503 with a Retry-After, not a 500 — the request is fine, the dependency is not.
    Explicitly no-store so a transient outage is never cached as a page."""
    return HTMLResponse(
        _render("error.html", code=503,
                detail="A data source is temporarily unavailable. "
                       "This page builds live from the reference databases — "
                       "please retry in a moment."),
        status_code=503,
        headers={"Retry-After": "30", "Cache-Control": "no-store"})


# ── routes ─────────────────────────────────────────────────────────────────────
_STATS_CACHE = None


def _corpus_stats():
    """Memoized on the index mtime. Previously a keyless permanent memo: if index.db was
    absent at the first request (deploy race, volume not yet mounted) it cached zeros
    FOREVER and the home page advertised "0 variants" until a restart — and it never
    picked up a reindex. The container HEALTHCHECK curls /, so that first request is
    guaranteed to be the risky one."""
    global _STATS_CACHE
    stamp = _index_lastmod()
    if _STATS_CACHE is None or _STATS_CACHE[0] != stamp:
        ix = _index()
        if not ix:
            return {"variants": 0, "genes": 0, "flagged": 0, "pathogenic": 0, "vus": 0}
        _STATS_CACHE = (stamp, IX.corpus_stats(ix))
    return _STATS_CACHE[1]


def _freq_ctx(rec):
    """(normalized gnomad dict, pops) for the population-frequency block.

    Every key is present: records on disk predate several fields, and in Jinja a missing
    key is Undefined where `Undefined is not none` is True — the trap that 500'd
    production (b598202). `pops` is (ancestry, af) sorted high→low."""
    g = {k: None for k in ("popmax", "ancestry", "af", "ac", "an", "ac_grpmax", "an_grpmax",
                           "faf", "faf95", "faf99", "faf95_exomes", "faf99_exomes",
                           "frequency", "band", "absent", "unknown", "is_common", "source",
                           "populations")}
    # Correct a stored block that claims an absence it cannot support, before anything
    # reads `absent` or `band` (pure; both bad shapes are visible in the stored dict).
    _refresh(rec)
    g.update(rec.get("gnomad") or {})
    if g.get("faf95") is None:
        g["faf95"] = g.get("faf")
    if g.get("af") is None and g.get("frequency") not in (None, ""):
        # dbSNP-fallback records cached before gnomad_for() started emitting `af` carry only
        # the raw `frequency` string. Without this they fail the "has a frequency" gate and
        # the whole block disappears — 942 of 56,544 sampled records.
        try:
            g["af"] = float(g["frequency"])
        except (TypeError, ValueError):
            pass
    g["band"] = g.get("band") or ""
    pops = sorted(((k, float(x)) for k, x in (g.get("populations") or {}).items() if x),
                  key=lambda kv: -kv[1])
    return g, pops


env.globals["anc_name"] = anc_name
env.globals["one_in"] = one_in
env.globals["grpmax_excluded"] = GRPMAX_EXCLUDED


_SET_SPLIT = re.compile(r"[\n;,]+")


@app.get("/")
def home(q: str = ""):
    q = q.strip()
    if q:
        # a comma / newline / semicolon list of 2+ items → a variant set
        parts = [p for p in (s.strip() for s in _SET_SPLIT.split(q)) if p]
        if len(parts) > 1:
            from urllib.parse import quote
            return RedirectResponse(f"{BASE}/set?q={quote(q)}", status_code=307)
        kind, term, hits = resolve_query(q)
        if not hits and _sv_lookup(q):
            return RedirectResponse(f"{BASE}/sv/{q}", status_code=307)
        return _resolution_response(kind, term, hits)
    s = _corpus_stats()
    return HTMLResponse(_render("home.html", n_variants=s["variants"],
                                n_genes=s["genes"], n_flagged=s["flagged"],
                                n_pathogenic=s.get("pathogenic"), n_vus=s.get("vus"),
                                nav="home"))


_SET_MAX = 60
_FLAG_SEVERITY = {"predictor_vs_clinvar": 3, "resolves_conflicting": 2,
                  "lof_resolves_conflicting": 2, "vus_predictors_lean": 2,
                  "predictors_split": 1}


@app.get("/set", response_class=HTMLResponse)
def variant_set(v: str = "", q: str = ""):
    ix = _index()
    # paste box → resolve each line via the shared resolver → canonical permalink
    if q.strip():
        slugs = []
        # Cap the INPUT, not just the output. resolve_query() can fall through to a live
        # per-gene build, so resolving every pasted line let one unauthenticated GET
        # trigger an unbounded number of them (_SET_MAX used to be applied only to the
        # resulting slug list, below).
        lines = [ln.strip() for ln in _SET_SPLIT.split(q) if ln.strip()][:_SET_MAX]
        for line in lines:
            for h in resolve_query(line)[2]:
                if h["canonical_slug"] not in slugs:
                    slugs.append(h["canonical_slug"])
        if slugs:
            return RedirectResponse(f"{BASE}/set?v=" + ",".join(slugs[:_SET_MAX]), status_code=307)
        miss = lines
        return HTMLResponse(_render("set.html", rows=[], missing=miss, stats=None,
                                    worklist=[], gene_groups=[], meta=CAT_META, nav="set", capped=False))
    reqs = [s.strip() for s in v.split(",") if s.strip()]
    capped = len(reqs) > _SET_MAX
    reqs = reqs[:_SET_MAX]
    if not ix or not reqs:
        return HTMLResponse(_render("set.html", rows=[], missing=[], stats=None,
                                    worklist=[], gene_groups=[], meta=CAT_META, nav="set", capped=False))
    found, missing = IX.set_rows(ix, reqs)
    import collections
    flagged = [r for r in found if r["flag"]]
    genes = collections.Counter(r["gene"] for r in found)
    conds = {r["primary_condition"] for r in found if r["primary_condition"]}
    sev = lambda r: _FLAG_SEVERITY.get(r["flag"], 0)
    worklist = sorted(flagged, key=lambda r: (-sev(r), r["gene"], r["slug"]))
    rows = sorted(found, key=lambda r: (0 if r["flag"] else 1, -sev(r), r["gene"], r["slug"]))
    gene_groups = [(g, [r for r in rows if r["gene"] == g]) for g, cnt in genes.most_common() if cnt > 1]
    stats = {"n": len(found), "genes": len(genes), "conditions": len(conds),
             "flagged": len(flagged), "high_sev": sum(1 for r in flagged if sev(r) >= 3),
             "cls": dict(collections.Counter(r["classification"] for r in found).most_common()),
             "tiers": dict(collections.Counter(tier_label(r["stars"]) for r in found).most_common()),
             "flag_ct": dict(collections.Counter(r["flag"] for r in flagged))}
    return HTMLResponse(_render("set.html", rows=rows, missing=missing, stats=stats,
                                worklist=worklist, gene_groups=gene_groups, meta=CAT_META,
                                severity=_FLAG_SEVERITY, nav="set", capped=capped))


@app.get("/suggest.json")
def suggest(q: str = ""):
    from fastapi.responses import JSONResponse
    ix = _index()
    items = IX.suggest(ix, q) if ix else []
    return JSONResponse(items, headers={"Cache-Control": "public, max-age=60"})


@app.get("/about", response_class=HTMLResponse)
def about():
    return _render("about.html", nav="about")


@app.get("/method", response_class=HTMLResponse)
def method():
    return _render("method.html", nav="method", preprint_body=PREPRINT_BODY)


@app.get("/genes", response_class=HTMLResponse)
def genes_directory():
    ix = _index()
    if not ix:
        raise StarletteHTTPException(503, "Gene directory needs the resolution index.")
    import collections
    groups = collections.OrderedDict()
    for gene, n in IX.directory_genes(ix):
        letter = gene[0].upper() if gene and gene[0].isalpha() else "#"
        groups.setdefault(letter, []).append({"gene": gene, "n": n})
    total = sum(len(v) for v in groups.values())
    return _render("genes.html", groups=groups, total=total, nav="genes")


_PUB = "https://sugi.bio" + (BASE or "/variant")
env.globals["pub"] = _PUB   # absolute public base, for canonical citation URLs


def _hit_small(r):
    gene = r.get("gene_symbol") or r.get("gene")
    slug = r.get("canonical_slug") or r.get("slug")
    return {
        "kind": "variant", "slug": slug, "gene": gene,
        "hgvs_p": r.get("hgvs_p"), "hgvs_c": r.get("hgvs_c"),
        "classification": r.get("classification"),
        "review_status": r.get("review_status"),
        "rsid": r.get("rsid"), "coordinate": r.get("coordinate"),
        "variation_id": r.get("variation_id") or r.get("vcv"),
        "canonical_url": f"{_PUB}/{slug}",
        "markdown_url": f"{_PUB}/{slug}.md",
        "cached": bool(gene and (_RECORDS_DIR / f"{gene.upper()}.pkl.gz").exists()),
    }


def _hit_sv(r):
    return {
        "kind": "structural_variant", "nssv": r.get("nssv"), "nsv": r.get("nsv"),
        "chromosome": r.get("chrom"), "start": r.get("start"), "end": r.get("end"),
        "svtype": r.get("svtype"), "significance": r.get("significance"),
        "canonical_url": f"{_PUB}/sv/{r.get('nssv')}",
        "markdown_url": f"{_PUB}/sv/{r.get('nssv')}.md",
        "cached": True,
    }


@app.get("/resolve.json")
def resolve_json(q: str = "", gene: str = "", slug: str = "", limit: int = 10):
    """Read-only resolver for agents/MCP: identifier -> canonical page(s).

    Small variants and structural variants share the response shape. An rsID,
    coordinate, GENE+HGVS, slug, NSSV/nsv, VCV or locus resolves; multiple hits
    are returned as `ambiguous` rather than silently choosing one.
    """
    limit = max(1, min(int(limit or 10), 50))
    q = (q or "").strip()
    key = (slug or q).strip()
    hits = []
    if key:
        if slug:
            small = (_index_hits(slug) or []) or _vcv_hits(slug)
            hits = [_hit_small(r) for r in small]
            if not hits:
                hits = [_hit_sv(r) for r in _sv_lookup(slug)]
        else:
            sv = _sv_lookup(q)
            if sv:
                hits = [_hit_sv(r) for r in sv]
            else:
                kind, _term, recs = resolve_query(q)
                if kind == "gene" and recs:
                    hits = [{"kind": "gene", "gene": h["gene"], "variants": h["variants"],
                             "canonical_url": f"{_PUB}/gene/{h['gene']}",
                             "markdown_url": None, "cached": True,
                             "note": "gene-level match; add an HGVS, rsID or coordinate "
                                     "for a variant page"}
                            for h in recs]
                else:
                    hits = [_hit_small(r) for r in recs]
                    if not hits:
                        hits = [_hit_small(r)
                                for r in (_vcv_hits(q) or _index_hits(q) or [])]
    if gene:
        g = gene.upper()
        hits = [h for h in hits
                if (h.get("gene") or "").upper() == g
                or (h.get("kind") == "structural_variant" and g in _sv_genes(h["nssv"]))]
    hits = hits[:limit]
    kinds = {h["kind"] for h in hits}
    return JSONResponse(
        {"query": q, "kind": (kinds.pop() if len(kinds) == 1 else ("mixed" if kinds else None)),
         "ambiguous": len(hits) > 1, "hits": hits, "not_found": not hits,
         "release": _RELEASE_VERSION if RELEASE_MODE else None},
        headers={"Cache-Control": "public, max-age=300"})



def _index_lastmod():
    """Build-level lastmod = the index db's mtime (ISO date). Honest: the page was
    (re)built when the index was."""
    try:
        import datetime
        return datetime.date.fromtimestamp(INDEX_DB.stat().st_mtime).isoformat()
    except Exception:
        return None


# Sitemap index memo: enumerating 6.5k genes takes ~5s, so cache the built XML and
# regenerate only when index.db changes (keyed on its lastmod). Per-worker + in-memory:
# first hit after a rebuild is slow, every hit after is instant. Auto-invalidates.
_SITEMAP_INDEX_CACHE = None  # (index_lastmod, sv_lastmod, xml)
_SV_CHUNK = 50000            # sitemaps.org cap per urlset


def _sv_lastmod():
    try:
        import datetime
        return datetime.date.fromtimestamp(SV_DB.stat().st_mtime).isoformat()
    except Exception:
        return None


@app.get("/sitemap.xml")
def sitemap():
    # Sitemap INDEX (not a flat urlset): a 600k-URL urlset breaches the sitemaps.org
    # 50k/50MB cap and Google silently drops it. One child sitemap per gene, plus
    # chunks for the structural-variant corpus.
    global _SITEMAP_INDEX_CACHE
    ix = _index()
    if not ix:
        raise StarletteHTTPException(503, "Sitemap needs the resolution index; run sugivariant.index.")
    lm = _index_lastmod()
    sv_lm = _sv_lastmod()
    if _SITEMAP_INDEX_CACHE and _SITEMAP_INDEX_CACHE[0] == lm and _SITEMAP_INDEX_CACHE[1] == sv_lm:
        xml = _SITEMAP_INDEX_CACHE[2]
    else:
        lm_tag = f"<lastmod>{lm}</lastmod>" if lm else ""
        children = "".join(
            f"<sitemap><loc>{_PUB}/sitemap-{g.lower()}.xml</loc>{lm_tag}</sitemap>"
            for g, _ in IX.sitemap_genes(ix))
        sv_conn = _sv_index()
        if sv_conn:
            n = sv_conn.execute("SELECT COUNT(*) FROM sv").fetchone()[0]
            sv_tag = f"<lastmod>{sv_lm}</lastmod>" if sv_lm else ""
            children += "".join(
                f"<sitemap><loc>{_PUB}/sitemap-sv-{i}.xml</loc>{sv_tag}</sitemap>"
                for i in range(1, (n + _SV_CHUNK - 1) // _SV_CHUNK + 1))
        xml = ('<?xml version="1.0" encoding="UTF-8"?>'
               '<sitemapindex xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
               f'{children}</sitemapindex>')
        _SITEMAP_INDEX_CACHE = (lm, sv_lm, xml)
    return Response(content=xml, media_type="application/xml",
                    headers={"Cache-Control": "public, max-age=3600"})


@app.get("/sitemap-sv-{chunk}.xml")
def sitemap_sv(chunk: int):
    conn = _sv_index()
    if not conn or chunk < 1:
        raise StarletteHTTPException(404, "No structural-variant sitemap.")
    rows = conn.execute("SELECT nssv FROM sv ORDER BY nssv LIMIT ? OFFSET ?",
                        (_SV_CHUNK, (chunk - 1) * _SV_CHUNK)).fetchall()
    if not rows:
        raise StarletteHTTPException(404, "No structural-variant sitemap.")
    lm = _sv_lastmod()
    lm_tag = f"<lastmod>{lm}</lastmod>" if lm else ""
    urls = "".join(f"<url><loc>{_PUB}/sv/{r[0]}</loc>{lm_tag}</url>" for r in rows)
    xml = ('<?xml version="1.0" encoding="UTF-8"?>'
           '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
           f'{urls}</urlset>')
    return Response(content=xml, media_type="application/xml",
                    headers={"Cache-Control": "public, max-age=3600"})


@app.get("/sitemap-{gene}.xml")
def sitemap_gene(gene: str):
    ix = _index()
    if not ix:
        raise StarletteHTTPException(503, "Sitemap needs the resolution index.")
    slugs = IX.sitemap_gene_slugs(ix, gene)
    if not slugs:
        raise StarletteHTTPException(404, f"No sitemap for “{gene}”.")
    lm = _index_lastmod()
    lm_tag = f"<lastmod>{lm}</lastmod>" if lm else ""
    urls = "".join(f"<url><loc>{_PUB}/{s}</loc>{lm_tag}</url>" for s in slugs)
    xml = ('<?xml version="1.0" encoding="UTF-8"?>'
           '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
           f'{urls}</urlset>')
    return Response(content=xml, media_type="application/xml",
                    headers={"Cache-Control": "public, max-age=3600"})


@app.get("/disagreements")
def disagreements_hub():
    # Folded into the per-gene hubs (each gene's QC section); browse via /genes.
    return RedirectResponse(f"{BASE}/genes", status_code=308)


_DZ_CAP = 150   # per-category display cap (best-reviewed first); avoids a 2000-row whale dump


@app.get("/disagreements/{gene}")
def disagreements_gene(gene: str):
    # The per-gene disagreements now live in a section on the gene hub; keep the old
    # URL working by redirecting to that anchor.
    return RedirectResponse(f"{BASE}/gene/{gene.upper().strip('/')}#disagreements",
                            status_code=308)


_CLS_ORDER = ["Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic",
              "Conflicting classifications of pathogenicity"]
_HUB_CAP = 300   # per-classification display cap (sitemap carries the full set)
_SV_CAP = 100    # structural-variant rows on a gene hub (DMD has 1,785)


@app.get("/gene/{gene}", response_class=HTMLResponse)
def gene_hub(gene: str):
    g = gene.upper().strip("/")
    ix = _index()
    rows = IX.gene_rows(ix, g) if ix else None
    if not rows:                                   # not indexed yet → live build fallback
        # The index is authoritative once built. If it already holds a gene_meta row for
        # this symbol saying zero, we have enumerated it and nothing was filed under it —
        # do NOT fall through to a live build. For a phantom symbol (TTN-AS1 claimed 6,410
        # variants it had no rows for) that build is a ~80,000-call enumeration on the
        # request path, with no deadline, triggerable by anyone, and advertised by
        # type-ahead. It also short-circuits the repeated re-enumeration of genes we
        # already know are empty, which previously cost a full search on every cold hit.
        if ix is not None and IX.gene_known_empty(ix, g):
            recs = []
        else:
            recs = _records(g)
        if not recs:
            # A real gene with no ClinVar variants in our corpus (e.g. Sugi Atlas links here
            # for tRNA/Y-chromosome/ORF genes we don't build) → a graceful cross-link back to
            # Atlas instead of a bare 404. Garbage symbols fall through to the generic 404.
            try:
                from sugivariant.build import _hgnc_id
                real_gene = bool(_hgnc_id(g))
            except Exception:
                real_gene = False
            if real_gene:
                return HTMLResponse(_render("error.html", code=404,
                                            detail=f"No variant page for “{g}”.", gene=g),
                                    status_code=404)
            raise StarletteHTTPException(404, f"No variants built for “{g}”.")
        for x in recs:
            _refresh(x)               # pickled readouts may predate a fix to them
        rows = sorted(({"slug": x["canonical_slug"], "hgvs_p": x.get("hgvs_p"),
                        "hgvs_c": x.get("hgvs_c"), "classification": x["classification"],
                        "stars": review_stars(x.get("review_status")),
                        "primary_condition": (x.get("conditions") or [{}])[0].get("name"),
                        "flag": (disagreement_flag(x) or {}).get("category")} for x in recs),
                      key=lambda r: (-r["stars"], r["slug"]))
    # group by classification, cap the display per group (full set is in the sitemap)
    groups = []
    seen = {r["classification"] for r in rows}
    for cls in _CLS_ORDER + sorted(seen - set(_CLS_ORDER)):
        items = [r for r in rows if r["classification"] == cls]
        if items:
            groups.append({"cls": cls, "total": len(items), "shown": items[:_HUB_CAP]})
    flagged = sum(1 for r in rows if r["flag"])
    # evidence-disagreement (QC) subset, folded into the gene hub (its own section)
    dfull = {cat: [r for r in rows if r["flag"] == cat] for cat in DISAGREEMENT_CATEGORIES}
    dcounts = {k: len(v) for k, v in dfull.items()}
    dgroups = {cat: items[:_DZ_CAP] for cat, items in dfull.items()}
    sv_conn = _sv_index()
    sv_total = SV.gene_sv_count(sv_conn, g) if sv_conn else 0
    sv_rows = SV.gene_svs(sv_conn, g, limit=_SV_CAP) if sv_conn else []
    return _render("gene_hub.html", gene=g, groups=groups, total=len(rows),
                   flagged=flagged, cap=_HUB_CAP, dgroups=dgroups, dcounts=dcounts,
                   dcap=_DZ_CAP, cats=DISAGREEMENT_CATEGORIES, meta=CAT_META,
                   sv_rows=sv_rows, sv_total=sv_total, sv_cap=_SV_CAP)


@app.get("/sv", response_class=HTMLResponse)
def sv_home():
    """Landing page for the structural-variant corpus."""
    conn = _sv_index()
    if not conn:
        raise StarletteHTTPException(404, "Structural variants are not in this release.")
    total = conn.execute("SELECT COUNT(*) FROM sv").fetchone()[0]
    sig = [{"significance": r[0] or "not provided", "n": r[1]} for r in conn.execute(
        "SELECT significance, COUNT(*) n FROM sv GROUP BY significance ORDER BY n DESC")]
    genes = [{"gene": r[0], "n": r[1]} for r in conn.execute(
        "SELECT gene, COUNT(*) n FROM sv_gene GROUP BY gene ORDER BY n DESC LIMIT 30")]
    return HTMLResponse(_render("sv_index.html", total=total, sig=sig,
                                genes=genes, nav="variant"))


@app.get("/sv/{key}.md", response_class=PlainTextResponse)
def sv_md(key: str, section: str = ""):
    """Markdown twin of a dbVar structural-variant page (NSSV / nsv / VCV / locus).

    `?section=` works as it does on a variant page."""
    rows = _sv_lookup(key)
    if not rows:
        return PlainTextResponse("This structural variant is not in the current release.\n",
                                 status_code=404)
    headers = {"Cache-Control": "public, max-age=3600",
               "Access-Control-Allow-Origin": "*",
               "Content-Type": "text/markdown; charset=utf-8"}
    release_line = release_line_md()
    if len(rows) > 1:
        md = "# Structural variants\n\n"
        if release_line:
            md += release_line + "\n\n"
        md += "\n".join(
            f"- {r['nssv']} — {r.get('svtype')} at chr{r['chrom']}:{r['start']}-{r['end']}"
            + (f", {r['significance']}" if r.get("significance") else "")
            for r in rows)
        return PlainTextResponse(md + "\n", headers=headers)
    rec = rows[0]
    md = SV.render_md(rec, _sv_genes(rec["nssv"]), release_line=release_line)
    headings = _md_headings(md)
    headers["X-Sugi-Sections"] = ", ".join(headings)
    if section:
        block = _md_section(md, section)
        if block is None:
            return PlainTextResponse(
                f"No section {section!r}. Available: {', '.join(headings)}\n",
                status_code=404, headers=headers)
        md = block
    return PlainTextResponse(md, headers=headers)


@app.get("/sv/{key}", response_class=HTMLResponse)
def sv_page(key: str):
    """dbVar structural-variant page. Read-only; never builds."""
    rows = _sv_lookup(key)
    if not rows:
        if RELEASE_MODE:
            raise StarletteHTTPException(404, f"“{key}” is not in the current release.")
        raise StarletteHTTPException(404, f"No structural variant for “{key}”.")
    if len(rows) == 1:
        rec, genes, hits = rows[0], _sv_genes(rows[0]["nssv"]), None
    else:
        rec, genes, hits = None, None, rows
    return HTMLResponse(_render("sv.html", v=rec, genes=genes, hits=hits, nav="variant"))


def _md_headings(md):
    """Level-2 headings of a rendered .md page, anchor and spacing stripped."""
    out = []
    for line in md.splitlines():
        if line.startswith("## "):
            name = line[3:].split("{#", 1)[0].strip()
            if name:
                out.append(name)
    return out


def _md_section(md, want):
    """The `## <want>` block of a .md page, or None. Case-insensitive, anchor-free.

    Agents fetch pages to cite them; a batch of 30 pages is ~230k characters, so the
    caller needs to ask for one section rather than receive everything. This is the
    server side of that, so the client stays a thin fetch.
    """
    lines = md.splitlines()
    key = (want or "").strip().lower().rstrip("s")
    start = None
    for i, line in enumerate(lines):
        if line.startswith("## "):
            name = line[3:].split("{#", 1)[0].strip().lower().rstrip("s")
            if name == key:
                start = i
                break
    if start is None:
        return None
    end = next((j for j in range(start + 1, len(lines)) if lines[j].startswith("## ")),
               len(lines))
    return "\n".join(lines[start:end]).strip() + "\n"


@app.get("/{slug}.md", response_class=PlainTextResponse)
def variant_md(slug: str, section: str = ""):
    """Markdown twin of a variant page — the same prose the page renders, as plain text
    for LLMs and AI agents (which ingest Markdown better than HTML), at the same URL as
    the page + '.md'. Deterministic; GET + cached; CORS-open.

    `?section=Identity` returns just that block (case-insensitive). Without it the whole
    page is returned; `X-Sugi-Sections` names the blocks a caller can ask for."""
    rec = _resolve_page(slug.lower().strip("/"))
    if not rec:
        if RELEASE_MODE:
            headers = {"X-Sugi-Variant-Release": str(_RELEASE_ID)} if _RELEASE_ID else {}
            return PlainTextResponse("This variant is not in the current release.\n",
                                     status_code=404, headers=headers)
        raise StarletteHTTPException(404, f"No variant page for “{slug}”.")
    from sugivariant.render import render_body
    title = rec.get("name") or f"{rec.get('gene_symbol', '')} {rec.get('hgvs_p') or rec.get('hgvs_c') or ''}".strip()
    release_line = (release_line_md() + "\n\n") if release_line_md() else ""
    md = (f"# {title}\n\n{_PUB}/{rec['canonical_slug']}\n\n"
          f"{release_line}" + render_body(rec))
    headings = _md_headings(md)
    headers = {"Cache-Control": "public, max-age=3600",
               "Access-Control-Allow-Origin": "*",
               "Content-Type": "text/markdown; charset=utf-8",
               "X-Sugi-Sections": ", ".join(headings)}
    if section:
        block = _md_section(md, section)
        if block is None:
            return PlainTextResponse(
                f"No section {section!r}. Available: {', '.join(headings)}\n",
                status_code=404, headers=headers)
        md = block
    return PlainTextResponse(md, headers=headers)


@app.get("/{slug}", response_class=HTMLResponse)
def variant_page(slug: str, view: str = ""):
    slug = slug.lower().strip("/")
    if _RSID_RE.match(slug):                       # rsID URL → dbSNP reverse-map
        return _resolution_response("rsid", slug, resolve_rsid(slug))
    rec = _resolve_page(slug)
    if not rec:
        if RELEASE_MODE:
            raise StarletteHTTPException(
                404, f"“{slug}” is not in the current release.")
        raise StarletteHTTPException(404, f"No variant page for “{slug}”.")
    # Add a live discordance cue to the similar-variant cross-links (the flag isn't in
    # the sibling records at build time; the index carries it per slug — cheap lookup).
    sims = rec.get("similar")
    ix = _index()
    if sims and ix:
        flags = IX.slug_flags(ix, [s["slug"] for s in sims])
        for s in sims:
            s["flag"] = flags.get(s["slug"])
    # One canonical layout (the dashboard). A legacy ?view= param is accepted but ignored
    # so old bookmarks/links still resolve to the same page instead of 404ing.
    # schema.org structured data: jsonld.py was written and tested but never wired in, so
    # every page shipped without it. Defensive — a schema error must never break a page.
    freq_g, freq_pops = _freq_ctx(rec)
    # Re-derive the pickled readouts from the record's own fields before ANYTHING reads
    # them (pure, no lookups). Explicit rather than relying on disagreement_flag()'s side
    # effect, which happens to run first today.
    #
    # This must precede the JSON-LD. refresh_frequency() exists to strip a stored gnomAD
    # band that asserts an absence the data cannot support, and the JSON-LD used to be
    # built ABOVE it -- so 187,493 pages (84,456 of them P/LP) served HTML saying "no
    # allele count reported in gnomAD v4.1 (AN 1,383,496)" beside structured data saying
    # "is absent from gnomAD v4.1". The machine-readable surface is the one search engines
    # and language models quote, so it was manufacturing PM2-supporting rarity on exactly
    # the pages where the human-readable one had already withdrawn the claim.
    _refresh(rec)
    # schema.org structured data: jsonld.py was written and tested but never wired in, so
    # every page shipped without it. Defensive — a schema error must never break a page.
    try:
        from sugivariant.jsonld import as_script_tag
        jsonld = as_script_tag(rec, {"generated_at": _index_lastmod()})
    except Exception:
        jsonld = ""
    html = _render("variant.html", v=rec, canonical=rec["canonical_slug"], nav="variant",
                   jsonld=jsonld, freq_g=freq_g, freq_pops=freq_pops)
    if rec.get("_degraded"):
        # Built while an upstream fetch was failing → incomplete evidence. Serve it,
        # but keep it out of every cache so the next request re-builds it cleanly.
        return HTMLResponse(html, headers={"Cache-Control": "no-store"})
    return html
