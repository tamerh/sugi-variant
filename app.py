"""Sugi Variant — the web. A FastAPI + Jinja2 + ISR renderer over the self-contained
`sugivariant` package (decoupled from Sugi Atlas). Server-rendered
HTML (crawlable / AI-citable — the whole differentiator), disk-cached after first
hit. Sibling of the Sugi Predict app; reuses the atlas.css design system.

Run: uvicorn app:app --host 0.0.0.0 --port 8013
"""
import hashlib
import os
import pathlib
import re
import sys

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, Response, RedirectResponse
from fastapi.staticfiles import StaticFiles
from starlette.exceptions import HTTPException as StarletteHTTPException
from jinja2 import Environment, FileSystemLoader, select_autoescape

# The variant science lives in the local, self-contained sugivariant package.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sugibiobtree import map_all                          # noqa: E402
from sugivariant.build import enriched_records            # noqa: E402
from sugivariant.collect import collect                   # noqa: E402
from sugivariant.render import _label as variant_label    # noqa: E402
from sugivariant.slug import _norm, _norm_hgvs            # noqa: E402
from sugivariant.enrich import (review_stars, review_tier, tier_label,  # noqa: E402
                                missense_short, disagreement_flag,
                                DISAGREEMENT_CATEGORIES)
from sugivariant import index as IX                       # noqa: E402

ROOT = pathlib.Path(__file__).parent
BASE = os.environ.get("BASE_PATH", "").rstrip("/")
CACHE_DIR = pathlib.Path(os.environ.get("CACHE_DIR") or (ROOT / "cache"))

env = Environment(loader=FileSystemLoader(str(ROOT / "templates")),
                  autoescape=select_autoescape(["html"]))
env.globals["base"] = BASE
env.globals["css_v"] = lambda: str(int(max(
    (ROOT / "static" / "css" / f).stat().st_mtime for f in ("app.css", "atlas.css"))))

# ── classification → colour class + display + star rendering (template helpers) ─
_CLS_CLASS = {
    "pathogenic": "path", "likely pathogenic": "lpath",
    "pathogenic/likely pathogenic": "path",
    "conflicting classifications of pathogenicity": "conflict",
    "uncertain significance": "vus", "likely benign": "lben", "benign": "ben",
}


def cls_class(classification):
    return _CLS_CLASS.get((classification or "").strip().lower(), "vus")


def stars(review_status):
    # named review-status tier (not a ★ rating) — {n, label, cls}
    return review_tier(review_status)


from sugivariant.render import short_hgvs as _short_hgvs   # noqa: E402
env.globals.update(label=variant_label, cls_class=cls_class, stars=stars, hgvs_disp=_short_hgvs)
# Data-source attribution (AlphaMissense CC BY 4.0 + REVEL ODbL require it) — single
# source of truth in render.py so HTML and the markdown twin can't drift.
from sugivariant.render import data_provenance           # noqa: E402
env.globals["data_provenance"] = data_provenance
env.globals["tier_label"] = tier_label

# ClinGen dosage haploinsufficiency: 0–3 is the evidence scale; 30/40 are category codes
# (not points on the scale), so don't render them as "30/3".
_HI_CATEGORY = {"30": "gene assoc. with recessive phenotype", "40": "dosage-sensitivity unlikely"}
env.globals["hi_display"] = lambda h: _HI_CATEGORY.get(str(h).strip(), f"{h}/3")

# Inline markdown-bold in generated verdict text: **x** → <strong>x</strong>.
# `mdbold` returns safe HTML (rest escaped); `mdstrip` yields plain text (for
# headings / meta where we just want the words, no markers).
import markupsafe as _ms                                  # noqa: E402
_MD_BOLD = re.compile(r"\*\*(.+?)\*\*")
env.filters["mdbold"] = lambda s: _ms.Markup(
    _MD_BOLD.sub(r"<strong>\1</strong>", str(_ms.escape(s)))) if s else ""
env.filters["mdstrip"] = lambda s: _MD_BOLD.sub(r"\1", s or "")

app = FastAPI(title="Sugi Variant")
# Serve at ROOT (like Sugi Predict): nginx `proxy_pass …:8013/;` strips the
# /variant/ prefix, so routes + static are unprefixed. BASE_PATH is used ONLY
# in templates ({{ base }}) to generate the public /variant/… links.
app.mount("/static", StaticFiles(directory=str(ROOT / "static")), name="static")


# ── preprint (Method page) ─────────────────────────────────────────────────────
# The /method page inlines the make4ht HTML export of the preprint. Build flow
# (mirrors Sugi Predict): in the preprint repo run `make html`, then
# copy main.html / main.pdf / main*.svg into static/preprint/ and main.css into
# static/preprint/main.raw.css, and run scripts/build_preprint_css.py to scope the
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
            _, ev = _GENE_CACHE.popitem(last=False)    # drop the least-recently-used gene
            _GENE_CACHE_BYTES -= len(ev) * _BYTES_PER_VARIANT


def _gene_of(slug):
    m = _SLUG_SEP.search(slug)
    return slug[:m.start()] if m else None


import threading                                        # noqa: E402
_GENE_LOCKS, _LOCKS_LOCK = {}, threading.Lock()


def _gene_lock(g):
    with _LOCKS_LOCK:
        return _GENE_LOCKS.setdefault(g, threading.Lock())


def _records(gene):
    g = gene.upper()
    recs = _cache_get(g)
    if recs is not None:
        return recs
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
                _cache_put(g, recs)          # LRU insert (may evict cold genes; they stay on disk)
                return recs
            except Exception:
                pass
        recs = enriched_records(g) or []
        _cache_put(g, recs)
        if recs:                     # don't persist empty/invalid-gene results
            try:
                _RECORDS_DIR.mkdir(parents=True, exist_ok=True)
                tmp = path.with_suffix(".tmp")
                with gzip.open(tmp, "wb") as f:
                    pickle.dump(recs, f)
                tmp.replace(path)    # atomic
            except Exception:
                pass
    return recs


def _collect_one(vcv):
    """Build a single variant record straight from its VCV (index hit) — no per-gene
    build. Used as a resolution fallback and for the whale fast path."""
    from sugivariant.collect import attach_enrichment
    try:
        rec = collect(str(vcv))
        if rec:
            attach_enrichment(rec, {})              # single-variant enrichment (no per-gene caches)
            rec["_partial"] = True
            return rec
    except Exception:
        pass
    return None


def _resolve(slug):
    gene = _gene_of(slug)
    recs = _records(gene) if gene else []
    for r in recs:
        if r["canonical_slug"] == slug or slug in (r.get("slugs") or []):
            return r
    # Not in the gene's own records — resolve via the index alias table.
    ix = _index()
    if ix:
        hits = IX.lookup(ix, slug)
        if hits:
            vcv = str(hits[0]["vcv"])
            for r in recs:                          # legacy/alias-only slug (e.g. pre-cap delins URL)
                if str(r.get("variation_id")) == vcv:
                    return r
            # The VCV isn't under this gene at all — ClinVar's sort-gene differs from the
            # HGVS-derived gene (overlapping loci, e.g. MT-ATP8 filed under MT-ATP6). The
            # per-gene build will never contain it, so build the single variant directly.
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
        return None
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


def _render(tpl, **ctx):
    return env.get_template(tpl).render(**ctx)


# ── persistent resolution index (optional; falls back to the live path if absent) ─
INDEX_DB = pathlib.Path(os.environ.get("INDEX_DB") or (CACHE_DIR / "index.db"))
_tls = threading.local()


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


def _ix_rec(row):
    """Adapt an index row to the record shape templates/redirects expect."""
    return {"canonical_slug": row["slug"], "gene_symbol": row["gene"],
            "hgvs_p": row["hgvs_p"], "hgvs_c": row["hgvs_c"],
            "classification": row["classification"], "review_status": row["review_status"],
            "rsid": row["rsid"], "variation_id": row["vcv"]}


def _index_hits(key):
    """[adapted records] for a normalized key via the index, or None if no index."""
    ix = _index()
    return [_ix_rec(r) for r in IX.lookup(ix, key)] if ix else None


# ── query resolution (rsID reverse-map + gene-first HGVS match) ─────────────────
# Three arrival keys are supported (biobtree keys on identifiers, not free-text
# HGVS): an rsID resolves directly via dbSNP reverse-map; a GENE + HGVS/AA query
# resolves gene-first (enumerate + slug-match, all on our side). See the biobtree
# key-scheme note — raw HGVS is not searchable, so we normalize it to OUR slug.
_RSID_RE = re.compile(r"^rs\d+$", re.I)


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
        return RedirectResponse(f"{BASE}/{hits[0]['canonical_slug']}", status_code=307)
    return HTMLResponse(_render("results.html", kind=kind, term=term, hits=hits))


# ── disagreement browse view ────────────────────────────────────────────────────
# The one signal no competitor surfaces (benchmark 2026-07): variants where the
# evidence doesn't line up. Descriptive QC, not reclassification (HANDOVER §8).
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
        h["Cache-Control"] = "public, max-age=0, s-maxage=86400, stale-while-revalidate=604800, stale-if-error=86400"
        if request.headers.get("if-none-match") == etag:
            return Response(status_code=304, headers=h)
        return HTMLResponse(content=body, headers=h)
    return resp


@app.exception_handler(StarletteHTTPException)
async def _err(request, exc):
    return HTMLResponse(_render("error.html", code=exc.status_code,
                                detail=exc.detail), status_code=exc.status_code)


# ── routes ─────────────────────────────────────────────────────────────────────
_STATS_CACHE = None


def _corpus_stats():
    global _STATS_CACHE
    if _STATS_CACHE is None:
        ix = _index()
        _STATS_CACHE = (IX.corpus_stats(ix) if ix
                        else {"variants": 0, "genes": 0, "flagged": 0})
    return _STATS_CACHE


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
        return _resolution_response(*resolve_query(q))
    s = _corpus_stats()
    return HTMLResponse(_render("home.html", n_variants=s["variants"],
                                n_genes=s["genes"], n_flagged=s["flagged"], nav="home"))


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
        for line in _SET_SPLIT.split(q):
            line = line.strip()
            if not line:
                continue
            for h in resolve_query(line)[2]:
                if h["canonical_slug"] not in slugs:
                    slugs.append(h["canonical_slug"])
        if slugs:
            return RedirectResponse(f"{BASE}/set?v=" + ",".join(slugs[:_SET_MAX]), status_code=307)
        miss = [ln.strip() for ln in _SET_SPLIT.split(q) if ln.strip()]
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


def _index_lastmod():
    """Build-level lastmod = the index db's mtime (ISO date). Honest: the page was
    (re)built when the index was."""
    try:
        import datetime
        return datetime.date.fromtimestamp(INDEX_DB.stat().st_mtime).isoformat()
    except Exception:
        return None


@app.get("/sitemap.xml")
def sitemap():
    # Sitemap INDEX (not a flat urlset): a 600k-URL urlset breaches the sitemaps.org
    # 50k/50MB cap and Google silently drops it. One child sitemap per gene, each
    # well under the cap.
    ix = _index()
    if not ix:
        raise StarletteHTTPException(503, "Sitemap needs the resolution index; run sugivariant.index.")
    lm = _index_lastmod()
    lm_tag = f"<lastmod>{lm}</lastmod>" if lm else ""
    children = "".join(
        f"<sitemap><loc>{_PUB}/sitemap-{g.lower()}.xml</loc>{lm_tag}</sitemap>"
        for g, _ in IX.sitemap_genes(ix))
    xml = ('<?xml version="1.0" encoding="UTF-8"?>'
           '<sitemapindex xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
           f'{children}</sitemapindex>')
    return Response(content=xml, media_type="application/xml")


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
    return Response(content=xml, media_type="application/xml")


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


@app.get("/gene/{gene}", response_class=HTMLResponse)
def gene_hub(gene: str):
    g = gene.upper().strip("/")
    ix = _index()
    rows = IX.gene_rows(ix, g) if ix else None
    if not rows:                                   # not indexed yet → live build fallback
        recs = _records(g)
        if not recs:
            raise StarletteHTTPException(404, f"No variants built for “{g}”.")
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
    return _render("gene_hub.html", gene=g, groups=groups, total=len(rows),
                   flagged=flagged, cap=_HUB_CAP, dgroups=dgroups, dcounts=dcounts,
                   dcap=_DZ_CAP, cats=DISAGREEMENT_CATEGORIES, meta=CAT_META)


_VIEW_TEMPLATES = {"datasheet": "variant_v3.html", "classic": "variant_classic.html"}


@app.get("/{slug}", response_class=HTMLResponse)
def variant_page(slug: str, view: str = ""):
    slug = slug.lower().strip("/")
    if _RSID_RE.match(slug):                       # rsID URL → dbSNP reverse-map
        return _resolution_response("rsid", slug, resolve_rsid(slug))
    rec = _resolve_page(slug)
    if not rec:
        raise StarletteHTTPException(404, f"No variant page for “{slug}”.")
    # ?view= selects a layout preview (Default / Dashboard / Datasheet); the switcher
    # bar links between them. canonical stays the p-slug regardless.
    tpl = _VIEW_TEMPLATES.get(view, "variant.html")
    return _render(tpl, v=rec, canonical=rec["canonical_slug"], nav="variant",
                   view=view, disagreement=disagreement_flag(rec))
