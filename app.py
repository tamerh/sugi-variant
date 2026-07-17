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


env.globals.update(label=variant_label, cls_class=cls_class, stars=stars)
# Data-source attribution (AlphaMissense CC BY 4.0 + REVEL ODbL require it) — single
# source of truth in render.py so HTML and the markdown twin can't drift.
from sugivariant.render import data_provenance           # noqa: E402
env.globals["data_provenance"] = data_provenance
env.globals["tier_label"] = tier_label

app = FastAPI(title="Sugi Variant")
# Serve at ROOT (like Sugi Predict): nginx `proxy_pass …:8013/;` strips the
# /variant/ prefix, so routes + static are unprefixed. BASE_PATH is used ONLY
# in templates ({{ base }}) to generate the public /variant/… links.
app.mount("/static", StaticFiles(directory=str(ROOT / "static")), name="static")

# ── record resolution (in-memory + persistent per-gene record cache) ───────────
# A gene's enriched records are expensive to build (a whale like BRCA2 ≈ minutes),
# so once built they're persisted to disk (gzip pickle) and loaded on later hits —
# this survives restarts and stops the cold-build Cloudflare timeout.
_GENE_CACHE = {}
_RECORDS_DIR = CACHE_DIR / "records"
_SLUG_SEP = re.compile(r"-[pcnm]-")   # p./c./n.(ncRNA)/m.(mito)


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
    if g in _GENE_CACHE:
        return _GENE_CACHE[g]
    import gzip
    import pickle
    path = _RECORDS_DIR / f"{g}.pkl.gz"
    # per-gene lock: only ONE thread builds a given gene; concurrent hits wait and
    # then get the cached result (no duplicate multi-minute whale builds).
    with _gene_lock(g):
        if g in _GENE_CACHE:
            return _GENE_CACHE[g]
        if path.exists():
            try:
                with gzip.open(path, "rb") as f:
                    _GENE_CACHE[g] = pickle.load(f)
                return _GENE_CACHE[g]
            except Exception:
                pass
        recs = enriched_records(g) or []
        _GENE_CACHE[g] = recs
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


def _resolve(slug):
    gene = _gene_of(slug)
    if not gene:
        return None
    for r in _records(gene):
        if r["canonical_slug"] == slug or slug in (r.get("slugs") or []):
            return r
    return None


def _render(tpl, **ctx):
    return env.get_template(tpl).render(**ctx)


# ── persistent resolution index (optional; falls back to the live path if absent) ─
INDEX_DB = pathlib.Path(os.environ.get("INDEX_DB") or (CACHE_DIR / "index.db"))
_INDEX = None


def _index():
    """The SQLite resolution index, or None if it hasn't been built. Opened once,
    read-only across request threads."""
    global _INDEX
    if _INDEX is None and INDEX_DB.exists():
        _INDEX = IX.open_db(str(INDEX_DB), check_same_thread=False)
    return _INDEX


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
        "label": "Resolves a conflicting call", "tone": "info",
        "blurb": "ClinVar submitters conflict, but the independent predictors "
                 "unanimously agree — a resolving in-silico read."},
    "lof_resolves_conflicting": {
        "label": "Predicted LoF vs a conflicting call", "tone": "info",
        "blurb": "ClinVar submitters conflict, and it's a predicted loss-of-function "
                 "change in a gene where LoF causes disease — flagged for review, not "
                 "resolved (verify NMD-escape for C-terminal truncations). Non-missense."},
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


@app.get("/")
def home(q: str = ""):
    if q.strip():
        return _resolution_response(*resolve_query(q))
    s = _corpus_stats()
    return HTMLResponse(_render("home.html", n_variants=s["variants"],
                                n_genes=s["genes"], n_flagged=s["flagged"]))


_SET_MAX = 60
_FLAG_SEVERITY = {"predictor_vs_clinvar": 3, "resolves_conflicting": 2,
                  "lof_resolves_conflicting": 2, "predictors_split": 1}


@app.get("/set", response_class=HTMLResponse)
def variant_set(v: str = "", q: str = ""):
    ix = _index()
    # paste box → resolve each line via the shared resolver → canonical permalink
    if q.strip():
        slugs = []
        for line in re.split(r"[\n;]+", q):
            line = line.strip()
            if not line:
                continue
            for h in resolve_query(line)[2]:
                if h["canonical_slug"] not in slugs:
                    slugs.append(h["canonical_slug"])
        if slugs:
            return RedirectResponse(f"{BASE}/set?v=" + ",".join(slugs[:_SET_MAX]), status_code=307)
        miss = [ln.strip() for ln in re.split(r"[\n;]+", q) if ln.strip()]
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
    return _render("method.html", nav="method")


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


@app.get("/disagreements", response_class=HTMLResponse)
def disagreements_hub():
    # Counts only for genes already built (cold BRCA2 ≈ minutes) — link the rest.
    featured = []
    for g in FEATURED_GENES:
        counts = None
        if g in _GENE_CACHE:
            counts = {k: len(v) for k, v in _grouped_flags(_GENE_CACHE[g]).items()}
        featured.append({"gene": g, "counts": counts})
    return _render("disagreements_hub.html", featured=featured,
                   cats=DISAGREEMENT_CATEGORIES, meta=CAT_META)


@app.get("/disagreements/{gene}", response_class=HTMLResponse)
def disagreements_gene(gene: str):
    gene = gene.upper().strip("/")
    recs = _records(gene)
    if not recs:
        raise StarletteHTTPException(404, f"No variants built for “{gene}”.")
    groups = _grouped_flags(recs)
    counts = {k: len(v) for k, v in groups.items()}
    return _render("disagreements_gene.html", gene=gene, groups=groups, counts=counts,
                   flagged=sum(counts.values()), total=len(recs),
                   cats=DISAGREEMENT_CATEGORIES, meta=CAT_META)


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
    return _render("gene_hub.html", gene=g, groups=groups, total=len(rows),
                   flagged=flagged, cap=_HUB_CAP)


_VIEW_TEMPLATES = {"datasheet": "variant_v3.html", "classic": "variant_classic.html"}


@app.get("/{slug}", response_class=HTMLResponse)
def variant_page(slug: str, view: str = ""):
    slug = slug.lower().strip("/")
    if _RSID_RE.match(slug):                       # rsID URL → dbSNP reverse-map
        return _resolution_response("rsid", slug, resolve_rsid(slug))
    rec = _resolve(slug)
    if not rec:
        raise StarletteHTTPException(404, f"No variant page for “{slug}”.")
    # ?view= selects a layout preview (Default / Dashboard / Datasheet); the switcher
    # bar links between them. canonical stays the p-slug regardless.
    tpl = _VIEW_TEMPLATES.get(view, "variant.html")
    return _render(tpl, v=rec, canonical=rec["canonical_slug"], nav="variant",
                   view=view, disagreement=disagreement_flag(rec))
