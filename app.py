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
from sugivariant.render import DATA_SOURCES              # noqa: E402
env.globals["data_sources"] = DATA_SOURCES
env.globals["tier_label"] = tier_label

app = FastAPI(title="Sugi Variant")
# Serve at ROOT (like Sugi Predict): nginx `proxy_pass …:8013/;` strips the
# /variant/ prefix, so routes + static are unprefixed. BASE_PATH is used ONLY
# in templates ({{ base }}) to generate the public /variant/… links.
app.mount("/static", StaticFiles(directory=str(ROOT / "static")), name="static")

# ── record resolution (in-process per-gene cache; ISR html-cache on top) ───────
_GENE_CACHE = {}
_SLUG_SEP = re.compile(r"-[pc]-")


def _gene_of(slug):
    m = _SLUG_SEP.search(slug)
    return slug[:m.start()] if m else None


def _records(gene):
    g = gene.upper()
    if g not in _GENE_CACHE:
        _GENE_CACHE[g] = enriched_records(g) or []
    return _GENE_CACHE[g]


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
        "label": "LoF resolves a conflicting call", "tone": "info",
        "blurb": "ClinVar submitters conflict, but it's a predicted loss-of-function "
                 "change in a loss-of-function-intolerant gene — a mechanism-based "
                 "resolving signal (the non-missense analog)."},
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
        h["Cache-Control"] = "no-cache, stale-if-error=86400"
        if request.headers.get("if-none-match") == etag:
            return Response(status_code=304, headers=h)
        return HTMLResponse(content=body, headers=h)
    return resp


@app.exception_handler(StarletteHTTPException)
async def _err(request, exc):
    return HTMLResponse(_render("error.html", code=exc.status_code,
                                detail=exc.detail), status_code=exc.status_code)


# ── routes ─────────────────────────────────────────────────────────────────────
@app.get("/")
async def home(q: str = ""):
    if q.strip():
        return _resolution_response(*resolve_query(q))
    demos = ["pten-p-arg173cys", "acta1-p-pro309ala", "asxl1-p-gly646trp",
             "pten-p-arg130gln", "acta1-c-809-1g-t"]
    return HTMLResponse(_render("home.html", demos=demos))


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
async def sitemap():
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
async def sitemap_gene(gene: str):
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
async def disagreements_hub():
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
async def disagreements_gene(gene: str):
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
async def gene_hub(gene: str):
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


@app.get("/{slug}", response_class=HTMLResponse)
async def variant_page(slug: str):
    slug = slug.lower().strip("/")
    if _RSID_RE.match(slug):                       # rsID URL → dbSNP reverse-map
        return _resolution_response("rsid", slug, resolve_rsid(slug))
    rec = _resolve(slug)
    if not rec:
        raise StarletteHTTPException(404, f"No variant page for “{slug}”.")
    # canonical: if hit via an alias, the template sets rel=canonical to the p-slug
    return _render("variant.html", v=rec, canonical=rec["canonical_slug"], nav="variant",
                   disagreement=disagreement_flag(rec))
