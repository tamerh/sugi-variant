"""Persistent resolution index (SQLite).

A normalized lookup key → ClinVar VCV map, plus a light per-variant row. Built
OFFLINE in one pass over the genes (see build_index); the live app reads it for
instant resolution across all three supported arrival keys — rsID, GRCh38
coordinate (which has NO biobtree→ClinVar edge, so the index is the only path),
and GENE + HGVS/AA — and for a global disagreements list and the sitemap, without
rebuilding a gene just to resolve a query.

Scope: this is resolution + listing, NOT a page store. Rendering a page still
builds the gene (enrichment needs the per-gene caches); the index just tells the
app which VCV/slug a query means and whether a page exists at all.
"""
import os
import re
import sqlite3

from sugivariant.build import enriched_records
from sugivariant.slug import _norm, _norm_hgvs
from sugivariant.enrich import disagreement_flag, missense_short, review_stars

SCHEMA = """
CREATE TABLE IF NOT EXISTS variant(
  vcv INTEGER PRIMARY KEY, gene TEXT, slug TEXT, hgvs_p TEXT, hgvs_c TEXT,
  classification TEXT, review_status TEXT, stars INTEGER, rsid TEXT, coordinate TEXT,
  primary_condition TEXT, flag TEXT);
CREATE TABLE IF NOT EXISTS alias(key TEXT, vcv INTEGER, PRIMARY KEY(key, vcv));
CREATE INDEX IF NOT EXISTS ix_alias ON alias(key);
CREATE INDEX IF NOT EXISTS ix_gene  ON variant(gene);
CREATE INDEX IF NOT EXISTS ix_flag  ON variant(flag);
CREATE TABLE IF NOT EXISTS gene_meta(gene TEXT PRIMARY KEY, n INTEGER);
"""


def open_db(path, check_same_thread=True):
    conn = sqlite3.connect(path, check_same_thread=check_same_thread)
    conn.row_factory = sqlite3.Row
    conn.executescript(SCHEMA)
    return conn


def norm_coordinate(s):
    """Normalize a user coordinate to the stored 'chr:pos:ref:alt' key (drop a
    'chr' prefix, uppercase alleles, accept ':' or '-' separators)."""
    s = (s or "").strip()
    m = re.match(r"^(?:chr)?([0-9xymt]+)[:\-](\d+)[:\-]([acgt]+)[:\-]([acgt]+)$", s, re.I)
    return f"{m.group(1)}:{m.group(2)}:{m.group(3).upper()}:{m.group(4).upper()}" if m else None


def variant_keys(rec):
    """All normalized resolution keys for a record: slug aliases, rsID, GRCh38
    coordinate, 1-letter short (R130Q), and the bare (operator-stripped) p./c.
    forms — mirrors the live resolver's _hgvs_keys so both paths agree."""
    keys = set(rec.get("slugs") or [])
    if rec.get("rsid"):
        keys.add(rec["rsid"].lower())
    if rec.get("coordinate"):
        keys.add(rec["coordinate"])
    g = _norm(rec["gene_symbol"])
    short = missense_short(rec.get("hgvs_p"))
    if short:
        keys.add(f"{g}-{_norm_hgvs(short)}")
    for form in (rec.get("hgvs_p"), rec.get("hgvs_c")):
        if form:
            keys.add(f"{g}-{_norm_hgvs(re.sub(r'^[pc][.]', '', form))}")
    return {k for k in keys if k}


def index_gene(conn, gene, recs=None):
    """(Re)index one gene: replace its rows + aliases. Returns the record count."""
    recs = recs if recs is not None else (enriched_records(gene) or [])
    cur = conn.cursor()
    vcvs = [int(r["variation_id"]) for r in recs if r.get("variation_id")]
    if vcvs:
        qs = ",".join("?" * len(vcvs))
        cur.execute(f"DELETE FROM alias WHERE vcv IN ({qs})", vcvs)
        cur.execute(f"DELETE FROM variant WHERE vcv IN ({qs})", vcvs)
    for r in recs:
        vcv = int(r["variation_id"])
        flag = disagreement_flag(r)
        cond = (r.get("conditions") or [{}])[0].get("name")
        cur.execute(
            "INSERT OR REPLACE INTO variant VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
            (vcv, r["gene_symbol"], r["canonical_slug"], r.get("hgvs_p"), r.get("hgvs_c"),
             r.get("classification"), r.get("review_status"), review_stars(r.get("review_status")),
             r.get("rsid"), r.get("coordinate"), cond, flag["category"] if flag else None))
        cur.executemany("INSERT OR IGNORE INTO alias VALUES(?,?)",
                        [(k.lower(), vcv) for k in variant_keys(r)])  # lookup() lowercases
    cur.execute("INSERT OR REPLACE INTO gene_meta VALUES(?,?)", (gene.upper(), len(recs)))
    conn.commit()
    return len(recs)


def build_index(genes, db_path):
    """Full/partial build over a gene list. Idempotent per gene (re-runnable)."""
    conn = open_db(db_path)
    total = 0
    for g in genes:
        total += index_gene(conn, g.upper())
    conn.close()
    return total


# ── read side (used by the live app) ────────────────────────────────────────────
def lookup(conn, key):
    """Normalized key → [variant rows]. 0, 1 or many (an rsID is position-level)."""
    rows = conn.execute(
        "SELECT v.* FROM alias a JOIN variant v ON v.vcv=a.vcv WHERE a.key=? "
        "ORDER BY v.gene, v.slug", (key.lower() if key else "",)).fetchall()
    return [dict(r) for r in rows]


def all_slugs(conn):
    return [r["slug"] for r in conn.execute("SELECT slug FROM variant ORDER BY slug")]


def genes(conn):
    """[(gene, n)] over ALL built genes (from gene_meta)."""
    return [(r["gene"], r["n"]) for r in
            conn.execute("SELECT gene, n FROM gene_meta ORDER BY gene")]


def gene_slugs(conn, gene):
    return [r["slug"] for r in conn.execute(
        "SELECT slug FROM variant WHERE gene=? ORDER BY slug", (gene.upper(),))]


def gene_rows(conn, gene):
    """All light variant rows for a gene (for the /gene/{SYM} hub), best first
    (review stars desc, then slug)."""
    return [dict(r) for r in conn.execute(
        "SELECT slug, hgvs_p, hgvs_c, classification, stars, primary_condition, flag "
        "FROM variant WHERE gene=? ORDER BY stars DESC, slug", (gene.upper(),))]


# Sitemap inclusion gate: advertise to Google only pages with assertion criteria
# (>=1 star) AND a named condition. The thin 0-star / no-condition tail is still
# SERVED (the resolver answers any hit) — just not advertised, to avoid the
# scaled-thin-content risk on a young YMYL domain. Serving is NEVER gated here.
_SITEMAP_WHERE = "stars >= 1 AND primary_condition IS NOT NULL AND primary_condition != ''"


def sitemap_genes(conn):
    """[(gene, eligible_n)] — only genes with >=1 sitemap-eligible page."""
    return [(r["gene"], r["n"]) for r in conn.execute(
        f"SELECT gene, COUNT(*) n FROM variant WHERE {_SITEMAP_WHERE} GROUP BY gene ORDER BY gene")]


def sitemap_gene_slugs(conn, gene):
    return [r["slug"] for r in conn.execute(
        f"SELECT slug FROM variant WHERE gene=? AND {_SITEMAP_WHERE} ORDER BY slug", (gene.upper(),))]


def disagreements(conn, category=None, limit=None):
    q = "SELECT * FROM variant WHERE flag IS NOT NULL"
    args = []
    if category:
        q += " AND flag=?"
        args.append(category)
    q += " ORDER BY gene, slug"
    if limit:
        q += " LIMIT ?"
        args.append(limit)
    return [dict(r) for r in conn.execute(q, args)]


def stats(conn):
    n = conn.execute("SELECT COUNT(*) FROM variant").fetchone()[0]
    g = conn.execute("SELECT COUNT(*) FROM gene_meta").fetchone()[0]
    f = conn.execute("SELECT COUNT(*) FROM variant WHERE flag IS NOT NULL").fetchone()[0]
    return {"variants": n, "genes": g, "flagged": f}


if __name__ == "__main__":
    import sys
    db = os.environ.get("INDEX_DB", "cache/index.db")
    genes = sys.argv[1:]
    print(f"building {db} over {len(genes)} genes…")
    n = build_index(genes, db)
    print(f"indexed {n} variants -> {db} ({os.path.getsize(db)/1024/1024:.1f} MB)")
    print(stats(open_db(db)))
