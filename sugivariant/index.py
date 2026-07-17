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
            keys.add(f"{g}-{_norm_hgvs(re.sub(r'^[pc][.]', '', form))}")  # operator-stripped
            keys.add(f"{g}-{_norm_hgvs(form)}")   # full form incl. c./p. — the pre-cap slug, so an
                                                  # over-long variant's old (uncapped) URL still resolves
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


def build_index(genes, db_path, skip_done=True, log_every=200):
    """Full/partial build over a gene list. RESUMABLE: genes already in gene_meta are
    skipped (idempotent + survives interruption). Per-gene failures are recorded as
    n=-1 and skipped, not fatal — re-run after deleting the -1 rows to retry them."""
    import time
    conn = open_db(db_path)
    done = ({r["gene"] for r in conn.execute("SELECT gene FROM gene_meta")} if skip_done else set())
    todo = [g.upper() for g in genes if g.upper() not in done]
    total, failed, t0 = 0, 0, time.time()
    print(f"  {len(done)} already done, {len(todo)} to build", flush=True)
    for i, g in enumerate(todo, 1):
        try:
            total += index_gene(conn, g)
        except Exception as e:
            failed += 1
            conn.execute("INSERT OR REPLACE INTO gene_meta VALUES(?,?)", (g, -1))
            conn.commit()
            print(f"  ERR {g}: {str(e)[:80]}", flush=True)
        if i % log_every == 0:
            dt = time.time() - t0
            print(f"  {i}/{len(todo)} genes · {total} variants · {failed} failed · "
                  f"{dt:.0f}s · {i/dt:.2f} genes/s", flush=True)
    conn.close()
    print(f"  built {total} variants, {failed} genes failed", flush=True)
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


def directory_genes(conn):
    """[(gene, variant_count)] for every gene that has at least one built variant —
    the /genes A-Z directory."""
    return [(r["gene"], r["n"]) for r in conn.execute(
        "SELECT gene, COUNT(*) n FROM variant GROUP BY gene ORDER BY gene")]


def suggest(conn, q, limit=8):
    """Type-ahead suggestions. A single token → matching GENES (prefix). A gene +
    partial change (e.g. 'PTEN R130') → matching VARIANTS via the alias keys. Returns
    [{kind, label, sub, url}] (url is relative to the app base)."""
    q = (q or "").strip()
    if len(q) < 2:
        return []
    toks = [t for t in re.split(r"[\s,]+", q) if t]
    gene_tok = next((t for t in toks if "." not in t and ">" not in t), None)
    rest = " ".join(t for t in toks if t != gene_tok).strip() if gene_tok else ""
    out = []
    # gene + partial change → variant suggestions (match the normalized alias keys)
    if gene_tok and rest:
        key = f"{_norm(gene_tok)}-%{_norm_hgvs(rest)}%"
        seen = set()
        # best-first: higher review tier, then stronger classification, then slug
        for r in conn.execute(
                "SELECT v.slug, v.gene, v.hgvs_p, v.hgvs_c, v.classification "
                "FROM alias a JOIN variant v ON v.vcv=a.vcv WHERE a.key LIKE ? "
                "ORDER BY v.stars DESC, "
                "  CASE WHEN v.classification='Pathogenic' THEN 0 "
                "       WHEN v.classification LIKE 'Pathogenic/%' THEN 1 "
                "       WHEN v.classification LIKE 'Likely%' THEN 2 ELSE 3 END, v.slug "
                "LIMIT ?",
                (key, limit * 3)):
            if r["slug"] in seen:
                continue
            seen.add(r["slug"])
            out.append({"kind": "variant", "url": r["slug"], "sub": r["classification"],
                        "label": f"{r['gene']} {r['hgvs_p'] or r['hgvs_c']}"})
            if len(out) >= limit:
                break
        if out:
            return out
    # otherwise → gene prefix suggestions (most-populated first)
    pref = (gene_tok or toks[0]).upper() + "%"
    for r in conn.execute(
            "SELECT gene, n FROM gene_meta WHERE gene LIKE ? AND n>0 ORDER BY n DESC, gene LIMIT ?",
            (pref, limit)):
        out.append({"kind": "gene", "url": f"gene/{r['gene']}", "label": r["gene"],
                    "sub": f"{r['n']} variants"})
    return out


def set_rows(conn, slugs):
    """Resolve a list of requested slugs (canonical or alias) to light variant rows
    via the alias table — no per-gene build. Returns (found_rows, missing_slugs),
    order-preserving and deduped by VCV."""
    found, missing, seen = [], [], set()
    for s in slugs:
        r = conn.execute(
            "SELECT v.* FROM alias a JOIN variant v ON v.vcv=a.vcv WHERE a.key=? LIMIT 1",
            ((s or "").lower(),)).fetchone()
        if r and r["vcv"] not in seen:
            seen.add(r["vcv"])
            found.append(dict(r))
        elif not r:
            missing.append(s)
    return found, missing


def gene_count(conn, gene):
    r = conn.execute("SELECT n FROM gene_meta WHERE gene=?", (gene.upper(),)).fetchone()
    return (r["n"] if r else 0) or 0


def corpus_stats(conn):
    """Headline counts for the home page."""
    return {
        "variants": conn.execute("SELECT COUNT(*) FROM variant").fetchone()[0],
        "genes": conn.execute("SELECT COUNT(DISTINCT gene) FROM variant").fetchone()[0],
        "flagged": conn.execute("SELECT COUNT(*) FROM variant WHERE flag IS NOT NULL").fetchone()[0],
    }


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
    import argparse
    ap = argparse.ArgumentParser(description="Build the resolution index (resumable).")
    ap.add_argument("--genes-file", help="one gene symbol per line")
    ap.add_argument("--db", default=os.environ.get("INDEX_DB", "cache/index.db"))
    ap.add_argument("genes", nargs="*", help="gene symbols (or use --genes-file)")
    a = ap.parse_args()
    genes = a.genes
    if a.genes_file:
        genes = [ln.strip() for ln in open(a.genes_file) if ln.strip() and not ln.startswith("#")]
    print(f"building {a.db} over {len(genes)} genes (resumable)…", flush=True)
    build_index(genes, a.db)
    sz = os.path.getsize(a.db) / 1024 / 1024
    print(f"-> {a.db} ({sz:.1f} MB) | {stats(open_db(a.db))}")
