"""dbVar structural-variant corpus (prototype).

A separate corpus from the small-variant records. The entry key is the dbVar
NSSV accession, the location is an interval, and the source is biobtree dataset
826. dbVar carries no gene assignment; genes are found by overlapping the SV
interval with Ensembl gene coordinates.

Storage is one SQLite file, `sv.db`:

    sv            one row per NSSV (the record)
    sv_id         NSSV -> integer id (R-tree keys are integers)
    sv_interval   R-tree over start/end
    sv_alias      NSSV / nsv / VCV / chr:start-end -> NSSV
    sv_gene       gene -> NSSV overlaps (for gene hubs)

Read-only for serving; `build` is the only writer. This is deliberately separate
from records/<GENE>.pkl.gz and from the small-variant index: a CNV is not a VCV
missense record, and its predictors are not the SNV predictors.

CLI:
    python -m sugivariant.dbvar build --db /tmp/sv.db \
        --intervals 1:10001,22:19700000 --limit 200 --genes TMC1,CFTR --apply
    python -m sugivariant.dbvar lookup --db /tmp/sv.db --key NSSV16255736
    python -m sugivariant.dbvar overlap --db /tmp/sv.db --chrom 9 --start 72500000 --end 72900000
    python -m sugivariant.dbvar show --db /tmp/sv.db --nssv NSSV16255736
"""
from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from urllib.parse import unquote

SCHEMA = "sugi-variant-sv/1"

DB_SCHEMA = """
CREATE TABLE IF NOT EXISTS sv(
  nssv TEXT PRIMARY KEY, nsv TEXT, chrom TEXT, start INTEGER, end INTEGER,
  svtype TEXT, svlen TEXT, significance TEXT, origin TEXT, somatic INTEGER,
  description TEXT, phenotypes TEXT, clinvar TEXT);
CREATE TABLE IF NOT EXISTS sv_id(nssv TEXT PRIMARY KEY, id INTEGER UNIQUE);
CREATE VIRTUAL TABLE IF NOT EXISTS sv_interval USING rtree(id, start, end);
CREATE TABLE IF NOT EXISTS sv_alias(key TEXT, nssv TEXT, PRIMARY KEY(key, nssv));
CREATE INDEX IF NOT EXISTS ix_sv_chrom ON sv(chrom);
CREATE TABLE IF NOT EXISTS sv_gene(gene TEXT, nssv TEXT, PRIMARY KEY(gene, nssv));
"""


def _txt(x):
    return unquote(str(x)).strip() if x is not None else ""


def _list(x):
    if isinstance(x, list):
        return [unquote(str(v)).strip() for v in x if str(v).strip()]
    if isinstance(x, str) and x.strip():
        return [unquote(v).strip() for v in x.split(",") if v.strip()]
    return []


def record_from_attrs(nssv, attrs):
    """Build an SV record from a biobtree dbvar entry's Attributes."""
    d = (attrs or {}).get("Dbvar") or (attrs or {})
    try:
        start, end = int(d.get("start")), int(d.get("end"))
    except (TypeError, ValueError):
        return None
    chrom = _txt(d.get("chromosome")) or _txt(d.get("seq_region"))
    if not chrom:
        return None
    return {
        "schema": SCHEMA,
        "nssv": nssv,
        "nsv": _txt(d.get("region_id")) or None,
        "chrom": chrom,
        "start": min(start, end),
        "end": max(start, end),
        "svtype": _txt(d.get("svtype")) or None,
        "svlen": _txt(d.get("svlen")) or None,
        "significance": _txt(d.get("clinical_significance")) or None,
        "origin": _txt(d.get("origin")) or None,
        "somatic": 1 if d.get("somatic") in (True, 1, "1", "true", "True") else 0,
        "description": _txt(d.get("description")) or None,
        "phenotypes": _list(d.get("phenotype_ids")),
        "clinvar": _list(d.get("clinvar_accessions")),
    }


def collect(nssv):
    """Fetch one NSSV's full dbVar record, or None."""
    from sugibiobtree import entry
    attrs = (entry(nssv, "dbvar") or {}).get("Attributes") or {}
    return record_from_attrs(nssv, attrs)


def aliases(rec):
    """Every way a reader may ask for this SV."""
    keys = {rec["nssv"].lower()}
    if rec.get("nsv"):
        keys.add(rec["nsv"].lower())
    keys.add(f"{rec['chrom']}:{rec['start']}-{rec['end']}".lower())
    keys.add(f"{rec['chrom']}:{rec['start']}-{rec['end']}:{(rec.get('svtype') or '').lower()}"
             .lower())
    for a in rec.get("clinvar") or []:
        keys.add(a.lower())
        if a.upper().startswith("VCV"):
            num = a[3:].split(".")[0].lstrip("0") or "0"
            keys.add(num.lower())
    return {k for k in keys if k}


def _insert(conn, rec):
    cur = conn.cursor()
    row = cur.execute("SELECT id FROM sv_id WHERE nssv=?", (rec["nssv"],)).fetchone()
    if row:
        ident = row[0]
        cur.execute("DELETE FROM sv_interval WHERE id=?", (ident,))
        cur.execute("DELETE FROM sv_alias WHERE nssv=?", (rec["nssv"],))
        cur.execute("DELETE FROM sv_gene WHERE nssv=?", (rec["nssv"],))
    else:
        ident = (cur.execute("SELECT COALESCE(MAX(id),0)+1 FROM sv_id").fetchone()[0])
        cur.execute("INSERT INTO sv_id VALUES(?,?)", (rec["nssv"], ident))
    cur.execute(
        "INSERT OR REPLACE INTO sv VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (rec["nssv"], rec.get("nsv"), rec["chrom"], rec["start"], rec["end"],
         rec.get("svtype"), rec.get("svlen"), rec.get("significance"),
         rec.get("origin"), rec.get("somatic"), rec.get("description"),
         json.dumps(rec.get("phenotypes") or []),
         json.dumps(rec.get("clinvar") or [])))
    cur.execute("INSERT INTO sv_interval VALUES(?,?,?)", (ident, rec["start"], rec["end"]))
    cur.executemany("INSERT OR IGNORE INTO sv_alias VALUES(?,?)",
                    [(k, rec["nssv"]) for k in aliases(rec)])
    return ident


def open_db(path):
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.executescript(DB_SCHEMA)
    return conn


def open_ro(path):
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def lookup(conn, key):
    """Rows for an NSSV / nsv / VCV / chr:start-end alias."""
    return [dict(r) for r in conn.execute(
        "SELECT s.* FROM sv_alias a JOIN sv s ON s.nssv=a.nssv WHERE a.key=? "
        "ORDER BY s.nssv", ((key or "").lower(),))]


def overlapping(conn, chrom, start, end, limit=None):
    """SVs whose interval overlaps [start, end] on a chromosome."""
    q = ("SELECT s.* FROM sv_interval i JOIN sv_id d ON d.id=i.id "
         "JOIN sv s ON s.nssv=d.nssv "
         "WHERE i.start<=? AND i.end>=? AND s.chrom=? ORDER BY s.start, s.nssv")
    args = [int(end), int(start), str(chrom).lstrip("chr")]
    if limit:
        q += " LIMIT ?"
        args.append(int(limit))
    return [dict(r) for r in conn.execute(q, args)]


def genes_for_sv(conn, nssv):
    return [r[0] for r in conn.execute(
        "SELECT gene FROM sv_gene WHERE nssv=? ORDER BY gene", (nssv,))]


def gene_coords(symbol):
    """(chrom, start, end) for a gene symbol via HGNC -> Ensembl, or None."""
    from sugibiobtree import entry, map_all
    from sugivariant.build import _hgnc_id
    h = _hgnc_id(symbol)
    if not h:
        return None
    rows = list(map_all(h, ">>hgnc>>ensembl", cap=1))
    if not rows:
        return None
    attrs = (entry(rows[0].get("id"), "ensembl") or {}).get("Attributes") or {}
    e = attrs.get("Ensembl") or {}
    try:
        return (str(e.get("seq_region")), int(e.get("start")), int(e.get("end")))
    except (TypeError, ValueError):
        return None


def annotate_genes(conn, genes):
    """Insert sv_gene rows for the given gene symbols. Returns (gene, n_svs) list."""
    out = []
    for g in genes:
        c = gene_coords(g)
        if not c:
            out.append((g, 0))
            continue
        chrom, start, end = c
        n = 0
        for rec in overlapping(conn, chrom, start, end):
            conn.execute("INSERT OR IGNORE INTO sv_gene VALUES(?,?)", (g, rec["nssv"]))
            n += 1
        out.append((g, n))
    conn.commit()
    return out


def build(db_path, intervals, limit=200, genes=None, apply=True, verbose=True):
    """Build a prototype SV corpus from dbVar intervals (or explicit NSSVs)."""
    from sugibiobtree import map_all
    nssvs = []
    # Spread the sample across the supplied intervals instead of filling the limit
    # from the first one: a 1p36 interval alone can return hundreds of large CNVs.
    per = max(1, limit // max(1, len(intervals))) if limit else None
    for coord in intervals:
        got = 0
        for r in map_all(coord, ">>dbvar", cap=None):
            nssv = r.get("id")
            if nssv and nssv not in nssvs:
                nssvs.append(nssv)
                got += 1
            if per and got >= per:
                break
    nssvs = nssvs[:limit] if limit else nssvs
    conn = open_db(db_path) if apply else None
    built = 0
    for nssv in nssvs:
        rec = collect(nssv)
        if not rec:
            continue
        built += 1
        if verbose and (built % 50 == 0):
            print(f"  collected {built}/{len(nssvs)}", flush=True)
        if conn:
            _insert(conn, rec)
    if conn:
        conn.commit()
        if genes:
            for g, n in annotate_genes(conn, genes):
                if verbose:
                    print(f"  {g}: {n} overlapping SVs")
        conn.close()
    if verbose:
        print(f"{'built' if apply else 'would build'} {built} SVs from "
              f"{len(intervals)} interval(s) -> {db_path}")
    return built


def render_md(rec, genes=None):
    size = ""
    try:
        size = f"{abs(int(rec['svlen'])):,} bp"
    except (TypeError, ValueError, KeyError):
        pass
    L = [f"# dbVar {rec['nssv']} — {rec.get('svtype') or 'structural variant'}", "",
         f"**{rec.get('svtype') or 'SV'}** at "
         f"chr{rec['chrom']}:{rec['start']}-{rec['end']}"
         + (f" ({size})" if size else "") + ".", ""]
    rows = [
        ("Type", rec.get("svtype")),
        ("Size", size),
        ("Location", f"chr{rec['chrom']}:{rec['start']}-{rec['end']} (GRCh38)"),
        ("Clinical significance", rec.get("significance")),
        ("Origin", rec.get("origin")),
        ("Somatic", "yes" if rec.get("somatic") else None),
        ("dbVar region", rec.get("nsv")),
        ("ClinVar", ", ".join(rec.get("clinvar") or []) or None),
        ("Phenotypes", ", ".join(rec.get("phenotypes") or []) or None),
    ]
    L += ["## Identity {#identity}", "", "| Field | Value |", "| --- | --- |"]
    L += [f"| {k} | {v} |" for k, v in rows if v]
    L += ["", "## Overlapping genes {#genes}", ""]
    if genes:
        L.append("Genes whose Ensembl span overlaps this interval: "
                 + ", ".join(f"[{g}](https://sugi.bio/atlas/gene/{g}/)" for g in genes) + ".")
    else:
        L.append("No gene overlap annotated in this prototype build.")
    if rec.get("description"):
        L += ["", "## dbVar description", "", rec["description"]]
    L += ["", "## Interpretation", "",
          "This page reports the dbVar record and its ClinVar accession. The "
          "small-variant in-silico predictors (AlphaMissense, REVEL, SaProt, "
          "SpliceAI) do not apply to a structural variant; no ACMG code is "
          "assigned here.", ""]
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--db", required=True)
    b.add_argument("--intervals", required=True,
                   help="comma-separated chr:pos points, e.g. 1:10001,22:19700000")
    b.add_argument("--limit", type=int, default=200)
    b.add_argument("--genes", default="", help="comma-separated symbols to annotate")
    b.add_argument("--apply", action="store_true")
    for name in ("lookup", "overlap", "show"):
        p = sub.add_parser(name)
        p.add_argument("--db", required=True)
    s = sub.choices
    s["lookup"].add_argument("--key", required=True)
    s["overlap"].add_argument("--chrom", required=True)
    s["overlap"].add_argument("--start", type=int, required=True)
    s["overlap"].add_argument("--end", type=int, required=True)
    s["show"].add_argument("--nssv", required=True)
    a = ap.parse_args(argv)

    if a.cmd == "build":
        genes = [g.strip().upper() for g in a.genes.split(",") if g.strip()]
        intervals = [x.strip() for x in a.intervals.split(",") if x.strip()]
        return build(a.db, intervals, a.limit, genes, a.apply)

    conn = open_ro(a.db)
    if a.cmd == "lookup":
        rows = lookup(conn, a.key)
    elif a.cmd == "overlap":
        rows = overlapping(conn, a.chrom, a.start, a.end)
    else:
        rows = [dict(r) for r in conn.execute("SELECT * FROM sv WHERE nssv=?", (a.nssv,))]
    if not rows:
        print("no results")
        return 1
    if a.cmd == "show" and len(rows) == 1:
        rec = rows[0]
        for k in ("phenotypes", "clinvar"):
            try:
                rec[k] = json.loads(rec.get(k) or "[]")
            except Exception:
                rec[k] = []
        print(render_md(rec, genes_for_sv(conn, rec["nssv"])))
    else:
        print(json.dumps(rows, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
