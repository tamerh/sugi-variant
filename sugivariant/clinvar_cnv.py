"""ClinVar CNVs: copy-number variants whose clinical record is ClinVar, not dbVar.

The dbVar corpus covers what was deposited in the archive. A large share of
clinical CNVs never were: they are ClinVar submissions with a coordinate-bearing
name (`GRCh38/hg38 7q11.23(chr7:74031409-74104712)x3`) and a germline
classification, and no dbVar NSSV points at them. Before this pass those variants
had no page anywhere — the small-variant corpus drops them because a CNV has no
`c.` HGVS to make a slug from, and the SV corpus is dbVar-derived.

The rule that keeps the two sources honest:

    every clinical CNV gets exactly one page.

If a VCV is already reachable through a dbVar NSSV, that page stays and the VCV is
aliased to it. If it is not, a VCV-keyed page is written from ClinVar's own
coordinates and classification, with `source='clinvar'` and no dbVar id claimed.

CLI:
    python -m sugivariant.clinvar_cnv build --db cache/sv.db \
        --records cache/records --limit 50 --apply
    python -m sugivariant.clinvar_cnv build --db cache/sv.db --records cache/records \
        --shard 0/6 --apply
"""
from __future__ import annotations

import argparse
import json
import re
import sys

from sugibiobtree import entry

from sugivariant import collect as VC
from sugivariant import dbvar as SV
from sugivariant.build import _hgnc_id

SCHEMA = "sugi-variant-sv/1"

# `GRCh38/hg38 7q11.23(chr7:74031409-74104712)x3`, and the hg19 equivalent we skip.
_CNV_NAME = re.compile(
    r"^(?P<build>GRCh38/hg38|GRCh37/hg19)\s+[^(]*"
    r"\(chr(?P<chrom>[0-9XYM]+):(?P<start>\d+)-(?P<end>\d+)\)", re.I)

# ClinVar's type -> the type vocabulary the SV pages already use.
_TYPE = {
    "copy number loss": "DEL",
    "copy number gain": "DUP",
    "copy number variation": "CNV",
}


def vcv_id(variation_id):
    return f"VCV{int(variation_id):09d}"


def cnv_from_attrs(variation_id, attrs, gene=None):
    """A corpus record for one ClinVar CNV, or None if it cannot be placed.

    Coordinates come from the structured start/stop when ClinVar carries them,
    otherwise from the name. Only GRCh38 is accepted: the page states its assembly,
    and silently placing an hg19 interval on hg38 would be wrong by megabases.
    """
    d = (attrs or {}).get("Clinvar") or (attrs or {})
    if not isinstance(d, dict):
        return None
    kind = (d.get("type") or "").strip().lower()
    if kind not in _TYPE:
        return None
    name = (d.get("name") or "").strip()
    m = _CNV_NAME.match(name)
    chrom = start = end = None
    if m and m.group("build").lower().startswith("grch38"):
        chrom, start, end = m.group("chrom"), int(m.group("start")), int(m.group("end"))
    if start is None and d.get("start") and d.get("stop"):
        chrom = str(d.get("chromosome") or "").strip()
        try:
            start, end = int(d["start"]), int(d["stop"])
        except (TypeError, ValueError):
            start = end = None
    if not (chrom and start and end):
        return None
    ident = vcv_id(variation_id)
    rec = {
        "schema": SCHEMA,
        "nssv": ident,
        "nsv": None,
        "chrom": chrom,
        "start": min(start, end),
        "end": max(start, end),
        "svtype": _TYPE[kind],
        "svlen": str(abs(end - start)),
        "significance": (d.get("germline_classification") or "").strip() or None,
        "origin": "germline",
        "somatic": 0,
        "description": name or None,
        "phenotypes": [p for p in (d.get("phenotype_ids") or []) if p],
        "clinvar": [ident],
        "source": "clinvar",
    }
    if not SV.should_build(rec):
        return None
    return rec


def _covered_by_dbvar(conn, rec):
    """The NSSV page that already serves this VCV, if any."""
    keys = [k for k in SV.aliases(rec) if k]
    if not keys:
        return None
    qs = ",".join("?" * len(keys))
    rows = conn.execute(
        f"SELECT DISTINCT nssv FROM sv_alias WHERE key IN ({qs})", keys).fetchall()
    for r in rows:
        nssv = r[0]
        if nssv and not nssv.upper().startswith("VCV"):
            return nssv
    return None


def _cnv_ids(rows):
    """ClinVar ids whose name is a placeable CNV. The gene comes from the caller."""
    return [str(vid) for vid, _cls, name in rows
            if name and _CNV_NAME.match(name.strip())]


def build(db_path, records_dir=None, genes=None, limit=0, apply=False, shard=None):
    import os
    if not genes:
        if not records_dir:
            raise SystemExit("need --records or --genes")
        genes = sorted(p[:-len(".pkl.gz")] for p in os.listdir(records_dir)
                       if p.endswith(".pkl.gz"))
    if shard:
        k, n = (int(x) for x in shard.split("/"))
        genes = [g for i, g in enumerate(genes) if i % n == k]
    if limit:
        genes = genes[:limit]

    conn = SV.open_db(db_path) if apply else SV.open_ro(db_path)
    seen, added, linked, skipped = set(), 0, 0, 0
    for i, gene in enumerate(genes, 1):
        hgnc = _hgnc_id(gene)
        if not hgnc:
            continue
        try:
            rows = VC.enumerate_gene(hgnc) or []
        except Exception:
            continue
        for vid in _cnv_ids(rows):
            if vid in seen:
                continue
            seen.add(vid)
            try:
                attrs = (entry(int(vid), "clinvar") or {}).get("Attributes") or {}
            except Exception:
                skipped += 1
                continue
            rec = cnv_from_attrs(vid, attrs)
            if not rec:
                skipped += 1
                continue
            existing = _covered_by_dbvar(conn, rec)
            if existing:
                linked += 1
                if apply:
                    conn.execute("INSERT OR IGNORE INTO sv_alias VALUES(?,?)",
                                 (rec["nssv"].lower(), existing))
                    conn.execute("INSERT OR IGNORE INTO sv_gene VALUES(?,?)",
                                 (gene, existing))
                continue
            added += 1
            if apply:
                SV._insert(conn, rec)
                conn.execute("INSERT OR IGNORE INTO sv_gene VALUES(?,?)",
                             (gene, rec["nssv"]))
                if added % 200 == 0:
                    conn.commit()
        if apply and i % 25 == 0:
            conn.commit()
        if i % 500 == 0:
            print(f"  {i}/{len(genes)} genes · added {added} · linked {linked} · skipped {skipped}",
                  flush=True)
    if apply:
        conn.commit()
    if genes and not seen:
        print("WARNING: no gene resolved to an HGNC id -- check --records", flush=True)
    print(f"{'APPLY' if apply else 'DRY RUN'} -- genes {len(genes)} · CNV pages {added} · "
          f"already in dbVar {linked} · skipped {skipped}", flush=True)
    return {"added": added, "linked": linked, "skipped": skipped}


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--db", required=True)
    b.add_argument("--records", default=None)
    b.add_argument("--genes", default=None, help="comma-separated symbols")
    b.add_argument("--shard", default=None, help="K/N")
    b.add_argument("--limit", type=int, default=0)
    b.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "build":
        genes = [g.strip() for g in a.genes.split(",") if g.strip()] if a.genes else None
        build(a.db, a.records, genes, a.limit, a.apply, a.shard)
    return 0


if __name__ == "__main__":
    sys.exit(main())
