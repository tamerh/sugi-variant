"""Snapshot manifest for the records corpus.

READ-ONLY by construction: this module only *describes* what is on disk (the
per-gene pickles and their ``.dbver.json`` stamps) and in the resolution index.
It never builds a gene, writes a record, or touches biobtree. That is what lets
it be run safely against the live cache before any serving path changes.

Two jobs:

    build    scan records/ + index.db -> snapshot.json + an audit report
    verify   re-check a manifest against the files it describes

The manifest is deliberately multi-corpus from the start (``small_variant``
today, ``structural_variant`` when dbVar lands) so adding a corpus does not
require a schema change. Only the small-variant corpus is populated here.

Usage
-----
    python -m sugivariant.snapshot build \
        --records /data2/sugi-variant-cache-dev/records \
        --index   /data2/sugi-variant-cache-dev/index.db \
        --out     /tmp/snapshot.json \
        --report  /tmp/snapshot-audit.md

    python -m sugivariant.snapshot verify \
        --manifest /tmp/snapshot.json \
        --records  /data2/sugi-variant-cache-dev/records \
        --index    /data2/sugi-variant-cache-dev/index.db

Add ``--count-records`` to load every pickle and compare its true length with
``gene_meta.n`` (slower: one gzip read per gene). ``--count-sample N`` does that
for a deterministic sample only. ``--current-versions versions.json`` compares
each gene's stamp with the current dataset versions and reports staleness.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import pickle
import random
import sqlite3
import sys
from datetime import datetime, timezone

SCHEMA = "sugi-variant-snapshot/1"
RECORD_SUFFIX = ".pkl.gz"
STAMP_SUFFIX = ".dbver.json"
CORPUS = "small_variant"

# Status priority: the first matching state wins as the primary `status`; every
# problem found is also recorded in `issues`.
#
# `known_empty` is NOT a problem: gene_meta says the gene was enumerated and
# nothing was filed under it, so there is no record file by design. Counting
# those as "missing" (they are 10.9k of the 29.8k gene_meta rows) hid the 218
# genes that have index rows and no record file -- the real gap a snapshot
# would serve as a 404.
_STATUS_ORDER = ("missing", "corrupt", "mismatch", "records_only", "orphan",
                 "stale", "unstamped", "empty", "known_empty", "ok")


def _utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _stamp_path(records_path):
    d, base = os.path.split(str(records_path))
    return os.path.join(d, base[:-len(RECORD_SUFFIX)] + STAMP_SUFFIX)


def read_stamp(records_path):
    """The sidecar stamp for a gene, or None when it is absent/unreadable."""
    try:
        with open(_stamp_path(records_path), encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else None
    except Exception:
        return None


def _count_records(records_path):
    """(count, error). Counts a pickle without importing the app."""
    try:
        with gzip.open(records_path, "rb") as f:
            recs = pickle.load(f)
        return (len(recs), None) if isinstance(recs, list) else (None, "not-a-list")
    except Exception as e:                       # corrupt / truncated / not a pickle
        return None, type(e).__name__


def scan_records(records_dir, count_records=False, count_sample=None):
    """{GENE: {path, size, mtime, n, count_error, stamp}} for every *.pkl.gz.

    File existence is always scanned for every gene; ``n`` is populated only when
    ``count_records`` is set, or for a deterministic sample with ``count_sample``.
    """
    found = {}
    try:
        entries = list(os.scandir(records_dir))
    except FileNotFoundError:
        raise SystemExit(f"records dir not found: {records_dir}")
    for e in entries:
        if not e.name.endswith(RECORD_SUFFIX):
            continue
        gene = e.name[:-len(RECORD_SUFFIX)]
        try:
            st = e.stat()
        except OSError:
            continue
        found[gene] = {
            "path": e.name,
            "size": st.st_size,
            "mtime": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(timespec="seconds"),
            "n": None,
            "count_error": None,
            "stamp": read_stamp(e.path),
        }
    if count_records:
        targets = sorted(found)
    elif count_sample:
        targets = sorted(found)
        if len(targets) > count_sample:
            targets = random.Random(0).sample(targets, count_sample)
    else:
        targets = []
    for gene in targets:
        n, err = _count_records(os.path.join(records_dir, found[gene]["path"]))
        found[gene]["n"] = n
        found[gene]["count_error"] = err
    return found


def read_index(index_db):
    """Read-only view of the index: gene_meta, per-gene variant counts, total."""
    conn = sqlite3.connect(f"file:{index_db}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        gene_meta = {r["gene"]: r["n"] for r in conn.execute("SELECT gene, n FROM gene_meta")}
        variant_counts = {
            r["gene"]: r["n"]
            for r in conn.execute("SELECT gene, COUNT(*) AS n FROM variant GROUP BY gene")
        }
        total = conn.execute("SELECT COUNT(*) FROM variant").fetchone()[0]
        flagged = conn.execute("SELECT COUNT(*) FROM variant WHERE flag IS NOT NULL").fetchone()[0]
    finally:
        conn.close()
    return {
        "gene_meta": gene_meta,
        "variant_counts": variant_counts,
        "total": total,
        "flagged": flagged,
    }


def _stale_from(stamp, current):
    if not current or not isinstance(stamp, dict):
        return []
    return sorted(d for d, v in current.items() if stamp.get(d) != v)


def build_manifest(records_dir, index_db, current_versions=None,
                   count_records=False, count_sample=None):
    """(manifest, report). Pure read; no writes, no biotree."""
    records = scan_records(records_dir, count_records=count_records,
                           count_sample=count_sample)
    index = read_index(index_db)
    index_genes = set(index["gene_meta"]) | set(index["variant_counts"])
    record_genes = set(records)
    genes = sorted(index_genes | record_genes)

    gene_state = {}
    missing, known_empty, empty, corrupt, mismatch, records_only, orphan = (
        [], [], [], [], [], [], [])
    stale, unstamped = [], []
    for g in genes:
        r = records.get(g)
        meta_n = index["gene_meta"].get(g)
        idx_n = index["variant_counts"].get(g, 0)
        issues = []
        if r is None:
            # A gene with no record file is only a gap if the index has rows for it.
            # gene_meta.n == 0 means "enumerated, nothing filed here" -- an expected
            # empty, not a missing page.
            if idx_n > 0:
                issues.append("missing")
                missing.append(g)
                status = "missing"
            else:
                issues.append("known_empty")
                known_empty.append(g)
                status = "known_empty"
            state = {"path": None, "size": None, "mtime": None, "n": None,
                     "stamp": None, "status": status, "issues": issues}
        else:
            n = r["n"]
            if idx_n == 0:
                # A record file the index has no rows for. Distinguish a legacy file
                # under an old symbol (gene_meta absent or 0: not part of the current
                # corpus) from a current gene_meta entry with no rows (real drift).
                if meta_n in (None, 0):
                    issues.append("orphan")
                    orphan.append(g)
                else:
                    issues.append("records_only")
                    records_only.append(g)
            if r["count_error"]:
                issues.append("corrupt")
                corrupt.append({"gene": g, "error": r["count_error"]})
            if r["size"] == 0 or (n is not None and n == 0):
                issues.append("empty")
                empty.append(g)
            expected = meta_n if meta_n not in (None, -1) else idx_n
            if n is not None and expected is not None and n != expected:
                issues.append("mismatch")
                mismatch.append({"gene": g, "records": n, "gene_meta": meta_n,
                                 "index_rows": idx_n})
            stale_on = _stale_from(r["stamp"], current_versions)
            if stale_on:
                issues.append("stale")
                stale.append({"gene": g, "datasets": stale_on})
            if r["stamp"] is None:
                issues.append("unstamped")
                unstamped.append(g)
            status = next((s for s in _STATUS_ORDER if s in issues), "ok")
            state = {"path": r["path"], "size": r["size"], "mtime": r["mtime"],
                     "n": n, "stamp": r["stamp"], "status": status,
                     "issues": issues}
        gene_state[g] = state

    counted = [s for s in gene_state.values() if s["n"] is not None]
    all_counted = bool(gene_state) and len(counted) == len(gene_state)
    # The SERVED corpus is the genes the index has rows for AND a record file exists
    # for. Orphan files (legacy symbols) and missing genes are reported, but they are
    # not part of the snapshot's addressable set.
    served = {g: s for g, s in gene_state.items()
              if s["path"] is not None and "orphan" not in s["issues"]}
    served_counted = [s for s in served.values() if s["n"] is not None]
    served_all_counted = bool(served) and len(served_counted) == len(served)
    record_total = sum(s["n"] for s in served_counted) if served_counted else None
    status = "ok" if not (missing or corrupt or mismatch or records_only) else "partial"

    manifest = {
        "schema": SCHEMA,
        "snapshot_id": _utc_now(),
        "created_at": _utc_now(),
        "corpora": {
            CORPUS: {
                "format": "per-gene-pickle-v1",
                "records_dir": os.path.abspath(records_dir),
                "index_db": os.path.abspath(index_db),
                "genes": len(served),
                "index_genes": len(index_genes),
                "record_files": len(record_genes),
                "orphan_files": len(orphan),
                "known_empty": len(known_empty),
                "missing": len(missing),
                "records": record_total if served_all_counted else None,
                "index_records": index["total"],
                "counted": served_all_counted,
                "datasets": dict(current_versions or {}),
                "scope": "all-cached-genes",
                "status": status,
                "gene_state": gene_state,
            }
        },
    }
    report = {
        "snapshot_id": manifest["snapshot_id"],
        "records_dir": os.path.abspath(records_dir),
        "index_db": os.path.abspath(index_db),
        "totals": {
            "record_genes": len(record_genes),
            "served_genes": len(served),
            "index_genes": len(index_genes),
            "index_records": index["total"],
            "flagged": index["flagged"],
            "counted_genes": len(served_counted),
            "record_total": record_total,
        },
        "index_only": sorted(index_genes - record_genes),
        "records_only": records_only,
        "orphan": orphan,
        "missing": missing,
        "known_empty": known_empty,
        "empty": empty,
        "corrupt": corrupt,
        "mismatch": mismatch,
        "stale": stale,
        "unstamped": unstamped,
    }
    return manifest, report


def verify_manifest(manifest, records_dir, index_db, count_records=False, sample=None):
    """Re-check a manifest against the files/index it describes. Returns a report."""
    corpus = ((manifest or {}).get("corpora") or {}).get(CORPUS) or {}
    state = corpus.get("gene_state") or {}
    genes = sorted(state)
    if sample and len(genes) > sample:
        genes = random.Random(0).sample(genes, sample)
    index = read_index(index_db)
    errors, checked = [], 0
    for g in genes:
        s = state[g]
        checked += 1
        if s.get("status") == "missing":
            continue
        path = os.path.join(records_dir, s.get("path") or f"{g}{RECORD_SUFFIX}")
        if not os.path.exists(path):
            errors.append({"gene": g, "error": "file-missing", "path": path})
            continue
        if s.get("size") is not None and os.path.getsize(path) != s["size"]:
            errors.append({"gene": g, "error": "size-changed",
                           "expected": s["size"], "actual": os.path.getsize(path)})
        if count_records:
            n, err = _count_records(path)
            if err:
                errors.append({"gene": g, "error": f"corrupt:{err}"})
            elif s.get("n") is not None and n != s["n"]:
                errors.append({"gene": g, "error": "count-changed",
                               "expected": s["n"], "actual": n})
        meta_n = index["gene_meta"].get(g)
        if s.get("n") is not None and meta_n not in (None, -1) and s["n"] != meta_n:
            errors.append({"gene": g, "error": "index-count-changed",
                           "expected": s["n"], "actual": meta_n})
    return {"ok": not errors, "checked": checked, "errors": errors,
            "snapshot_id": manifest.get("snapshot_id")}


def _md_report(report, manifest):
    t = report["totals"]

    def block(title, items, fmt=lambda x: str(x), cap=100):
        out = [f"## {title} ({len(items)})", ""]
        if not items:
            return out + ["None.", ""]
        out += [fmt(x) for x in items[:cap]]
        if len(items) > cap:
            out.append(f"… and {len(items) - cap} more")
        return out + [""]

    L = ["# Snapshot audit", "",
         f"- snapshot: `{report['snapshot_id']}`",
         f"- records: `{report['records_dir']}`",
         f"- index: `{report['index_db']}`",
         f"- schema: `{manifest['schema']}`", "",
         "## Totals", "",
         f"- record files: **{t['record_genes']:,}**",
         f"- served genes: **{t['served_genes']:,}**",
         f"- index genes: **{t['index_genes']:,}**",
         f"- index rows: **{t['index_records']:,}**",
         f"- flagged: **{t['flagged']:,}**",
         f"- counted genes: **{t['counted_genes']:,}**",
         f"- counted records: **{t['record_total']:,}**" if t["record_total"] is not None
         else "- counted records: not counted",
         ""]
    L += block("Index genes with rows but NO record file (real gap)", report["missing"])
    L += block("Record files with no current index rows (orphan/legacy)", report["orphan"])
    L += block("Record files whose gene_meta says rows but index has none", report["records_only"])
    L += block("Empty record files", report["empty"])
    L += block("Known-empty genes (no page expected)", report["known_empty"])
    L += block("Corrupt record files", report["corrupt"],
               lambda x: f"- {x['gene']}: {x['error']}")
    L += block("Count mismatches (pickle vs gene_meta)", report["mismatch"],
               lambda x: f"- {x['gene']}: pickle {x['records']}, gene_meta {x['gene_meta']}, "
                         f"index rows {x['index_rows']}")
    L += block("Stale by dataset version", report["stale"],
               lambda x: f"- {x['gene']}: {', '.join(x['datasets'])}")
    L += block("Unstamped (no .dbver.json)", report["unstamped"])
    return "\n".join(L).rstrip() + "\n"


def _write_atomic(path, text):
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="scan the cache and write a manifest + audit")
    b.add_argument("--records", required=True, help="records/ directory")
    b.add_argument("--index", required=True, help="index.db")
    b.add_argument("--out", default="snapshot.json", help="manifest path")
    b.add_argument("--report", help="optional Markdown audit path")
    b.add_argument("--current-versions", help="JSON {dataset: version} to check staleness")
    b.add_argument("--count-records", action="store_true",
                   help="load every pickle and compare with gene_meta.n")
    b.add_argument("--count-sample", type=int, metavar="N",
                   help="count only a deterministic sample of N genes")
    b.add_argument("--strict", action="store_true",
                   help="exit 1 if any missing/empty/corrupt/count-mismatch is found")

    v = sub.add_parser("verify", help="re-check a manifest against the cache")
    v.add_argument("--manifest", required=True)
    v.add_argument("--records", required=True)
    v.add_argument("--index", required=True)
    v.add_argument("--count-records", action="store_true")
    v.add_argument("--sample", type=int, metavar="N")
    v.add_argument("--strict", action="store_true")

    a = ap.parse_args(argv)

    if a.cmd == "build":
        current = None
        if a.current_versions:
            with open(a.current_versions, encoding="utf-8") as f:
                current = json.load(f)
        manifest, report = build_manifest(
            a.records, a.index, current_versions=current,
            count_records=a.count_records, count_sample=a.count_sample)
        _write_atomic(a.out, json.dumps(manifest, indent=2, sort_keys=False) + "\n")
        md = _md_report(report, manifest)
        if a.report:
            _write_atomic(a.report, md)
        t = report["totals"]
        print(f"snapshot {report['snapshot_id']}")
        print(f"  record files : {t['record_genes']:,}")
        print(f"  served genes : {t['served_genes']:,}")
        print(f"  index genes  : {t['index_genes']:,}")
        print(f"  index rows   : {t['index_records']:,}")
        print(f"  missing      : {len(report['missing']):,}  (index rows, no record)")
        print(f"  orphan       : {len(report['orphan']):,}  (records, no current index rows)")
        print(f"  records-only : {len(report['records_only']):,}")
        print(f"  known-empty  : {len(report['known_empty']):,}  (no page expected)")
        print(f"  empty        : {len(report['empty']):,}")
        print(f"  corrupt      : {len(report['corrupt']):,}")
        print(f"  mismatches   : {len(report['mismatch']):,}")
        print(f"  stale        : {len(report['stale']):,}")
        print(f"  unstamped    : {len(report['unstamped']):,}")
        if t["record_total"] is not None:
            print(f"  counted      : {t['record_total']:,} records")
        print(f"  manifest     : {a.out}")
        if a.report:
            print(f"  report       : {a.report}")
        if a.strict and (report["missing"] or report["empty"] or
                         report["corrupt"] or report["mismatch"]):
            return 1
        return 0

    with open(a.manifest, encoding="utf-8") as f:
        manifest = json.load(f)
    res = verify_manifest(manifest, a.records, a.index,
                          count_records=a.count_records, sample=a.sample)
    print(f"verify {res['snapshot_id']}  checked {res['checked']:,}  "
          f"{'OK' if res['ok'] else str(len(res['errors'])) + ' errors'}")
    for e in res["errors"][:50]:
        print("  ", e)
    if len(res["errors"]) > 50:
        print(f"  … and {len(res['errors']) - 50} more")
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
