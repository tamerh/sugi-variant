"""Run a rebuild-queue.json: rebuild each gene, write its pickle, reindex it.

Resumable: a progress file records completed and failed genes, so an interrupted
run continues where it stopped. This is the worker that closes the page-level
queue; it is not part of the request path.

    queue_rebuild.py --queue cache/rebuild-queue.json \
        --records cache/records --index cache/index.db \
        --progress cache/queue-progress.json --apply [--limit N]

The queue entries are the page-level gaps from `snapshot.py build`:
`queue` (no_record / pages_behind) plus `reindex` (stale slugs). `alias_only`
entries are not rebuilt.
"""
import argparse
import gzip
import json
import os
import pickle
import sys
import time

from sugivariant.build import enriched_records
from sugivariant.enrich import (
    records_collapsed, stored_record_count, write_data_stamp,
)
from sugivariant.index import index_gene, open_db


def _load_json(path, default):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default


def _save_json(path, obj):
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
    os.replace(tmp, path)


def _write(records_path, recs):
    tmp = f"{records_path}.{os.getpid()}.tmp"
    with gzip.open(tmp, "wb") as f:
        pickle.dump(recs, f)
    os.replace(tmp, records_path)
    return write_data_stamp(records_path)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--queue", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--progress", required=True)
    ap.add_argument("--limit", type=int, default=0, help="0 means all")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)

    q = _load_json(a.queue, {})
    entries = list(q.get("queue") or []) + list(q.get("reindex") or [])
    # One gene once, preserving queue order (largest gaps first).
    genes = list(dict.fromkeys(e["gene"] for e in entries if e.get("gene")))
    prog = _load_json(a.progress, {"done": [], "failed": {}})
    done = set(prog.get("done") or [])
    todo = [g for g in genes if g not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"{'APPLY' if a.apply else 'DRY RUN'} -- {len(todo)} of {len(genes)} genes "
          f"({len(done)} already done)", flush=True)

    conn = open_db(a.index) if a.apply else None
    t0 = time.time()
    for i, g in enumerate(todo, 1):
        path = os.path.join(a.records, f"{g}.pkl.gz")
        before = stored_record_count(path)
        t = time.time()
        try:
            recs = enriched_records(g) or []
        except Exception as e:
            print(f"  {g:<12} ERROR {type(e).__name__}: {str(e)[:100]}", flush=True)
            prog.setdefault("failed", {})[g] = f"{type(e).__name__}: {e}"
            _save_json(a.progress, prog)
            continue
        if not recs:
            print(f"  {g:<12} 0 records -- skipped", flush=True)
            prog.setdefault("failed", {})[g] = "0 records"
            _save_json(a.progress, prog)
            continue
        collapsed = records_collapsed(path, len(recs))
        if collapsed:
            print(f"  {g:<12} REFUSED {len(recs)} would replace {collapsed}", flush=True)
            prog.setdefault("failed", {})[g] = f"collapse {len(recs)}<{collapsed}"
            _save_json(a.progress, prog)
            continue
        print(f"  {g:<12} {str(before):>7} -> {len(recs):>7} records  {time.time()-t:6.1f}s",
              flush=True)
        if a.apply:
            _write(path, recs)
            index_gene(conn, g, recs=recs)
            done.add(g)
            prog["done"] = sorted(done)
            prog.setdefault("failed", {}).pop(g, None)
            _save_json(a.progress, prog)
            if i % 10 == 0:
                print(f"  ... {i}/{len(todo)}  {time.time()-t0:.0f}s", flush=True)
    if conn:
        conn.commit()
        conn.close()
    print(f"finished {len(todo)} genes in {time.time()-t0:.0f}s; "
          f"{len(done)} done, {len(prog.get('failed') or {})} failed", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
