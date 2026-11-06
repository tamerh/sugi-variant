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


def _stored_vcvs(path):
    """(variation_id, gene_symbol, hgvs_p) for every record in a stored pickle."""
    try:
        with gzip.open(path, "rb") as f:
            recs = pickle.load(f)
    except Exception:
        return None
    out = []
    for r in recs or []:
        v = r.get("variation_id") or r.get("vcv")
        if v is not None:
            out.append((str(v), (r.get("gene_symbol") or "").upper(),
                        (r.get("hgvs_p") or "").strip()))
    return out


def carrier_reshard(records_dir, gene, recs):
    """Is a shrinking write a legitimate reshard, or a real collapse?

    A gene's pickle can legitimately carry another gene's records: ClinVar files
    that overlap the locus come back under the host gene, and `gene_symbol` then
    comes from the HGVS transcript (the "carrier-gene defect"). When those carrier
    rows are dropped from the host gene, the record count falls, but nothing is
    lost *if* the other gene's own pickle already serves every dropped variant.

    A variant counts as served when the other gene's pickle carries its VCV **or**
    a record with the same protein change (`hgvs_p`): two ClinVar VCVs that differ
    only in the nucleotide change (c.135G>A vs c.135G>C, both p.Pro45=) are one
    page by the build's own slug dedupe, so dropping the second VCV loses no page.

    Returns (safe, dropped, uncovered) where `uncovered` names the carrier genes
    whose own pickle does not yet serve the dropped variants. The caller must
    still refuse when `safe` is False.
    """
    path = os.path.join(records_dir, f"{gene}.pkl.gz")
    stored = _stored_vcvs(path)
    if stored is None:
        return False, 0, {}
    new = {str(r.get("variation_id") or r.get("vcv")) for r in recs
           if (r.get("variation_id") or r.get("vcv")) is not None}
    dropped = [row for row in stored if row[0] not in new]
    if not dropped:
        return True, 0, {}
    by_gene = {}
    host = (gene or "").upper()
    own_rows = 0
    for v, g, p in dropped:
        if not g or g == host:
            own_rows += 1        # the gene's own lost variants: never a reshard
            continue
        by_gene.setdefault(g, set()).add((v, p))
    uncovered = {gene: own_rows} if own_rows else {}
    for g, want in by_gene.items():
        own = _stored_vcvs(os.path.join(records_dir, f"{g}.pkl.gz"))
        have_v = {v for v, _, _ in own} if own else set()
        have_p = {f"{gg}|{pp}" for _, gg, pp in own if pp} if own else set()
        missing = {v for v, p in want
                   if v not in have_v and not (p and f"{g}|{p}" in have_p)}
        if missing:
            uncovered[g] = len(missing)
    return not uncovered, len(dropped), uncovered


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
            safe, dropped, uncovered = carrier_reshard(a.records, g, recs)
            if not safe:
                print(f"  {g:<12} REFUSED {len(recs)} would replace {collapsed} "
                      f"(uncovered carrier {uncovered})", flush=True)
                prog.setdefault("failed", {})[g] = (
                    f"collapse {len(recs)}<{collapsed} uncovered {uncovered}")
                _save_json(a.progress, prog)
                continue
            print(f"  {g:<12} RESHARD  {collapsed:>7} -> {len(recs):>7} records "
                  f"({dropped} carrier rows covered by their own genes)  {time.time()-t:6.1f}s",
                  flush=True)
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
