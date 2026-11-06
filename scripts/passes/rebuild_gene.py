"""Re-derive a gene's records with the current build pipeline.

A pass for `variant.sh preview` / `variant.sh pass`: it does not edit one field
in each record, it runs `build.enriched_records(gene)` -- the single record
assembly path the app itself uses -- and writes the result atomically with a
fresh data stamp.

    rebuild_gene.py <records-dir> --genes=TTN,FLG --apply

Without --apply it builds but writes nothing. It refuses a write that would
collapse an existing gene (records_collapsed), and it skips a gene that comes
back empty rather than replacing records with nothing.
"""
import argparse
import gzip
import os
import pickle
import sys
import time

from sugivariant.build import enriched_records
from sugivariant.enrich import (
    records_collapsed, stored_record_count, write_data_stamp,
)


def _write(records_path, recs):
    """Atomic write + stamp. The temp file lives beside the target, so a dev
    symlink is replaced (not written through) and prod's inode is untouched."""
    tmp = f"{records_path}.{os.getpid()}.tmp"
    with gzip.open(tmp, "wb") as f:
        pickle.dump(recs, f)
    os.replace(tmp, records_path)
    return write_data_stamp(records_path)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("records", help="records directory")
    ap.add_argument("--genes", default="", help="comma-separated symbols (default: all)")
    ap.add_argument("--apply", action="store_true", help="write the records")
    ap.add_argument("--limit", type=int, help="only the first N genes")
    a = ap.parse_args(argv)

    genes = [g.strip().upper() for g in a.genes.split(",") if g.strip()]
    if not genes:
        genes = sorted(f[:-len(".pkl.gz")] for f in os.listdir(a.records)
                       if f.endswith(".pkl.gz"))
    if a.limit:
        genes = genes[:a.limit]

    print(f"{'APPLY' if a.apply else 'DRY RUN'} -- {len(genes)} gene(s) against {a.records}")
    total_before = total_after = 0
    rebuilt = skipped = refused = 0
    t0 = time.time()
    for g in genes:
        path = os.path.join(a.records, f"{g}.pkl.gz")
        before = stored_record_count(path)
        t = time.time()
        try:
            recs = enriched_records(g) or []
        except Exception as e:
            print(f"  {g:<12} ERROR {type(e).__name__}: {str(e)[:90]}")
            skipped += 1
            continue
        dt = time.time() - t
        if not recs:
            print(f"  {g:<12} 0 records -- skipped (nothing written)")
            skipped += 1
            continue
        collapsed = records_collapsed(path, len(recs))
        if collapsed:
            print(f"  {g:<12} REFUSED {len(recs)} records would replace {collapsed}")
            refused += 1
            continue
        print(f"  {g:<12} {str(before):>7} -> {len(recs):>7} records  {dt:6.1f}s")
        if a.apply:
            stamped = _write(path, recs)
            if not stamped:
                print(f"  {g:<12} NOTE no data stamp (biobtree meta unreachable?)")
        total_before += before or 0
        total_after += len(recs)
        rebuilt += 1
    print(f"\n{rebuilt} rebuilt, {skipped} skipped, {refused} refused; "
          f"{total_before:,} -> {total_after:,} records; {time.time() - t0:.1f}s"
          + ("" if a.apply else "  (nothing written)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
