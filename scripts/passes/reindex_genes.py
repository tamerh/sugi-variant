"""Reindex genes from their EXISTING pickles.

After a record rebuild the index can still carry the old canonical slug for a VCV
that is now in the pickle. `index_gene()` replaces exactly those rows and their
aliases, so the new slugs land without rebuilding a record.

Rows for VCVs that are NOT in the pickle are deliberately left alone. They are the
alias-only / stale-slug tail: deleting them would drop a URL that may still be
linked, and merging them into the page that does exist is a separate alias pass.
This pass therefore never makes a page less reachable.

    reindex_genes.py <records-dir> --genes=ATM,FLNC [--index PATH] [--apply]

No record is built and no pickle is written. Without --apply it reports what it
would touch. Safe to run twice.
"""
import argparse
import gzip
import os
import pickle
import sys

from sugivariant.index import index_gene, open_db


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("records", help="records directory")
    ap.add_argument("--genes", required=True, help="comma-separated symbols")
    ap.add_argument("--index", help="index.db (default: <records>/../index.db)")
    ap.add_argument("--apply", action="store_true", help="write the index")
    a = ap.parse_args(argv)

    index_path = a.index or os.path.join(os.path.dirname(os.path.abspath(a.records)),
                                         "index.db")
    genes = [g.strip().upper() for g in a.genes.split(",") if g.strip()]
    print(f"{'APPLY' if a.apply else 'DRY RUN'} -- {len(genes)} gene(s) -> {index_path}")

    conn = open_db(index_path) if a.apply else None
    for g in genes:
        path = os.path.join(a.records, f"{g}.pkl.gz")
        try:
            with gzip.open(path, "rb") as f:
                recs = pickle.load(f)
        except Exception as e:
            print(f"  {g:<12} SKIP cannot read pickle: {type(e).__name__}")
            continue
        if not isinstance(recs, list):
            print(f"  {g:<12} SKIP pickle is not a list")
            continue
        if conn is None:
            import sqlite3
            ro = sqlite3.connect(f"file:{index_path}?mode=ro", uri=True)
            old = ro.execute("SELECT COUNT(*) FROM variant WHERE gene=?", (g,)).fetchone()[0]
            ro.close()
            print(f"  {g:<12} would replace {old:>6} VCV rows, index {len(recs):>6} records")
            continue

        old = conn.execute("SELECT COUNT(*) FROM variant WHERE gene=?", (g,)).fetchone()[0]
        n = index_gene(conn, g, recs=recs)      # replaces the pickle's VCV rows + aliases
        print(f"  {g:<12} replaced rows for {n:>6} records (gene had {old:>6} rows)")
    if conn is not None:
        conn.commit()
        conn.close()
    print(("applied" if a.apply else "nothing written"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
