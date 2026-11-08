"""Reconcile gene_meta.n with the rows the variant table actually holds.

The invariant: `gene_meta[G].n == COUNT(variant WHERE gene=G)`. `index_gene`
maintains it for the symbols it writes, but a rebuild can move rows between
symbols (the carrier-gene defect) or drop a duplicate VCV whose row was filed under
another symbol, and then the count drifts. The live-index test asserts this
invariant, and a drifted count also feeds "known empty" detection and the gene
hubs, so a refresh has to leave it exact.

    reconcile_gene_meta.py --index cache/index.db [--apply]
"""
import argparse
import sqlite3
import sys


def reconcile(conn, apply=False):
    """Return (drifted, fixed) lists of (gene, claimed, actual)."""
    actual = dict(conn.execute(
        "SELECT gene, COUNT(*) FROM variant WHERE gene IS NOT NULL GROUP BY gene").fetchall())
    claimed = dict(conn.execute("SELECT gene, n FROM gene_meta").fetchall())
    drifted = []
    for gene in set(actual) | set(claimed):
        a, c = actual.get(gene, 0), claimed.get(gene)
        if c is None or c != a:
            drifted.append((gene, c, a))
    if apply:
        conn.executemany("UPDATE gene_meta SET n=? WHERE gene=?",
                         [(a, g) for g, _c, a in drifted])
        conn.commit()
    return drifted


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)
    conn = (sqlite3.connect(a.index) if a.apply
            else sqlite3.connect(f"file:{a.index}?mode=ro", uri=True))
    conn.execute("PRAGMA busy_timeout=120000")
    drifted = reconcile(conn, apply=a.apply)
    for gene, claimed, actual in sorted(drifted, key=lambda r: -(abs((r[2] or 0) - (r[1] or 0))))[:20]:
        print(f"  {gene:<14} gene_meta {claimed!s:>7} -> {actual:>7}")
    print(f"{'repaired' if a.apply else 'would repair'} {len(drifted)} gene_meta rows")
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
