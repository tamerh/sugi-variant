"""Merge index rows of variants that no longer have a page into the page that does.

The build collapses ClinVar duplicates: two VCVs with the same protein change (or
the same coding change) become one page, and only the surviving VCV is written to
the pickle. ``index_gene`` replaces rows for the VCVs it writes, so the *dropped*
VCV keeps its old row and aliases. That row is a dead URL: it resolves to a gene,
the gene's pickle has no record for it, and the page 404s — even though the same
variant is served one slug away.

This pass re-points such a row's aliases at the surviving VCV and deletes the dead
row. A URL that used to 404 (for example an old coordinate slug) then resolves to
the real page. A row whose slug is served by nothing is left alone: that is a
genuine gap, not a duplicate, and deleting it would hide it.

    merge_stale_slugs.py --queue cache/rebuild-queue.json \
        --records cache/records --index cache/index.db [--apply]

Only genes the queue lists as ``reindex`` (index slugs missing from the pickle) are
touched, so this is a bounded post-rebuild pass, not a corpus scan.
"""
import argparse
import gzip
import json
import pickle
import sqlite3
import sys

from sugivariant.index import open_db


def _load(path):
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def _live_maps(recs):
    """Slug / protein / coding / coordinate -> the VCV that serves that variant."""
    by_slug, by_protein, by_c, by_coord = {}, {}, {}, {}
    for r in recs or []:
        v = r.get("variation_id") or r.get("vcv")
        if v is None:
            continue
        v = int(v)
        for s in (r.get("slugs") or [r.get("canonical_slug")]):
            if s:
                by_slug[s.lower()] = v
        if r.get("hgvs_p"):
            by_protein[r["hgvs_p"].strip()] = v
        if r.get("hgvs_c"):
            by_c[r["hgvs_c"].strip()] = v
        if r.get("coordinate"):
            by_coord[r["coordinate"].strip().upper()] = v
    return by_slug, by_protein, by_c, by_coord


def merge_gene(conn, gene, recs, apply=False):
    """Re-point and delete dead rows for one gene. Returns (merged, kept).

    A dropped VCV's row can carry an old coordinate slug, so the URL alone does not
    identify the survivor. The variant does: the build collapses duplicates on the
    protein change (the canonical slug is gene + p. form), so the row is matched to
    the live record sharing its ``hgvs_p`` first, then its coding change, then its
    allele-level coordinate. That is the same identity the build used to decide they
    were one page — matching on anything looser would merge two real variants.
    """
    by_slug, by_protein, by_c, by_coord = _live_maps(recs)
    if not by_slug:
        return [], []
    live_vcvs = set(by_slug.values())
    merged, kept = [], []
    rows = conn.execute(
        "SELECT vcv, slug, hgvs_p, hgvs_c, coordinate FROM variant WHERE gene=?",
        (gene,)).fetchall()
    for vcv, slug, hgvs_p, hgvs_c, coord in rows:
        if vcv in live_vcvs:
            continue
        new = (by_slug.get((slug or "").lower())
               or (by_protein.get((hgvs_p or "").strip()) if hgvs_p else None)
               or (by_c.get((hgvs_c or "").strip()) if hgvs_c else None)
               or (by_coord.get((coord or "").strip().upper()) if coord else None))
        if new is None:
            kept.append({"vcv": vcv, "slug": slug})
            continue
        merged.append({"vcv": vcv, "slug": slug, "into": new})
        if apply:
            # Every key that used to reach the dead VCV now reaches the survivor.
            conn.execute("UPDATE OR IGNORE alias SET vcv=? WHERE vcv=?", (new, vcv))
            conn.execute("DELETE FROM alias WHERE vcv=?", (vcv,))
            conn.execute("DELETE FROM variant WHERE vcv=?", (vcv,))
    if apply and merged:
        conn.execute(
            "UPDATE gene_meta SET n=(SELECT COUNT(*) FROM variant WHERE gene=?) "
            "WHERE gene=?", (gene, gene))
    return merged, kept


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--queue", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)

    q = json.load(open(a.queue, encoding="utf-8"))
    genes = [e["gene"] for e in (q.get("reindex") or []) if e.get("gene")]
    if a.limit:
        genes = genes[:a.limit]
    print(f"{'APPLY' if a.apply else 'DRY RUN'} -- {len(genes)} genes", flush=True)

    conn = open_db(a.index) if a.apply else sqlite3.connect(
        f"file:{a.index}?mode=ro", uri=True)
    total_merged = total_kept = 0
    for g in genes:
        try:
            recs = _load(f"{a.records}/{g}.pkl.gz")
        except Exception as e:
            print(f"  {g:<12} SKIP {type(e).__name__}", flush=True)
            continue
        merged, kept = merge_gene(conn, g, recs, apply=a.apply)
        total_merged += len(merged)
        total_kept += len(kept)
        if merged or kept:
            print(f"  {g:<12} merged {len(merged):>3}  left {len(kept):>3}", flush=True)
    if a.apply:
        conn.commit()
        conn.close()
    print(f"merged {total_merged} dead rows, left {total_kept} genuine gaps", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
