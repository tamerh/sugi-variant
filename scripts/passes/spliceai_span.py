"""Upgrade SpliceAI readouts stored as `absent` that are really "evaluated, below 0.2".

WHY. SpliceAI's precomputed release scores every possible single-nucleotide substitution
inside a gene's body, then keeps only the alleles whose largest delta reaches 0.2. So the
region it covers for a gene is ONE contiguous interval, sampled sparsely by that filter.
A position lying between the lowest and the highest position the release holds any record
for in the gene is therefore inside the covered interval: SpliceAI did evaluate it, and
the absence of all three alternate alleles means every delta of each fell under 0.2 --
negative evidence, not missing data. spliceai_for() now returns `below_threshold` for
that case; this pass applies the same test to records already on disk.

SAFETY. Only `absent` records are ever touched, and only ever upgraded. A gene whose
SpliceAI index comes back empty (biobtree hiccup, or the gene genuinely has none) yields
an empty span, every test fails, and nothing changes -- a transport failure can only cost
an upgrade, never cause a wrong one. Files are rewritten atomically. --apply is required
to write; without it this only counts.

Usage:  spliceai_span.py <records-dir> [--apply] [--workers N] [--limit N]
"""
import concurrent.futures as cf
import glob
import gzip
import json
import os
import pickle
import sys
import tempfile
import threading
import time

# The repo root, two levels up, so the pass runs from any checkout rather than
# naming one machine's layout.
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
from sugivariant.enrich import gene_spliceai, spliceai_for   # noqa: E402

RECORDS = sys.argv[1]
APPLY = "--apply" in sys.argv
WORKERS = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--workers=")), 4))
LIMIT = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--limit=")), 0))
CKPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "spliceai_span.done")

lock = threading.Lock()
tot = {"genes": 0, "files_written": 0, "absent": 0, "upgraded": 0,
       "still_absent": 0, "no_index": 0, "errors": 0, "no_coord": 0, "would_score": 0}
t0 = time.time()


def done_set():
    if not os.path.exists(CKPT):
        return set()
    with open(CKPT) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def one(path):
    out = {"genes": 1, "absent": 0, "upgraded": 0, "still_absent": 0,
           "no_index": 0, "errors": 0, "no_coord": 0, "files_written": 0, "would_score": 0}
    try:
        with gzip.open(path, "rb") as fh:
            recs = pickle.load(fh)
    except Exception:
        out["errors"] = 1
        return path, out
    absents = [r for r in recs if (r.get("spliceai") or {}).get("status") == "absent"]
    out["absent"] = len(absents)
    if not absents:
        return path, out                      # no biobtree call for a gene with nothing to fix
    # Bracket each record against ITS OWN gene's index, not the file's. A record file is
    # keyed by one gene but carries the overlapping genes' variants too (a WRAP53 missense
    # is filed under TP53), and the span that decides this question is the span of the gene
    # the variant is actually in. gene_spliceai() is not memoised -- the caller holds the
    # result -- so keep one cache per file; in practice that is one entry.
    caches, changed = {}, 0
    for r in absents:
        hgnc, coord = r.get("hgnc_id"), r.get("coordinate")
        if not coord:
            out["no_coord"] += 1
            continue
        if not hgnc:
            out["no_index"] += 1
            continue
        if hgnc not in caches:
            try:
                caches[hgnc] = gene_spliceai(hgnc)
            except Exception:
                caches[hgnc] = None
                out["errors"] += 1
        cache = caches[hgnc]
        if not (cache or {}).get("span"):
            out["no_index"] += 1              # nothing to bracket against: leave as-is
            continue
        new = spliceai_for(coord, cache)
        st = (new or {}).get("status")
        if st == "below_threshold":
            r["spliceai"] = new
            changed += 1
        else:
            # Deliberately narrow: this pass only ever turns `absent` into
            # `below_threshold`. Against the record's own gene index a few of these come
            # back `scored` or `other_allele` -- i.e. the original build used the FILE
            # gene's index and missed a real score. That is the separate carrier-gene
            # defect; fixing it here would quietly widen this change, so they are counted
            # and left alone.
            out["still_absent"] += 1
            if st in ("scored", "other_allele"):
                out["would_score"] = out.get("would_score", 0) + 1
    out["upgraded"] = changed
    if changed and APPLY:
        d = os.path.dirname(path)
        fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
        os.close(fd)
        try:
            with gzip.open(tmp, "wb") as fh:
                pickle.dump(recs, fh, protocol=pickle.HIGHEST_PROTOCOL)
            os.chmod(tmp, 0o644)      # mkstemp gives 0600; records must stay world-readable
            os.replace(tmp, path)
            out["files_written"] = 1
        except Exception:
            out["errors"] = 1
            if os.path.exists(tmp):
                os.unlink(tmp)
    return path, out


def main():
    only = next((a.split("=")[1] for a in sys.argv if a.startswith("--genes=")), "")
    if only:
        files = [os.path.join(RECORDS, g.strip() + ".pkl.gz") for g in only.split(",")]
        files = [f for f in files if os.path.exists(f)]
    else:
        files = sorted(glob.glob(os.path.join(RECORDS, "*.pkl.gz")))
    # NOTE on the dev tree: every <GENE>.pkl.gz there is a SYMLINK into the prod records
    # directory, so dev and prod read the same bytes. os.replace() below swaps the dev
    # DIRECTORY ENTRY, which breaks that symlink and leaves prod's file untouched -- which
    # is what makes a dev-first preview possible at all. Restoring is `ln -sfn` back.
    skip = done_set() if (APPLY and not only) else set()
    files = [f for f in files if f not in skip]
    # --limit takes the head of an alphabetical list, which is a biased timing sample
    # (gene size is not uniform across the alphabet). --sample takes a seeded random one.
    SAMPLE = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--sample=")), 0))
    if SAMPLE:
        import random
        random.Random(17).shuffle(files)
        files = files[:SAMPLE]
    elif LIMIT:
        files = files[:LIMIT]
    print(f"{'APPLY' if APPLY else 'DRY RUN'}  {len(files):,} gene files in {RECORDS}"
          f"  ({len(skip):,} already done)  workers={WORKERS}", flush=True)
    ck = open(CKPT, "a") if (APPLY and not only) else None
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for i, (path, out) in enumerate(ex.map(one, files), 1):
            with lock:
                for k, v in out.items():
                    tot[k] += v
                if ck:
                    ck.write(path + "\n")
                    if i % 200 == 0:
                        ck.flush()
            if i % 250 == 0 or i == len(files):
                el = time.time() - t0
                print(f"{i:6,}/{len(files):,}  {el/60:6.1f}m  "
                      f"absent {tot['absent']:,}  upgraded {tot['upgraded']:,} "
                      f"({100*tot['upgraded']/max(tot['absent'],1):.1f}%)  "
                      f"still absent {tot['still_absent']:,}  no index {tot['no_index']:,}  "
                      f"written {tot['files_written']:,}  err {tot['errors']}", flush=True)
    if ck:
        ck.close()
    print("---")
    print(json.dumps(tot, indent=2))


main()
