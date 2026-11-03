"""Give a misfiled record its OWN gene's structural context instead of the file gene's.

WHY. A record file is keyed by one gene and carries other genes' variants too, because a
build enumerates a gene's ClinVar set and biobtree's symbol search crosses previous-symbol
collisions. collect.py then does `rec["structure"] = ctx["structure"]` for every record in
the build, so one gene's UniProt features, PDB entries and AlphaFold readout get stamped
onto all of them.

The result is a confident, sourced-looking, wrong statement. Live examples before this pass:

  USE1.pkl.gz      36 USE1 records (HGNC:30882, 4 intervals) + 10 UBE2Z (HGNC:25847, 25).
                   UBE2Z's features were on all 46, so USE1 variants read
                   "the ubc core domain" -- a UBE2Z domain.
  SLC19A1.pkl.gz   236 SLC19A1 + 189 RFC1 + 96 COL18A1. The RFC1 records are genuine
                   replication-factor-C variants (NM_002913, chr4; SLC19A1 is chr21),
                   correctly labelled and simply misfiled.

Nothing looks broken, because the residue number lands inside a real interval. It just
describes a different protein. structure_context.py counted these at 31,052 records and
deliberately left them alone, because it had a rule never to remove a feature. Removing is
the whole point here.

THE SAFETY RULE. `hgnc_id` cannot be trusted on its own. In SLC19A1.pkl.gz the 96 COL18A1
records carry HGNC:10937, which is SLC19A1's id, not COL18A1's -- so reading features at
that id would swap one wrong protein for another. A record is therefore only rewritten when
its `hgnc_id`'s APPROVED SYMBOL matches its own `gene_symbol`. Over a 1,200-gene sample all
5,697 removal candidates passed that check, so the guard costs nothing and bounds the
blast radius if an inconsistent record ever becomes one.

A gene whose fetch fails or returns no intervals is skipped, so a transport failure can
never strip a record.

WHAT IT WRITES. Two things, both per record and both from the record's own gene:
  · `structural` -- recomputed, and allowed to LOSE features, which is the fix.
  · `structure`  -- in a MIXED file only, each gene group gets its own dict, so the
    Structure row stops showing a neighbour's pLDDT and PDB entries. One dict per gene
    per file, not per record, so sharing and file size are preserved.

Usage:  carrier_structure.py <records-dir> [--apply] [--workers N] [--genes=A,B] [--sample N]
"""
import concurrent.futures as cf
import glob
import gzip
import json
import os
import pickle
import re
import sys
import tempfile
import threading
import time

# The repo root, two levels up, so the pass runs from any checkout rather than
# naming one machine's layout.
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
from sugivariant.build import _approved_symbols                      # noqa: E402
from sugivariant.enrich import gene_structure, structural_context    # noqa: E402

RECORDS = sys.argv[1]
APPLY = "--apply" in sys.argv
WORKERS = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--workers=")), 4))
CKPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "carrier_structure.done")

_AM_SHORT_POS = re.compile(r"^[A-Z](\d+)[A-Z]$")

lock = threading.Lock()
tot = {"genes": 0, "files_written": 0, "records": 0, "fixed": 0, "features_removed": 0,
       "structure_rehomed": 0, "unsafe_hgnc": 0, "no_hgnc": 0, "no_intervals": 0,
       "errors": 0, "already_own_gene": 0}
t0 = time.time()

# One approved-symbol lookup per hgnc id, shared across worker threads. dict writes are
# atomic under the GIL and a duplicate fetch is harmless, so no lock is needed.
_SYMS = {}


def approved(hgnc):
    if hgnc not in _SYMS:
        try:
            _SYMS[hgnc] = {s.upper() for s in (_approved_symbols(hgnc) or [])}
        except Exception:
            _SYMS[hgnc] = set()
    return _SYMS[hgnc]


def trustworthy(rec):
    """True when the record's hgnc_id and its own gene_symbol name the same gene.

    This is the guard that keeps the fix from swapping one wrong protein for another:
    SLC19A1.pkl.gz holds COL18A1-labelled records carrying SLC19A1's HGNC id."""
    hg, sym = rec.get("hgnc_id"), rec.get("gene_symbol")
    return bool(hg and sym and sym.upper() in approved(hg))


def is_foreign(rec, gs):
    """True when the record's stored structure block describes a DIFFERENT protein.

    The accession settles it when both sides have one. pLDDT is only a fallback, because
    it is not identity: C7's stored block carries the same accession as its own gene
    (P10643) but a pLDDT of 89.67 against the current 80.11, which is upstream drift
    since the build, not a misfiled record. Treating that as foreign pulled 664 C7
    records into a pass about carriers. Refreshing drifted values is a separate job.
    """
    cur = (rec.get("structure") or {}).get("alphafold") or {}
    own = (gs.get("alphafold") or {}) or {}
    if cur.get("acc") and own.get("acc"):
        return cur["acc"] != own["acc"]
    return cur.get("plddt") != own.get("plddt")


def done_set():
    if not os.path.exists(CKPT):
        return set()
    with open(CKPT) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def _canonical_pos(r):
    """Mirrors collect.py: an AlphaMissense position counts as canonical only when a
    SaProt score proves it, otherwise structural_context is handed None and abstains."""
    if not (r.get("am_isoform_mismatch") and r.get("saprot")):
        return None
    m = _AM_SHORT_POS.match(((r.get("alphamissense") or {}).get("short")) or "")
    return int(m.group(1)) if m else None


def one(path):
    out = {"genes": 1, "records": 0, "fixed": 0, "features_removed": 0,
           "structure_rehomed": 0, "unsafe_hgnc": 0, "no_hgnc": 0, "no_intervals": 0,
           "errors": 0, "files_written": 0, "already_own_gene": 0}
    try:
        with gzip.open(path, "rb") as fh:
            recs = pickle.load(fh)
    except Exception:
        out["errors"] = 1
        return path, out
    if not recs:
        return path, out

    caches, changed = {}, 0
    for r in recs:
        out["records"] += 1
        hg = r.get("hgnc_id")
        if not hg:
            out["no_hgnc"] += 1
            continue
        if not trustworthy(r):
            out["unsafe_hgnc"] += 1
            continue
        if hg not in caches:
            try:
                caches[hg] = gene_structure(hg)
            except Exception:
                caches[hg] = None
                out["errors"] += 1
        gs = caches[hg]
        if not gs:
            continue
        intervals = gs.get("intervals") or []
        if not intervals:
            out["no_intervals"] += 1       # nothing to read: leave the record alone
            continue
        if not is_foreign(r, gs):
            # Its stored context already belongs to its own gene. Any difference from a
            # fresh read is upstream drift, which this pass deliberately does not chase.
            out["already_own_gene"] += 1
            continue

        before = r.get("structural")
        new_ctx = structural_context(r.get("hgvs_p"), intervals,
                                     r.get("am_isoform_mismatch"), _canonical_pos(r))
        lost = set((before or {}).get("features") or []) - \
               set((new_ctx or {}).get("features") or [])
        if new_ctx != before:
            r["structural"] = new_ctx
            out["fixed"] += 1
            out["features_removed"] += len(lost)
            changed += 1
        # and the block itself, so the Structure row stops reporting a neighbour's pLDDT
        # and PDB entries. One dict per gene per file, not per record, so pickle's sharing
        # and the file size are preserved.
        r["structure"] = gs
        out["structure_rehomed"] += 1
        changed += 1

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
    # On the dev tree each <GENE>.pkl.gz is a SYMLINK into prod; os.replace() swaps the dev
    # directory entry and leaves prod's file untouched, which is what makes a preview work.
    skip = done_set() if (APPLY and not only) else set()
    files = [f for f in files if f not in skip]
    SAMPLE = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--sample=")), 0))
    if SAMPLE:
        import random
        random.Random(17).shuffle(files)
        files = files[:SAMPLE]
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
                      f"records {tot['records']:,}  fixed {tot['fixed']:,}  "
                      f"features removed {tot['features_removed']:,}  "
                      f"structure rehomed {tot['structure_rehomed']:,}  "
                      f"own-gene {tot['already_own_gene']:,}  unsafe {tot['unsafe_hgnc']:,}  "
                      f"written {tot['files_written']:,}  err {tot['errors']}", flush=True)
    if ck:
        ck.close()
    print("---")
    print(json.dumps(tot, indent=2))


main()
