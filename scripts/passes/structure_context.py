"""Add membrane topology to the stored structural context, and the AlphaFold accession.

WHY. `structural_context()` reports the UniProt features overlapping a variant's residue,
but `_UF_KEEP` never listed UniProt's membrane-topology types, so the two facts that
matter most for a channel or transporter were dropped at fetch time: which transmembrane
helix a residue sits in, and which side of the membrane it faces. Probed on staging
2026-11-03, TMC1 (HGNC:16513) returns 10 `transmembrane region` and 11 `topological
domain` features with clean intervals and cytoplasmic/extracellular descriptions. Over 40
sampled genes, 25% carry them and 7.5% of residue-bearing records gain a structural line
they did not have, taking coverage from 49.5% to 57.0%.

It also stores `structure.alphafold.acc`, the UniProt accession, which the alphafold row
already carries as its own id. Without it a page can only link the AlphaFold model for a
MISSENSE variant, because the accession otherwise lives only on the AlphaMissense block.

The phrase deliberately does not number the helices. UniProt describes each one only as
"helical", so calling one TM4 would be this pass counting intervals, and literature TM
numbering does not reliably follow UniProt's order. The interval carries the detail.

SAFETY. The recompute replicates collect.py's call EXACTLY, including the isoform guard:
where ClinVar's numbering disagrees with UniProt's, `structural_context` reads the
AlphaMissense residue only when a SaProt score proves that position is canonical, and
otherwise reports nothing. 5.6% of records are flagged that way, so getting this wrong
would put features on the wrong residue at scale.

Each record is bracketed against ITS OWN gene, never the file's. A record file is keyed
by one gene and carries the overlapping genes' variants too: USE1.pkl.gz holds 36 USE1
records (HGNC:30882, 4 UniProt intervals) and 10 UBE2Z records (HGNC:25847, 25), so the
file gene's features describe the wrong protein for part of the file.

That also means this pass REFUSES to remove a feature. A record whose own gene yields
fewer features than the stored value almost always means the stored value was computed
from the file gene and is describing the wrong protein. Real, and visible on prod today
(USE1 variants currently show "the ubc core domain", which is UBE2Z's), but it is the
separate carrier-gene defect. Fixing it here would widen this pass from "add membrane
topology" into "rewrite structural context", so those records are counted as `regressed`
and left exactly as they are, the same choice spliceai_span makes for its own carrier
cases. Measured at 2.1% of records over a 400-gene sample. A gene whose fetch fails or
returns nothing is likewise skipped, so a transport failure can only cost an addition.

A NOTE FOR WHOEVER COMES NEXT. This writes the DERIVED field (`structural`) and not the
raw `structure.intervals` the features were read from, so on a record written before this
pass the stored intervals still lack the two new types while the derived line reflects
them. Nothing reads stored intervals at serve time, so no page is affected, and a rebuilt
gene gets both. But do not treat stored intervals as current: fetch them, as this pass
does. Bringing them into line would rewrite nearly every file for no visible change, so it
is left for whenever a pass touches them for a reason of its own.

Files are rewritten atomically. --apply is required to write; without it this only counts.

Usage:  structure_context.py <records-dir> [--apply] [--workers N] [--genes=A,B] [--sample N]
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
from sugivariant.enrich import gene_structure, structural_context   # noqa: E402

RECORDS = sys.argv[1]
APPLY = "--apply" in sys.argv
WORKERS = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--workers=")), 4))
CKPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "structure_context.done")

# collect.py's own pattern for reading AlphaMissense's canonical residue out of "A462T".
_AM_SHORT_POS = re.compile(r"^[A-Z](\d+)[A-Z]$")

lock = threading.Lock()
tot = {"genes": 0, "files_written": 0, "records": 0, "gained": 0, "changed": 0,
       "regressed": 0, "no_hgnc": 0, "no_intervals": 0, "errors": 0, "acc_set": 0,
       "mixed_file": 0}
t0 = time.time()


def done_set():
    if not os.path.exists(CKPT):
        return set()
    with open(CKPT) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def _canonical_pos(r):
    """The residue to read features at when ClinVar's numbering is not UniProt's.

    Mirrors collect.py: only an AlphaMissense position BACKED BY A SAPROT SCORE counts as
    proof that the position resolved against the current canonical sequence. Without that
    proof structural_context() is handed None and reports nothing at all."""
    if not (r.get("am_isoform_mismatch") and r.get("saprot")):
        return None
    m = _AM_SHORT_POS.match(((r.get("alphamissense") or {}).get("short")) or "")
    return int(m.group(1)) if m else None


def one(path):
    out = {"genes": 1, "records": 0, "gained": 0, "changed": 0, "regressed": 0,
           "no_hgnc": 0, "no_intervals": 0, "errors": 0, "files_written": 0,
           "acc_set": 0, "mixed_file": 0}
    try:
        with gzip.open(path, "rb") as fh:
            recs = pickle.load(fh)
    except Exception:
        out["errors"] = 1
        return path, out
    if not recs:
        return path, out

    # Bracket each record against ITS OWN gene, not the file's. A record file is keyed by
    # one gene but carries the overlapping genes' variants too, and the features that
    # describe a residue are the features of the gene the variant is actually in.
    # USE1.pkl.gz holds 36 USE1 (HGNC:30882, 4 intervals) and 10 UBE2Z (HGNC:25847, 25
    # intervals) records; using the file gene for all of them is exactly the defect this
    # pass would otherwise spread. gene_structure() is not memoised, so keep one cache per
    # file -- in practice one or two entries.
    caches, changed = {}, 0
    for r in recs:
        out["records"] += 1
        hgnc = r.get("hgnc_id")
        if not hgnc:
            out["no_hgnc"] += 1
            continue
        if hgnc not in caches:
            try:
                caches[hgnc] = gene_structure(hgnc)
            except Exception:
                caches[hgnc] = None
                out["errors"] += 1
        gs = caches[hgnc]
        intervals = (gs or {}).get("intervals") or []
        if not intervals:
            out["no_intervals"] += 1     # nothing to bracket against: leave as-is
            continue

        before = r.get("structural")
        new_ctx = structural_context(r.get("hgvs_p"), intervals,
                                     r.get("am_isoform_mismatch"), _canonical_pos(r))
        bf = set((before or {}).get("features") or [])
        nf = set((new_ctx or {}).get("features") or [])
        if bf - nf:
            # This record would LOSE a feature. Against its own gene's intervals that
            # usually means the stored value was computed from the FILE gene's features
            # and is describing the wrong protein -- a real defect, but the separate
            # carrier-gene one. Fixing it here would quietly widen this pass from "add
            # membrane topology" into "rewrite structural context", so it is counted and
            # left alone, the same choice spliceai_span makes for its own carrier cases.
            out["regressed"] += 1
            continue
        if new_ctx != before:
            r["structural"] = new_ctx
            changed += 1
            out["changed"] += 1
            if not before:
                out["gained"] += 1

    # The accession, so a non-missense page can link the model too. collect.py does
    # `rec["structure"] = ctx["structure"]`, so every record in a gene file shares ONE
    # dict object (verified: 8,749 MSH6 records, 1 object) and pickle preserves the
    # sharing, so this costs ~20 bytes per FILE rather than per record. The flip side is
    # that a MIXED file cannot hold per-record accessions without un-sharing it, so the
    # write happens only where the whole file is one gene; elsewhere the AlphaMissense
    # block still carries the accession for missense records.
    hgncs = {r.get("hgnc_id") for r in recs if r.get("hgnc_id")}
    if len(hgncs) == 1:
        acc = ((caches.get(next(iter(hgncs))) or {}).get("alphafold") or {}).get("acc")
        st = recs[0].get("structure")
        if acc and isinstance(st, dict) and isinstance(st.get("alphafold"), dict) \
                and st["alphafold"].get("acc") != acc:
            st["alphafold"]["acc"] = acc
            changed += 1
            out["acc_set"] += 1
    else:
        out["mixed_file"] = 1

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
    # directory. os.replace() swaps the dev DIRECTORY ENTRY, breaking that symlink and
    # leaving prod's file untouched, which is what makes a dev-first preview possible.
    skip = done_set() if (APPLY and not only) else set()
    files = [f for f in files if f not in skip]
    SAMPLE = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--sample=")), 0))
    if SAMPLE:
        import random
        random.Random(17).shuffle(files)       # seeded: gene size is not uniform A-Z
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
                      f"records {tot['records']:,}  changed {tot['changed']:,}  "
                      f"gained {tot['gained']:,}  acc {tot['acc_set']:,}  "
                      f"written {tot['files_written']:,}  no-intervals {tot['no_intervals']:,}  "
                      f"carrier-skipped {tot['regressed']:,}  mixed {tot['mixed_file']:,}  "
                      f"err {tot['errors']}", flush=True)
    if ck:
        ck.close()
    print("---")
    print(json.dumps(tot, indent=2))


main()
