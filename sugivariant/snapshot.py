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
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone

SCHEMA = "sugi-variant-release/1"
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


def _read_pickle(records_path):
    """(n_records, n_distinct_slugs, error). Read once; no app import."""
    try:
        with gzip.open(records_path, "rb") as f:
            recs = pickle.load(f)
        if not isinstance(recs, list):
            return None, None, "not-a-list"
        slugs = {r.get("canonical_slug") for r in recs if r.get("canonical_slug")}
        return len(recs), len(slugs), None
    except Exception as e:                       # corrupt / truncated / not a pickle
        return None, None, type(e).__name__


def _count_records(records_path):
    """(count, error). Counts a pickle without importing the app."""
    n, _slugs, err = _read_pickle(records_path)
    return n, err


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
            "n_slugs": None,
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
        n, n_slugs, err = _read_pickle(os.path.join(records_dir, found[gene]["path"]))
        found[gene]["n"] = n
        found[gene]["n_slugs"] = n_slugs
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
        slug_counts = {
            r["gene"]: r["n"]
            for r in conn.execute(
                "SELECT gene, COUNT(DISTINCT slug) AS n FROM variant GROUP BY gene")
        }
        total = conn.execute("SELECT COUNT(*) FROM variant").fetchone()[0]
        flagged = conn.execute("SELECT COUNT(*) FROM variant WHERE flag IS NOT NULL").fetchone()[0]
    finally:
        conn.close()
    return {
        "gene_meta": gene_meta,
        "variant_counts": variant_counts,
        "slug_counts": slug_counts,
        "total": total,
        "flagged": flagged,
    }


def _stale_from(stamp, current):
    if not current or not isinstance(stamp, dict):
        return []
    return sorted(d for d, v in current.items() if stamp.get(d) != v)


def _page_diff(conn, records_dir, gene):
    """Compare the INDEX's distinct slugs with the PICKLE's canonical slugs.

    Row counts are not pages: two ClinVar VCVs can share one canonical slug, and an
    index built before a slug-shape fix can still carry a coordinate slug for a page
    that now has a gene-HGVS slug. This separates four things:

      missing_slugs  index pages with no pickle page and no alias pointing at one
                     (the real rebuild gap)
      stale_slugs    index slugs whose VCV's alias keys already name a pickle page
                     (the page exists; reindex, do not rebuild)
      extra_slugs    pickle pages the index does not list (a surplus)
      alias_only     index VCVs whose slug IS in the pickle but whose VCV is not
                     (duplicate VCVs resolving to an existing page)
    """
    p = os.path.join(records_dir, f"{gene}{RECORD_SUFFIX}")
    pkl_slugs, pkl_vcvs = set(), set()
    if os.path.exists(p):
        try:
            with gzip.open(p, "rb") as f:
                recs = pickle.load(f)
            if isinstance(recs, list):
                pkl_slugs = {r.get("canonical_slug") for r in recs if r.get("canonical_slug")}
                pkl_vcvs = {str(r.get("variation_id")) for r in recs
                            if r.get("variation_id") is not None}
        except Exception:
            pass
    rows = conn.execute("SELECT vcv, slug FROM variant WHERE gene=?", (gene,)).fetchall()
    idx_slugs = {r[1] for r in rows if r[1]}
    slug_to_vcv = {r[1]: r[0] for r in rows if r[1]}
    missing = sorted(idx_slugs - pkl_slugs)
    extra = sorted(pkl_slugs - idx_slugs)
    alias_only = sum(1 for slug, vcv in ((r[1], r[0]) for r in rows)
                     if slug in pkl_slugs and str(vcv) not in pkl_vcvs)
    stale, true_missing = [], []
    if missing and pkl_slugs:
        # One batched alias lookup for the missing slugs' VCVs: a key that is already
        # a pickle page means the index slug is stale, not that the page is absent.
        by_vcv = defaultdict(set)
        vcvs = [slug_to_vcv[s] for s in missing]
        qs = ",".join("?" * len(vcvs))
        for vcv, key in conn.execute(f"SELECT vcv, key FROM alias WHERE vcv IN ({qs})", vcvs):
            by_vcv[str(vcv)].add(key)
        for s in missing:
            v = str(slug_to_vcv.get(s))
            if by_vcv.get(v, set()) & pkl_slugs:
                stale.append(s)
            else:
                true_missing.append(s)
    else:
        true_missing = missing
    return {"index_slugs": len(idx_slugs), "missing_slugs": true_missing,
            "stale_slugs": stale, "extra_slugs": extra, "alias_only": alias_only}


def _git(*args):
    try:
        return subprocess.check_output(
            ["git", *args], cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


def _release_meta(release=None, released_at=None, variant_version=None,
                  variant_commit=None, biobtree_version=None, biobtree_commit=None):
    """The release identity stamped into the manifest and served pages.

    Precedence: explicit argument, VARIANT_BUILD_VERSION, git describe, 'dev'.
    The private release script sets the env vars once at build start, exactly as
    Atlas's stamp_version() does, so every artifact of one run shares the stamp.
    """
    rel = (release or os.environ.get("VARIANT_BUILD_VERSION")
           or _git("describe", "--tags", "--always") or "dev")
    return {
        "release": rel,
        "released_at": (released_at or _utc_now()[:10]),
        "variant_version": variant_version or rel,
        "variant_commit": (variant_commit or os.environ.get("VARIANT_BUILD_COMMIT")
                           or _git("rev-parse", "--short", "HEAD")),
        "biobtree_version": biobtree_version or os.environ.get("BIOBTREE_VERSION"),
        "biobtree_commit": biobtree_commit or os.environ.get("BIOBTREE_COMMIT"),
    }


def build_manifest(records_dir, index_db, current_versions=None,
                   count_records=False, count_sample=None,
                   release=None, released_at=None, variant_version=None,
                   variant_commit=None, biobtree_version=None, biobtree_commit=None):
    """(manifest, report). Pure read; no writes, no biotree."""
    records = scan_records(records_dir, count_records=count_records,
                           count_sample=count_sample)
    index = read_index(index_db)
    pd_conn = sqlite3.connect(f"file:{index_db}?mode=ro", uri=True)
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
                     "index_rows": idx_n, "stamp": None, "stale_on": [],
                     "status": status, "issues": issues}
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
                     "n": n, "index_rows": idx_n, "stamp": r["stamp"],
                     "stale_on": stale_on or [],
                     "status": status, "issues": issues}
        idx_slugs = index["slug_counts"].get(g, 0)
        state["index_slugs"] = idx_slugs
        countable = r is None or (r["count_error"] is None and r["n"] is not None)
        if idx_slugs and countable:
            diff = _page_diff(pd_conn, records_dir, g)
            state.update(diff)
            state["pages_missing"] = len(diff["missing_slugs"])
            state["pages_stale"] = len(diff["stale_slugs"])
            state["rows_missing"] = max(0, idx_n - (r["n"] if r else 0))
        else:
            state["missing_slugs"] = None
            state["stale_slugs"] = []
            state["extra_slugs"] = []
            state["alias_only"] = 0
            state["pages_missing"] = 0 if idx_slugs == 0 else None
            state["pages_stale"] = 0
            state["rows_missing"] = None
        gene_state[g] = state
    pd_conn.close()

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
    pages_missing_total = sum(s.get("pages_missing") or 0 for s in gene_state.values())
    pages_stale_total = sum(s.get("pages_stale") or 0 for s in gene_state.values())
    alias_only_rows = sum(s.get("alias_only") or 0 for s in gene_state.values())
    status = "ok" if not (missing or corrupt or mismatch or records_only) else "partial"

    meta = _release_meta(release, released_at, variant_version, variant_commit,
                         biobtree_version, biobtree_commit)
    manifest = {
        "schema": SCHEMA,
        "release": meta["release"],
        "released_at": meta["released_at"],
        "variant_version": meta["variant_version"],
        "variant_commit": meta["variant_commit"],
        "biobtree_version": meta["biobtree_version"],
        "biobtree_commit": meta["biobtree_commit"],
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
                "index_slugs": sum(index["slug_counts"].values()),
                "pages_missing": pages_missing_total,
                "pages_stale": pages_stale_total,
                "alias_only_rows": alias_only_rows,
                "counted": served_all_counted,
                "datasets": dict(current_versions or {}),
                "scope": "all-cached-genes",
                "status": status,
                "gene_state": gene_state,
            }
        },
    }
    report = {
        "release": manifest["release"],
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
            "pages_missing": pages_missing_total,
            "pages_stale": pages_stale_total,
            "alias_only_rows": alias_only_rows,
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


def rebuild_queue(manifest):
    """The rebuilds needed to make the snapshot complete, counted in PAGES.

    A page is a distinct canonical slug. Two ClinVar VCVs can share one slug, and a
    stale index can carry a coordinate slug for a page that now has a gene-HGVS
    slug, so a row-count deficit overstates the work. The queue separates:

      pages_behind   index slugs with no pickle record (the real gap)
      no_record      no pickle at all; every index slug is missing
      unbuildable    a rebuild was attempted and produced no buildable variant
                     (no enumerable ClinVar set). Listed separately: re-running can
                     never fix it, so it must not sit in the work queue forever.
      alias_only     index VCVs whose slug already exists in the pickle
                     (duplicates resolving to an existing page -- reindex, no rebuild)

    Freshness is NOT a rebuild reason for ``queue``: a dataset version bump marks
    almost every gene stale, and mixing that into the page gaps would drown them.
    It is reported separately in ``refresh`` -- the genes whose stamp is behind the
    current dataset versions, with the datasets that moved. That is a whole-corpus
    refresh (its own scheduled run), and it is the only way new variants in an
    existing gene appear: a gene whose pickle and index agree is invisible to the
    page diff even when its data is two ClinVar releases old.

    Requires counted records (``--count-records`` or ``--count-sample``).
    Uncounted genes with index rows are reported in ``uncounted`` rather than
    guessed at.
    """
    corpus = ((manifest or {}).get("corpora") or {}).get(CORPUS) or {}
    state = corpus.get("gene_state") or {}
    # Genes a rebuild proved unbuildable. The audit cannot know this without building
    # them, so the queue runner records it and the queue honours it -- otherwise those
    # genes are re-attempted on every refresh and inflate "missing pages" with pages
    # that cannot exist.
    unbuildable = set()
    _rdir = corpus.get("records_dir") or (manifest or {}).get("records_dir")
    if _rdir:
        try:
            with open(os.path.join(_rdir, "unbuildable.json"), encoding="utf-8") as f:
                unbuildable = set(json.load(f).get("genes") or [])
        except Exception:
            unbuildable = set()
    entries, alias_entries, reindex_entries, uncounted, unbuildable_entries = [], [], [], [], []
    for g, s in state.items():
        if s.get("status") in ("known_empty", "orphan"):
            continue
        idx_slugs = s.get("index_slugs") or 0
        n = s.get("n")
        pages = s.get("pages_missing")
        pages_stale = s.get("pages_stale") or 0
        rows = s.get("rows_missing")
        alias = s.get("alias_only") or 0
        stale = "stale" in (s.get("issues") or [])
        base = {"gene": g, "records": n, "index_rows": s.get("index_rows") or 0,
                "index_slugs": idx_slugs, "pages_stale": pages_stale, "stale": stale}
        if s.get("status") == "missing" and idx_slugs > 0 and g in unbuildable:
            unbuildable_entries.append({**base, "reason": "unbuildable", "records": 0,
                                        "pages": pages if pages is not None else idx_slugs,
                                        "rows": rows if rows is not None else (s.get("index_rows") or 0),
                                        "alias_only": 0})
        elif s.get("status") == "missing" and idx_slugs > 0:
            entries.append({**base, "reason": "no_record", "records": 0,
                            "pages": pages if pages is not None else idx_slugs,
                            "rows": rows if rows is not None else (s.get("index_rows") or 0),
                            "alias_only": 0})
        elif n is None or pages is None:
            if idx_slugs > 0:
                uncounted.append(g)
        elif pages > 0:
            entries.append({**base, "reason": "pages_behind", "pages": pages,
                            "rows": rows or 0, "alias_only": alias})
        elif pages_stale > 0:
            reindex_entries.append({**base, "reason": "stale_slug", "pages": 0,
                                    "rows": rows or 0, "alias_only": alias})
        elif (rows or 0) > 0 or alias > 0:
            alias_entries.append({**base, "reason": "alias_only", "pages": 0,
                                  "rows": rows or 0, "alias_only": alias})
    refresh_entries = []
    for g, s in state.items():
        if s.get("status") in ("known_empty", "orphan"):
            continue
        if "stale" not in (s.get("issues") or []):
            continue
        refresh_entries.append({
            "gene": g,
            "reason": "stale_dataset",
            "datasets": (s.get("stale_on") or []),
            "records": s.get("n"),
            "index_rows": s.get("index_rows") or 0,
        })
    # Smallest genes first: most of the corpus is refreshed early, and an
    # interrupted run leaves only the whales to resume.
    refresh_entries.sort(key=lambda e: (e["records"] if e["records"] is not None else 0, e["gene"]))
    entries.sort(key=lambda e: (-e["pages"], -e["rows"], e["gene"]))
    reindex_entries.sort(key=lambda e: (-e["pages_stale"], -e["rows"], e["gene"]))
    alias_entries.sort(key=lambda e: (-e["alias_only"], -e["rows"], e["gene"]))
    by_reason = {}
    for e in entries:
        by_reason[e["reason"]] = by_reason.get(e["reason"], 0) + 1
    return {
        "schema": "sugi-variant-rebuild-queue/2",
        "release": manifest.get("release"),
        "generated_at": _utc_now(),
        "totals": {
            "genes": len(entries),
            "pages_to_add": sum(e["pages"] for e in entries),
            "rows_to_add": sum(e["rows"] for e in entries),
            "variants_to_add": sum(e["rows"] for e in entries),   # old name, rows not pages
            "by_reason": by_reason,
            "reindex_genes": len(reindex_entries),
            "stale_pages_to_reindex": sum(e["pages_stale"] for e in reindex_entries),
            "alias_only_genes": len(alias_entries),
            "alias_only_rows": sum(e["alias_only"] for e in alias_entries),
            "uncounted_genes": len(uncounted),
            "refresh_genes": len(refresh_entries),
            "unbuildable_genes": len(unbuildable_entries),
            "unbuildable_pages": sum(e["pages"] for e in unbuildable_entries),
        },
        "queue": entries,
        "reindex": reindex_entries,
        "alias_only": alias_entries,
        "uncounted": sorted(uncounted),
        "refresh": refresh_entries,
        "unbuildable": unbuildable_entries,
    }


SERVABLE = ("ok", "stale", "mismatch")


def serve_manifest(manifest, sv_db=None):
    """Compact manifest for the read-only serving mode.

    Only genes with a record file and a servable status. Missing/known-empty/
    orphan genes are absent: the app answers them as "not in snapshot" without
    touching a record build. `stale` is served (the page carries its source dates)
    and `mismatch` is served (the pickle exists; its page is usable).

    If `sv_db` is given, the structural-variant corpus is named in the release
    alongside the small-variant corpus.
    """
    corpus = ((manifest or {}).get("corpora") or {}).get(CORPUS) or {}
    state = corpus.get("gene_state") or {}
    # Genes a rebuild proved unbuildable. The audit cannot know this without building
    # them, so the queue runner records it and the queue honours it -- otherwise those
    # genes are re-attempted on every refresh and inflate "missing pages" with pages
    # that cannot exist.
    unbuildable = set()
    _rdir = corpus.get("records_dir") or (manifest or {}).get("records_dir")
    if _rdir:
        try:
            with open(os.path.join(_rdir, "unbuildable.json"), encoding="utf-8") as f:
                unbuildable = set(json.load(f).get("genes") or [])
        except Exception:
            unbuildable = set()
    genes = {}
    for g, s in state.items():
        if s.get("path") and s.get("status") in SERVABLE:
            genes[g] = {"path": s["path"], "status": s["status"],
                        "stale": "stale" in (s.get("issues") or [])}
    out = {
        "schema": "sugi-variant-release-serve/1",
        "release": manifest.get("release"),
        "released_at": manifest.get("released_at"),
        "variant_version": manifest.get("variant_version"),
        "variant_commit": manifest.get("variant_commit"),
        "biobtree_version": manifest.get("biobtree_version"),
        "biobtree_commit": manifest.get("biobtree_commit"),
        "created_at": manifest.get("created_at"),
        "records_dir": corpus.get("records_dir"),
        "index_db": corpus.get("index_db"),
        "corpora": {"small_variant": {"genes": genes, "status": corpus.get("status"),
                                      "datasets": corpus.get("datasets")}},
    }
    if sv_db:
        n = None
        try:
            conn = sqlite3.connect(f"file:{sv_db}?mode=ro", uri=True)
            n = conn.execute("SELECT COUNT(*) FROM sv").fetchone()[0]
            conn.close()
        except Exception:
            pass
        out["corpora"]["structural_variant"] = {"db": os.path.abspath(sv_db),
                                                "records": n}
    return out


def verify_manifest(manifest, records_dir, index_db, count_records=False, sample=None):
    """Re-check a manifest against the files/index it describes. Returns a report."""
    corpus = ((manifest or {}).get("corpora") or {}).get(CORPUS) or {}
    state = corpus.get("gene_state") or {}
    # Genes a rebuild proved unbuildable. The audit cannot know this without building
    # them, so the queue runner records it and the queue honours it -- otherwise those
    # genes are re-attempted on every refresh and inflate "missing pages" with pages
    # that cannot exist.
    unbuildable = set()
    _rdir = corpus.get("records_dir") or (manifest or {}).get("records_dir")
    if _rdir:
        try:
            with open(os.path.join(_rdir, "unbuildable.json"), encoding="utf-8") as f:
                unbuildable = set(json.load(f).get("genes") or [])
        except Exception:
            unbuildable = set()
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
            "release": manifest.get("release")}


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
         f"- release: `{report['release']}`",
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
         f"- **pages missing** (distinct slugs): **{t['pages_missing']:,}**",
         f"- stale index slugs (page exists; reindex): **{t['pages_stale']:,}**",
         f"- alias-only rows (duplicate VCVs on an existing page): **{t['alias_only_rows']:,}**",
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
    b.add_argument("--queue-out", help="write the rebuild queue JSON here")
    b.add_argument("--serve-out", help="write the compact serve manifest here")
    b.add_argument("--sv-db", help="include the dbVar sv.db in the release manifest")
    b.add_argument("--current-versions", help="JSON {dataset: version} to check staleness")
    b.add_argument("--release", help="release id (vX.Y.Z); default VARIANT_BUILD_VERSION or git describe")
    b.add_argument("--released-at", help="release date (YYYY-MM-DD; default today)")
    b.add_argument("--variant-version", help="code version stamp (default: release id)")
    b.add_argument("--variant-commit", help="code commit stamp (default: git short HEAD)")
    b.add_argument("--biobtree-version", help="biobtree release (default: $BIOBTREE_VERSION)")
    b.add_argument("--biobtree-commit", help="biobtree commit (default: $BIOBTREE_COMMIT)")
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
            count_records=a.count_records, count_sample=a.count_sample,
            release=a.release, released_at=a.released_at,
            variant_version=a.variant_version, variant_commit=a.variant_commit,
            biobtree_version=a.biobtree_version, biobtree_commit=a.biobtree_commit)
        _write_atomic(a.out, json.dumps(manifest, indent=2, sort_keys=False) + "\n")
        md = _md_report(report, manifest)
        if a.report:
            _write_atomic(a.report, md)
        queue = rebuild_queue(manifest)
        if a.queue_out:
            _write_atomic(a.queue_out, json.dumps(queue, indent=2, sort_keys=False) + "\n")
        if a.serve_out:
            serve = serve_manifest(manifest, sv_db=a.sv_db)
            _write_atomic(a.serve_out,
                          json.dumps(serve, indent=2, sort_keys=False) + "\n")
        t = report["totals"]
        print(f"release {report['release']}")
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
        t = queue["totals"]
        print(f"  rebuild queue: {t['genes']:,} genes / {t['pages_to_add']:,} missing pages"
              f" / {t['rows_to_add']:,} rows"
              + (f"  ({t['reindex_genes']:,} reindex-only genes, "
                 f"{t['stale_pages_to_reindex']:,} stale slugs)" if t["reindex_genes"] else "")
              + (f"  ({t['alias_only_genes']:,} alias-only genes, "
                 f"{t['alias_only_rows']:,} rows)" if t["alias_only_genes"] else "")
              + (f"  ({t['uncounted_genes']:,} uncounted)" if t["uncounted_genes"] else ""))
        for e in queue["queue"][:10]:
            print(f"    {e['pages']:>6,} pages {e['rows']:>7,} rows  {e['gene']:<12} "
                  f"{e['reason']:<14} records {e['records'] or 0:,} / index {e['index_rows']:,}")
        if a.queue_out:
            print(f"  queue        : {a.queue_out}")
        if a.serve_out:
            print(f"  serve        : {a.serve_out} "
                  f"({len(serve['corpora']['small_variant']['genes']):,} genes)")
        if a.strict and (report["missing"] or report["empty"] or
                         report["corrupt"] or report["mismatch"]):
            return 1
        return 0

    with open(a.manifest, encoding="utf-8") as f:
        manifest = json.load(f)
    res = verify_manifest(manifest, a.records, a.index,
                          count_records=a.count_records, sample=a.sample)
    print(f"verify {res['release']}  checked {res['checked']:,}  "
          f"{'OK' if res['ok'] else str(len(res['errors'])) + ' errors'}")
    for e in res["errors"][:50]:
        print("  ", e)
    if len(res["errors"]) > 50:
        print(f"  … and {len(res['errors']) - 50} more")
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
