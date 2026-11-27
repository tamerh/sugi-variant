"""The Atlas entity corpus, as the release manifest describes it.

A release has to say *which* Atlas corpus it serves. Listing 54,351 entities in
`release.json` would add megabytes for no gain: the corpus is a directory of
per-entity pickles plus one SQLite index, and the index is already the authority on
membership. So the manifest names the index and anchors it with a sha256 and its
counts. At startup the app re-checks those, and if they disagree it refuses to serve
Atlas rather than serving a corpus the manifest does not describe.

What the anchors are for:

* `index_sha256` — the membership and the metadata (title, canonical, fragment sha).
* `root_index_sha256` / `search_sha256` — the site-level pages and the search data,
  which still come from the tree (about, the type indexes, `/atlas/search/`).
  Hashing the whole 5.7 GB tree on every start is not an option, so the two files
  that identify the tree are anchored instead.
* `redirect_rules` — the 301 table that makes old Atlas URLs resolve.

Nothing here writes, builds or renders. It reads.
"""
from __future__ import annotations

import hashlib
import json
import os
import sqlite3
from pathlib import Path

FORMAT = "atlas-entity-pickle-v1"
SCHEMA = "sugi-atlas-corpus/1"


def sha256_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _dir_bytes(root: Path) -> int:
    total = 0
    for base, _dirs, files in os.walk(root):
        for name in files:
            try:
                total += os.path.getsize(os.path.join(base, name))
            except OSError:
                pass
    return total


def _file_count(root: Path) -> int:
    return sum(len(files) for _base, _dirs, files in os.walk(root))


def stats(index_db) -> dict:
    """Counts from the entity index: entities, by type, fragment bytes, built_at."""
    conn = sqlite3.connect(f"file:{index_db}?mode=ro", uri=True)
    try:
        by_type = dict(conn.execute(
            "SELECT type, COUNT(*) FROM entity GROUP BY type").fetchall())
        total = conn.execute("SELECT COUNT(*) FROM entity").fetchone()[0]
        frag = conn.execute("SELECT COALESCE(SUM(bytes), 0) FROM entity").fetchone()[0]
        with_md = conn.execute(
            "SELECT COUNT(*) FROM entity WHERE md_sha IS NOT NULL").fetchone()[0]
        try:
            meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
        except sqlite3.Error:
            meta = {}
    finally:
        conn.close()
    return {"entities": total, "by_type": by_type, "fragment_bytes": frag,
            "with_markdown": with_md, "built_at": meta.get("built_at"),
            "index_schema": meta.get("schema")}


def describe(cache, index_db, root=None, redirects=None) -> dict:
    """The `corpora.atlas` block for a release manifest."""
    cache, index_db = Path(cache), Path(index_db)
    if not index_db.is_file():
        raise FileNotFoundError(f"no Atlas entity index at {index_db}")
    if not cache.is_dir():
        raise FileNotFoundError(f"no Atlas entity cache at {cache}")
    s = stats(index_db)
    desc = {
        "schema": SCHEMA,
        "format": FORMAT,
        "cache": os.path.abspath(cache),
        "index_db": os.path.abspath(index_db),
        "index_sha256": sha256_file(index_db),
        "entities": s["entities"],
        "by_type": s["by_type"],
        "fragment_bytes": s["fragment_bytes"],
        "corpus_bytes": _dir_bytes(cache),
        "corpus_files": _file_count(cache),
        "with_markdown": s["with_markdown"],
        "built_at": s["built_at"],
        "scope": "all-static-entities",
    }
    if root:
        root = Path(root)
        desc["root"] = os.path.abspath(root)
        if (root / "index.html").is_file():
            desc["root_index_sha256"] = sha256_file(root / "index.html")
        if (root / "manifest.json").is_file():
            desc["search"] = {"path": "manifest.json",
                              "bytes": (root / "manifest.json").stat().st_size,
                              "sha256": sha256_file(root / "manifest.json")}
    if redirects:
        redirects = Path(redirects)
        if redirects.is_file():
            desc["redirects"] = os.path.abspath(redirects)
            desc["redirects_sha256"] = sha256_file(redirects)
            try:
                desc["redirect_rules"] = len(
                    json.loads(redirects.read_text(encoding="utf-8")).get("rules") or [])
            except Exception:
                desc["redirect_rules"] = None
    return desc


def check(desc, cache=None, index_db=None, root=None, redirects=None) -> list[str]:
    """Re-check a manifest block against the files on disk.

    Paths in the manifest are the ones the build saw, which are host paths. A
    container mounts the same directory elsewhere, so the caller resolves the
    location and passes it in; the checks are on the resolved paths.
    """
    problems: list[str] = []
    cache = Path(cache or desc.get("cache") or "/nonexistent")
    index_db = Path(index_db or desc.get("index_db") or "/nonexistent")
    if not cache.is_dir():
        problems.append(f"cache missing: {cache}")
    if not index_db.is_file():
        problems.append(f"index missing: {index_db}")
        return problems
    want = desc.get("index_sha256")
    if want:
        got = sha256_file(index_db)
        if got != want:
            problems.append(f"index sha {got[:16]} != manifest {want[:16]}")
    s = stats(index_db)
    if desc.get("entities") is not None and s["entities"] != desc["entities"]:
        problems.append(f"entities {s['entities']} != manifest {desc['entities']}")
    if desc.get("by_type") and s["by_type"] != desc["by_type"]:
        problems.append(f"by_type {s['by_type']} != manifest {desc['by_type']}")
    has_md = desc["with_markdown"] if desc.get("with_markdown") is not None else None
    if has_md is not None and s["with_markdown"] != has_md:
        problems.append(f"with_markdown {s['with_markdown']} != manifest {has_md}")
    root = Path(root or desc.get("root") or "/nonexistent")
    if desc.get("root"):
        if not root.is_dir():
            problems.append(f"root missing: {root}")
        elif desc.get("root_index_sha256"):
            got = sha256_file(root / "index.html")
            if got != desc["root_index_sha256"]:
                problems.append(f"root index sha {got[:16]} != manifest")
        search = desc.get("search") or {}
        if search.get("sha256"):
            path = root / (search.get("path") or "manifest.json")
            if not path.is_file():
                problems.append(f"search data missing: {path}")
            elif sha256_file(path) != search["sha256"]:
                problems.append("search data sha differs from manifest")
    red = desc.get("redirects")
    if red:
        path = Path(redirects or red)
        if not path.is_file():
            problems.append(f"redirects missing: {path}")
        elif desc.get("redirects_sha256") and \
                sha256_file(path) != desc["redirects_sha256"]:
            problems.append("redirects sha differs from manifest")
    return problems
