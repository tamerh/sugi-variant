"""Atlas passthrough (Phase 1a).

The first step of moving Atlas into the app is **not** to re-render it. It is to make
the app answer the same URLs with the same bytes, under the same release headers, so
the URL semantics — 301, trailing slash, 404, `.md` twins, ETag — can be proven
independently of any renderer.

So this module serves the existing static tree verbatim:

    /atlas/gene/BRCA1/         → <root>/gene/BRCA1/index.html
    /atlas/gene/BRCA1          → 301 /atlas/gene/BRCA1/
    /atlas/gene/BRCA1/index.md → <root>/gene/BRCA1/index.md   (text/markdown)
    /atlas/<renamed-slug>/     → 301 to its successor, from the nginx table

Enabled only when `ATLAS_ROOT` is set, so a deployment that does not opt in behaves
exactly as before (Atlas then 404s from this app, as it does today).

Phase 1b replaces the file read with a pickle record and the fragment with the shared
shell; the routing and the gates built here stay.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("ATLAS_ROOT") or "/nonexistent")
REDIRECTS = Path(os.environ.get("ATLAS_REDIRECTS") or "/nonexistent")
ENABLED = ROOT.is_dir()

_TYPES = {"gene", "disease", "drug"}
_MD = {".md": "text/markdown; charset=utf-8"}

_sha_cache: dict[tuple[str, int, int], str] = {}
_redirects: dict[str, str] | None = None


def enabled() -> bool:
    return ENABLED


def redirects() -> dict[str, str]:
    """{old path: new path} from the Phase 0 table, loaded once."""
    global _redirects
    if _redirects is None:
        _redirects = {}
        try:
            with open(REDIRECTS, encoding="utf-8") as f:
                for r in json.load(f).get("rules") or []:
                    src, dst = r.get("from"), r.get("to")
                    if src and dst:
                        _redirects[src] = dst
        except Exception:
            _redirects = {}
    return _redirects


def _safe(rel: str) -> str | None:
    """A tree-relative path, or None if it tries to escape the root."""
    # Only the leading slash: the trailing one is meaningful here (Atlas canonical
    # URLs end in "/"), and stripping it made the canonical redirect loop.
    rel = (rel or "").lstrip("/")
    if not rel:
        return ""
    if ".." in rel.split("/") or rel.startswith("/") or "\x00" in rel:
        return None
    return rel


def target_for(rel: str) -> tuple[str, str] | None:
    """('file', path) | ('redirect', to) | ('404', reason) for a tree-relative path.

    Kept separate from I/O so the URL rules are unit-testable without a corpus.
    """
    rel = _safe(rel)
    if rel is None:
        return ("404", "unsafe")
    public = "/atlas/" + rel
    table = redirects()
    # The rename table is consulted first, exactly as nginx does, for both the
    # slash and non-slash forms it contains.
    for key in (public, public.rstrip("/"), public.rstrip("/") + "/"):
        if key in table:
            return ("redirect", table[key])
    if rel == "":
        return ("file", "index.html")
    clean = rel.rstrip("/")
    last = clean.rsplit("/", 1)[-1]
    if last == "index.md":
        return ("file", clean)
    if "." in last and not last.endswith(".html"):
        # any other file (assets) is served from the tree as-is
        return ("file", clean)
    if not public.endswith("/"):
        # Atlas canonical form is the trailing slash.
        return ("redirect", public + "/")
    if last == "index.html":
        return ("file", clean)
    return ("file", clean + "/index.html")


def etag_for(path: Path) -> str:
    st = path.stat()
    key = (str(path), st.st_mtime_ns, st.st_size)
    hit = _sha_cache.get(key)
    if hit:
        return hit
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    digest = '"' + h.hexdigest()[:32] + '"'
    if len(_sha_cache) > 20000:          # bounded; pages are re-hashed after eviction
        _sha_cache.clear()
    _sha_cache[key] = digest
    return digest


def read(rel: str):
    """(bytes, content_type, etag) for a tree-relative path, or None if absent."""
    rel = _safe(rel)
    if rel is None:
        return None
    path = ROOT / rel
    if not path.is_file():
        return None
    ctype = _MD.get(path.suffix.lower()) or "text/html; charset=utf-8"
    return path.read_bytes(), ctype, etag_for(path)
