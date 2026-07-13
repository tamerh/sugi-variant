"""Small vendored helpers (decoupled from atlas.render_common / atlas.pipeline):
markdown table, YAML escaping, version stamps, and the page-frontmatter builder."""
import html
import json
import os
import re
import subprocess
import urllib.request
from datetime import datetime, timezone

from sugivariant.biobtree import API

GENERATED_BY = "Sugi Variant"
_VERSION = "sugi-variant-0.1"
_CTRL = re.compile(r"[\x00-\x1f\x7f]")


def _cell(c):
    if c is None:
        return ""
    s = html.unescape(str(c))
    if s.strip().lower() == "nan":
        return ""
    return s.replace("|", "\\|")


def table(headers, rows):
    """GitHub-flavored markdown table; blank cells empty, pipes escaped, identical
    rows collapsed, header-only table renders as ''."""
    out = ["| " + " | ".join(headers) + " |",
           "| " + " | ".join("---" for _ in headers) + " |"]
    seen = set()
    for r in rows:
        cells = tuple(_cell(c) for c in r)
        if not any(cells) or cells in seen:
            continue
        seen.add(cells)
        out.append("| " + " | ".join(cells) + " |")
    return "" if len(out) == 2 else "\n".join(out)


def _yaml_escape(s):
    return _CTRL.sub("", str(s)).replace("\\", "\\\\").replace('"', '\\"')


_ATLAS_VERSION = None


def atlas_version():
    """The Sugi Variant build stamp (git describe of this repo, else the packaged
    version). Kept named atlas_version for frontmatter-key compatibility."""
    global _ATLAS_VERSION
    if _ATLAS_VERSION is None:
        try:
            _ATLAS_VERSION = subprocess.check_output(
                ["git", "describe", "--tags", "--always"],
                cwd=os.path.dirname(__file__), stderr=subprocess.DEVNULL).decode().strip() or _VERSION
        except Exception:
            _ATLAS_VERSION = _VERSION
    return _ATLAS_VERSION


_BB_META = None


def _biobtree_meta():
    global _BB_META
    if _BB_META is None:
        try:
            d = json.loads(urllib.request.urlopen(f"{API}/ws/meta", timeout=5).read())
            _BB_META = d.get("appparams") or {}
        except Exception:
            _BB_META = {}
    return _BB_META


def biobtree_version():
    return _biobtree_meta().get("biobtree_version") or "unknown"


def biobtree_commit():
    return _biobtree_meta().get("biobtree_commit") or "unknown"


def build_meta(entity_type, slug, title, datasets, generated_at=None):
    """Page-frontmatter meta (variant needs no bundle/evidence machinery)."""
    return {
        "title": title, "symbol": slug, "entity_type": entity_type,
        "generated_at": generated_at or datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "atlas_version": atlas_version(), "biobtree_version": biobtree_version(),
        "biobtree_commit": biobtree_commit(), "generated_by": GENERATED_BY,
        "datasets": datasets,
    }
