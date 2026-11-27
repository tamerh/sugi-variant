"""Phase 1b: turn the static Atlas tree into a per-entity pickle corpus.

Atlas's generation is not touched. This reads what the build already produced
(`index.html` + `index.md` per entity) and writes one record per entity:

    atlas-cache/gene/BRCA1.pkl.gz
      {type, slug, title, canonical, body_html, body_md, sections, sha256, bytes}

`body_html` is the page's `<article class=entity-page>` **fragment**, not the whole
page: the app supplies the shell. Keeping the fragment means the content is never
re-rendered, so there is nothing to drift.

`atlas_index.db` is the lookup the app needs before loading a record (hit vs 404,
sitemaps, resolver): (type, slug) -> path, title, fragment sha, md sha, bytes. It
also carries the alias table from the tree's `manifest.json` (names and synonyms,
normalised the way Atlas's own search normalises them), so the resolver can answer
"carasil syndrome" without loading a 13 MB JSON at serve time.

    python migration/build_atlas_corpus.py \
        --html /data2/atlas-migration/atlas-html \
        --out /data2/atlas-migration/atlas-cache \
        --index /data2/atlas-migration/atlas_index.db \
        --apply

Site-level pages (about, the type indexes, the home page) are not entities: they stay
files and the app falls back to the tree for them.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import pickle
import re
import sqlite3
import sys
import time

TYPES = ("gene", "disease", "drug")
_ART_OPEN = re.compile(r"<article\b[^>]*>", re.I)
_ART_CLOSE = "</article>"
_TITLE = re.compile(r"<title>(.*?)</title>", re.S | re.I)
_CANON = re.compile(r'<link[^>]*rel=["\']?canonical["\']?[^>]*href=["\']?([^"\'\s>]+)', re.I)
_CANON2 = re.compile(r'<link[^>]*href=["\']?([^"\'\s>]+)["\']?[^>]*rel=["\']?canonical', re.I)
_H2 = re.compile(r"<h2\b[^>]*\bid=[\"']?([^\"'\s>]+)[\"']?[^>]*>(.*?)</h2>", re.S | re.I)

SCHEMA = """
CREATE TABLE IF NOT EXISTS entity(
  type TEXT, slug TEXT, path TEXT, title TEXT, canonical TEXT,
  frag_sha TEXT, md_sha TEXT, bytes INTEGER, sections TEXT,
  PRIMARY KEY(type, slug));
CREATE TABLE IF NOT EXISTS alias(
  key TEXT, type TEXT, slug TEXT, name TEXT,
  PRIMARY KEY(key, type, slug));
CREATE INDEX IF NOT EXISTS alias_key ON alias(key);
CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT);
"""


def normalize_name(name):
    """Atlas's own search normalisation: lowercase, punctuation to one space."""
    return re.sub(r"[^a-z0-9]+", " ", (name or "").lower()).strip()


def aliases(root):
    """(type, name, slug) from manifest.json's name maps and synonym maps."""
    path = os.path.join(root, "manifest.json")
    try:
        with open(path, encoding="utf-8") as f:
            m = json.load(f)
    except Exception:
        return []
    out = []
    syn = m.get("syn") or {}
    for etype in TYPES:
        for source in (m.get(etype) or {}, syn.get(etype) or {}):
            for name, slug in source.items():
                if name and slug:
                    out.append((etype, name, slug))
    return out


def _text(html):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", html or "")).strip()


def fragment(html):
    m = _ART_OPEN.search(html)
    if not m:
        return None
    end = html.find(_ART_CLOSE, m.start())
    if end < 0:
        return None
    return html[m.start():end + len(_ART_CLOSE)]


def sections(frag):
    out = []
    for sid, label in _H2.findall(frag or ""):
        out.append({"id": sid, "title": _text(label)})
    return out


def entity(html_path, md_path, etype, slug):
    html = open(html_path, encoding="utf-8", errors="replace").read()
    frag = fragment(html)
    if not frag:
        return None
    title_m = _TITLE.search(html)
    title = _text(title_m.group(1)) if title_m else slug
    canon = None
    for rx in (_CANON, _CANON2):
        cm = rx.search(html)
        if cm:
            canon = cm.group(1)
            break
    body_md = ""
    if md_path and os.path.isfile(md_path):
        body_md = open(md_path, encoding="utf-8", errors="replace").read()
    return {
        "schema": "sugi-atlas-entity/1",
        "type": etype, "slug": slug, "title": title, "canonical": canon,
        "body_html": frag, "body_md": body_md, "sections": sections(frag),
        "frag_sha": hashlib.sha256(frag.encode()).hexdigest(),
        "md_sha": hashlib.sha256(body_md.encode()).hexdigest() if body_md else None,
        "bytes": len(frag.encode()),
    }


def _write(path, rec):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.{os.getpid()}.tmp"
    with gzip.open(tmp, "wb", compresslevel=6) as f:
        pickle.dump(rec, f, protocol=4)
    os.replace(tmp, path)


def scan(html_root):
    for etype in TYPES:
        base = os.path.join(html_root, etype)
        if not os.path.isdir(base):
            continue
        for slug in sorted(os.listdir(base)):
            d = os.path.join(base, slug)
            idx = os.path.join(d, "index.html")
            if not os.path.isfile(idx):
                continue
            yield etype, slug, idx, os.path.join(d, "index.md")


def build(html_root, out, index, apply, limit=0, shard=None):
    jobs = list(scan(html_root))
    if shard:
        k, n = (int(x) for x in shard.split("/"))
        jobs = [j for i, j in enumerate(jobs) if i % n == k]
    if limit:
        jobs = jobs[:limit]
    conn = sqlite3.connect(index) if apply else None
    full = apply and not limit and not shard
    if conn:
        conn.executescript(SCHEMA)
        if full:
            # A full build replaces the corpus; stale rows from an older tree would
            # otherwise stay resolvable after their page disappears.
            conn.execute("DELETE FROM entity")
            conn.execute("DELETE FROM alias")
        conn.commit()
    t0, done, skipped, written = time.time(), 0, 0, 0
    for etype, slug, html_path, md_path in jobs:
        rec = entity(html_path, md_path, etype, slug)
        done += 1
        if not rec:
            skipped += 1
            continue
        if apply:
            _write(os.path.join(out, etype, slug + ".pkl.gz"), rec)
            conn.execute(
                "INSERT OR REPLACE INTO entity VALUES(?,?,?,?,?,?,?,?,?)",
                (etype, slug, f"{etype}/{slug}.pkl.gz", rec["title"], rec["canonical"],
                 rec["frag_sha"], rec["md_sha"], rec["bytes"],
                 json.dumps(rec["sections"])))
            written += 1
            if written % 2000 == 0:
                conn.commit()
                print(f"  {done:,}/{len(jobs):,}  {time.time()-t0:.0f}s", flush=True)
        elif done % 5000 == 0:
            print(f"  {done:,}/{len(jobs):,}", flush=True)
    if conn and full:
        # Names and synonyms, keyed the way the Atlas search normalises them.
        canon = {}
        try:
            with open(os.path.join(html_root, "manifest.json"), encoding="utf-8") as f:
                canon = (json.load(f).get("canon") or {})
        except Exception:
            canon = {}
        rows = []
        known = set((t, s) for t, s in conn.execute("SELECT type, slug FROM entity"))
        for etype, name, slug in aliases(html_root):
            if (etype, slug) not in known:
                continue
            key = normalize_name(name)
            if not key:
                continue
            display = (canon.get(etype) or {}).get(slug) or name
            rows.append((key, etype, slug, display))
        conn.executemany(
            "INSERT OR REPLACE INTO alias(key, type, slug, name) VALUES(?,?,?,?)", rows)
        n_alias = conn.execute("SELECT COUNT(*) FROM alias").fetchone()[0]
        conn.execute("INSERT OR REPLACE INTO meta VALUES('aliases',?)", (str(n_alias),))
        print(f"  aliases      : {n_alias:,} (from {len(rows):,} names)", flush=True)
    if conn:
        conn.execute("INSERT OR REPLACE INTO meta VALUES('schema','sugi-atlas-index/2')")
        conn.execute("INSERT OR REPLACE INTO meta VALUES('built_at',?)",
                     (time.strftime("%Y-%m-%dT%H:%M:%S"),))
        conn.commit()
        conn.close()
    print(f"{'APPLY' if apply else 'DRY RUN'} -- {done:,} entities scanned · "
          f"{written:,} written · {skipped} without a fragment · {time.time()-t0:.0f}s")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--shard", default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)
    return build(a.html, a.out, a.index, a.apply, a.limit, a.shard)


if __name__ == "__main__":
    sys.exit(main())
