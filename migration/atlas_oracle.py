"""Phase 0 for the Atlas migration: freeze a page-for-page oracle of the live site.

The Atlas site is 54k static pages. Before any of it is re-served by an application,
we need a baseline that says exactly what those pages are, so every later phase can
prove parity instead of asserting it.

    python migration/atlas_oracle.py snapshot --root /var/www/sugi.bio/atlas \
        --out /data2/atlas-migration/oracle-20261127.jsonl
    python migration/atlas_oracle.py verify --root /var/www/sugi.bio/atlas \
        --snapshot /data2/atlas-migration/oracle-20261127.jsonl

A snapshot line is one page:

    {"url": "/atlas/gene/BRCA1/", "type": "gene", "slug": "BRCA1",
     "file": "gene/BRCA1/index.html", "bytes": 256349, "sha256": "…", "mtime": "…"}

`url` is the public path the browser asks for; `file` is the path under the root.
Markdown twins (`index.md`) are recorded too, so the `.md` surface is covered by the
same gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

EXTS = (".html", ".md")


def _sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _record(root, path):
    rel = os.path.relpath(path, root).replace(os.sep, "/")
    parts = rel.split("/")
    if parts[-1] not in ("index.html", "index.md"):
        return None
    entity = parts[:-1]
    if not entity:
        return None
    url = "/atlas/" + "/".join(entity) + "/" + ("" if parts[-1] == "index.html" else "index.md")
    st = os.stat(path)
    return {
        "url": url,
        "type": entity[0] if len(entity) > 1 else "",
        "slug": entity[1] if len(entity) > 1 else entity[0],
        "file": rel,
        "bytes": st.st_size,
        "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)),
        "sha256": _sha256(path),
    }


def snapshot(root, out):
    rows = []
    for dirpath, _dirs, files in os.walk(root):
        for name in files:
            if not name.endswith(EXTS) or name not in ("index.html", "index.md"):
                continue
            rec = _record(root, os.path.join(dirpath, name))
            if rec:
                rows.append(rec)
    rows.sort(key=lambda r: r["url"])
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        for r in rows:
            f.write(json.dumps(r, sort_keys=True) + "\n")
    by_type = {}
    for r in rows:
        by_type[r["type"]] = by_type.get(r["type"], 0) + 1
    print(f"snapshot: {len(rows):,} pages → {out}")
    for t, n in sorted(by_type.items()):
        print(f"  {t:<10} {n:>7,}")
    return rows


def load(path):
    with open(path, encoding="utf-8") as f:
        return {json.loads(line)["url"]: json.loads(line) for line in f if line.strip()}


def verify(root, snap_path):
    old = load(snap_path)
    new = {r["url"]: r for r in ( _record(root, os.path.join(dp, fn))
                                  for dp, _d, fs in os.walk(root)
                                  for fn in fs if fn in ("index.html", "index.md") )
           if r}
    added = sorted(set(new) - set(old))
    removed = sorted(set(old) - set(new))
    changed = sorted(u for u in set(old) & set(new) if old[u]["sha256"] != new[u]["sha256"])
    print(f"oracle verify: {len(old):,} → {len(new):,} pages")
    print(f"  added {len(added):,} · removed {len(removed):,} · changed {len(changed):,}")
    for label, rows in (("added", added[:10]), ("removed", removed[:10]), ("changed", changed[:10])):
        for u in rows:
            print(f"  {label:<7} {u}")
    return 0 if not (added or removed or changed) else 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("snapshot"); s.add_argument("--root", required=True); s.add_argument("--out", required=True)
    v = sub.add_parser("verify"); v.add_argument("--root", required=True); v.add_argument("--snapshot", required=True)
    a = ap.parse_args(argv)
    return snapshot(a.root, a.out) and 0 if a.cmd == "snapshot" else verify(a.root, a.snapshot)


if __name__ == "__main__":
    sys.exit(main())
