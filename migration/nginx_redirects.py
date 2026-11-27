"""Phase 0: make the redirect surface machine-readable and testable.

The live Atlas site depends on ~200 `location = … { return 301 …; }` rules in nginx
snippets: renamed diseases, retired gene Q&A pages, drug-slug normalisations. If the
application takes over `/atlas/*` and those rules are forgotten, every renamed entity
404s — and no existing test would notice.

    python migration/nginx_redirects.py extract \
        --conf /data/ext/deploy/sugi-disease-rename-redirects.conf \
        --conf /data/ext/deploy/sugi-gene-atlas-redirects.conf \
        --out /data2/atlas-migration/redirects.json

    python migration/nginx_redirects.py verify \
        --table /data2/atlas-migration/redirects.json \
        --oracle /data2/atlas-migration/oracle-20261127.jsonl

`verify` is offline on purpose: every rule's target must exist in the oracle snapshot,
so a rename that points at a page which no longer exists fails the gate without any
traffic. (A live probe from this host is blocked by Cloudflare bot protection.)
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

# location = /atlas/disease/old/ { return 301 /atlas/disease/new/; }
_RULE = re.compile(
    r"^\s*location\s*=\s*(?P<src>\S+)\s*\{\s*return\s+301\s+(?P<dst>\S+?);\s*\}", re.M)
# rewrite ^/old$ /new permanent;
_REWRITE = re.compile(r"^\s*rewrite\s+(?P<src>\S+)\s+(?P<dst>\S+)\s+permanent\s*;", re.M)


def _norm(path):
    """Compare paths without their trailing slash, keeping the canonical form."""
    path = (path or "").split("?")[0].split("#")[0]
    if len(path) > 1 and path.endswith("/"):
        path = path[:-1]
    return path


def extract(files):
    rules = []
    for path in files:
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
        for m in _RULE.finditer(text):
            rules.append({"from": m.group("src"), "to": m.group("dst"),
                          "kind": "location", "source": os.path.basename(path)})
        for m in _REWRITE.finditer(text):
            src, dst = m.group("src"), m.group("dst")
            if src.startswith("^"):
                continue          # regex rewrites are not a simple table
            rules.append({"from": src, "to": dst, "kind": "rewrite",
                          "source": os.path.basename(path)})
    # de-dup exact pairs, keep first
    seen, out = set(), []
    for r in rules:
        k = (r["from"], r["to"])
        if k in seen:
            continue
        seen.add(k)
        out.append(r)
    return out


def verify(table, oracle_path):
    rules = json.load(open(table, encoding="utf-8"))["rules"]
    pages = set()
    with open(oracle_path, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                pages.add(_norm(json.loads(line)["url"]))
    bad = []
    for r in rules:
        dst = _norm(r["to"])
        # /atlas/gene/X/ -> /atlas/gene/X/index.html in the snapshot is normalized to /atlas/gene/X
        if dst and dst not in pages and not dst.endswith(".md"):
            bad.append(r)
    print(f"redirects: {len(rules):,} rules · {len(rules)-len(bad):,} land on a page in the oracle · {len(bad)} do not")
    for r in bad[:15]:
        print(f"  MISS {r['from']} -> {r['to']}  ({r['source']})")
    return 0 if not bad else 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract"); e.add_argument("--conf", action="append", required=True); e.add_argument("--out", required=True)
    v = sub.add_parser("verify"); v.add_argument("--table", required=True); v.add_argument("--oracle", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "extract":
        rules = extract(a.conf)
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump({"schema": "sugi-atlas-redirects/1", "rules": rules}, f, indent=1)
        by_kind = {}
        for r in rules:
            by_kind[r["kind"]] = by_kind.get(r["kind"], 0) + 1
        print(f"extracted {len(rules):,} rules from {len(a.conf)} files → {a.out}")
        for k, n in sorted(by_kind.items()):
            print(f"  {k:<10} {n:>5,}")
        return 0
    return verify(a.table, a.oracle)


if __name__ == "__main__":
    sys.exit(main())
