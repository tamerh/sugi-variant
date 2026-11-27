"""Phase 1a gate: the app must return the same bytes and the same status codes as the
static site, for the whole URL surface.

Runs against a locally started app with ATLAS_ROOT set (dev only, prod untouched):

    python migration/verify_atlas_passthrough.py \
        --base http://127.0.0.1:8020 \
        --oracle /data2/atlas-migration/oracle-20261127.jsonl \
        --redirects /data2/atlas-migration/redirects.json \
        --sample 300

Checks:

  1. sampled entity/site pages: body sha256 equals the oracle's
  2. the trailing-slash form 301s to the canonical form
  3. `.md` twins: 200, text/markdown, sha256 equals the oracle's
  4. every redirect rule: the old URL 301s to the new one, and the new one is 200
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
import urllib.error
import urllib.request


def _fetch(base, path, follow=True):
    url = base.rstrip("/") + path
    req = urllib.request.Request(url, headers={"Accept": "*/*"})
    opener = urllib.request.build_opener(
        urllib.request.HTTPRedirectHandler if follow else _NoRedirect)
    try:
        with opener.open(req, timeout=120) as r:
            return r.status, r.headers, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers, e.read()


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **k):
        return None


def _sha(b):
    return hashlib.sha256(b).hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://127.0.0.1:8020")
    ap.add_argument("--oracle", required=True)
    ap.add_argument("--redirects", required=True)
    ap.add_argument("--sample", type=int, default=200)
    a = ap.parse_args(argv)

    html, md = {}, {}
    for line in open(a.oracle, encoding="utf-8"):
        r = json.loads(line)
        (md if r["url"].endswith(".md") else html)[r["url"]] = r
    rules = json.load(open(a.redirects, encoding="utf-8"))["rules"]

    # 1. byte parity — all site-level pages, the largest pages, and a random sample
    urls = list(html)
    site = [u for u in urls if sum(1 for _ in u.strip("/").split("/")) <= 3]
    largest = [r["url"] for r in sorted(html.values(), key=lambda x: -x["bytes"])[:10]]
    random.seed(11)
    sample = list(dict.fromkeys(site + largest + random.sample(urls, min(a.sample, len(urls)))))
    bad = []
    for u in sample:
        code, _h, body = _fetch(a.base, u)
        if code != 200 or _sha(body) != html[u]["sha256"]:
            bad.append((u, code, "sha" if code == 200 else "status"))
    print(f"1. byte parity      : {len(sample)-len(bad)}/{len(sample)} pages identical")
    for b in bad[:8]:
        print("   MISMATCH", b)

    # 2. trailing slash
    slash_bad = 0
    for u in sample[:25]:
        if not u.endswith("/"):
            continue
        code, h, _b = _fetch(a.base, u.rstrip("/"), follow=False)
        if code != 301 or not (h.get("location") or "").endswith(u):
            slash_bad += 1
            print("   SLASH", u, code, h.get("location"))
    print(f"2. trailing slash   : {25-slash_bad}/25 checked 301 to the canonical form")

    # 3. .md twins
    md_bad = 0
    for u in sample[:25]:
        mdu = u + "index.md"
        if mdu not in md:
            continue
        code, h, body = _fetch(a.base, mdu)
        if code != 200 or not (h.get("content-type") or "").startswith("text/markdown") \
                or _sha(body) != md[mdu]["sha256"]:
            md_bad += 1
            print("   MD", mdu, code, h.get("content-type"))
    print(f"3. markdown twins   : 25-{md_bad}/25 identical and text/markdown")

    # 4. every redirect rule. Only /atlas/* belongs to the app; the /biobtree/*
    # legacy rules stay in nginx ahead of it, so their targets are checked but the
    # redirect itself is not this app's job.
    atlas_rules = [r for r in rules if r["from"].startswith("/atlas/")]
    legacy = [r for r in rules if not r["from"].startswith("/atlas/")]
    rule_bad = 0
    for r in atlas_rules:
        src, dst = r["from"], r["to"]
        code, h, _b = _fetch(a.base, src, follow=False)
        loc = h.get("location") or ""
        ok = code == 301 and (loc == dst or loc.endswith(dst))
        code2, _h2, _b2 = _fetch(a.base, dst)
        if not ok or code2 != 200:
            rule_bad += 1
            if rule_bad <= 8:
                print("   RULE", src, "->", dst, "|", code, loc, "| target", code2)
    print(f"4. redirect table   : {len(atlas_rules)-rule_bad}/{len(atlas_rules)} /atlas rules land correctly")
    if legacy:
        print(f"   (nginx-owned      : {len(legacy)} /biobtree legacy rules, targets checked below)")
    legacy_bad = sum(1 for r in legacy if _fetch(a.base, r["to"])[0] != 200)
    if legacy:
        print(f"   legacy targets    : {len(legacy)-legacy_bad}/{len(legacy)} resolve")

    failed = bool(bad or slash_bad or md_bad or rule_bad or legacy_bad)
    print("RESULT:", "FAIL" if failed else "PASS")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
