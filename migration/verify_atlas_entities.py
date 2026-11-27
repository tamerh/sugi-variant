"""Phase 1b gate: every entity page must serve the fragment Atlas built, unchanged.

The app now renders the entity content from a pickle through the shared shell. The
shell is expected to differ from the static page; the *entity content* is not. So this
gate compares the `<article>` block the app serves against the `frag_sha` recorded when
the corpus was built from the static tree — plus the `.md` twin by sha.

    python migration/verify_atlas_entities.py \
        --base http://127.0.0.1:8111 \
        --index /data2/atlas-migration/atlas_index.db \
        --sample 2000 --threads 8
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sqlite3
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

_ART = re.compile(r"<article\b[^>]*>.*?</article>", re.S | re.I)


def _get(base, path):
    req = urllib.request.Request(base.rstrip("/") + path, headers={"Accept": "*/*"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def _sha(b):
    return hashlib.sha256(b).hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://127.0.0.1:8111")
    ap.add_argument("--index", required=True)
    ap.add_argument("--sample", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)

    conn = sqlite3.connect(f"file:{a.index}?mode=ro", uri=True)
    rows = conn.execute(
        "SELECT type, slug, frag_sha, md_sha, bytes FROM entity "
        "ORDER BY bytes DESC").fetchall()
    biggest = rows[:20]
    rest = rows[20:]
    step = max(1, len(rest) // max(1, a.sample))
    sample = biggest + rest[::step]
    print(f"checking {len(sample):,} of {len(rows):,} entities "
          f"({len(biggest)} largest + every {step}th)")

    def check(row):
        etype, slug, frag_sha, md_sha, _b = row
        code, body = _get(a.base, f"/atlas/{etype}/{slug}/")
        if code != 200:
            return (slug, f"page {code}")
        m = _ART.search(body.decode("utf-8", "replace"))
        if not m:
            return (slug, "no <article> in response")
        if _sha(m.group(0).encode()) != frag_sha:
            return (slug, "fragment bytes differ")
        mcode, mbody = _get(a.base, f"/atlas/{etype}/{slug}/index.md")
        if md_sha and (mcode != 200 or _sha(mbody) != md_sha):
            return (slug, f"markdown {mcode}")
        return None

    bad = []
    with ThreadPoolExecutor(max_workers=a.threads) as pool:
        for res in pool.map(check, sample):
            if res:
                bad.append(res)
    print(f"entity content: {len(sample)-len(bad):,}/{len(sample):,} identical")
    for slug, why in bad[:10]:
        print(f"  {slug}: {why}")
    print("RESULT:", "FAIL" if bad else "PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
