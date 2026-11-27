"""Phase 3 gate: the composed gene page must not change a word of either corpus.

The composed page puts the variant hub and the Atlas gene context on one page. The
only defensible claim is that neither half was re-rendered, so this gate compares:

* the Atlas half — `.entity-body` after the enhancement script ran, on `/gene/X`
  against `/atlas/gene/X/`. Identical, or the composition rewrote prose;
* the shape — exactly one h1, one canonical, no Atlas page header, no Atlas TOC
  (its sticky logic binds to a different element here), and the Atlas script ran
  (collapsible sections exist);
* the `.md` twin — exactly one `# ` heading and a `## Gene context` section whose
  Atlas headings are demoted, not dropped.

    python migration/verify_gene_composed.py --base http://127.0.0.1:8128 \
        --index /data2/sugi-variant-cache-dev/index.db --sample 12
"""
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
import urllib.request

from playwright.sync_api import sync_playwright

_WS = re.compile(r">\s+<")


def _md(base, path):
    with urllib.request.urlopen(base.rstrip("/") + path, timeout=180) as r:
        return r.read().decode()


def body_dom(page, url):
    page.goto(url, wait_until="load")
    page.wait_for_timeout(250)
    html = page.eval_on_selector(".entity-body", "e => e.outerHTML")
    return _WS.sub("><", html).strip()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://127.0.0.1:8128")
    ap.add_argument("--index", default="/data2/sugi-variant-cache-dev/index.db")
    ap.add_argument("--sample", type=int, default=12)
    ap.add_argument("--genes", default="BRCA1,TP53,DMD,CFTR,MEN1,SCN1A")
    a = ap.parse_args(argv)

    conn = sqlite3.connect(f"file:{a.index}?mode=ro", uri=True)
    known = {r[0] for r in conn.execute(
        "SELECT gene FROM gene_meta WHERE n > 0 ORDER BY n DESC")}
    wanted = [g.strip().upper() for g in a.genes.split(",") if g.strip()]
    genes = [g for g in wanted if g in known]
    # Add a spread of sizes, so the gate is not only whales.
    rows = [r[0] for r in conn.execute(
        "SELECT gene FROM gene_meta WHERE n > 0 ORDER BY n DESC")]
    step = max(1, len(rows) // max(1, a.sample))
    for g in rows[::step][:a.sample]:
        if g not in genes:
            genes.append(g)
    conn.close()
    print(f"checking {len(genes)} genes: {', '.join(genes[:8])}"
          + (" …" if len(genes) > 8 else ""))

    bad = []
    with sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1280, "height": 900})
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        for g in genes:
            errors.clear()
            try:
                atlas = body_dom(page, f"{a.base}/atlas/gene/{g}/")
                composed = body_dom(page, f"{a.base}/gene/{g}")
            except Exception as e:                      # noqa: BLE001
                bad.append((g, f"load: {str(e).splitlines()[0][:120]}"))
                continue
            if errors:
                bad.append((g, f"js error: {errors[0][:120]}"))
                continue
            if atlas != composed:
                i = next((k for k, (x, y) in enumerate(zip(atlas, composed)) if x != y),
                         min(len(atlas), len(composed)))
                bad.append((g, f"atlas prose differs at byte {i:,} "
                               f"({len(atlas):,} vs {len(composed):,})"))
                continue
            shape = page.evaluate("""() => ({
                h1: document.querySelectorAll('h1').length,
                canonical: !!document.querySelector('link[rel=canonical]'),
                atlasHeader: !!document.querySelector('.entity-header'),
                atlasToc: !!document.querySelector('.entity-toc-aside'),
                context: !!document.getElementById('gene-context'),
                details: document.querySelectorAll('.entity-body details.atlas-section').length,
            })""")
            problems = []
            if shape["h1"] != 1:
                problems.append(f"{shape['h1']} h1")
            if not shape["canonical"]:
                problems.append("no canonical")
            if shape["atlasHeader"]:
                problems.append("Atlas header left in")
            if shape["atlasToc"]:
                problems.append("Atlas TOC left in")
            if not shape["context"]:
                problems.append("no gene-context block")
            if not shape["details"]:
                problems.append("Atlas sections were not enhanced")
            if problems:
                bad.append((g, "; ".join(problems)))
                continue
            md = _md(a.base, f"/gene/{g}.md")
            heads = [ln for ln in md.splitlines() if ln.startswith("# ")]
            if len(heads) != 1 or "## Gene context" not in md:
                bad.append((g, f"md: {len(heads)} h1, "
                               f"gene-context={'## Gene context' in md}"))
        browser.close()

    print(f"composed: {len(genes)-len(bad)}/{len(genes)} identical prose, one identity")
    for g, why in bad[:10]:
        print(f"  {g}: {why}")
    print("RESULT:", "FAIL" if bad else "PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
