"""Phase 1b render gate: the app page must behave like the static Atlas page.

The fragment comparison (`verify_atlas_entities.py`) proves the *bytes* match. It
cannot see that Atlas enhances the entity body at runtime, in the page shell: H2
sections are collapsed into `<details>`, the TOC gains chevrons and a filter, tables
gain a Copy CSV button. That script is now carried by the app too, so this gate loads
both pages in a real browser and compares the article DOM after the script ran.

    python migration/verify_atlas_render.py \
        --base http://127.0.0.1:8115 \
        --root /data2/atlas-migration/atlas-html \
        --sample 25

Static page: file://<root>/<type>/<slug>/index.html
App page:    <base>/atlas/<type>/<slug>/   (run with BASE_PATH="" so links match)

The TOC's sticky `top` is a layout measurement of the surrounding shell, which the app
deliberately owns, so that one style attribute is normalised away.
"""
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

_TOP = re.compile(r'(class="?entity-toc-(?:aside|reopen)"?[^>]*?)style="top:[^"]*"')


def article_dom(page, url):
    page.goto(url, wait_until="load")
    page.wait_for_timeout(250)
    html = page.eval_on_selector("article.entity-page", "e => e.outerHTML")
    html = _TOP.sub(r"\1", html)
    # Collapse whitespace between tags so indentation in the shell cannot matter.
    return re.sub(r">\s+<", "><", html).strip()


def stats(page):
    return page.evaluate(
        """() => {
            const a = document.querySelector('article.entity-page');
            const d = [...a.querySelectorAll('details.atlas-section')];
            return {
              sections: d.length,
              open: d.filter(x => x.open).length,
              subtables: a.querySelectorAll('.table-wrapper').length,
              copy: a.querySelectorAll('.table-copy-btn').length,
              chevrons: a.querySelectorAll('.toc-chevron').length,
              filter: a.querySelectorAll('.toc-filter input').length,
              expand: a.querySelectorAll('.atlas-section-toggle').length,
            };
        }"""
    )


def check_search(page, base, query="imatinib"):
    """The header search must reach Atlas entities, and opening one must navigate.

    A JS change is invisible to the HTML gates: this types into the real search box,
    picks the first row with the keyboard, and checks where it lands.
    """
    page.goto(base.rstrip("/") + "/", wait_until="load")
    box = page.locator("textarea.hero-search-input, input[type=search][name=q]").first
    if box.count() == 0:
        print("search: no search box on this deployment (skipped)")
        return []
    box.click()
    box.type(query, delay=30)
    page.wait_for_timeout(900)
    rows = page.eval_on_selector_all(
        ".ac-item", "els => els.map(e => e.innerText.replace(/\\s+/g, ' ').trim())")
    if not any("atlas" in r.lower() for r in rows):
        return [f"search: no Atlas entity for {query!r} (rows: {rows[:3]})"]
    page.keyboard.press("ArrowDown")
    page.keyboard.press("Enter")
    page.wait_for_load_state("load")
    page.wait_for_timeout(300)
    if f"/atlas/" not in page.url:
        return [f"search: opened {page.url}, expected an /atlas/ page"]
    print(f"search: {query!r} -> {page.url.split('/', 3)[-1]}")
    return []


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://127.0.0.1:8115")
    ap.add_argument("--root", default="/data2/atlas-migration/atlas-html")
    ap.add_argument("--index", default="/data2/atlas-migration/atlas_index.db")
    ap.add_argument("--sample", type=int, default=25)
    ap.add_argument("--search", default="imatinib",
                    help="query for the search-box check (empty string skips it)")
    a = ap.parse_args(argv)

    root = Path(a.root)
    rows = sqlite3.connect(f"file:{a.index}?mode=ro", uri=True).execute(
        "SELECT type, slug FROM entity ORDER BY bytes DESC").fetchall()
    biggest = rows[:5]
    rest = rows[5:]
    step = max(1, len(rest) // max(1, a.sample - len(biggest)))
    sample = biggest + rest[::step]
    print(f"checking {len(sample):,} of {len(rows):,} entities "
          f"({len(biggest)} largest + every {step}th)")

    bad = []
    with sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1280, "height": 900})
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        for etype, slug in sample:
            static_url = (root / etype / slug / "index.html").resolve().as_uri()
            app_url = f"{a.base.rstrip('/')}/atlas/{etype}/{slug}/"
            try:
                errors.clear()
                want = article_dom(page, static_url)
                want_s = stats(page)
                got = article_dom(page, app_url)
                got_s = stats(page)
            except Exception as e:  # noqa: BLE001 - a broken page is a finding
                bad.append((slug, f"load: {str(e).splitlines()[0][:120]}"))
                continue
            if errors:
                bad.append((slug, f"js error: {errors[0][:120]}"))
                continue
            if want != got:
                first = next((i for i, (x, y) in enumerate(zip(want, got)) if x != y),
                             min(len(want), len(got)))
                bad.append((slug, f"dom differs at byte {first:,} "
                                  f"({len(want):,} vs {len(got):,})"))
                continue
            if want_s != got_s:
                bad.append((slug, f"enhancement stats differ {want_s} vs {got_s}"))
        if a.search:
            bad += [(f"search:{a.search}", why) for why in check_search(page, a.base, a.search)]
        browser.close()

    print(f"rendered article: {len(sample)-len(bad):,}/{len(sample):,} identical")
    for slug, why in bad[:10]:
        print(f"  {slug}: {why}")
    print("RESULT:", "FAIL" if bad else "PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
