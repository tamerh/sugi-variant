"""Phase 3: the composed gene page.

Two corpora, one page. The Variant gene hub supplies the evidence (counts, the variant
list, the disagreements view, overlapping structural variants); Atlas supplies the gene
context (function, structure, pathways, disease and drug links). This module does the
one thing that has to be mechanical: take the Atlas fragment the build produced and
remove the furniture the composed page supplies itself.

What is removed, and why:

* the breadcrumb and the `<h1>`/kind header — the composed page has its own gene
  header, and a second h1 would tell a crawler the page is two pages;
* the on-page navigation (the sticky aside, the mobile `<details>` TOC and the reopen
  button) — it lists Atlas's sections only, so on a page whose first half is the
  variant list it would be a table of contents for the second half, positioned above
  neither. The section list survives as the collapsible headings;
* the `container` wrapper — otherwise the Atlas block would nest one max-width inside
  another.

What is NOT touched: `.entity-body`, the sections, the tables, the About/Sources
footer. That is the whole point — `migration/verify_gene_composed.py` compares the
`.entity-body` DOM of `/gene/X` against `/atlas/gene/X/` after the enhancement script
has run, and it must be identical. If someone starts "improving" the prose here, that
gate fails, which is what it is for.

Nothing here re-derives or re-formats a clinical statement. Ordering and chrome only.
"""
from __future__ import annotations

import re

from sugivariant import atlas_web as ATW

_BREADCRUMB = re.compile(r"<p class=[\"']?breadcrumb[\"']?>.*?</p>", re.S)
_HEADER = re.compile(r"<header class=[\"']?entity-header[\"']?>.*?</header>", re.S)
_ASIDE = re.compile(r"<aside class=[\"']?entity-toc-aside[\"']?.*?</aside>", re.S)
_INLINE_TOC = re.compile(r"<details class=[\"']?entity-toc-inline[\"']?>.*?</details>", re.S)
_REOPEN = re.compile(r"<button[^>]*id=[\"']?entity-toc-reopen[\"']?[^>]*>.*?</button>", re.S)
_CONTAINER = '<div class="container container-wide">'


def strip_page_furniture(html: str) -> str:
    """The Atlas article without the parts the composed page provides."""
    for rx in (_BREADCRUMB, _HEADER, _ASIDE, _INLINE_TOC, _REOPEN):
        html = rx.sub("", html, count=1)
    return html.replace(_CONTAINER, '<div class="composed-atlas">', 1)


def demote_markdown(md: str) -> str:
    """Atlas markdown as a section of a bigger document.

    The `# Gene` heading goes (the composed page has its own), and every other heading
    drops one level so the Atlas sections nest under the composed section. The Atlas
    provenance lines (`> Gene · canonical`, `> Atlas vX · generated`) stay: they name
    the build, which is the one thing a reader quoting the prose needs.
    """
    lines = (md or "").splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)
    if lines and lines[0].startswith("# "):
        lines.pop(0)
    out = []
    for line in lines:
        if re.match(r"^#{2,6} ", line):
            line = "#" + line
        out.append(line)
    return "\n".join(out).strip("\n")


def gene_block(gene: str, base: str = "") -> dict | None:
    """The Atlas half of a composed gene page, or None when Atlas has no such gene.

    None is not an error: a gene can exist in the variant corpus and not in Atlas
    (different HGNC snapshots), and then the page is simply the variant hub.
    """
    if not (ATW.enabled() and ATW.cache_enabled()):
        return None
    if not ATW.has_entity("gene", gene):
        return None
    rec = ATW.entity("gene", gene)
    if not rec:
        return None
    return {
        "html": strip_page_furniture(ATW.rewrite_links(rec.get("body_html") or "", base)),
        "md": demote_markdown(rec.get("body_md") or ""),
        "slug": gene,
        "title": rec.get("title"),
        "built": (rec.get("canonical") or ""),
    }
