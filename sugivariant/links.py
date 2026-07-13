"""Cross-link URL helpers — standalone shim (decoupled from atlas.page.links).
Gene/disease links point at the public Sugi Atlas site (they exist there); no
manifest gating is needed, so maybe_link always links and load() is a no-op."""
import re

_ATLAS = "https://sugi.bio/atlas"


def gene_url(symbol=None, hgnc_id=None):
    return f"{_ATLAS}/gene/{symbol}/" if symbol else None


def disease_url(mondo_id=None, name=None):
    if not name:
        return None
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
    return f"{_ATLAS}/disease/{slug}/"


def maybe_link(text, url):
    if not text:
        return text or ""
    return f"[{text}]({url})" if url else str(text)


def load(dist_root=None):
    """No-op: the standalone links to the external atlas, no local manifest."""
    return None
