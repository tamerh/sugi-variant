"""Deterministic slugs for variant pages, from GENE + HGVS.

biobtree can't resolve HGVS free-text (`search("p.Pro334Arg")` → 0), so variant
pages are enumerated gene-first and the URL must be derivable from `GENE` + the
normalized HGVS — in BOTH the protein (p.) and coding (c.) forms, since the
search demand comes as both `acta1 "pro309ala"` and `"c.925c>g" acta1`.
Canonical = the p. form when present (the more-searched shape); the c. form is a
same-page alias.
"""
import hashlib
import re

_SLUG_MAX = 72

# ClinVar `name` looks like: "NM_001100.4(ACTA1):c.983_985del (p.Lys328del)".
# The nucleotide-level form is coding c. for protein genes, n. for non-coding RNA
# genes (RMRP, RNU4-2, TERC…), or m. for mitochondrial genes (MT-TL1 m.3243A>G).
# The negative lookbehind stops c/n/m matching inside a token (e.g. the M in NM_).
_C_RE = re.compile(r"(?<![A-Za-z])((?:c|n|m)\.[^\s)]+)")
_P_RE = re.compile(r"\(p\.([^)]+)\)")


_NAME_GENE_RE = re.compile(r"\(([^()]+)\):[cnmg]\.")


def name_gene(name):
    """The gene from a ClinVar NAME's transcript — 'NM_…(TTN):c.…' -> 'TTN'. This is
    the AUTHORITATIVE gene (the transcript the variant is annotated on). Preferred
    over ClinVar's primary gene_symbol, which for a genomic overlap can point at an
    antisense/lncRNA gene (e.g. TTN-AS1) and mis-file the host gene's variant."""
    m = _NAME_GENE_RE.search(name or "")
    return m.group(1) if m else None


def parse_hgvs(name):
    """(nucleotide_form, p_form) from a ClinVar `name`; either may be None. The
    nucleotide form is c./n./m. — so ncRNA and mitochondrial variants (which have no
    p. change) still yield a slug instead of being dropped."""
    c = _C_RE.search(name or "")
    p = _P_RE.search(name or "")
    return (c.group(1) if c else None), (f"p.{p.group(1)}" if p else None)


def _norm(s):
    """Slug-normalize a plain token (e.g. gene symbol): lowercase, non-alnum → '-'.
    A gene like NKX2-1 stays 'nkx2-1'."""
    s = (s or "").lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-")


def _norm_hgvs(s):
    """Slug-normalize an HGVS form, preserving the position operators that are
    MEANINGFUL and collide when dropped (audit P1): `*` (3'UTR) — `c.*667A>T` vs
    `c.667A>T`; and `+`/`-` (intronic direction) — `c.616+4` vs `c.616-4`. Map
    them to distinct word tokens before the non-alnum collapse. (In coding HGVS a
    bare '-' only ever means the intronic-upstream operator.)"""
    s = (s or "").lower()
    s = s.replace("*", " star ").replace("+", " plus ").replace("-", " minus ")
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-")


def _cap(slug):
    """Bound a slug's length. A large delins/ins insertion sequence (satellite DNA,
    etc.) would otherwise yield a 300+ char URL. Keep a readable prefix and append a
    short deterministic hash of the full form — stable across builds and collision-
    safe. Normal variants (well under the cap) are returned unchanged."""
    if len(slug) <= _SLUG_MAX:
        return slug
    h = hashlib.sha1(slug.encode()).hexdigest()[:8]
    return slug[:_SLUG_MAX - 9].rstrip("-") + "-" + h


def variant_slugs(gene, c_form, p_form):
    """(canonical_slug, [all_slugs]) — canonical is the p. form when present,
    else c.; every distinct form is emitted as an alias so either query shape
    lands the same page. None,[] if neither HGVS parses. Over-long slugs (huge
    delins insertions) are capped to a prefix + hash so URLs stay sane."""
    g = _norm(gene)
    if not g:
        return None, []
    out = []
    for form in (p_form, c_form):        # p. first → canonical
        if form:
            slug = _cap(f"{g}-{_norm_hgvs(form)}")
            if slug not in out:
                out.append(slug)
    return (out[0] if out else None), out
