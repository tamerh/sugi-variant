"""schema.org JSON-LD for variant pages — a machine-readable structured record
+ a FAQPage keyed to the literal 'Is GENE HGVS pathogenic?' question, so an AI
assistant can lift the answer verbatim and attribute it. Every fact traces to a
source; the block carries provenance + a dateModified.

schema.org has no genetic-variant type, so the variant is modelled as a
MedicalEntity carrying its identifiers (HGVS / rsID / ClinVar VCV) + sameAs to
the authoritative databases; the FAQPage is the citation payload.
"""
import json

_HOST = "https://sugi.bio"


def _url(rec):
    return f"{_HOST}/atlas/variant/{rec['canonical_slug']}/"


def _label(rec):
    from sugivariant.render import _label as lbl
    return lbl(rec)


def _identifiers(rec):
    out = []
    if rec.get("rsid"):
        out.append({"@type": "PropertyValue", "propertyID": "dbSNP", "value": rec["rsid"]})
    out.append({"@type": "PropertyValue", "propertyID": "ClinVar",
                "value": f"VCV{rec['variation_id']}"})
    for form, sys in ((rec.get("hgvs_p"), "HGVS.p"), (rec.get("hgvs_c"), "HGVS.c")):
        if form:
            out.append({"@type": "PropertyValue", "propertyID": sys, "value": form})
    return out


def _same_as(rec):
    urls = [f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{rec['variation_id']}/"]
    if rec.get("rsid"):
        urls.append(f"https://www.ncbi.nlm.nih.gov/snp/{rec['rsid']}/")
    return urls


def _faq(rec):
    """Deterministic Q&A pairs — only those we can actually answer from the
    record. The lead 'is it pathogenic?' question mirrors the machine-query shape."""
    from sugivariant.render import declarative_plain
    label = _label(rec)
    qas = [(f"Is {label} pathogenic?", declarative_plain(rec))]

    conds = [c.get("name") for c in (rec.get("conditions") or []) if c.get("name")]
    if conds:
        qas.append((f"What condition is associated with {label}?",
                    f"{label} is associated with {', '.join(conds[:3])} (ClinVar)."))

    g = rec.get("gnomad")
    if g:
        qas.append((f"How common is {label} in the general population?",
                    f"{label} is {g['band']} (gnomAD)."))

    dg = rec.get("digest") or {}
    inh = dg.get("inheritance") or []
    if inh:
        cond = dg.get("name") or "this condition"
        qas.append((f"How is {cond} inherited?",
                    f"Typically {', '.join(inh[:2])} (Orphanet). A genetic "
                    "counselor can assess individual risk for a family."))

    am = rec.get("alphamissense")
    if am and am.get("class"):
        qas.append((f"What do computational predictors say about {label}?",
                    f"AlphaMissense predicts {am['class'].replace('_', ' ')} "
                    f"(score {am['score']}) — an in-silico prediction, not a clinical call."))

    return {
        "@type": "FAQPage",
        "@id": _url(rec) + "#faq",
        "mainEntity": [
            {"@type": "Question", "name": q,
             "acceptedAnswer": {"@type": "Answer", "text": a}}
            for q, a in qas],
    }


def build_jsonld(rec, meta=None):
    """Full @graph: the variant MedicalEntity + the FAQPage citation payload."""
    url = _url(rec)
    entity = {
        "@type": "MedicalEntity",
        "@id": url + "#variant",
        "name": _label(rec),
        "description": _plain_lead(rec),
        "identifier": _identifiers(rec),
        "sameAs": _same_as(rec),
        "url": url,
    }
    alt = [x for x in (rec.get("hgvs_p"), rec.get("hgvs_c")) if x] + (rec.get("hgvs_expressions") or [])
    if alt:
        entity["alternateName"] = alt
    gene = rec.get("gene_symbol")
    if gene:
        from sugivariant import links
        entity["isPartOf"] = {"@type": "Gene", "name": gene,
                              "url": _HOST + (links.gene_url(symbol=gene, hgnc_id=rec.get("hgnc_id")) or f"/atlas/gene/{gene}/")}
    if meta and meta.get("generated_at"):
        entity["dateModified"] = meta["generated_at"]
    # provenance
    cites = ["https://www.ncbi.nlm.nih.gov/clinvar/"]
    if rec.get("alphamissense"):
        cites.append("https://alphamissense.hegelab.org/")
    if rec.get("gnomad"):
        cites.append("https://gnomad.broadinstitute.org/")
    entity["citation"] = cites

    return {"@context": "https://schema.org", "@graph": [entity, _faq(rec)]}


def _plain_lead(rec):
    from sugivariant.render import declarative_plain
    return declarative_plain(rec)


def as_script_tag(rec, meta=None):
    """Inline <script type=application/ld+json> block for the variant page body."""
    body = json.dumps(build_jsonld(rec, meta), indent=2)
    return f'<script type="application/ld+json">\n{body}\n</script>'
