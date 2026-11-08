# Resolver API (`/resolve.json`)

A small, read-only JSON endpoint for agents, MCP tools and scripts. It turns any
supported identifier into the canonical page URL(s), without scraping HTML.

```
GET https://sugi.bio/variant/resolve.json?q=<query>&gene=<gene>&slug=<slug>&limit=<n>
```

All parameters are optional, but at least one of `q`/`slug` should be given.
`limit` is capped at 50 (default 10).

## Response

```json
{
  "query": "PTEN R130Q",
  "kind": "variant",
  "ambiguous": false,
  "not_found": false,
  "release": "v1.0.0",
  "hits": [
    {
      "kind": "variant",
      "slug": "pten-p-arg130gln",
      "gene": "PTEN",
      "hgvs_p": "p.Arg130Gln",
      "hgvs_c": "c.389G>A",
      "classification": "Pathogenic",
      "review_status": "reviewed by expert panel",
      "rsid": "rs121909229",
      "coordinate": "10:87933148:G:A",
      "variation_id": 7829,
      "canonical_url": "https://sugi.bio/variant/pten-p-arg130gln",
      "markdown_url": "https://sugi.bio/variant/pten-p-arg130gln.md",
      "cached": true
    }
  ]
}
```

- `kind` — `variant`, `structural_variant`, `mixed`, or `null`.
- `ambiguous` — more than one hit. The API never silently chooses one.
- `not_found` — no hit. Distinct from a failed request.
- `release` — the release serving this answer, or `null` outside release mode.
- `cached` — whether the gene or SV is present in this deployment. In a release
  every hit is servable.

Structural hits carry `nssv`, `nsv`, `chromosome`, `start`, `end`, `svtype` and
`significance` instead of the small-variant fields, with the same
`canonical_url` / `markdown_url` pair.

## Gene queries

A bare gene symbol is a valid query and resolves to the gene hub rather than
`not_found`:

```json
{"kind": "gene", "not_found": false,
 "hits": [{"kind": "gene", "gene": "BRCA1", "variants": 11591,
           "canonical_url": "https://sugi.bio/variant/gene/BRCA1",
           "note": "gene-level match; add an HGVS, rsID or coordinate for a variant page"}]}
```

It is one page, so `ambiguous` is false — but the query is under-specified, and the
note says so. Do not treat a gene hit as a variant answer.

**A gene hit is fetchable too.** `/gene/<SYMBOL>.md` is the hub's Markdown twin, so a
client can treat every hit the same way: fetch `markdown_url`. The hit also carries
`variants` (the exact count), which answers "how many variants does BRCA1 have"
without fetching anything.

```
GET /variant/gene/BRCA1.md              # 64 KB, capped lists, counts exact
GET /variant/gene/BRCA1.md?section=Summary   # 523 B
```

Sections: `Summary`, `Variants`, `Evidence disagreements`, `Structural variants`,
`Sources` — listed in `X-Sugi-Sections`. Lists are capped and say how many they
omitted; the counts above them are exact.

## Fetching one section

`.md` pages accept `?section=<name>` (case-insensitive), which returns just that
block. This is what a batch needs: a full variant page is ~16k characters, while
`?section=Identity` is ~1.1k.

```
GET /variant/brca1-p-leu1570pro.md?section=Identity
```

Every `.md` response carries `X-Sugi-Sections` listing the available blocks. An
unknown section is a `404` whose body lists them, so a client can fall back
deliberately instead of silently receiving the whole page.

## A query can resolve to a different HGVS

Two ClinVar VCVs can share one page — most often the same protein change on
slightly different nucleotides, or the same variant annotated on different
transcripts. The page keeps one of them; the other's HGVS resolves to that page
through an alias. `BRCA1 c.68_69del` therefore returns `brca1-p-glu23fs`, whose
record carries `c.67del`. That is the same variant, not a wrong match: report the
page's own HGVS, not the one in the query.

## Structural variants: clinical corpus only

`/sv` covers dbVar entries with a clinical assertion, gated to the germline ACMG
classes. A valid-looking dbVar id outside that set (for example `nsv483076`, study
nstd41, BAC aCGH, no clinical assertion) is `not_found` — it is out of scope, not
missing.

## Keys

| Key | Example |
|---|---|
| Page slug | `pten-p-arg130gln` |
| rsID | `rs121909229` |
| GRCh38 coordinate | `10:87933148:G:A` |
| Gene + HGVS | `PTEN p.Arg130Gln`, `PTEN c.389G>A` |
| ClinVar accession | `VCV000007829` (small), `VCV000059108` (structural) |
| dbVar accession | `NSSV15120701`, `nsv3923526` |
| SV coordinate | `9:72224348-77332127` |
| Locus | `9:72500000-72900000` (all overlapping SVs) |

An rsID is position-level and can map to several variants, so it is a list, never
a single pick. A locus returns every SV whose span covers it.

## Recommended client behaviour

1. Call the resolver once per identifier.
2. If `not_found`, report that — do not guess or search a different identifier.
3. If `ambiguous`, show the list (gene, HGVS, classification, URL) and ask which
   one is meant. Do not summarise all candidates as if they were one variant.
4. Otherwise fetch `markdown_url`, which is the same page as `canonical_url` in
   Markdown for citation.

The `.md` twin and every response carry the release: the footer and the
`X-Sugi-Variant-Release`, `X-Sugi-Variant-Version` and
`X-Sugi-Variant-Released-At` headers.

## Scope

The resolver answers *where a variant's page is*, not *what it means*. The page
is the authoritative evidence view; the API does not add interpretation, scores
or a classification of its own.
