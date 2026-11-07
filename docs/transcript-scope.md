# Scope: transcript-explicit variant pages

## Why

MANE Select is a *reporting convention*, not a claim that other isoforms do not
matter. Our own corpus shows how often that distinction is real: in a 250-gene
sample (52,188 records) **100%** of records carry more than one HGVS expression
and **68.4%** carry two or more distinct `NM_`/`NR_` accessions. A second sample
of 120 genes (34,804 records) puts the multi-accession share at **73.1%**.

Today the page flattens that into one line: c./p. from ClinVar's reported
transcript, MANE Select named, and the rest dumped as raw text. That is not
wrong, but it hides which transcript each statement belongs to, and it is
exactly where a reader with a different isoform loses confidence.

## What we have, and what we do not

`hgvs_expressions` is present on every record. Its shape (34,804-record sample):

| Prefix | Count | Meaning |
|---|---|---|
| `NM_` | 104,721 | RefSeq mRNA |
| `NC_` | 69,578 | genomic (chromosome / mitochondrial) |
| `NG_` | 37,927 | genomic (reference region) |
| `LRG_` | 35,246 | Locus Reference Genomic |
| `NR_` | 11,032 | RefSeq non-coding |

**Zero expressions carry a protein (`p.`) form.** The record has exactly one
p. change, on ClinVar's reported transcript. So we can table the per-transcript
**c./n.** change from data in hand, but a per-transcript **protein** change would
require translation we cannot do from a c. string alone (isoforms differ in CDS
and frame). That is an upstream ask, not a page decision.

## Design

One section, present on every page that has at least one `NM_`/`NR_` expression:

```
Transcripts
  NM_001754.5   c.952T>C    MANE Select · ClinVar-reported
  NM_001001890.3 c.871T>C   alternative
  NR_...        n.123A>G    non-coding transcript
  (genomic)     NC_000021.9:g.34799316A>G · NG_011402.2:g.1190396T>C · LRG_482:g....
```

Rules:

1. **One row per accession**, highest version present. When two versions of the
   same accession appear (`NM_001754.4`, `NM_001754.5`), keep the higher version
   and say so once — a version mismatch is itself information, not noise.
2. **Roles, not ranks**: `MANE Select` (accession equals `mane.mrna`),
   `ClinVar-reported` (`mane.clinvar_transcript`), `MANE Plus Clinical` (when the
   upstream flag lands). A row can carry more than one role.
3. **Never derive a protein change per transcript.** Show the c./n. change only.
   The page's single p. change stays where it is, explicitly attributed to the
   ClinVar-reported transcript.
4. **Genomic expressions go in one trailing row**, not as transcript rows: they
   describe the locus, not an isoform.
5. **State the boundary in one sentence**: numbering is transcript-specific; the
   consequence/exon context on this page is computed on the named transcript; we
   do not choose an interpretation between isoforms.
6. **Absent is stated, not silent**: a gene with no MANE Select (some
   non-coding genes) says so rather than showing an empty table.

## Where it goes

- HTML: a "Transcripts" block in the identity section, after the HGVS rows.
- `.md`: the same table as a Markdown table, replacing the current
  `**All HGVS expressions:** …` line.
- `jsonld.py`: `alternateName` already carries the expressions; leave as is.

Parity is guarded by `tests/test_surface_parity.py`; the block must render or be
absent in both surfaces together.

## Edge cases to cover

| Case | Behaviour |
|---|---|
| No `NM_`/`NR_`, only `NC_`/`NG_`/`LRG` | show the genomic row only, no transcript table |
| Only one transcript | show it (one row is still useful: it names the transcript) |
| Mitochondrial (`NC_012920`, `m.`) | genomic row, `m.` change |
| Non-coding transcript (`NR_`, `n.`) | transcript row with `n.` change |
| CNV / no c. form | no transcript table; current single line remains |
| Two versions of one accession | one row, higher version, "also submitted as NM_…​.4" |
| Accession in `hgvs_expressions` but not in `mane` | role `alternative`, never guessed |
| Gene with no MANE Select | explicit "no MANE Select for this gene" |

## Files touched

- `sugivariant/enrich.py` — `transcript_expressions(rec)`: parse, group, dedupe,
  role-tag. Pure, no network.
- `sugivariant/present.py` — the section state (present/absent, row list) so the
  template decides nothing.
- `templates/variant.html` + `sugivariant/render.py` (`render_body`) — render the
  same rows in both surfaces.
- `tests/test_transcripts.py` — parser + roles + version dedupe + edge cases.
- `tests/test_surface_parity.py` — extend for the new block.

## Effort

Half a day: the data is already in every record, the parsing is string work, and
no network or build change is involved. No record rebuild is needed — only a
renderer change, so it ships as a code-only release (`release` → `deploy`), no
corpus refresh.

## Out of scope here (see the biotree request)

- Per-transcript protein change and consequence (needs upstream annotation).
- Per-transcript exon numbering and NMD-escape (needs per-transcript exons).
- MANE Plus Clinical flag (needs the flag on the transcript relation).
- Any claim that one isoform is "the disease transcript": we use exactly three
  defensible anchors — MANE Select, MANE Plus Clinical, ClinVar-reported.
