"""Deterministic enrichment for variant pages — the layers that turn a ClinVar
echo into an integrative reference. All computed from primary data:

  - AlphaMissense in-silico pathogenicity (per missense, joined by protein_variant)
  - gnomAD population frequency + a rarity band (joined by rsID)
  - cross-source CONCORDANCE verdict (ClinVar vs AlphaMissense vs gnomAD)
  - submitter-consensus statistics (from the ClinVar submissions[])
  - same-residue hotspot context (from the gene's own ClinVar enumeration)

Nothing here is a clinical call: concordance/consensus are transparent
descriptions of independent evidence, each shown with its source.
"""
import re
from collections import Counter

from sugibiobtree import map_all

_AA3to1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
           "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
           "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
           "Tyr": "Y", "Val": "V", "Ter": "*"}
_MISSENSE_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def missense_short(hgvs_p):
    """'p.Pro309Ala' → 'P309A' (AlphaMissense key). None for non-missense
    (frameshift/del/dup/nonsense) — AlphaMissense only models missense SNVs."""
    m = _MISSENSE_RE.match((hgvs_p or "").strip())
    if not m:
        return None
    a1, a2 = _AA3to1.get(m.group(1)), _AA3to1.get(m.group(3))
    return f"{a1}{m.group(2)}{a2}" if a1 and a2 else None


def protein_position(hgvs_p):
    """Integer residue position from any p. HGVS ('p.Pro309Ala'→309), or None."""
    m = re.search(r"p\.[A-Za-z]{3}(\d+)", hgvs_p or "")
    return int(m.group(1)) if m else None


# ── per-gene caches (fetch once, reuse across all the gene's variants) ────────
def gene_alphamissense(hgnc_id):
    """{protein_short: (am_class, am_pathogenicity)} for the whole gene."""
    out = {}
    # uncapped: the percentile stats must see the WHOLE distribution — a big
    # protein (e.g. TTN, ~100k substitutions) far exceeds the old 60-page/6k cap,
    # which silently skewed its am_percentile. Runs once per gene (cached).
    for r in map_all(hgnc_id, ">>hgnc>>uniprot>>alphamissense", cap=None):
        pv = (r.get("protein_variant") or "").strip()
        if pv:
            out[pv] = (r.get("am_class"), r.get("am_pathogenicity"))
    return out


_GENOMIC_SNV = re.compile(r"NC_0*\d+\.\d+:g\.(\d+)([ACGT]+)>([ACGT]+)$")


def variant_coordinate(rec):
    """GRCh38 gnomAD-style key 'chr:pos:ref:alt' (no chr prefix). SNVs only.

    The variant carries genomic HGVS in BOTH assemblies (NC_..10 = GRCh37,
    NC_..11 = GRCh38); to pick the GRCh38 one unambiguously we require the HGVS
    position to equal the ClinVar `start` (which ClinVar reports on GRCh38).
    Using the wrong assembly would key gnomAD wrong and read as a false absence.
    """
    start = rec.get("start")
    chrom = str(rec.get("chromosome") or "").strip()
    if not (start and chrom):
        return None
    for e in rec.get("hgvs_expressions") or []:
        m = _GENOMIC_SNV.match(e)
        if m and m.group(1) == str(start):
            return f"{chrom}:{start}:{m.group(2)}:{m.group(3)}"
    return None


def _coord_entry(coord, dataset):
    """Attributes dict for a coordinate-keyed dataset via entry() — the ONLY
    working access for gnomad_variant / alphamissense / conservation (they are
    NOT map-chainable: `map(coord, '>>gnomad_variant')` returns 0). Audit Tier 1/2."""
    if not coord:
        return None
    from sugibiobtree import entry
    try:
        a = (entry(coord, dataset) or {}).get("Attributes") or {}
    except Exception:
        return None
    if not a:
        return None
    # single-key Attributes wrapper (e.g. {"GnomadVariant": {...}}); biobtree
    # returns {"Empty": true} for a key with no data — guard to dict-only.
    v = next(iter(a.values())) if len(a) == 1 else a
    return v if isinstance(v, dict) else None


def gnomad_frequency(rec):
    """gnomAD v4.1 per-variant frequency — the ACMG BA1/BS1/PM2 layer. Looked up
    by coordinate via entry() (NOT map — that returns 0; audit Tier 1 false-Absent
    bug). Surfaces popmax + faf + per-ancestry. Falls back to the dbSNP global
    frequency only when no genomic coordinate can be parsed (e.g. indels)."""
    coord = variant_coordinate(rec)
    if coord:
        g = _coord_entry(coord, "gnomad_variant")
        if g:
            popmax = _f(g.get("af_grpmax"))
            anc = g.get("grpmax_ancestry")
            pops = {k[3:]: g[k] for k in g if k.startswith("af_") and k not in ("af_grpmax", "af_exomes", "af_genomes") and g.get(k)}
            # Correct the AF denominator: the combined `af` adds a callset's AN even when the
            # variant is absent from it, so exome-only / genome-only variants under-report
            # (COL4A1 P484A: combined 3.72e-6 vs gnomAD's exome-only 4.10e-6). Recompute from the
            # per-callset AC/AN, summing only callsets that actually observed the variant (ac>0).
            ace, ane = _f(g.get("ac_exomes")), _f(g.get("an_exomes"))
            acg, ang = _f(g.get("ac_genomes")), _f(g.get("an_genomes"))
            obs_ac = obs_an = 0.0
            if ace:
                obs_ac += ace; obs_an += (ane or 0)
            if acg:
                obs_ac += acg; obs_an += (ang or 0)
            if obs_an:                       # per-callset available → corrected AF + honest AC/AN
                af, ac, an = obs_ac / obs_an, int(obs_ac), int(obs_an)
            else:                            # fall back to the combined values
                af = _f(g.get("af"))
                ac = int(_f(g.get("ac"))) if _f(g.get("ac")) is not None else None
                an = int(_f(g.get("an"))) if _f(g.get("an")) is not None else None
            return {"af": af, "popmax": popmax, "ancestry": anc,
                    "faf": g.get("faf"), "faf99": _f(g.get("faf99")),
                    "ac": ac, "an": an,
                    "ac_grpmax": _f(g.get("ac_grpmax")), "an_grpmax": _f(g.get("an_grpmax")),
                    "populations": pops,
                    "absent": False, "is_common": (popmax or 0) >= 0.05,
                    "band": _gnomad_band(af, popmax, anc), "source": "gnomAD v4.1"}
        return {"absent": True, "is_common": False, "popmax": None,
                "band": "absent from gnomAD v4.1", "source": "gnomAD v4.1"}
    # fallback — dbSNP inline global frequency (no coordinate to key gnomAD v4.1)
    return gnomad_for(rec.get("rsid"))


def alphamissense_for(coord):
    """AlphaMissense for the variant, looked up BY COORDINATE (audit Tier 2: the
    per-gene protein-keyed map uses a different isoform's numbering, so ~43% of
    ASXL1 missense scores were missed). {class, score, short} or None."""
    a = _coord_entry(coord, "alphamissense")
    if not a or a.get("am_pathogenicity") is None:
        return None
    return {"class": a.get("am_class"), "score": str(a.get("am_pathogenicity")),
            "short": a.get("protein_variant"), "uniprot": a.get("uniprot_id")}


_DATASET_VERSIONS = None


def dataset_versions():
    """{biobtree group → last_built date YYYY-MM-DD} from the live meta — the HONEST
    per-dataset build date. HANDOVER §7: the single biobtree 'dev' version stamp lies,
    so pin per-dataset last_built instead. Cached; {} if meta is unreachable."""
    global _DATASET_VERSIONS
    if _DATASET_VERSIONS:                      # cache only a SUCCESSFUL (non-empty) fetch —
        return _DATASET_VERSIONS               # a transient failure must not permanently disable this
    out = {}
    try:
        import urllib.request
        import json as _json
        import os
        # BIOBTREE_WS overrides; otherwise use the SAME endpoint the rest of the app
        # reaches biobtree on (ATLAS_BIOBTREE). Without this the container queries
        # localhost:9291 (nothing there) and every version stamp goes blank in prod.
        base = (os.environ.get("BIOBTREE_WS")
                or os.environ.get("ATLAS_BIOBTREE", "http://localhost:9291"))
        with urllib.request.urlopen(base + "/ws/meta", timeout=5) as r:
            meta = _json.load(r)
        for v in (meta.get("datasets") or {}).values():
            g, lb = v.get("group"), v.get("last_built")
            if g and lb and g not in out:
                out[g] = lb[:10]
    except Exception:
        return {}                              # leave the cache empty → retried on the next call
    _DATASET_VERSIONS = out
    return out


def molecular_consequence(rec):
    """Deterministic molecular consequence from HGVS grammar — DESCRIPTIVE typing,
    never an applied PVS1/ACMG code (§8). Predicted-LoF classes (nonsense /
    frameshift / canonical splice ±1,2 / start-loss) set lof=True. {type, label, lof}
    or None. Mainly for the non-missense ~51% where the predictor panel is dark."""
    p = rec.get("hgvs_p") or ""
    c = rec.get("hgvs_c") or ""
    if "fs" in p:
        return {"type": "frameshift", "label": "frameshift", "lof": True}
    if re.match(r"p\.Met1(\?|=|[A-Za-z]{3})", p):
        return {"type": "start_loss", "label": "start-loss", "lof": True}
    if re.match(r"p\.Ter\d", p) or "ext" in p:
        return {"type": "stop_loss", "label": "stop-loss", "lof": True}
    if "Ter" in p or re.search(r"\*\d*$", p):
        return {"type": "nonsense", "label": "nonsense (stop-gain)", "lof": True}
    if re.search(r"\d[+-][12](?![0-9])", c):
        return {"type": "splice", "label": "canonical splice-site (±1/±2)", "lof": True}
    if "delins" in p:
        return {"type": "inframe_delins", "label": "in-frame delins", "lof": False}
    if any(k in p for k in ("del", "dup", "ins")) and "fs" not in p:
        return {"type": "inframe_indel", "label": "in-frame indel", "lof": False}
    return None


# Somatic/myeloid disease contexts where a germline dosage-mechanism narrative would
# be §8-inappropriate (the ASXL1/CHIP audit case). Deliberately narrow — clear somatic
# drivers only; hereditary-cancer conditions are germline and must NOT match.
_SOMATIC_COND_RE = re.compile(r"leukemi|myelodysplas|myeloproliferat|myeloid|clonal h", re.I)


def lof_context(rec):
    """Descriptive loss-of-function mechanism read: a predicted-LoF consequence
    paired with the evidence that LoF causes disease for this gene. The PRIMARY
    signal is ClinGen dosage haploinsufficiency (curated: score 3 = sufficient
    evidence); gnomAD constraint (LOEUF/pLI) is supporting population evidence.
    NOT an applied PVS1 code (§8) — it never assigns ACMG evidence. Returns
    {label, loeuf, pli, haplo, haploinsufficient, constrained, lof_disease_gene}
    or None."""
    cons = rec.get("consequence")
    if not cons or not cons.get("lof"):
        return None
    gc = rec.get("gene_context") or {}
    con, dos = gc.get("constraint") or {}, gc.get("dosage") or {}
    loeuf, pli, haplo = con.get("loeuf"), con.get("pli"), dos.get("haplo")

    def _int(x):
        try:
            return int(float(x))
        except (TypeError, ValueError):
            return None

    haploinsufficient = _int(haplo) == 3          # ClinGen: sufficient evidence
    try:
        constrained = ((loeuf is not None and float(loeuf) < 0.35)
                       or (pli is not None and float(pli) >= 0.9))
    except (TypeError, ValueError):
        constrained = False
    # §8: a germline dosage / LoF-disease-mechanism narrative must NOT attach to a
    # SOMATIC-condition variant (the ASXL1/CHIP hazard). Suppress when the variant's
    # primary condition is a somatic/myeloid context; a germline Orphanet digest always
    # overrides. (Gating on the digest alone was too strict — it dropped germline genes
    # like LDLR whose condition has no Orphanet Disease entry.)
    primary = ((rec.get("conditions") or [{}])[0].get("name") or "")
    germline = not _SOMATIC_COND_RE.search(primary)   # the DISPLAYED (primary) condition governs
    # NMD: for a truncating change, LoF impact is position-dependent (C-terminal /
    # last-exon truncations may escape nonsense-mediated decay). We have no exon model,
    # so flag the caveat rather than assert clean LoF.
    nmd_caveat = cons["type"] in ("nonsense", "frameshift")
    return {"label": cons["label"], "loeuf": loeuf, "pli": pli, "haplo": haplo,
            "haploinsufficient": haploinsufficient, "constrained": constrained,
            "germline": germline, "nmd_caveat": nmd_caveat,
            "lof_disease_gene": germline and (haploinsufficient or constrained)}


def _fix_acmg_modifier(code):
    """Repair an ACMG strength modifier truncated in the biobtree staging data: 'PM3_Very'
    is not a valid code (the SVI closed set is _Supporting/_Moderate/_Strong/_VeryStrong), so
    it's a cut of 'PM3_VeryStrong'. Only the clearly-invalid trailing '_Very' is repaired;
    valid codes pass through unchanged. (Temporary guard pending an upstream fix.)"""
    return code[:-5] + "_VeryStrong" if isinstance(code, str) and code.endswith("_Very") else code


_ACMG_FULL = re.compile(r"\b([PB][MVSP]\d)(_[A-Za-z]+)\b")


def _reconcile_codes(codes, summary):
    """Repair truncated/fully-stripped ACMG modifiers by reconciling each applied code against
    the rationale summary, which carries the full form (benchmark 2026-07: chip 'PS4' shown
    while the verbatim rationale says 'PS4_VeryStrong'). Fixes the '_Very' cut too."""
    full = {b: b + m for b, m in _ACMG_FULL.findall(summary or "")}
    out = []
    for code in codes:
        code = _fix_acmg_modifier(code)
        if isinstance(code, str) and "_" not in code and code in full:
            code = full[code]                        # bare code → full form from the rationale
        out.append(code)
    return out


def clingen_criteria(ca_id):
    """Full ClinGen VCEP variant-pathogenicity record by allele-registry (CA) id:
    the applied ACMG codes, per-criterion rationale, and provenance — the authority
    tier. We REPORT the expert panel's applied codes verbatim; we NEVER re-tally or
    combine them into our own score (§8, ClinGen SVI). None if unavailable."""
    if not ca_id:
        return None
    from sugibiobtree import entry
    try:
        a = (entry(ca_id, "clingen_variant") or {}).get("Attributes") or {}
    except Exception:
        return None
    c = a.get("ClingenVariant") or {}
    if not c.get("evidence_codes_met"):
        return None
    return {"codes_met": _reconcile_codes(c.get("evidence_codes_met") or [], c.get("summary")),
            "codes_not_met": _reconcile_codes(c.get("evidence_codes_not_met") or [], c.get("summary")),
            "summary": c.get("summary"), "moi": c.get("moi"),
            "guideline": c.get("guideline"), "approval_date": c.get("approval_date"),
            "published_date": c.get("published_date"),
            "erepo": c.get("evidence_repo_link")}


def am_isoform_mismatch(am, hgvs_p):
    """AlphaMissense is looked up by genomic coordinate and returns ITS transcript's
    protein change. If that residue disagrees with the ClinVar p.HGVS, AM is numbered
    on a DIFFERENT isoform than the reported variant (the known Tier-2 hazard) — a QC
    flag so the AM read isn't silently trusted against the wrong residue. Never
    suppresses the datum; only flags the disagreement. {am, clinvar} or None."""
    if not am or not hgvs_p:
        return None
    am_short = (am.get("short") or "").upper()
    cv_short = (missense_short(hgvs_p) or "").upper()
    if am_short and cv_short and am_short != cv_short:
        return {"am": am.get("short"), "clinvar": missense_short(hgvs_p)}
    return None


_NM_ACCESSION_RE = re.compile(r"^\s*([NX][MR]_\d+)")


def clinvar_transcript(name):
    """Base RefSeq transcript accession (version stripped) from a ClinVar variant name
    like 'NM_000546.6(TP53):c.328C>G' -> 'NM_000546'. None if the name has no NM_/NR_
    prefix (e.g. genomic-only names)."""
    m = _NM_ACCESSION_RE.match(name or "")
    return m.group(1) if m else None


def mane_select(hgnc_id):
    """MANE Select transcript for a gene, read the way sugi-atlas reads it — the
    is_mane_select flag lives on the RefSeq relation, not the Ensembl-transcript one:
    map_all(hgnc, '>>hgnc>>ensembl>>refseq[is_mane_select==true]'). Returns
    {mrna, protein} (RefSeq accessions, no version) or None. MANE Select is the single
    NCBI+EMBL-EBI-agreed reference transcript for clinical reporting."""
    if not hgnc_id:
        return None
    rows = map_all(hgnc_id, ">>hgnc>>ensembl>>refseq[is_mane_select==true]") or []
    mrna = next((t["id"] for t in rows if t.get("type") == "mRNA"), None)
    if not mrna:
        return None
    prot = next((t["id"] for t in rows if t.get("type") == "protein"), None)
    return {"mrna": mrna, "protein": prot}


# ClinVar-calibrated SaProt damaging divider (Youden-optimal; tools/eval/saprot_calibration.py,
# 2026-07-20). Was a -7.5 heuristic; -10.0 lifts specificity 0.69 -> 0.94.
_SAPROT_DAMAGING = -10.0


def saprot_for(uniprot, protein_variant):
    """SaProt-650M structure-aware protein-language-model variant effect (Su et
    al. 2023; computed in-house from MIT weights → redistributable), keyed by
    uniprot:protein_variant (both from AlphaMissense). LLR ≤ 0, more negative =
    more damaging; NO calibrated ACMG threshold — surface the raw LLR and let the
    concordance panel do the work. The 'damaging vs tolerated' divider is used ONLY
    for the agreement read, and is now DATA-CALIBRATED against ClinVar (not a guess):
    on 476 ClinVar P/LP-vs-B/LB missense variants SaProt LLR has AUC 0.917; the
    Youden-optimal cut is -10.0 (sens 0.80, spec 0.94) vs the old heuristic -7.5
    (sens 0.91, spec 0.69 — over-called benign). See tools/eval/saprot_calibration.py.
    Unsupervised → a genuinely orthogonal predictor next to supervised AlphaMissense/
    REVEL. Covers ~98.6% of the proteome; the ~1.4% without an AlphaFold structure
    fall back to AlphaMissense (SaProt simply returns None here)."""
    if not (uniprot and protein_variant):
        return None
    a = _coord_entry(f"{uniprot}:{protein_variant}", "saprot")
    llr = _f((a or {}).get("saprot_llr"))
    if llr is None:
        return None
    return {"llr": llr, "damaging": llr <= _SAPROT_DAMAGING}


def conservation_for(coord):
    """Per-position evolutionary conservation (phyloP / GERP++ / phastCons) via
    entry() on the chr:pos (ref/alt-agnostic) key. The ONLY computational signal
    for the non-missense/splice/intronic pages. GERP may be missing where a
    position isn't in the mammalian alignment → None, not zero."""
    if not coord:
        return None
    pos = ":".join(coord.split(":")[:2])          # chr:pos:ref:alt → chr:pos
    a = _coord_entry(pos, "conservation")
    if not a:
        return None
    phylop = _f(a.get("phylop"))
    phastcons = _f(a.get("phastcons"))
    gerp = _f(a.get("gerp"))
    if phylop is None and phastcons is None and gerp is None:
        return None
    return {"phylop": phylop, "phastcons": phastcons, "gerp": gerp}


# REVEL → descriptive strength BAND (NOT an ACMG code). The thresholds are the
# Pejaver-2022 calibration points, but §8 forbids presenting REVEL as an ACMG
# PP3/BP4 vote (it is an agreement signal, not additive evidence). So we describe
# the strength of REVEL's own leaning, never emit a PP3/BP4 code.
def _revel_band(s):
    if s >= 0.932:
        return "strongly pathogenic-leaning"
    if s >= 0.773:
        return "moderately pathogenic-leaning"
    if s >= 0.644:
        return "pathogenic-leaning"
    if s > 0.290:
        return "indeterminate"
    if s >= 0.183:
        return "benign-leaning"
    if s > 0.016:
        return "moderately benign-leaning"
    return "strongly benign-leaning"


def revel_for(coord):
    """REVEL ensemble missense pathogenicity (0–1) by chr:pos:ref:alt, with a
    descriptive strength band (Pejaver-2022 thresholds) + direction. Presented as
    an AGREEMENT signal only — never an additive ACMG PP3/BP4 code (§8, ClinGen SVI)."""
    a = _coord_entry(coord, "revel")
    s = _f((a or {}).get("revel"))
    if s is None:
        return None
    direction = "pathogenic" if s >= 0.644 else "benign" if s <= 0.290 else "indeterminate"
    return {"score": s, "band": _revel_band(s), "direction": direction}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _gnomad_band(af, popmax, ancestry):
    """ACMG-flavoured frequency band from gnomAD popmax (the ACMG-standard
    metric). Descriptive, not a clinical criterion call."""
    p = popmax if popmax is not None else af
    if not p:
        return "absent from gnomAD v4.1"
    # ONE unit for the popmax value everywhere (audit P2a): a percentage, 3 s.f.
    pct = f"popmax {p * 100:.3g}%" + (f", {ancestry}" if ancestry else "")
    if p >= 0.05:
        return f"common ({pct}) — too common for a highly-penetrant pathogenic allele"
    if p >= 0.01:
        return f"low-frequency ({pct})"
    if p >= 1e-3:
        return f"rare ({pct})"
    if p >= 1e-4:
        return f"ultra-rare ({pct})"
    return f"very rare ({pct})"


def _sci(x):
    return f"{x:.1e}"


def gnomad_for(rsid):
    """Population-frequency read from the dbSNP inline global frequency (fallback
    when no coordinate is available). Absence is itself a signal."""
    if not rsid:
        return None
    d = map_all(rsid, ">>dbsnp")
    if not d:
        return None
    freq = (d[0].get("gnomad_frequency") or "").strip()
    absent = freq in ("", "0", "0.0")
    # Return the SAME key set as the coordinate path (gnomad_frequency) — a divergent
    # shape here crashed the variant template on every fallback record with a frequency
    # (Jinja: a missing key is Undefined, and `Undefined is not none` is True).
    return {"frequency": freq, "af": _f(freq), "absent": absent,
            "popmax": None, "ancestry": None, "faf": None, "faf99": None,
            "ac": None, "an": None, "ac_grpmax": None, "an_grpmax": None,
            "populations": {},
            "is_common": (d[0].get("is_common") == "true"),
            "band": _freq_band(freq, absent), "source": "dbSNP/gnomAD"}


def _freq_band(freq, absent):
    if absent:
        return "absent from gnomAD"
    try:
        f = float(freq)
    except ValueError:
        return "present in gnomAD"
    if f < 1e-4:
        return f"ultra-rare (gnomAD MAF {freq})"
    if f < 1e-3:
        return f"rare (gnomAD MAF {freq})"
    if f < 1e-2:
        return f"low-frequency (gnomAD MAF {freq})"
    return f"common (gnomAD MAF {freq})"


# ── derived analyses ─────────────────────────────────────────────────────────
def submitter_consensus(submissions):
    """Agreement among ClinVar submitters — from submissions[]. Returns a small
    dict {n, breakdown, verdict} or None."""
    calls = [(s.get("classification") or "").strip() for s in submissions if s.get("classification")]
    if not calls:
        return None
    bd = Counter(calls)
    top, topn = bd.most_common(1)[0]
    if len(bd) == 1:
        verdict = f"unanimous ({topn}/{len(calls)} {top})"
    else:
        parts = ", ".join(f"{n} {c}" for c, n in bd.most_common())
        verdict = f"conflicting ({parts})"
    return {"n": len(calls), "breakdown": dict(bd), "verdict": verdict}


def concordance(classification, am, gnomad, spliceai=None, conservation=None, revel=None, saprot=None):
    """Cross-source concordance readout: do the computational predictors agree with
    ClinVar? Returns {lines, verdict, flags}. ClinGen-SVI framing: in-silico
    predictors are NOT independent, so exactly ONE calibrated tool carries the
    ACMG weight — AlphaMissense for missense, conservation for non-missense — and
    REVEL (+ conservation, on missense) is shown as an AGREEMENT signal, never an
    additive vote. Population rarity is shown but never counted (audit P2). Never a
    clinical determination."""
    cls = (classification or "").lower()
    clinvar_path = "pathogenic" in cls and "conflict" not in cls
    clinvar_vus = "uncertain" in cls
    lines, flags, agree, total, consensus = [], [], 0, 0, None

    if am:
        amclass, amscore = am.get("class"), am.get("score")
        am_path = amclass == "likely_pathogenic"
        total += 1
        pretty = (amclass or "").replace("_", " ")
        if clinvar_path and am_path:
            agree += 1
            lines.append(f"AlphaMissense **{pretty}** ({amscore}) — concordant with the pathogenic call")
        elif clinvar_path and amclass == "likely_benign":
            flags.append("AlphaMissense predicts likely-benign")
            lines.append(f"AlphaMissense **{pretty}** ({amscore}) — ⚠ discordant with the pathogenic call")
        else:
            lines.append(f"AlphaMissense **{pretty}** ({amscore})")
        # REVEL — an AGREEMENT check on AlphaMissense (NOT a second vote; ClinGen SVI)
        if revel:
            if revel["direction"] == "indeterminate":
                rel = "indeterminate"
            else:
                rel = "agrees with" if (revel["direction"] == "pathogenic") == am_path else "differs from"
            lines.append(f"REVEL {revel['score']} ({revel['band']}) — {rel} AlphaMissense")
        # SaProt — orthogonal (ClinVar-independent) protein-language-model opinion,
        # also an agreement signal (not additive). Raw LLR surfaced; the ClinVar-
        # calibrated -10.0 cut divides the damaging/tolerated *read* only.
        if saprot:
            rel = "agrees with" if saprot["damaging"] == am_path else "differs from"
            lines.append(f"SaProt LLR {saprot['llr']:g} "
                         f"({'damaging' if saprot['damaging'] else 'tolerated'}) — {rel} AlphaMissense")
        # In-silico CONSENSUS — descriptive agreement of the independent missense
        # predictors (do they concur?). NOT additive ACMG evidence (ClinGen SVI).
        preds = [("AlphaMissense", am_path)]
        if revel and revel["direction"] != "indeterminate":
            preds.append(("REVEL", revel["direction"] == "pathogenic"))
        if saprot:
            preds.append(("SaProt", saprot["damaging"]))
        if len(preds) >= 2:
            n_dmg = sum(1 for _, d in preds if d)
            names = ", ".join(n for n, _ in preds)
            alln = "both" if len(preds) == 2 else f"all {len(preds)}"
            if n_dmg == len(preds):
                summary = f"{alln} predictors call this damaging"
            elif n_dmg == 0:
                summary = f"{alln} predictors call this tolerated"
            else:
                summary = f"{n_dmg}/{len(preds)} damaging — mixed"
            consensus = {"n_damaging": n_dmg, "total": len(preds),
                         "names": names, "summary": summary,
                         "unanimous": n_dmg in (0, len(preds))}
            lines.append(f"→ **In-silico consensus:** {summary} ({names})")
    else:
        lines.append("**Missense predictors (AlphaMissense/REVEL/SaProt) do not apply** to this "
                     "variant type — they score amino-acid substitutions only. For non-missense "
                     "variants the computational signal comes from conservation"
                     + (" and, for splice-region changes, SpliceAI" if spliceai else "")
                     + " (below), not from a missing predictor score.")

    if spliceai:
        # SpliceAI is a descriptive AGREEMENT signal for splice-region variants, NOT an additive
        # vote. For a non-missense variant the single calibrated tool is conservation (below), so
        # SpliceAI does not increment the concordance counter. (It is only surfaced at delta >= 0.2,
        # i.e. it always predicts a splice-altering effect; counting it both double-counted with
        # conservation on splice variants and let it only ever confirm — never dissent from — a
        # pathogenic call. It is now described and its agreement noted, without being counted.)
        effect = (spliceai.get('effect') or '').replace('_', ' ')
        if clinvar_path:
            note = "consistent with the pathogenic call"
        elif clinvar_vus:
            note = "a splice-altering prediction on an uncertain variant — a signal for review"
        else:
            note = "a splice-altering prediction the current classification does not reflect"
        lines.append(f"SpliceAI predicts **{effect}** (Δ {spliceai.get('score')}) — {note}")

    if conservation and conservation.get("phylop") is not None:
        p, gerp, pc = conservation["phylop"], conservation.get("gerp"), conservation.get("phastcons")
        conserved = p >= 2
        detail = ("phyloP %g" % p) + (", GERP %g" % gerp if gerp is not None else "") \
            + (", phastCons %g" % pc if pc is not None else "")
        if not am:
            # non-missense/splice/intronic: conservation IS the primary computational
            # signal (AlphaMissense/REVEL are missense-only) → it carries the weight.
            total += 1
            if conserved and clinvar_path:
                agree += 1
            note = ("highly conserved position (PP3-type support)" if conserved
                    else "fast-evolving / non-conserved position" if p <= -2
                    else "moderately conserved position")
            lines.append(f"Conservation: {detail} — {note}")
        else:
            # missense: conservation is an agreement check on AlphaMissense (not counted)
            lines.append(f"Conservation: {detail} — {'conserved (agrees)' if conserved else 'low conservation'}")

    # Population frequency — shown, but NEVER counted as concordant evidence.
    if gnomad:
        if gnomad.get("is_common"):
            flags.append(gnomad.get("band", "common in gnomAD"))
            lines.append(f"**{gnomad['band']}** ⚠")
        elif gnomad.get("absent"):
            lines.append("Absent from gnomAD v4.1 (very rare — ACMG PM2-supporting only)")
        else:
            faf = gnomad.get("faf")
            faf_note = (f"; filtering AF (faf95) {faf}" if faf
                        else "; faf/AC-AN not in this data projection")
            lines.append(gnomad["band"][0].upper() + gnomad["band"][1:]
                         + f" (gnomAD v4.1{faf_note} — BA1/BS1 are disease-specific thresholds, "
                         "not a blanket 5%)")

    if flags:
        verdict = "Evidence sources **disagree** — " + "; ".join(flags)
    elif total == 0:
        verdict = None
    elif agree == total and clinvar_path:
        verdict = (f"{agree} independent predictor{'s' if agree != 1 else ''} "
                   "**concordant** with the ClinVar classification")
    elif clinvar_vus and consensus and consensus.get("unanimous"):
        # VUS triage (§8): ClinVar has no definitive call, so there is nothing to be
        # concordant WITH — instead surface where the independent predictors unanimously
        # lean. Explicitly NOT a reclassification; a signal for prioritising review.
        lean = "damaging" if consensus["n_damaging"] else "tolerated"
        verdict = (f"ClinVar classifies this **uncertain**; the computational predictors "
                   f"unanimously lean **{lean}** — a triage signal, not a reclassification")
    else:
        verdict = "Mixed / partial computational evidence (see below)"
    return {"lines": lines, "verdict": verdict, "flags": flags, "consensus": consensus}


# Category order for the /disagreements browse view (higher = surfaced first).
DISAGREEMENT_CATEGORIES = ("predictor_vs_clinvar", "resolves_conflicting",
                           "lof_resolves_conflicting", "vus_predictors_lean",
                           "predictors_split")


def disagreement_flag(rec):
    """Classify a record's evidence disagreement for the /disagreements browse view,
    from the ALREADY-computed concordance (no new biobtree lookups). Returns
    {category, label, reason, severity} or None.

    The three disagreements no competitor surfaces (benchmark 2026-07):
      predictor_vs_clinvar     — a pathogenic/LP ClinVar call the predictors lean against
                                 (AlphaMissense — the ACMG-weighted tool — calls benign, or
                                 the in-silico consensus is majority-tolerant)
      resolves_conflicting     — ClinVar 'conflicting', but the missense predictors agree
      lof_resolves_conflicting — ClinVar 'conflicting', but it's a predicted loss-of-function
                                 change in a LoF-intolerant gene (the non-missense analog)
      predictors_split         — the independent predictors disagree with each other

    Descriptive only — a QC signal, NOT a reclassification (ClinGen SVI, HANDOVER §8)."""
    c = rec.get("concordance") or {}
    cons = c.get("consensus") or {}
    total = cons.get("total") or 0
    n_dmg = cons.get("n_damaging") or 0
    cls = (rec.get("classification") or "").lower()
    is_path = "pathogenic" in cls and "conflict" not in cls
    is_conf = "conflict" in cls
    is_vus = "uncertain" in cls
    am_benign = any("likely-benign" in f for f in (c.get("flags") or []))

    if is_path and (am_benign or (total and n_dmg * 2 < total)):
        detail = cons.get("summary") if total else "AlphaMissense predicts likely-benign"
        return {"category": "predictor_vs_clinvar", "label": "Predictors vs ClinVar",
                "reason": f"ClinVar {rec.get('classification')}, but {detail}", "severity": 3}
    if is_conf and total >= 2 and cons.get("unanimous"):
        return {"category": "resolves_conflicting", "label": "Predictors agree, ClinVar conflicts",
                "reason": f"ClinVar conflicting; {cons.get('summary')} — a QC flag, not a reclassification",
                "severity": 2}
    # non-missense: a predicted-LoF change in a LoF-intolerant gene is a mechanism-based
    # resolving signal on a conflicting call (brings the non-missense half into the QC view).
    lof = rec.get("lof_context")
    if is_conf and lof and lof.get("lof_disease_gene"):
        why = ("a gene with sufficient ClinGen haploinsufficiency evidence"
               if lof.get("haploinsufficient") else "a LoF-intolerant gene")
        nmd = " (verify NMD-escape for C-terminal truncations)" if lof.get("nmd_caveat") else ""
        return {"category": "lof_resolves_conflicting", "label": "Predicted LoF vs a conflicting call",
                "reason": f"ClinVar conflicting; predicted {lof['label']} in {why}{nmd} — flagged for review",
                "severity": 2}
    # VUS with a unanimous predictor lean — the highest-value new triage class (Phase 2).
    # ClinVar has no definitive call; the independent predictors all point one way. A
    # prioritisation signal for review, explicitly NOT a reclassification (ClinGen SVI, §8).
    if is_vus and total >= 2 and cons.get("unanimous"):
        lean = "damaging" if n_dmg else "tolerated"
        return {"category": "vus_predictors_lean", "label": "Predictors lean, ClinVar uncertain",
                "reason": f"ClinVar uncertain significance; {cons.get('summary')} ({lean}) "
                          "— a triage signal, not a reclassification", "severity": 2}
    if total >= 2 and not cons.get("unanimous"):
        return {"category": "predictors_split", "label": "Predictors split",
                "reason": cons.get("summary"), "severity": 1}
    return None


# ── Batch 3: per-gene context (fetched once per gene, cached in ctx) ─────────
def gene_context(hgnc_id, gene_symbol=None):
    """ACMG-adjacent gene block: gnomAD constraint + ClinGen dosage + gene-disease
    validity + inheritance (GenCC/validity). Fetched once per gene."""
    con = map_all(hgnc_id, ">>hgnc>>gnomad_constraint")
    constraint = None
    if con:
        c = con[0]
        # Guard: the >>gnomad_constraint edge can return a DIFFERENT gene's row — e.g. RMRP
        # (an ncRNA with no gnomAD constraint) yields NME1's LOEUF/pLI (benchmark 2026-07).
        # gnomAD constraint is protein-coding only; only trust an exact gene-symbol match.
        row_gene = (c.get("gene_symbol") or "").strip().upper()
        if not gene_symbol or (row_gene and row_gene == gene_symbol.strip().upper()):
            constraint = {"pli": c.get("pli"), "loeuf": c.get("loeuf"), "mis_z": c.get("mis_z")}
    dos = map_all(hgnc_id, ">>hgnc>>clingen_dosage")
    dosage = ({"haplo": dos[0].get("haplo_score"), "triplo": dos[0].get("triplo_score")}
              if dos else None)
    validity = [{"disease": v.get("disease_label"), "moi": v.get("moi"),
                 "classification": v.get("classification")}
                for v in map_all(hgnc_id, ">>hgnc>>clingen_gene_validity") if v.get("disease_label")]
    inh = []
    for g in map_all(hgnc_id, ">>hgnc>>gencc"):
        t = (g.get("moi_title") or "").strip()
        if t and t not in inh:
            inh.append(t)
    return {"constraint": constraint, "dosage": dosage, "validity": validity[:4],
            "inheritance": inh}


# Meaningful UniProt feature types → a human phrase ('{d}' = the description).
_UF_KEEP = {
    "domain": "the {d} domain", "region of interest": "the '{d}' region",
    "binding site": "a binding site", "active site": "an active site",
    "modified residue": "a modified residue ({d})", "helix": "an α-helix",
    "strand": "a β-strand", "turn": "a turn", "motif": "the {d} motif",
    "zinc finger": "a zinc-finger region", "metal binding": "a metal-binding site",
    "dna-binding region": "a DNA-binding region", "site": "a functional site ({d})",
    "nucleotide binding region": "a nucleotide-binding region",
    "cross-link": "a cross-link site", "disulfide bond": "a disulfide bond",
    "sequence variant": "a UniProt-annotated variant site ({d})",
}


def gene_structure(hgnc_id):
    """PDB structures + AlphaFold confidence + parsed UniProt feature intervals
    (for per-residue structural context). Fetched once per gene."""
    pdb = [{"id": p.get("id"), "title": p.get("title"), "method": p.get("method"),
            "resolution": p.get("resolution")}
           for p in map_all(hgnc_id, ">>hgnc>>uniprot>>pdb")]
    af = map_all(hgnc_id, ">>hgnc>>uniprot>>alphafold")
    alphafold = ({"plddt": af[0].get("global_metric"),
                  "frac_high": af[0].get("fraction_plddt_very_high")} if af else None)
    intervals = []
    for f in map_all(hgnc_id, ">>hgnc>>uniprot>>ufeature", cap=60):
        t = (f.get("type") or "").strip().lower()
        if t not in _UF_KEEP:
            continue
        try:
            b, e = int(f.get("location_begin")), int(f.get("location_end"))
        except (TypeError, ValueError):
            continue
        intervals.append({"type": t, "desc": (f.get("description") or "").strip(),
                          "begin": b, "end": e})
    return {"pdb": pdb, "alphafold": alphafold, "intervals": intervals}


def gene_mavedb(hgnc_id):
    """{hgvs_pro: [{score, score_set, title, license}]} of the gene's MaveDB
    multiplexed functional-assay measurements, fetched once per gene. Empty for
    most genes (only MAVE-assayed genes carry scores). score_set_title is looked
    up once per score set (it's not in the lite map projection)."""
    from collections import defaultdict
    from sugibiobtree import entry
    by_hgvs, rep = defaultdict(list), {}
    for x in map_all(hgnc_id, ">>hgnc>>mavedb", cap=200):
        hp, ss = x.get("hgvs_pro"), x.get("score_set")
        if hp and ss and x.get("score") is not None:
            by_hgvs[hp].append({"score": x["score"], "score_set": ss,
                                "license": x.get("license")})
            rep.setdefault(ss, x.get("id"))
    titles = {}
    for ss, rid in rep.items():
        if rid:
            a = (entry(rid, "mavedb") or {}).get("Attributes") or {}
            m = a.get("Mavedb") or (next(iter(a.values()), {}) if a else {})
            titles[ss] = (m or {}).get("score_set_title")
    for lst in by_hgvs.values():
        for r in lst:
            r["title"] = titles.get(r["score_set"])
    return dict(by_hgvs)


def mavedb_for(hgvs_p, cache):
    """Deduped MaveDB assay measurements for a variant's protein change (one row
    per score set), or None."""
    rows_ = (cache or {}).get(hgvs_p)
    if not rows_:
        return None
    seen, out = set(), []
    for r in rows_:
        if r["score_set"] not in seen:
            seen.add(r["score_set"])
            out.append(r)
    return out[:6]


def gene_spliceai(hgnc_id):
    """{coordinate: {effect, score}} of the gene's SpliceAI splice-impact
    predictions (chr:pos:ref:alt keys), fetched once per gene."""
    out = {}
    # uncapped: every variant's coordinate must resolve its SpliceAI score — a
    # 60-page cap would drop splice annotations on big genes. Once per gene (cached).
    for r in map_all(hgnc_id, ">>hgnc>>spliceai", cap=None):
        cid = r.get("id")
        if cid and r.get("score"):
            out[cid] = {"effect": r.get("effect"), "score": r.get("score")}
    return out


def spliceai_for(coord, cache):
    """SpliceAI prediction for a variant's coordinate, if it has a meaningful
    (>=0.2) delta score — SpliceAI only annotates splice-relevant positions."""
    if not coord or not cache:
        return None
    hit = cache.get(coord)
    if not hit:
        return None
    try:
        if float(hit["score"]) < 0.2:
            return None
    except (TypeError, ValueError):
        return None
    return hit


# GO experimental-evidence ECO set (mirrors gene §7 s07_pathways).
_GO_EXPERIMENTAL = {"ECO:0000314", "ECO:0000315", "ECO:0000316", "ECO:0000353",
                    "ECO:0000270", "ECO:0000269", "ECO:0006056", "ECO:0007005",
                    "ECO:0007007", "ECO:0007003"}
_REACTOME_EV = {"TAS": "curated", "IEA": "electronic"}


def gene_pathways(hgnc_id):
    """Reactome pathways (disease-flagged, evidence-merged) + experimental GO
    terms for the gene, fetched once. GENE-level context — 'the gene this variant
    disrupts acts in these pathways', never a per-variant claim."""
    rx = {}
    for t in map_all(hgnc_id, ">>hgnc>>uniprot>>reactome"):
        pid = t.get("id")
        if not pid:
            continue
        ev = _REACTOME_EV.get((t.get("evidence") or "").strip(), "other")
        cur = rx.get(pid)
        best = "curated" if ev == "curated" or (cur and cur["evidence"] == "curated") else ev
        rx[pid] = {"id": pid, "name": t.get("name") or (cur or {}).get("name"),
                   "evidence": best,
                   "is_disease": (t.get("is_disease_pathway") == "true") or bool(cur and cur["is_disease"])}
    # disease pathways first, then curated, then by name
    pathways = sorted(rx.values(), key=lambda p: (not p["is_disease"],
                                                  p["evidence"] != "curated", p["name"] or ""))
    go_rows = list(map_all(hgnc_id, ">>hgnc>>uniprot>>go"))
    has_eco = any(r.get("evidence", "").startswith("ECO:") for r in go_rows)

    def _pick(ns):
        # experimental terms when ECO present; else top terms w/o the exp. claim
        rows = [r for r in go_rows if r.get("type") == ns
                and (not has_eco or r.get("evidence") in _GO_EXPERIMENTAL)]
        return [{"id": r.get("id"), "name": r.get("name")} for r in rows if r.get("name")]

    return {"pathways": pathways,
            "disease_pathways": [p for p in pathways if p["is_disease"]],
            "go_mf": _pick("molecular_function"),
            "go_bp": _pick("biological_process"),
            "go_experimental": has_eco,
            "top_function": (_pick("molecular_function")[:1] or [{}])[0].get("name")}


def mechanism_narrative(rec, pathways):
    """Deterministic variant→gene→pathway→disease sentence (NOT an LLM). Only
    fires with a real anchor (a Reactome disease pathway OR an experimental GO
    molecular-function term). Gene-level, framed 'thought to act' — never 'causes'."""
    disease_pw = (pathways or {}).get("disease_pathways") or []
    top_fn = (pathways or {}).get("top_function")
    if not (disease_pw or top_fn):
        return None
    from sugivariant.render import short_hgvs   # local: render has no top-level enrich dep
    gene = rec.get("gene_symbol")
    pchange = short_hgvs(rec.get("hgvs_p") or rec.get("hgvs_c"))
    st = rec.get("structural") or {}
    dom = next((f for f in (st.get("features") or []) if "domain" in f or "region" in f), None)
    if dom:
        dom = dom.replace("'", "")           # drop stray quotes from the UniProt feature label
    # Gene function first (this section is gene-level), then the variant's position as a
    # plain locational fact — never an asserted effect.
    parts = []
    if top_fn:
        parts.append(f"{gene}'s established molecular role includes {top_fn}.")
    if dom:
        parts.append(f"{' ' if parts else ''}{pchange} falls in {dom}.")
    elif not top_fn:
        parts.append(f"{gene} participates in disease-associated pathways.")
    # (per-variant AlphaMissense/SpliceAI predictions live in the Computational card — this
    #  gene-function narrative stays gene-level; benchmark 2026-07 flagged a SpliceAI overclaim.)
    if disease_pw:
        # Report the Reactome annotation as a FACT — do NOT assert an arbitrary disease
        # pathway is "the mechanism" of this variant's condition (benchmark 2026-07: that
        # glued the first/alphabetical disease pathway to conditions[0], fabricating causal
        # links — e.g. a somatic-cancer pathway as the mechanism of germline Loeys-Dietz).
        names = ", ".join(f"**{p['name']}**" for p in disease_pw[:2] if p.get("name"))
        if names:
            parts.append(f" In Reactome, {gene} is annotated in disease pathway(s) including {names}.")
    return "".join(parts)


def gene_pharmgkb(hgnc_id):
    """{rsID: {annotations, clinical}} pharmacogenomics for the gene, fetched once.
    Empty for the vast majority of genes (only pharmacogenes carry PGx)."""
    out = {}
    for a in map_all(hgnc_id, ">>hgnc>>pharmgkb_var_annotation", cap=20):
        rs = a.get("variant")
        if rs and a.get("drugs"):
            out.setdefault(rs, {"annotations": [], "clinical": []})["annotations"].append(
                {"drugs": a.get("drugs"), "category": a.get("phenotype_category"),
                 "significance": a.get("significance"), "sentence": a.get("sentence"),
                 "pmid": a.get("pmid")})
    for c in map_all(hgnc_id, ">>hgnc>>pharmgkb_clinical", cap=20):
        rs = c.get("variant")
        if rs and c.get("chemicals"):
            out.setdefault(rs, {"annotations": [], "clinical": []})["clinical"].append(
                {"chemicals": c.get("chemicals"), "level": c.get("level_of_evidence"),
                 "type": c.get("type"), "phenotypes": c.get("phenotypes")})
    return out


def pharmgkb_for(rsid, cache):
    """PGx for a variant's rsID from the per-gene cache, deduped."""
    if not rsid or rsid not in (cache or {}):
        return None
    e = cache[rsid]
    drugs = sorted({a["drugs"] for a in e["annotations"] if a.get("drugs")})
    return {"drugs": drugs[:8], "n": len(e["annotations"]),
            "clinical": sorted(e["clinical"], key=lambda c: (c.get("level") or "9"))[:5]}


def gene_has_civic(hgnc_id):
    """Whether the gene has any CIViC curation — gates the per-variant CIViC
    lookup so non-oncology genes cost one probe, not one call per variant."""
    return bool(map_all(hgnc_id, ">>hgnc>>civic", cap=1))


def civic_for(variation_id, has_civic):
    """CIViC predictive/prognostic evidence for a variant (cancer genes). Joined
    via the clinvar↔civic_variant xref. None unless the gene has CIViC data."""
    if not has_civic:
        return None
    cv = map_all(variation_id, ">>clinvar>>civic_variant")
    if not cv:
        return None
    ev = map_all(cv[0]["id"], ">>civic_variant>>civic_evidence")
    if not ev:
        return None
    return {"name": cv[0].get("name"),
            "evidence": [{"disease": e.get("disease"), "therapies": e.get("therapies"),
                          "type": e.get("evidence_type"), "level": e.get("evidence_level"),
                          "significance": e.get("significance")} for e in ev[:8]]}


def gene_variant_landscape(recs):
    """Aggregate profile of the gene's built variants — for the per-gene index.
    Zero new calls: everything from the in-memory record list."""
    from collections import Counter
    by_cls = Counter(r.get("classification") for r in recs)
    by_type = Counter(r.get("variant_type") for r in recs if r.get("variant_type"))
    pos = Counter()
    for r in recs:
        p = protein_position(r.get("hgvs_p"))
        if p is not None and "pathogenic" in (r.get("classification") or "").lower():
            pos[p] += 1
    recurrent = [(p, n) for p, n in pos.most_common(8) if n > 1]
    span = None
    ps = [p for p in pos]
    if ps:
        span = (min(ps), max(ps))
    return {"by_class": dict(by_cls), "by_type": dict(by_type.most_common(6)),
            "recurrent": recurrent, "n": len(recs), "residue_span": span}


def condition_digest(conditions, cache):
    """Patient digest for the VARIANT's OWN condition (audit P1: not a gene-level
    constant). Routes each of the variant's ClinVar-linked conditions through
    Orphanet and returns the first that yields a germline Disease entry (climbing
    one MONDO parent level for thinly-annotated subtypes). Crucially, SOMATIC/
    acquired conditions (leukemia, mastocytosis…) have no Orphanet germline Disease
    → they yield no digest → no germline inheritance/onset/HPO framing is projected
    onto them. Cached per MONDO id."""
    # Route ONLY the PRIMARY condition (conditions[0], anchored to ClinVar's germline trait)
    # — never a sibling. Falling through to the next condition projected a DIFFERENT disease's
    # inheritance onto the headline (benchmark 2026-07: ACTA1 recessive alpha-actinopathy showed
    # a sibling's 'Autosomal dominant'). Better to show nothing than the wrong mode.
    c = (conditions or [None])[0]
    mid = c.get("mondo_id") if c else None
    if not mid:
        return None
    if mid not in cache:
        cache[mid] = _digest_via_mondo(mid) or _digest_via_parent(mid)
    return cache[mid]


def _digest_via_mondo(mondo_id):
    orphas = [o for o in map_all(mondo_id, ">>mondo>>orphanet")
              if (o.get("disorder_type") or "") == "Disease"]
    if not orphas:
        return None
    best = max(orphas, key=lambda o: int(o.get("phenotype_count") or 0))
    return _build_orphanet_digest(best.get("id"), best.get("name"))


def _digest_via_parent(mondo_id):
    for par in map_all(mondo_id, ">>mondo>>mondoparent")[:2]:
        d = _digest_via_mondo(par.get("id"))
        if d:
            return d
    return None


_ORPHA_PREFIX = re.compile(r"^[A-Z][A-Z0-9 ,/'()\-]{2,}:\s*")


def _clean_disorder_name(name):
    """Strip Orphanet classification-of-rarity prefixes (e.g. 'NON RARE IN EUROPE: ')
    that leak into the disorder name. Only removes a leading ALL-CAPS 'PREFIX: ' run,
    so real names like 'Hemochromatosis type 1' are untouched."""
    if not name:
        return name
    return _ORPHA_PREFIX.sub("", name).strip() or name


def _build_orphanet_digest(oid, fallback_name=None):
    o = _orphanet_entry(oid)
    if not o:
        return None
    phen = sorted((o.get("phenotypes") or []),
                  key=lambda p: -(p.get("frequency_value") or 0))
    prev = (o.get("prevalences") or [{}])[0]
    pc = prev.get("prevalence_class")
    return {
        "name": _clean_disorder_name(o.get("name")) or fallback_name,
        "inheritance": o.get("inheritance") or [],
        "onset": o.get("onset") or [],
        "prevalence": (f"{pc} ({prev.get('geographic')})"
                       if pc and pc.lower() != "unknown" else None),
        "phenotypes": [{"term": p.get("hpo_term"), "freq": p.get("frequency")}
                       for p in phen[:8] if p.get("hpo_term")],
    }


def _orphanet_entry(oid):
    from sugibiobtree import entry
    a = (entry(oid, "orphanet") or {}).get("Attributes") or {}
    return a.get("Orphanet") or (next(iter(a.values()), {}) if a else {})


def gene_panelapp(hgnc_id):
    """Green (diagnostic-grade) Genomics England panels the gene is on."""
    return [p.get("panel_name") for p in map_all(hgnc_id, ">>hgnc>>panelapp_gene")
            if (p.get("confidence") or "").lower() == "green" and p.get("panel_name")]


def condition_links(mondo_id, cache):
    """GARD registry + clinical-trial count for the variant's EXACT condition.
    Cached per MONDO id. Audit P2c: we no longer climb to a MONDO parent for
    trials — climbing to a broad umbrella term (e.g. autism spectrum disorder)
    inflated the count into the thousands and mis-attributed unrelated trials."""
    if not mondo_id:
        return {}
    if mondo_id in cache:
        return cache[mondo_id]
    gard = [g.get("id") for g in map_all(mondo_id, ">>mondo>>gard") if g.get("id")]
    trials = map_all(mondo_id, ">>mondo>>clinical_trials")
    out = {"gard": gard[0] if gard else None, "trial_count": len(trials)}
    cache[mondo_id] = out
    return out


# ── Batch 3: per-variant derived (mostly in-memory, no new calls) ────────────
def structural_context(hgvs_p, intervals):
    """UniProt features overlapping the variant's residue → human phrases."""
    pos = protein_position(hgvs_p)
    if pos is None or not intervals:
        return None
    out, seen = [], set()
    for f in intervals:
        if f["begin"] <= pos <= f["end"]:
            # Strip a foreign dbSNP id from the feature description — it's the rsID of the
            # UniProt-catalogued variant at this residue, not OUR variant (benchmark 2026-07:
            # ACADVL R459Q surfaced R459W's rs766742117).
            desc = re.sub(r"[;,]?\s*dbsnp:rs\d+\.?", "", f["desc"] or "", flags=re.I).strip(" ;.")
            phrase = _UF_KEEP[f["type"]].replace("{d}", desc or f["type"])
            if phrase not in seen:                      # dedup UniProt's duplicate curations
                seen.add(phrase)
                out.append(phrase)
    return {"position": pos, "features": out} if out else None


def am_percentile(score, am_map):
    """The variant's AlphaMissense percentile within the gene's full modeled set
    ('top 3% most-pathogenic-predicted of N'), or None."""
    try:
        s = float(score)
    except (TypeError, ValueError):
        return None
    vals = []
    for _, sc in am_map.values():
        try:
            vals.append(float(sc))
        except (TypeError, ValueError):
            pass
    if len(vals) < 20:
        return None
    n_below = sum(1 for v in vals if v < s)
    pct = round(100 * (len(vals) - n_below) / len(vals))
    return {"top_pct": max(pct, 1), "n": len(vals)}


def similar_variants(rec, recs):
    """Up to 6 sibling variant pages most related to this one — same residue,
    then same condition, then same variant type. Internal-mesh aid for both
    personas + crawlers."""
    pos = protein_position(rec.get("hgvs_p"))
    conds = {c.get("mondo_id") for c in (rec.get("conditions") or [])}
    vtype = rec.get("variant_type")
    me = rec["canonical_slug"]
    scored = []
    for r in recs:
        if r["canonical_slug"] == me:
            continue
        s = 0
        if pos is not None and protein_position(r.get("hgvs_p")) == pos:
            s += 3
        if conds & {c.get("mondo_id") for c in (r.get("conditions") or [])}:
            s += 2
        if vtype and r.get("variant_type") == vtype:
            s += 1
        if s:
            scored.append((s, r))
    # Rank within relation-score by what's most worth a click: more review stars, then
    # the more clinically interesting call (P/LP or VUS/conflicting over benign), then
    # slug for determinism. Keeps the alphabetical-tie problem from burying the useful
    # same-condition neighbors under intronic/splice siblings.
    def _interest(cls):
        cl = (cls or "").lower()
        if "pathogenic" in cl and "conflict" not in cl:
            return 3
        if "uncertain" in cl or "conflict" in cl:
            return 2
        if "benign" in cl:
            return 1
        return 0
    scored.sort(key=lambda x: (
        -x[0],
        -_STAR_N.get((x[1].get("review_status") or "").strip().lower(), 0),
        -_interest(x[1].get("classification")),
        x[1]["canonical_slug"]))
    out = []
    for _, r in scored[:8]:
        # the dominant relation (why it's similar) — for a scannable "why click" cue
        if pos is not None and protein_position(r.get("hgvs_p")) == pos:
            rel = "same residue"
        elif conds & {c.get("mondo_id") for c in (r.get("conditions") or [])}:
            rel = "same condition"
        else:
            rel = "same type"
        out.append({"label": f"{r['gene_symbol']} {r.get('hgvs_p') or r.get('hgvs_c')}",
                    "slug": r["canonical_slug"],
                    "classification": r.get("classification"),
                    "stars": _STAR_N.get((r.get("review_status") or "").strip().lower(), 0),
                    "relation": rel})
    return out


def submission_timeline(submissions):
    """First/last submission year + whether classifications diverged over time."""
    dated = [s for s in submissions if s.get("date")]
    years = sorted(int(s["date"][:4]) for s in dated if (s.get("date") or "")[:4].isdigit())
    if not years:
        return None
    calls = {(s.get("classification") or "").split("/")[0].strip().lower() for s in submissions}
    return {"first": years[0], "last": years[-1], "n": len(submissions),
            "stable": len(calls) <= 1}


_STAR_N = {"practice guideline": 4, "reviewed by expert panel": 3,
           "criteria provided, multiple submitters, no conflicts": 2,
           "criteria provided, single submitter": 1,
           "criteria provided, conflicting classifications": 1,
           "no assertion criteria provided": 0, "no classification provided": 0}


def review_stars(review_status):
    return _STAR_N.get((review_status or "").strip().lower(), 0)


# Named review-status tier (replaces the ★ rating metaphor — this is who reviewed
# the variant and how, not a quality score). label = short human name; cls = tier
# key for colour (rv4 = most authoritative … rv0 = unreviewed).
_REVIEW_LABEL = {
    "practice guideline": "Practice guideline",
    "reviewed by expert panel": "Expert panel",
    "criteria provided, multiple submitters, no conflicts": "Multiple submitters",
    "criteria provided, single submitter": "Single submitter",
    "criteria provided, conflicting classifications": "Conflicting submitters",
    "no assertion criteria provided": "No assertion criteria",
    "no classification provided": "Not classified",
}
_TIER_FALLBACK = {4: "Practice guideline", 3: "Expert panel", 2: "Multiple submitters",
                  1: "Single submitter", 0: "No assertion criteria"}


def review_tier(review_status):
    n = review_stars(review_status)
    label = _REVIEW_LABEL.get((review_status or "").strip().lower()) or _TIER_FALLBACK.get(n, "Unreviewed")
    return {"n": n, "label": label, "cls": f"rv{n}"}


def tier_label(n):
    """Short review-tier label from the numeric tier alone (for index-backed rows
    that carry only the tier, not the full review-status string)."""
    return _TIER_FALLBACK.get(n, "Unreviewed")


def plain_summary(rec):
    """Deterministic plain-language one-liner for patients (NOT an LLM). The
    confidence phrasing is CALIBRATED to classification strength + review stars +
    in-silico concordance (audit P2: don't flatten Likely-pathogenic, 1-star, and
    computationally-discordant calls into a bare 'disease-causing'). The condition
    is the VARIANT's own (audit P1: reconcile with the Summary)."""
    cls = (rec.get("classification") or "").lower()
    stars = review_stars(rec.get("review_status"))
    discordant = bool((rec.get("concordance") or {}).get("flags"))
    if "conflict" in cls:
        meaning = "of uncertain or conflicting significance (submitters disagree on it)"
    elif "likely pathogenic" in cls:
        meaning = "considered likely disease-causing"
    elif "pathogenic" in cls:
        if discordant:
            # 'flags' holds the single calibrated/weighted signal that dissents (e.g. AlphaMissense),
            # not the whole predictor set — REVEL/SaProt may still agree. Don't overstate as
            # 'predictors disagree' (benchmark 2026-07: CFTR R117H is 2/3 damaging).
            meaning = "reported as disease-causing, though the computational evidence is not fully concordant"
        elif stars >= 2:
            meaning = "considered disease-causing"
        else:
            meaning = "reported as disease-causing, but on limited review"
    elif "uncertain" in cls:
        # §8: never reclassify a VUS. State the uncertainty; the computational lean lives
        # in the concordance card, not here in the headline meaning.
        meaning = ("currently of uncertain significance — a definitive interpretation "
                   "has not yet been established")
    elif "benign" in cls:
        meaning = f"classified as {(rec.get('classification') or '').lower()} — not considered disease-causing"
    else:
        meaning = f"classified as {rec.get('classification')}"
    gene = rec.get("gene_symbol")
    n = rec.get("submitter_count") or 0
    who = (f", submitted by {n} submitter" + ("s" if n != 1 else "")) if n else ""
    cond = (rec.get("conditions") or [{}])[0].get("name")
    link = f", and is linked to {cond}" if cond else ""
    return f"This is a change in the {gene} gene that is {meaning}{who}{link}."


def residue_hotspot(hgvs_p, position_index):
    """Other pathogenic ClinVar variants at the same residue, from a prebuilt
    {position: [labels]} index over the gene's P/LP set. Returns a dict or None."""
    pos = protein_position(hgvs_p)
    if pos is None or not position_index:
        return None
    here = [x for x in position_index.get(pos, []) if x.get("hgvs_p") != hgvs_p]
    if not here:
        return None
    return {"position": pos, "others": here}
