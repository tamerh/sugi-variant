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
    # cap high enough for full per-gene coverage (~2.5k substitutions/small gene)
    # so the percentile stats see the whole distribution.
    for r in map_all(hgnc_id, ">>hgnc>>uniprot>>alphamissense", cap=60):
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
            af, popmax = _f(g.get("af")), _f(g.get("af_grpmax"))
            anc = g.get("grpmax_ancestry")
            pops = {k[3:]: g[k] for k in g if k.startswith("af_") and k != "af_grpmax" and g.get(k)}
            return {"af": g.get("af"), "popmax": popmax, "ancestry": anc,
                    "faf": g.get("faf"), "populations": pops,
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


def esm1b_for(uniprot, protein_variant):
    """ESM1b protein-language-model variant effect (Brandes 2023), keyed by
    uniprot:protein_variant (both from AlphaMissense). LLR ≤ 0, more negative =
    more damaging; NO fixed ACMG threshold — surface the raw LLR and let the
    concordance panel do the work (~-7.5 is the paper's rough divider, used only
    for the 'damaging vs tolerated' agreement read). Training is ClinVar-
    independent → a genuinely orthogonal predictor next to AlphaMissense/REVEL."""
    if not (uniprot and protein_variant):
        return None
    a = _coord_entry(f"{uniprot}:{protein_variant}", "esm1b")
    llr = _f((a or {}).get("esm1b_llr"))
    if llr is None:
        return None
    return {"llr": llr, "damaging": llr <= -7.5}


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


# REVEL → ACMG evidence tier (ClinGen SVI / Pejaver 2022 calibrated thresholds).
def _revel_tier(s):
    if s >= 0.932:
        return "PP3_Strong"
    if s >= 0.773:
        return "PP3_Moderate"
    if s >= 0.644:
        return "PP3_Supporting"
    if s > 0.290:
        return "indeterminate"
    if s >= 0.183:
        return "BP4_Supporting"
    if s > 0.016:
        return "BP4_Moderate"
    return "BP4_Strong"


def revel_for(coord):
    """REVEL ensemble missense pathogenicity (0–1) by chr:pos:ref:alt, with its
    ClinGen-calibrated ACMG tier + direction. Presented as an AGREEMENT signal
    (not additive PP3 — see concordance)."""
    a = _coord_entry(coord, "revel")
    s = _f((a or {}).get("revel"))
    if s is None:
        return None
    direction = "pathogenic" if s >= 0.644 else "benign" if s <= 0.290 else "indeterminate"
    return {"score": s, "tier": _revel_tier(s), "direction": direction}


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
    return {"frequency": freq, "absent": absent,
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


def concordance(classification, am, gnomad, spliceai=None, conservation=None, revel=None, esm1b=None):
    """Cross-source concordance readout: do the computational predictors agree with
    ClinVar? Returns {lines, verdict, flags}. ClinGen-SVI framing: in-silico
    predictors are NOT independent, so exactly ONE calibrated tool carries the
    ACMG weight — AlphaMissense for missense, conservation for non-missense — and
    REVEL (+ conservation, on missense) is shown as an AGREEMENT signal, never an
    additive vote. Population rarity is shown but never counted (audit P2). Never a
    clinical determination."""
    cls = (classification or "").lower()
    clinvar_path = "pathogenic" in cls and "conflict" not in cls
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
            lines.append(f"REVEL {revel['score']} ({revel['tier']}) — {rel} AlphaMissense")
        # ESM1b — orthogonal (ClinVar-independent) protein-language-model opinion,
        # also an agreement signal (not additive). Raw LLR surfaced; ~-7.5 divides
        # the damaging/tolerated *read* only.
        if esm1b:
            rel = "agrees with" if esm1b["damaging"] == am_path else "differs from"
            lines.append(f"ESM1b LLR {esm1b['llr']:g} "
                         f"({'damaging' if esm1b['damaging'] else 'tolerated'}) — {rel} AlphaMissense")
        # In-silico CONSENSUS — descriptive agreement of the independent missense
        # predictors (do they concur?). NOT additive ACMG evidence (ClinGen SVI).
        preds = [("AlphaMissense", am_path)]
        if revel and revel["direction"] != "indeterminate":
            preds.append(("REVEL", revel["direction"] == "pathogenic"))
        if esm1b:
            preds.append(("ESM1b", esm1b["damaging"]))
        if len(preds) >= 2:
            n_dmg = sum(1 for _, d in preds if d)
            names = ", ".join(n for n, _ in preds)
            if n_dmg == len(preds):
                summary = f"all {len(preds)} predictors call this damaging"
            elif n_dmg == 0:
                summary = f"all {len(preds)} predictors call this tolerated"
            else:
                summary = f"{n_dmg}/{len(preds)} damaging — mixed"
            consensus = {"n_damaging": n_dmg, "total": len(preds),
                         "names": names, "summary": summary,
                         "unanimous": n_dmg in (0, len(preds))}
            lines.append(f"→ **In-silico consensus:** {summary} ({names})")
    else:
        lines.append("AlphaMissense — not scored (not a missense SNV)")

    if spliceai:
        total += 1
        if clinvar_path:
            agree += 1
        lines.append(f"SpliceAI predicts **{(spliceai.get('effect') or '').replace('_', ' ')}** "
                     f"(Δ {spliceai.get('score')}) — a splice-altering effect consistent with pathogenicity")

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
            lines.append(gnomad["band"][0].upper() + gnomad["band"][1:]
                         + " (gnomAD v4.1 — PM2/BS1 context)")

    if flags:
        verdict = "Evidence sources **disagree** — " + "; ".join(flags)
    elif total == 0:
        verdict = None
    elif agree == total and clinvar_path:
        verdict = (f"{agree} independent predictor{'s' if agree != 1 else ''} "
                   "**concordant** with the ClinVar classification")
    else:
        verdict = "Mixed / partial computational evidence (see below)"
    return {"lines": lines, "verdict": verdict, "flags": flags, "consensus": consensus}


# ── Batch 3: per-gene context (fetched once per gene, cached in ctx) ─────────
def gene_context(hgnc_id):
    """ACMG-adjacent gene block: gnomAD constraint + ClinGen dosage + gene-disease
    validity + inheritance (GenCC/validity). Fetched once per gene."""
    con = map_all(hgnc_id, ">>hgnc>>gnomad_constraint")
    constraint = None
    if con:
        c = con[0]
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
    for r in map_all(hgnc_id, ">>hgnc>>spliceai", cap=60):
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
    gene = rec.get("gene_symbol")
    pchange = rec.get("hgvs_p") or rec.get("hgvs_c")
    vtype = (rec.get("variant_type") or "variant").lower()
    lead = f"This {vtype} alters {gene} at {pchange}"
    st = rec.get("structural") or {}
    dom = next((f for f in (st.get("features") or []) if "domain" in f or "region" in f), None)
    if dom:
        lead += f", in {dom}"
    parts = [lead + "."]
    if top_fn:
        parts.append(f" {gene}'s established role includes {top_fn}.")
    am = rec.get("alphamissense")
    if am and am.get("class") == "likely_pathogenic":
        parts.append(f" AlphaMissense predicts this substitution is likely pathogenic ({am['score']}).")
    if rec.get("spliceai"):
        parts.append(" SpliceAI predicts a splice-altering effect.")
    if disease_pw:
        cond = (rec.get("conditions") or [{}])[0].get("name")
        parts.append(f" Loss or alteration of {gene} function is curated by Reactome as acting "
                     f"through **{disease_pw[0]['name']}**"
                     + (f", the mechanism linked to {cond}" if cond else "") + ".")
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
    for c in conditions or []:
        mid = c.get("mondo_id")
        if not mid:
            continue
        if mid not in cache:
            cache[mid] = _digest_via_mondo(mid) or _digest_via_parent(mid)
        if cache[mid]:
            return cache[mid]
    return None


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


def _build_orphanet_digest(oid, fallback_name=None):
    o = _orphanet_entry(oid)
    if not o:
        return None
    phen = sorted((o.get("phenotypes") or []),
                  key=lambda p: -(p.get("frequency_value") or 0))
    prev = (o.get("prevalences") or [{}])[0]
    pc = prev.get("prevalence_class")
    return {
        "name": o.get("name") or fallback_name,
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
    out = []
    for f in intervals:
        if f["begin"] <= pos <= f["end"]:
            phrase = _UF_KEEP[f["type"]].replace("{d}", f["desc"] or f["type"])
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
    scored.sort(key=lambda x: (-x[0], x[1]["canonical_slug"]))
    return [{"label": f"{r['gene_symbol']} {r.get('hgvs_p') or r.get('hgvs_c')}",
             "slug": r["canonical_slug"]} for _, r in scored[:6]]


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
        meaning = "of uncertain or conflicting significance (clinical labs disagree on it)"
    elif "likely pathogenic" in cls:
        meaning = "considered likely disease-causing"
    elif "pathogenic" in cls:
        if discordant:
            meaning = "reported as disease-causing, though computational predictors disagree"
        elif stars >= 2:
            meaning = "considered disease-causing"
        else:
            meaning = "reported as disease-causing, but on limited review"
    else:
        meaning = f"classified as {rec.get('classification')}"
    gene = rec.get("gene_symbol")
    n = rec.get("submitter_count") or 0
    who = (f", submitted by {n} clinical lab" + ("s" if n != 1 else "")) if n else ""
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
