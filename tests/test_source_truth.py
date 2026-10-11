"""Cross-checks against the UPSTREAM sources, not against biobtree.

Every other test in this repo verifies that we render our own data correctly. These
verify that the data itself matches what the authoritative source says — ClinVar's
own API for the VCF triplet, gnomAD's for frequency and absence, CADD's for PHRED,
UCSC's for repeat tracks, Ensembl VEP's for the SpliceAI deltas, and Pangolin's
published archive for its scores.

They are OPT-IN because they hit live third-party APIs: slow, rate-limited, and able
to fail for reasons that have nothing to do with this codebase. A red build must mean
our bug, so they are skipped unless asked for:

    SUGI_SOURCE_TRUTH=1 python3 -m pytest tests/test_source_truth.py -v

The expected values below were verified by hand on 2026-10-11. Where one disagrees,
the disagreement itself is documented rather than asserted away — see the SpliceAI
case, which is a genuine raw-vs-masked convention difference and not an error.
"""
import json
import os
import time
import urllib.error
import urllib.request

import pytest

pytestmark = pytest.mark.skipif(
    not os.environ.get("SUGI_SOURCE_TRUTH"),
    reason="set SUGI_SOURCE_TRUTH=1 to cross-check against live upstream APIs")

UA = {"User-Agent": "sugi-variant source-truth check", "Content-Type": "application/json"}


def _get(url, timeout=60, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA["User-Agent"]})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if i == tries - 1:
                pytest.skip(f"upstream unreachable: {url}")
            time.sleep(3)


def _graphql(url, query, timeout=60, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=json.dumps({"query": query}).encode(),
                                         headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code not in (429, 502, 503) or i == tries - 1:
                pytest.fail(f"gnomAD API {e.code}: {e.read().decode()[:200]}")
            time.sleep(5)
        except Exception:
            if i == tries - 1:
                pytest.skip("gnomAD API unreachable")
            time.sleep(3)


# (variation_id, our coordinate, ClinVar's start) — start differs from positionVCF only
# on the deletion, which is the whole reason variant_coordinate uses positionVCF.
CLINVAR_CASES = [
    (1002925, "X:154030684:GTGGCACGGGGGCCTT:G", "154030685"),   # MECP2 15 bp del
    (1174899, "X:154097616:T:TCGC", "154097616"),               # MECP2 dup
    (1002697, "10:87960885:T:TC", "87960885"),                  # PTEN ins
    (1048544, "10:87965368:A:AGT", "87965368"),                 # PTEN dup
    (1003503, "10:87894028:T:G", "87894028"),                   # PTEN SNV
]


@pytest.mark.parametrize("vid,coord,start", CLINVAR_CASES)
def test_coordinate_matches_clinvars_own_vcf_fields(vid, coord, start):
    """The key we build must equal ClinVar's own positionVCF/ref/alt, and we must not
    have silently used `start` — which differs on the deletion."""
    import re
    x = _get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
             f"?db=clinvar&rettype=vcv&is_variationid&id={vid}")
    hit = None
    for m in re.finditer(r'<SequenceLocation([^>]*positionVCF[^>]*)/>', x):
        a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        if a.get("Assembly") == "GRCh38":
            hit = a
            break
    assert hit, f"ClinVar returned no GRCh38 positionVCF for {vid}"
    theirs = (f"{hit['Chr']}:{hit['positionVCF']}:{hit['referenceAlleleVCF']}"
              f":{hit['alternateAlleleVCF']}")
    assert theirs == coord, f"VCV{vid}: ClinVar says {theirs}, we build {coord}"
    assert hit["start"] == start, f"VCV{vid}: ClinVar start moved to {hit['start']}"


# (coordinate, expected gnomAD v4.1 AF or None when genuinely absent)
GNOMAD_CASES = [
    ("X-154097616-T-TCGC", 3.668244011591651e-05),   # an INDEL with a real frequency
    ("X-154032491-C-T", 9.111351650794692e-07),      # exome-only -> denominator fix
    ("X-154032202-C-T", None),
    ("10-87894028-T-G", None),
    ("10-87894058-C-A", None),
    ("10-87965368-A-AGT", None),                     # indel, genuinely absent
]


@pytest.mark.parametrize("vid,expect_af", GNOMAD_CASES)
def test_frequency_and_absence_match_gnomad(vid, expect_af):
    """Our AF must equal gnomAD's, and — more important — anything we report as ABSENT
    must genuinely be absent there, because that renders as an ACMG rarity claim."""
    q = ('{ variant(variantId:"%s", dataset:gnomad_r4){ variant_id'
         ' exome{ac an af} genome{ac an af} } }') % vid
    d = _graphql("https://gnomad.broadinstitute.org/api", q)
    v = (d.get("data") or {}).get("variant")
    if expect_af is None:
        assert v is None, f"{vid}: we claim absence but gnomAD has it: {v}"
        return
    assert v, f"{vid}: gnomAD no longer has a variant we report a frequency for"
    obs_ac = obs_an = 0.0
    for s in (v.get("exome"), v.get("genome")):
        if s and s.get("ac"):
            obs_ac += s["ac"]
            obs_an += s["an"]
    got = obs_ac / obs_an if obs_an else None
    assert got == pytest.approx(expect_af, rel=1e-9), f"{vid}: gnomAD now {got}"


CADD_CASES = [("X", 154032202, "C", "T", 34.0), ("X", 154032491, "C", "T", 18.21),
              ("10", 87894028, "T", "G", 32.0), ("10", 87894058, "C", "A", 31.0)]


@pytest.mark.parametrize("ch,pos,ref,alt,phred", CADD_CASES)
def test_cadd_phred_matches_cadd(ch, pos, ref, alt, phred):
    rows = json.loads(_get(f"https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{ch}:{pos}"))
    hit = [r for r in rows if r.get("Ref") == ref and r.get("Alt") == alt]
    assert hit, f"CADD has no {ref}>{alt} at {ch}:{pos}"
    assert float(hit[0]["PHRED"]) == pytest.approx(phred, abs=0.005)


# (chrom, pos, track, a name we display)
UCSC_CASES = [("chrX", 154097616, "rmsk", "(CTC)n"),
              ("chr10", 87960885, "rmsk", "(T)n"),
              ("chr10", 87965368, "genomicSuperDups", "chr9:33674359")]


@pytest.mark.parametrize("chrom,pos,track,name", UCSC_CASES)
def test_repeat_annotation_matches_ucsc(chrom, pos, track, name):
    """We present repeat context as a confidence qualifier, so the track and the repeat
    name both have to be the ones UCSC actually reports at that position."""
    d = json.loads(_get("https://api.genome.ucsc.edu/getData/track"
                        f"?genome=hg38;track={track};chrom={chrom};start={pos - 1};end={pos}"))
    rows = d.get(track) or []
    if isinstance(rows, dict):
        rows = sum(rows.values(), [])
    names = {r.get("repName") or r.get("name") for r in rows}
    assert name in names, f"UCSC {track} at {chrom}:{pos} now reports {names}"


SPLICEAI_CASES = [
    # (chrom, pos, alt, our four deltas AG/AL/DG/DL, whether VEP's masked AL/DL agree)
    ("X", 154032202, "T", (0.0, 0.0, 0.78, 0.52), True),
    ("X", 154032491, "T", (0.0, 0.0, 0.37, 0.0), True),
    ("10", 87894058, "T", (0.47, 0.0, 0.0, 0.0), True),
    # PTEN c.83T>G: we carry RAW acceptor-loss 0.21 with dp_al -10; VEP serves MASKED
    # and reports 0, because 87894018 is not an annotated acceptor (exon 2 starts at
    # 87894025). Every other field of the two records is identical. Documented, not
    # asserted away — the gain deltas must still agree exactly.
    ("10", 87894028, "G", (0.01, 0.21, 0.0, 0.0), False),
]


@pytest.mark.parametrize("ch,pos,alt,ours,loss_agrees", SPLICEAI_CASES)
def test_spliceai_deltas_match_ensembl_vep(ch, pos, alt, ours, loss_agrees):
    d = json.loads(_get(f"https://rest.ensembl.org/vep/human/region/{ch}:{pos}-{pos}:1/{alt}"
                        "?content-type=application/json&SpliceAI=1"))
    r0 = d[0] if isinstance(d, list) else d
    sp = next((tc["spliceai"] for tc in (r0.get("transcript_consequences") or [])
               if tc.get("spliceai")), None)
    assert sp, f"VEP returned no spliceai field for {ch}:{pos}>{alt}"
    ag, al, dg, dl = ours
    # Gains are a NEW site by definition, so masking never suppresses them: these must
    # always agree, and they are the check that we are reading the same computation.
    assert float(sp["DS_AG"]) == pytest.approx(ag, abs=0.005), "acceptor gain diverged"
    assert float(sp["DS_DG"]) == pytest.approx(dg, abs=0.005), "donor gain diverged"
    if loss_agrees:
        assert float(sp["DS_AL"]) == pytest.approx(al, abs=0.005), "acceptor loss diverged"
        assert float(sp["DS_DL"]) == pytest.approx(dl, abs=0.005), "donor loss diverged"
    else:
        assert float(sp["DS_AL"]) == 0.0 and al > 0.0, (
            "the known raw-vs-masked case changed shape: VEP now reports "
            f"AL={sp['DS_AL']} where we carry {al}. Re-verify which convention each "
            "side serves before trusting the rendered text.")


PANGOLIN_CASES = [("ENSG00000169057", "154032491", "C", "T", 0.23, -0.0),
                  ("ENSG00000169057", "154032202", "C", "T", 0.32, -0.64),
                  ("ENSG00000171862", "87894028", "T", "G", 0.1, -0.25),
                  # both deltas under 0.2 -> correctly ABSENT from the filtered release
                  ("ENSG00000171862", "87894058", "C", "A", 0.03, -0.03)]


@pytest.mark.parametrize("ensg,pos,ref,alt,gain,loss", PANGOLIN_CASES)
def test_pangolin_scores_match_the_published_archive(ensg, pos, ref, alt, gain, loss):
    """Verified against Pangolin's own Zenodo archive, which is the source biobtree
    ingests. The last case is the threshold semantics: both deltas below 0.2, so the
    variant is absent from the filtered release and we must show no score."""
    import gzip
    import io as _io
    import zipfile
    Z = "/data/biobtree/raw_data/pangolin/Pangolin_hg38_snvs_masked.zip"
    if not os.path.exists(Z):
        pytest.skip("Pangolin source archive not present on this host")
    z = zipfile.ZipFile(Z)
    member = f"Pangolin_hg38_snvs_masked/{ensg}.tsv.gz"
    if member not in z.namelist():
        pytest.skip(f"{ensg} not in archive")
    with z.open(member) as fh:
        data = gzip.GzipFile(fileobj=_io.BytesIO(fh.read()))
        for line in data:
            s = line.decode().rstrip("\n").split("\t")
            if len(s) >= 7 and s[1] == pos and s[2] == ref and s[3] == alt:
                assert float(s[4]) == pytest.approx(gain, abs=0.005)
                assert float(s[6]) == pytest.approx(loss, abs=0.005)
                if max(abs(float(s[4])), abs(float(s[6]))) < 0.2:
                    from sugivariant.enrich import _PANGOLIN_THRESHOLD
                    assert _PANGOLIN_THRESHOLD == 0.2
                return
    pytest.fail(f"{ensg} {pos}{ref}>{alt} not found in the Pangolin archive")
