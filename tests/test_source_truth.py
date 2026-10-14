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
    # Exome-only. gnomAD's own AF uses the JOINT denominator: the genome callset called
    # 114,170 alleles at this site and found none, which is evidence the estimate needs,
    # not padding. 8.25e-07, not the observing-callset 9.11e-07.
    ("X-154032491-C-T", 8.252854249642239e-07),
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
    # The JOINT denominator, which is what gnomAD's own page reports and what our page
    # must match: every callset's AN counts, including one that observed no copies,
    # because gnomAD's AN is alleles CALLED at this site. Summing only the observing
    # callsets (which this test used to do) over-reports an exome-only variant and
    # leaves the headline figure on a different denominator from grpmax and FAF95.
    tot_ac = tot_an = 0.0
    for s in (v.get("exome"), v.get("genome")):
        if s and s.get("an"):
            tot_ac += s.get("ac") or 0
            tot_an += s["an"]
    got = tot_ac / tot_an if tot_an else None
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


# Checked against the Ensembl VCF that biobtree actually ingests, via a remote tabix
# query. That is strictly better than going through Ensembl VEP: VEP serves the MASKED
# release and we serve the RAW one, so a VEP comparison tests two conventions at once
# and cannot distinguish our bug from their convention.
#
# Measured 2026-10-11 over the 318 scored variants in the MECP2+PTEN caches: raw and
# masked differ on at least one delta for 65.7% of them, and for 19.8% the masked
# scores all fall under 0.2 — so for about one in five variants where we show a score
# and a PP3-supporting band, a masked lookup shows no splice signal at all. Masking
# zeroes a LOSS scored against a position that is not an annotated site AND a GAIN that
# would merely strengthen an existing one, so neither direction is safe to assume.
SPLICEAI_SRC = ("https://ftp.ensembl.org/pub/current_variation/SpliceAI/"
                "spliceai_scores.{kind}.snv.ensembl_mane_v1.4.grch38.vcf.gz")

# (chrom, pos, ref, alt, our four deltas AG/AL/DG/DL as we display them)
SPLICEAI_CASES = [
    ("X", 154032202, "C", "T", (0.0, 0.0, 0.78, 0.52)),
    ("X", 154032491, "C", "T", (0.0, 0.0, 0.37, 0.0)),
    # the other allele at pten-p-pro38his, whose score the below-threshold row quotes.
    # Note AL is 0.34 in the RAW release; the 0.0 a VEP lookup shows there is masking.
    ("10", 87894058, "C", "T", (0.47, 0.34, 0.0, 0.0)),
    ("10", 87894028, "T", "G", (0.01, 0.21, 0.0, 0.0)),
]


def _tabix(kind, region):
    """Remote tabix against the published VCF. Skips unless tabix is installed."""
    import shutil
    import subprocess
    exe = shutil.which("tabix") or "/data/miniconda3/envs/biobtree/bin/tabix"
    if not os.path.exists(exe):
        pytest.skip("tabix not available for a remote VCF query")
    try:
        out = subprocess.run([exe, SPLICEAI_SRC.format(kind=kind), region],
                             capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        pytest.skip("remote tabix query timed out")
    if out.returncode != 0:
        pytest.skip(f"remote tabix failed: {out.stderr[:160]}")
    return out.stdout


def _deltas(vcf_text, ref, alt):
    """DS_AG/AL/DG/DL for one allele from SpliceAI's INFO field.
    Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"""
    for line in vcf_text.splitlines():
        f = line.split("\t")
        if len(f) < 8 or f[3] != ref or f[4] != alt or "SpliceAI=" not in f[7]:
            continue
        parts = f[7].split("SpliceAI=")[1].split(",")[0].split("|")
        if len(parts) >= 10:
            return tuple(float(x) for x in parts[2:6])
    return None


@pytest.mark.parametrize("ch,pos,ref,alt,ours", SPLICEAI_CASES)
def test_spliceai_deltas_match_the_raw_release_we_ingest(ch, pos, ref, alt, ours):
    got = _deltas(_tabix("raw", f"{ch}:{pos}-{pos}"), ref, alt)
    assert got is not None, f"raw release has no {ref}>{alt} at {ch}:{pos}"
    for name, g, o in zip(("DS_AG", "DS_AL", "DS_DG", "DS_DL"), got, ours):
        assert g == pytest.approx(o, abs=0.005), (
            f"{ch}:{pos}{ref}>{alt} {name}: raw release has {g}, we display {o}")


def test_the_below_threshold_claim_holds_in_the_raw_release():
    """pten-p-pro38his (C>A) renders "all four deltas below 0.2". Check that against
    the release itself: if any delta ever reaches 0.2 the sentence becomes false, and
    it is the sentence that turns an absent score into negative evidence."""
    got = _deltas(_tabix("raw", "10:87894058-87894058"), "C", "A")
    assert got is not None, "raw release has no C>A at 10:87894058"
    assert max(got) < 0.2, f"a delta now reaches the threshold: {got}"
    # and the allele we quote as the scored one must still be the larger
    other = _deltas(_tabix("raw", "10:87894058-87894058"), "C", "T")
    assert max(other) == pytest.approx(0.47, abs=0.005), other


def test_the_raw_masked_gap_is_real_and_tracked():
    """PTEN c.83T>G is the case that exposed the convention difference: we show raw
    acceptor-loss 0.21, Ensembl VEP shows 0. Neither is wrong — biobtree serves raw
    faithfully. This pins the gap so that if upstream ever switches which release it
    publishes, the page copy about it stops being true and a test says so."""
    raw = _deltas(_tabix("raw", "10:87894028-87894028"), "T", "G")
    msk = _deltas(_tabix("masked", "10:87894028-87894028"), "T", "G")
    assert raw and msk, "could not read both releases"
    assert raw[1] == pytest.approx(0.21, abs=0.005), f"raw DS_AL moved: {raw}"
    assert msk[1] == 0.0, f"masked DS_AL is no longer zeroed: {msk}"
    # Only the loss differs here; everything else is identical between the releases.
    assert raw[0] == msk[0] and raw[2] == msk[2] and raw[3] == msk[3], (raw, msk)


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


# CADD's public API serves the whole-genome SNV file only, so an indel score cannot be
# checked through it. The indel file is published with a .tbi, so check it there — this
# is the slice the ClinVar ref/alt unlock actually opened up, and the only way to
# confirm it is against the file CADD distributes.
CADD_INDEL_CASES = [
    # HFE p.His341fs, a 14 bp deletion. Had no coordinate at all before the unlock, so
    # it carried no CADD; now live on prod at 34.0.
    ("6", 26094198, "GGCACTACGTCTTA", "G", 34.0),
]


@pytest.mark.parametrize("ch,pos,ref,alt,phred", CADD_INDEL_CASES)
def test_cadd_indel_score_matches_the_published_indel_file(ch, pos, ref, alt, phred):
    import shutil
    import subprocess
    exe = shutil.which("tabix") or "/data/miniconda3/envs/biobtree/bin/tabix"
    if not os.path.exists(exe):
        pytest.skip("tabix not available for a remote lookup")
    url = ("https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/"
           "gnomad.genomes.r4.0.indel.tsv.gz")
    try:
        out = subprocess.run([exe, url, f"{ch}:{pos}-{pos}"],
                             capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        pytest.skip("remote tabix query timed out")
    if out.returncode != 0:
        pytest.skip(f"remote tabix failed: {out.stderr[:160]}")
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 6 and f[2] == ref and f[3] == alt:
            assert float(f[5]) == pytest.approx(phred, abs=0.005), (
                f"CADD now scores {ch}:{pos}{ref}>{alt} at {f[5]}, we display {phred}")
            return
    pytest.fail(f"{ch}:{pos} {ref}>{alt} is not in CADD's indel file")
