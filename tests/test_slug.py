"""Variant slug parsing/generation — the make-or-break piece (HGVS isn't
biobtree-searchable, so the URL must be derivable from GENE + HGVS in BOTH the
p. and c. forms). Pins the exact GSC demand queries → their slugs."""
from sugivariant.slug import parse_hgvs, variant_slugs


def test_parse_hgvs_from_clinvar_name():
    c, p = parse_hgvs("NM_001100.4(ACTA1):c.925C>G (p.Pro309Ala)")
    assert c == "c.925C>G"
    assert p == "p.Pro309Ala"


def test_parse_hgvs_partial():
    assert parse_hgvs("NM_001100.4(ACTA1):c.616-4C>G") == ("c.616-4C>G", None)
    assert parse_hgvs("") == (None, None)


def test_demand_queries_map_to_slugs():
    # The real GSC demand: `acta1 "pro309ala"` AND `"c.925c>g" acta1` must both
    # land the SAME variant page.
    c, p = parse_hgvs("NM_001100.4(ACTA1):c.925C>G (p.Pro309Ala)")
    canonical, slugs = variant_slugs("ACTA1", c, p)
    assert canonical == "acta1-p-pro309ala"        # p. form is canonical (more-searched)
    assert "acta1-p-pro309ala" in slugs
    assert "acta1-c-925c-g" in slugs               # c. form alias → same page


def test_variant_slugs_dedup_and_empty():
    # p.-only (no c.) → single slug; no HGVS → nothing.
    canonical, slugs = variant_slugs("TP53", None, "p.Arg175His")
    assert canonical == "tp53-p-arg175his" and slugs == ["tp53-p-arg175his"]
    assert variant_slugs("ACTA1", None, None) == (None, [])


def test_hgvs_operators_dont_collide():
    # audit P1: *, +, - are meaningful and must produce DISTINCT slugs
    from sugivariant.slug import variant_slugs
    utr = variant_slugs("PTEN", "c.*667A>T", None)[0]      # 3'UTR
    cod = variant_slugs("PTEN", "c.667A>T", None)[0]       # coding
    assert utr != cod and "star" in utr
    intron_after = variant_slugs("ACTA1", "c.616+4C>G", None)[0]
    intron_before = variant_slugs("ACTA1", "c.616-4C>G", None)[0]
    assert intron_after != intron_before                    # opposite introns
    assert "plus" in intron_after and "minus" in intron_before


def test_ncrna_and_mito_hgvs_parse_and_slug():
    # ncRNA (n.) and mitochondrial (m.) variants have no p. change — they must still
    # parse + slug (else they're dropped: RMRP, RNU4-2, TERC, MT-TL1…).
    c, p = parse_hgvs("NR_003051.4(RMRP):n.71A>G")
    assert c == "n.71A>G" and p is None
    assert variant_slugs("RMRP", c, p)[0] == "rmrp-n-71a-g"
    c2, _ = parse_hgvs("NC_012920.1(MT-TL1):m.3243A>G")
    assert c2 == "m.3243A>G"
    assert variant_slugs("MT-TL1", c2, None)[0] == "mt-tl1-m-3243a-g"
    # the c/n/m match must not fire inside a transcript token (the M in NM_)
    assert parse_hgvs("NM_000546.6(TP53):c.743G>A")[0] == "c.743G>A"


def test_gene_symbol_with_dash_preserved():
    # a gene like NKX2-1 must not get mangled by the HGVS operator mapping
    from sugivariant.slug import variant_slugs
    canonical, _ = variant_slugs("NKX2-1", None, "p.Arg100His")
    assert canonical == "nkx2-1-p-arg100his"


def test_name_gene_prefers_transcript_gene():
    from sugivariant.slug import name_gene
    # authoritative gene = the transcript in the name (not ClinVar's overlap gene_symbol)
    assert name_gene("NM_001267550.2(TTN):c.70690_70691dup (p.Thr23565fs)") == "TTN"
    assert name_gene("NR_003051.4(RMRP):n.71A>G") == "RMRP"
    assert name_gene("NC_012920.1(MT-TL1):m.3243A>G") == "MT-TL1"
    assert name_gene("NM_004006.3(NKX2-1):c.100A>G") == "NKX2-1"   # dash in gene preserved
    assert name_gene("no gene here") is None
