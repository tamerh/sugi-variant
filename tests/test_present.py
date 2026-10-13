"""Presentation states.

The rule this module exists to enforce: a visual claim is driven by an explicit
named state computed in Python, never by truthiness, a substring of prose, or an
internal measure that merely correlates with the reader's question. Each test
below names the real defect it guards against.
"""
import pathlib
import re

import pytest

from sugivariant.present import (TONE, UnknownState, alphamissense_view,
                                 consensus_view, evidence_view, frequency_view,
                                 tone, verdict_view)

TEMPLATES = pathlib.Path(__file__).resolve().parent.parent / "templates"
SEMANTIC = ("is-agree", "is-flag", "is-ok")


# ── the contract ─────────────────────────────────────────────────────────────
def test_every_state_any_view_can_emit_has_a_tone():
    """A state with no TONE entry must raise, not render unstyled: an unstyled
    element is indistinguishable from a deliberate neutral, so a typo'd or new state
    would silently drop a colour the reader needs."""
    with pytest.raises(UnknownState):
        tone("not_a_real_state")
    seen = set()
    for g in ({}, None, {"unknown": True, "band": "b"}, {"absent": True, "band": "b"},
              {"is_common": True, "band": "b"}, {"af": 1e-5, "band": "b"}):
        seen.add(frequency_view(g)["state"])
    for am in (None, {"class": "ambiguous"}, {"class": "likely_pathogenic"},
               {"class": "likely_benign"}):
        for cls in ("Pathogenic", "Benign", "Uncertain significance", None):
            seen.add(alphamissense_view(am, cls)["state"])
    for st in ("concordant", "discordant", "neutral", "mixed"):
        seen.add(consensus_view({"consensus": {"stance": st, "n_damaging": 1,
                                               "total": 2}})["state"])
    for k in (None, "triage"):
        seen.add(verdict_view({"verdict": "x", "verdict_kind": k})["state"])
    missing = seen - set(TONE)
    assert not missing, f"states with no TONE entry: {missing}"


def test_neutral_is_an_explicit_choice_not_a_fallthrough():
    """`None` in TONE means "deliberately no colour". It must be spelled out, so that
    adding a state and forgetting its tone fails loudly instead of rendering plain."""
    for state, css in TONE.items():
        assert css is None or isinstance(css, str) and css.startswith("is-"), (state, css)


# ── frequency: absence and missing data are different claims ─────────────────
def test_absent_and_unknown_never_collapse():
    """Absence from gnomAD is ACMG PM2-supporting evidence; a missing number is not.
    Two code paths used to turn the second into the first on 22.6% of pathogenic/LP
    records."""
    assert frequency_view({"absent": True, "band": "absent from gnomAD v4.1"})["label"] == "Absent"
    unk = frequency_view({"unknown": True, "band": "no allele count reported"})
    assert unk["label"] == "Not reported"
    assert unk["state"] != frequency_view({"absent": True, "band": ""})["state"]
    assert frequency_view(None)["label"] == "—"


def test_only_a_common_frequency_is_coloured():
    """Rare and absent are ordinary findings. Colouring them would imply a judgement
    the page does not make."""
    assert tone(frequency_view({"is_common": True, "band": ""})["state"])
    for g in ({"absent": True, "band": ""}, {"af": 1e-6, "band": ""},
              {"unknown": True, "band": ""}, None):
        assert tone(frequency_view(g)["state"]) == ""


# ── AlphaMissense: three-way, and an abstention earns no badge ───────────────
def test_an_ambiguous_call_gets_no_badge_in_either_direction():
    v = alphamissense_view({"class": "ambiguous"}, "Likely pathogenic")
    assert v["label"] is None and v["tone"] == ""


def test_the_badge_only_appears_where_there_is_a_settled_call_to_compare_with():
    assert alphamissense_view({"class": "likely_pathogenic"}, "Uncertain significance")["label"] is None
    assert alphamissense_view({"class": "likely_pathogenic"},
                              "Conflicting classifications of pathogenicity")["label"] is None
    assert "concordant" in alphamissense_view({"class": "likely_pathogenic"}, "Pathogenic")["label"]
    assert "discordant" in alphamissense_view({"class": "likely_benign"}, "Likely pathogenic")["label"]


# ── consensus: colour is the stance to ClinVar, not internal unanimity ───────
def test_unanimity_is_reported_as_text_and_never_as_colour():
    """"both scoring predictors call this tolerated" rendered GREEN beneath a Likely
    pathogenic classification, because the class came from `unanimous` — which says
    only that the predictors agree with each other."""
    disc = consensus_view({"consensus": {"stance": "discordant", "unanimous": True,
                                         "n_damaging": 0, "total": 2,
                                         "abstained": ["AlphaMissense"]}})
    assert disc["tone"] == "is-flag"
    assert "unanimous" in disc["detail"], "unanimity still belongs in the text"
    assert "AlphaMissense ambiguous" in disc["detail"], "the abstention must be disclosed"
    conc = consensus_view({"consensus": {"stance": "concordant", "unanimous": True,
                                         "n_damaging": 2, "total": 2}})
    assert conc["tone"] == "is-agree"
    for st in ("neutral", "mixed"):
        assert consensus_view({"consensus": {"stance": st, "n_damaging": 1,
                                             "total": 2}})["tone"] == ""


def test_no_consensus_renders_nothing_rather_than_an_empty_tile():
    assert consensus_view({"consensus": None}) is None
    assert consensus_view({}) is None
    assert consensus_view(None) is None


# ── verdict: styled by the branch that wrote it, not by its own prose ────────
def test_the_verdict_tone_survives_a_rewording():
    """The template tested `'triage' in verdict` — presentation keyed on an incidental
    substring of the sentence it was styling."""
    reworded = verdict_view({"verdict": "ClinVar is uncertain; predictors lean damaging "
                                        "— prioritise for review", "verdict_kind": "triage"})
    assert reworded["tone"] == "is-flag", "kind, not wording, decides the colour"
    assert verdict_view({"verdict": "a triage signal", "verdict_kind": None})["tone"] == ""
    assert verdict_view({"verdict": None}) is None


# ── structural: templates must carry no semantic judgement ───────────────────
def test_no_template_decides_a_semantic_class_for_itself():
    """Every defect in this audit lived in template logic. A semantic class may be
    interpolated from a precomputed tone, never chosen by an inline conditional."""
    offenders = []
    for f in sorted(TEMPLATES.glob("*.html")):
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            for cls in SEMANTIC:
                # a literal semantic class inside a Jinja conditional expression
                for m in re.finditer(r"\{\{[^}]*'" + cls + r"'[^}]*\}\}", line):
                    offenders.append(f"{f.name}:{i} chooses {cls} inline: {m.group(0)[:70]}")
    assert not offenders, "\n".join(offenders) + "\n→ compute the state in present.py"


def test_templates_do_not_test_prose_for_meaning():
    """`'triage' in c.verdict` styled a line by grepping its own sentence."""
    bad = []
    for f in sorted(TEMPLATES.glob("*.html")):
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if re.search(r"\{%[^%]*'[a-z ]{4,}'\s+in\s+\w+\.(verdict|summary|band|reason)", line):
                bad.append(f"{f.name}:{i}: {line.strip()[:90]}")
    assert not bad, "\n".join(bad)


def test_the_page_cannot_be_rendered_without_its_presentation_states():
    """_render derives `ui` itself, so a new call site cannot produce a page whose
    styling silently differs — the footgun that let index.py read a stale readout."""
    import os
    os.environ.setdefault("BASE_PATH", "/variant")
    import app
    rec = {"canonical_slug": "x-p-y", "gene_symbol": "X", "classification": "Pathogenic",
           "alphamissense": {"class": "likely_pathogenic", "score": "0.9"},
           "concordance": {"verdict": None, "lines": [], "flags": [], "consensus": None}}
    html = app._render("variant.html", v=rec, canonical="x-p-y", nav="variant", jsonld="",
                       freq_g=None, freq_pops=[])
    assert "concordant with ClinVar" in html


def test_evidence_view_adds_no_lookups(monkeypatch):
    """It must stay pure: it runs on every page render, including over whole genes."""
    import sugivariant.enrich as EN
    def boom(*a, **k):
        raise AssertionError("evidence_view performed a lookup")
    monkeypatch.setattr(EN, "map_all", boom, raising=False)
    monkeypatch.setattr(EN, "_coord_entry", boom, raising=False)
    evidence_view({"classification": "Pathogenic",
                   "alphamissense": {"class": "ambiguous", "score": "0.5"},
                   "gnomad": {"absent": True, "band": "absent from gnomAD v4.1"},
                   "concordance": {"consensus": {"stance": "neutral", "n_damaging": 0,
                                                 "total": 2}, "verdict": "x"}})


# ── evidence-card grouping ───────────────────────────────────────────────────
def test_no_group_heading_renders_without_rows_under_it():
    """evidence_groups() mirrors the template's row conditions, so the two can drift.
    An empty heading is worse than no grouping: it implies evidence that is not there."""
    import glob
    import gzip
    import os
    import pickle
    import random
    import re as _re
    os.environ.setdefault("BASE_PATH", "/variant")
    import app
    from sugivariant.enrich import refresh_record
    from sugivariant.present import EVIDENCE_GROUPS
    titles = {t for _, t in EVIDENCE_GROUPS}
    cache = os.environ.get("RECORDS_DIR") or "/data2/sugi-variant-cache/records"
    files = sorted(glob.glob(os.path.join(cache, "*.pkl.gz")))
    if not files:
        pytest.skip("no cached records")
    random.Random(5).shuffle(files)
    checked = 0
    for f in files[:12]:
        try:
            recs = pickle.load(gzip.open(f, "rb"))
        except Exception:
            continue
        for rec in recs[:40]:
            refresh_record(rec)
            g, pops = app._freq_ctx(rec)
            html = app._render("variant.html", v=rec, canonical=rec.get("canonical_slug"),
                               nav="variant", jsonld="", freq_g=g, freq_pops=pops)
            body = _re.search(r'<table class="vd-preds"><tbody>(.*?)</tbody>', html, _re.S)
            assert body, "evidence table not rendered"
            rows = _re.findall(r'<tr[^>]*class="(vd-grp|vd-note)"[^>]*>|<tr(?![^>]*vd-)',
                               body.group(1))
            # walk the rows: a vd-grp must be followed by at least one non-group row
            seq = _re.findall(r'<tr([^>]*)>', body.group(1))
            for i, attrs in enumerate(seq):
                if "vd-grp" not in attrs:
                    continue
                rest = seq[i + 1:]
                assert rest and "vd-grp" not in rest[0], (
                    f"{rec.get('canonical_slug')}: an empty group heading")
            checked += 1
    assert checked, "no records rendered — the sweep proved nothing"


def test_the_consensus_sits_inside_the_group_it_describes():
    """It used to be a paragraph below the whole 13-row card, which is why "both
    scoring predictors" was hard to place: the sentence was nowhere near its subjects."""
    import os
    os.environ.setdefault("BASE_PATH", "/variant")
    import app
    import re as _re
    rec = {"canonical_slug": "x-p-y", "gene_symbol": "X", "classification": "Likely pathogenic",
           "alphamissense": {"class": "ambiguous", "score": "0.44"},
           "revel": {"score": 0.2, "band": "benign-leaning", "direction": "benign"},
           "saprot": {"llr": -5.5, "damaging": False},
           "cadd": {"phred": 24.2}}
    from sugivariant.enrich import refresh_record
    refresh_record(rec)
    html = app._render("variant.html", v=rec, canonical="x-p-y", nav="variant", jsonld="",
                       freq_g=None, freq_pops=[])
    body = _re.search(r'<table class="vd-preds"><tbody>(.*?)</tbody>', html, _re.S).group(1)
    assert "vd-note" in body, "the consensus must render inside the table"
    # it must come after SaProt and before the next group heading
    i_sap, i_note = body.index("SaProt"), body.index("vd-note")
    assert i_sap < i_note, "the consensus must follow the predictors it summarises"
    nxt = body.find('class="vd-grp"', i_note)
    assert nxt == -1 or nxt > i_note
    assert "REVEL and SaProt" in body, "the summary must name its subjects"
