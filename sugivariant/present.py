"""Presentation state for the variant page.

Every bug in the 2026-10 evidence audit broke the same rule: a visual claim was
driven by something that merely correlated with the reader's question rather than
by the question itself.

  · the in-silico consensus took the agreement colour from `unanimous`, which says
    the predictors agree WITH EACH OTHER — so "both predictors call this tolerated"
    rendered green beneath a Likely pathogenic classification;
  · the triage verdict's colour came from `'triage' in verdict`, a substring of its
    own prose, so rewording the sentence would have silently restyled it;
  · `{% if v.spliceai %}` was true for {"status": "absent"}, so a page claimed
    SpliceAI "carries the signal" where it had never been assessed;
  · the gnomAD tile chose its word through a four-branch if/elif chain inlined in
    the template, which is where the "Rare" / "absent from gnomAD" contradiction
    could survive unseen.

So: this module computes a NAMED state for each element, and TONE maps state to a
CSS class in exactly one place. Templates read `ui.<element>.state`, `.label`,
`.detail` and `.tone` and make no judgement of their own. A state with no tone
entry raises rather than rendering unstyled, and the tests assert that templates
contain no semantic-class logic.

Pure: it reads the record and adds no lookups.
"""

# Tone is the ONLY place a state becomes a colour. `None` means "no colour" —
# which is the right answer far more often than a palette suggests, and must stay
# an explicit choice rather than a fallthrough.
TONE = {
    # in-silico consensus, relative to the ClinVar classification
    "concordant": "is-agree",
    "discordant": "is-flag",
    "neutral": None,
    "mixed": None,
    # AlphaMissense, relative to a pathogenic ClinVar call
    "am_concordant": "is-ok",
    "am_discordant": "is-flag",
    "am_no_call": None,
    "am_not_comparable": None,
    # population frequency
    "freq_common": "is-flag",
    "freq_rare": None,
    "freq_absent": None,
    "freq_unknown": None,
    "freq_none": None,
    # verdict line
    "triage": "is-flag",
    "plain": None,
}


class UnknownState(KeyError):
    """A state with no TONE entry. Raised rather than silently rendering unstyled,
    because an unstyled element is indistinguishable from a deliberate neutral."""


def tone(state):
    if state not in TONE:
        raise UnknownState(f"{state!r} has no TONE entry — add one deliberately")
    return TONE[state] or ""


def _el(state, label=None, detail=None, **extra):
    return {"state": state, "label": label, "detail": detail, "tone": tone(state), **extra}


def frequency_view(g):
    """The gnomAD tile. `absent` and `unknown` are different claims and must never
    collapse: absence is ACMG PM2-supporting evidence, a missing number is not."""
    if not g:
        return _el("freq_none", "—", "no coordinate frequency")
    band = g.get("band") or ""
    if g.get("unknown"):
        return _el("freq_unknown", "Not reported", band)
    if g.get("absent"):
        return _el("freq_absent", "Absent", band)
    if g.get("is_common"):
        return _el("freq_common", "Common", band)
    return _el("freq_rare", "Rare", band)


def alphamissense_view(am, classification):
    """The ✓/⚠ badge on the AlphaMissense row. AlphaMissense makes a THREE-way call
    and `ambiguous` is an abstention, so it earns no badge in either direction."""
    cls = (classification or "").lower()
    cv_path = "pathogenic" in cls and "conflict" not in cls
    amclass = (am or {}).get("class")
    if not amclass:
        return _el("am_not_comparable")
    if amclass == "ambiguous":
        return _el("am_no_call", None, "no call either way")
    if not cv_path:
        # nothing settled to be concordant WITH; the score still renders, unbadged
        return _el("am_not_comparable")
    if amclass == "likely_pathogenic":
        return _el("am_concordant", "✓ concordant with ClinVar")
    return _el("am_discordant", "⚠ discordant with ClinVar")


def consensus_view(c):
    """The consensus tile and line. `stance` is the relation to the clinical call;
    `unanimous` describes only the predictors' agreement with each other and is
    reported as text, never as colour."""
    cons = (c or {}).get("consensus")
    if not cons:
        return None
    abstained = cons.get("abstained") or []
    sub = "unanimous" if cons.get("unanimous") else "mixed"
    if abstained:
        sub += " · " + ", ".join(abstained) + " ambiguous"
    return _el(cons.get("stance") or "neutral",
               f"{cons.get('n_damaging')}/{cons.get('total')} damaging", sub,
               summary=cons.get("summary"), abstained=abstained)


def verdict_view(c):
    """The verdict sentence. Its kind is set by the branch that writes it, so the
    styling cannot drift when the wording changes."""
    c = c or {}
    if not c.get("verdict"):
        return None
    return _el(c.get("verdict_kind") or "plain", c["verdict"])


def evidence_view(rec, freq_g=None):
    """Every semantic state the variant page needs, in one call."""
    c = rec.get("concordance") or {}
    return {"frequency": frequency_view(freq_g if freq_g is not None else rec.get("gnomad")),
            "alphamissense": alphamissense_view(rec.get("alphamissense"),
                                                rec.get("classification")),
            "consensus": consensus_view(c),
            "verdict": verdict_view(c)}
