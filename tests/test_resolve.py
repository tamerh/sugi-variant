"""Resolver match-keys: a GENE + HGVS/AA query in any shape must normalize to the
same key the target record exposes. Pure (no biobtree / no server)."""
from sugivariant.slug import _norm_hgvs
from app import _hgvs_keys, _RSID_RE


PTEN_R130Q = {"hgvs_p": "p.Arg130Gln", "hgvs_c": "c.389G>A"}
ACTA1_P309A = {"hgvs_p": "p.Pro309Ala", "hgvs_c": "c.925C>G"}


def _matches(rec, query_hgvs):
    return _norm_hgvs(query_hgvs) in _hgvs_keys(rec)


def test_all_query_shapes_hit_the_same_record():
    # protein: 3-letter (with/without p.), 1-letter short
    assert _matches(PTEN_R130Q, "p.Arg130Gln")
    assert _matches(PTEN_R130Q, "Arg130Gln")
    assert _matches(PTEN_R130Q, "R130Q")
    # coding: with/without c.
    assert _matches(PTEN_R130Q, "c.389G>A")
    assert _matches(PTEN_R130Q, "389G>A")
    # the bare 3-letter demand shape
    assert _matches(ACTA1_P309A, "pro309ala")
    assert _matches(ACTA1_P309A, "P309A")


def test_no_false_match():
    assert not _matches(PTEN_R130Q, "R130W")     # different residue
    assert not _matches(PTEN_R130Q, "c.388G>A")  # different position


def test_rsid_pattern():
    assert _RSID_RE.match("rs121909229") and _RSID_RE.match("RS80357906")
    assert not _RSID_RE.match("rs") and not _RSID_RE.match("pten-p-arg130gln")


def test_non_missense_has_no_short_key():
    rec = {"hgvs_p": "p.Lys328del", "hgvs_c": "c.983_985del"}
    keys = _hgvs_keys(rec)
    assert _norm_hgvs("c.983_985del") in keys
    assert _norm_hgvs("p.Lys328del") in keys  # still matchable by the p. form


def test_multi_variant_query_routes_to_set():
    """A comma/newline/semicolon list of 2+ items is a set, not a single lookup —
    the home handler must redirect it to /set (pure: no biobtree call)."""
    from starlette.testclient import TestClient
    import app
    c = TestClient(app.app, follow_redirects=False)
    r = c.get("/", params={"q": "PTEN R130Q, BRCA2 c.1000C>T, rs121909229"})
    assert r.status_code == 307 and "/set?q=" in r.headers["location"]
    # a single variant still goes through single-variant resolution (not /set)
    assert app._SET_SPLIT.split("PTEN R130Q") == ["PTEN R130Q"]
    assert len([p for p in app._SET_SPLIT.split("A;B\nC") if p]) == 3


def test_a_wrong_verb_on_a_page_lands_on_the_page_not_on_an_error():
    """Every route in app.py is @app.get, no template has a POST form and the JS does no
    fetch, so a 405 is always a readable URL reached with the wrong verb. It happened for
    real: mobile Firefox POSTed /variant/jak3-c-2351-minus-15a-g on a browser refresh
    (2026-10-15, referer the page itself) and the visitor got a styled "Method Not
    Allowed" page for a variant that renders fine. A reload re-issues the method of the
    history entry, and nothing we serve creates a POST entry; the Cloudflare managed
    challenge sitting on /variant/ does, since it is served at the requested URL.

    303 and not 302: only 303 obliges the browser to redo the request as GET."""
    from starlette.testclient import TestClient
    import app
    c = TestClient(app.app, follow_redirects=False)
    r = c.post("/jak3-c-2351-minus-15a-g")
    assert r.status_code == 303, (
        "a POST on a page route still shows an error page instead of the page")
    assert r.headers["location"].endswith("/jak3-c-2351-minus-15a-g")
    # App-level intent only: the sugi.bio nginx vhost hides Cache-Control on /variant/
    # and stamps its own, so this never reaches Cloudflare. Asserted so the app stays
    # correct on its own terms, and correct if that proxy rule is ever lifted.
    assert "no-store" in r.headers.get("cache-control", "")
    # Only the method is forgiven. A genuinely unknown slug is still a 404 error page,
    # or this handler would paper over every miss with a redirect loop.
    assert c.get("/definitely-not-a-real-variant-slug").status_code == 404
