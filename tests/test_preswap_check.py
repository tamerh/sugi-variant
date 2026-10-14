"""The pre-swap gates must each catch the failure they are named for.

Every case below is a reconstruction of something that really happened while rebuilding
the index, and in every one the totals looked fine -- which is why these gates are
slug- and gene-level rather than count-level."""
import importlib.util
import pathlib
import sqlite3

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "preswap_check", pathlib.Path(__file__).resolve().parent.parent / "scripts" / "preswap_check.py")
PS = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(PS)

SCHEMA = """
CREATE TABLE variant(vcv INTEGER PRIMARY KEY, gene TEXT, slug TEXT, hgvs_p TEXT,
  hgvs_c TEXT, classification TEXT, review_status TEXT, stars INTEGER, rsid TEXT,
  coordinate TEXT, primary_condition TEXT, flag TEXT);
CREATE TABLE alias(key TEXT, vcv INTEGER, PRIMARY KEY(key, vcv));
CREATE TABLE gene_meta(gene TEXT PRIMARY KEY, n INTEGER);
"""


def _db(path, variants):
    """variants: [(vcv, gene, slug, classification)] -> a DB with aliases + gene_meta."""
    c = sqlite3.connect(path)
    c.executescript(SCHEMA)
    for vcv, gene, slug, cls in variants:
        c.execute("INSERT INTO variant VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                  (vcv, gene, slug, None, None, cls, None, 1, None, None, None, None))
        c.execute("INSERT OR IGNORE INTO alias VALUES(?,?)", (slug.lower(), vcv))
    for (g,) in c.execute("SELECT DISTINCT gene FROM variant").fetchall():
        n = c.execute("SELECT COUNT(*) FROM variant WHERE gene=?", (g,)).fetchone()[0]
        c.execute("INSERT OR REPLACE INTO gene_meta VALUES(?,?)", (g, n))
    c.commit()
    c.close()
    return str(path)


def _gate(tmp_path, old, new, fn):
    a = _db(tmp_path / "cur.db", old)
    b = _db(tmp_path / "new.db", new)
    c = PS._open(a, b)
    return fn(c, _Quiet())


class _Quiet:
    def line(self, *a):
        pass

    def table(self, *a, **k):
        pass


def _gene(name, n, cls="Uncertain significance", start=0):
    return [(start + i, name, f"{name.lower()}-v{i}", cls) for i in range(n)]


# --- the carrier bug: 36,374 variants gone while the TOTAL went UP ---------------
def test_lost_slug_catches_loss_hidden_by_growth(tmp_path):
    """TTN lost 5,899 pages to an unenumerated TTN-AS1 while a newer ClinVar more than
    made up the count. Totals rose; pages 404'd."""
    old = _gene("TTN", 40)
    new = _gene("TTN", 20) + _gene("NEWGENE", 60, start=500)
    assert len(new) > len(old)                      # totals say everything is fine
    ok, why = _gate(tmp_path, old, new, PS.gate_no_lost_slug)
    assert not ok and "404" in why


def test_not_shrunk_is_fooled_by_the_same_case(tmp_path):
    """Why the slug gate has to exist: the count-level gate PASSES this."""
    old = _gene("TTN", 40)
    new = _gene("TTN", 20) + _gene("NEWGENE", 60, start=500)
    ok, _ = _gate(tmp_path, old, new, PS.gate_not_shrunk)
    assert ok


def test_lost_slug_passes_when_nothing_is_lost(tmp_path):
    old = _gene("TTN", 20)
    ok, _ = _gate(tmp_path, old, old + _gene("NEW", 5, start=900), PS.gate_no_lost_slug)
    assert ok


def test_lost_slug_accepts_a_slug_kept_as_an_alias(tmp_path):
    """A renamed canonical slug is fine if the old URL still resolves."""
    old = [(1, "TP53", "tp53-old-slug", "Pathogenic")]
    a = _db(tmp_path / "cur.db", old)
    b = _db(tmp_path / "new.db", [(1, "TP53", "tp53-new-slug", "Pathogenic")])
    conn = sqlite3.connect(b)
    conn.execute("INSERT INTO alias VALUES(?,?)", ("tp53-old-slug", 1))
    conn.commit()
    conn.close()
    c = PS._open(a, b)
    ok, _ = PS.gate_no_lost_slug(c, _Quiet())
    assert ok


# --- the FLG bug: wrong gene resolved, 2,147 -> 201 -----------------------------
def test_gene_collapse_catches_the_flg_shape(tmp_path):
    old = _gene("FLG", 40)
    new = _gene("FLG", 4) + _gene("FGFR1", 80, start=500)
    ok, why = _gate(tmp_path, old, new, PS.gate_no_gene_collapse)
    assert not ok and "FLG" in why


def test_gene_collapse_allows_a_plausible_shrink(tmp_path):
    ok, _ = _gate(tmp_path, _gene("X", 40), _gene("X", 35), PS.gate_no_gene_collapse)
    assert ok


def test_gene_collapse_ignores_small_genes(tmp_path):
    """Under COLLAPSE_MIN the ratio is noise."""
    ok, _ = _gate(tmp_path, _gene("TINY", 6), _gene("TINY", 1), PS.gate_no_gene_collapse)
    assert ok


# --- clinically actionable pages --------------------------------------------------
def test_clinical_loss_catches_pathogenic_drop_under_growing_totals(tmp_path):
    old = _gene("A", 10, "Pathogenic") + _gene("B", 10, "Likely pathogenic", start=100)
    new = _gene("A", 2, "Pathogenic") + _gene("C", 90, "Benign", start=300)
    assert len(new) > len(old)
    ok, why = _gate(tmp_path, old, new, PS.gate_no_clinical_loss)
    assert not ok and "fell" in why


def test_clinical_loss_passes_when_pl_p_grows(tmp_path):
    old = _gene("A", 10, "Pathogenic")
    ok, _ = _gate(tmp_path, old, _gene("A", 12, "Pathogenic"), PS.gate_no_clinical_loss)
    assert ok


# --- a build that "succeeded" with genes missing -----------------------------------
def test_failed_genes_gate(tmp_path):
    a = _db(tmp_path / "cur.db", _gene("A", 5))
    b = _db(tmp_path / "new.db", _gene("A", 5))
    conn = sqlite3.connect(b)
    conn.execute("INSERT OR REPLACE INTO gene_meta VALUES('BROKEN', -1)")
    conn.commit()
    conn.close()
    ok, why = PS.gate_no_failed_genes(PS._open(a, b), _Quiet())
    assert not ok and "failed" in why


# --- a row is not a page ----------------------------------------------------------
def test_resolves_catches_a_row_with_no_alias(tmp_path):
    a = _db(tmp_path / "cur.db", _gene("A", 5))
    b = _db(tmp_path / "new.db", _gene("A", 5))
    conn = sqlite3.connect(b)
    conn.execute("DELETE FROM alias")               # rows present, nothing resolvable
    conn.commit()
    conn.close()
    ok, why = PS.gate_resolves(PS._open(a, b), _Quiet(), sample=5)
    assert not ok and "no alias" in why


def test_resolves_passes_on_a_sane_index(tmp_path):
    a = _db(tmp_path / "cur.db", _gene("A", 20))
    b = _db(tmp_path / "new.db", _gene("A", 20))
    ok, _ = PS.gate_resolves(PS._open(a, b), _Quiet(), sample=10)
    assert ok


def test_every_gate_is_registered():
    """A gate that is written but not in GATES never runs."""
    registered = {fn for _, fn in PS.GATES}
    defined = {v for k, v in vars(PS).items() if k.startswith("gate_") and callable(v)}
    assert defined == registered, f"unregistered: {defined - registered}"


# --- a gene RENAME is not a collapse ---------------------------------------------
def test_gene_collapse_accepts_a_rename_with_aliases(tmp_path):
    """HGNC renamed AATK to LMTK1, so ClinVar's HGVS says LMTK1 and all 314 rows move
    there. With the aatk-* keys still in the alias table, no URL breaks."""
    old = [(i, "AATK", f"aatk-v{i}", "Uncertain significance") for i in range(20)]
    new = [(i, "LMTK1", f"lmtk1-v{i}", "Uncertain significance") for i in range(20)]
    a = _db(tmp_path / "cur.db", old)
    b = _db(tmp_path / "new.db", new)
    conn = sqlite3.connect(b)
    for i in range(20):                               # old URLs kept resolvable
        conn.execute("INSERT OR IGNORE INTO alias VALUES(?,?)", (f"aatk-v{i}", i))
    conn.commit()
    conn.close()
    c = PS._open(a, b)
    ok, _ = PS.gate_no_gene_collapse(c, _Quiet())
    assert ok
    assert PS._renamed_to(c, "AATK") == ("LMTK1", 20)


def test_gene_collapse_rejects_a_rename_that_breaks_urls(tmp_path):
    """Same move, but WITHOUT the old slugs aliased -- every old URL 404s, so it is a
    loss however it is labelled."""
    old = [(i, "AATK", f"aatk-v{i}", "Uncertain significance") for i in range(20)]
    new = [(i, "LMTK1", f"lmtk1-v{i}", "Uncertain significance") for i in range(20)]
    c = PS._open(_db(tmp_path / "cur.db", old), _db(tmp_path / "new.db", new))
    assert PS._renamed_to(c, "AATK") is None
    ok, why = PS.gate_no_gene_collapse(c, _Quiet())
    assert not ok and "AATK" in why


def test_rename_detection_needs_a_single_destination(tmp_path):
    """Rows scattered across several genes is not a rename -- do not excuse it."""
    old = [(i, "X", f"x-v{i}", "Uncertain significance") for i in range(20)]
    new = ([(i, "A", f"a-v{i}", "Uncertain significance") for i in range(10)] +
           [(i, "B", f"b-v{i}", "Uncertain significance") for i in range(10, 20)])
    c = PS._open(_db(tmp_path / "cur.db", old), _db(tmp_path / "new.db", new))
    assert PS._renamed_to(c, "X") is None
