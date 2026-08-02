"""The products chapter: events render, dims drift-gate, full record (Plan 35).

`quiverlab.trace.products.products_chapter(A, kind, obj)` turns a Task-1 product
result object into the typed trace-event stream the HTML/JSON renderers turn into a
homework-grade chapter: a prose definition (StepNote), the authoritative HH
dimensions (ResultDims, which also DRIFT-GATES the narrated table dims against a
fresh Hochschild (co)homology), and one ProductStep per bidegree -- the nonzero
structure-constant equations (cup/cap/bracket) or the induced Connes matrix
(connes_b). `spec.run` backs the worked-steps bundle with the first product kind
when no HH/module trace claimed it.
"""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_selfcert]


def test_chapter_events_render_and_carry_tables(tmp_path):
    # Marco 2026-08-01: the whole family renders as ONE big degree-major CAYLEY table --
    # rows/cols over every cohomology class, each cell the product in the target basis,
    # zeros SHOWN, a beyond-window cell an em dash (not computed, not zero).
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    from tests.trace._matrix_grid import cayley_cells, cayley_headers
    A = ql.truncated_polynomial(2, field=ql.GF(2))       # k[x]/(x^2), the canonical case
    events = products_chapter(A, "cup", A.cup_products(2))
    assert events, "chapter must not be empty"
    html = render_html(events, title="cup", algebra=A)
    assert "cup product" in html.lower() and "HH" in html
    assert "ql-cayley" in html
    tables = cayley_cells(html)
    assert len(tables) == 1                               # ONE big table for the family
    big = tables[0]
    # degree-major axis: alpha^0_1, alpha^0_2, alpha^1_1, alpha^1_2, alpha^2_1, alpha^2_2
    cols, rows = cayley_headers(html)
    assert rows == [r"\alpha^{0}_{1}", r"\alpha^{0}_{2}", r"\alpha^{1}_{1}",
                    r"\alpha^{1}_{2}", r"\alpha^{2}_{1}", r"\alpha^{2}_{2}"]
    # a known nonzero product (alpha^1_1 cup alpha^1_1 = alpha^2_1), a SHOWN computed
    # zero (alpha^1_2 cup alpha^1_2 = 0), and a BEYOND-WINDOW em dash (target 3 > 2).
    assert big[2][2] == r"\alpha^{2}_{1}"
    assert big[3][3] == "0"
    assert big[2][4] == "—"
    assert "beyond the computed window" in html


def test_chapter_carries_notation_legend_before_the_tables():
    # Plan-35 follow-up (Marco): a legend defining alpha/beta/gamma/z/w and naming
    # the concrete recorded basis appears BEFORE the structure-constant tables.
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    hp = A.cup_products(2)
    html = render_html(products_chapter(A, "cup", hp), title="cup", algebra=A)
    assert "recorded basis" in html
    assert hp.basis in html                       # the concrete basis, e.g. bar/GF(7)
    # the legend precedes the first table's map label (\otimes) and equations.
    assert html.index("recorded basis") < html.index("otimes")


def test_drift_gate_fires_on_dim_mismatch():
    from quiverlab.trace.products import products_chapter
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    hp = A.cup_products(1)
    # sabotage a dim -- the builder must refuse to narrate inconsistent data
    hp.tables[(0, 0)] = hp.tables[(0, 0)].__class__(
        kind="cup", degrees=(0, 0), out_degree=0, dims=(99, 1, 1),
        constants=hp.tables[(0, 0)].constants)
    with pytest.raises(Exception):
        products_chapter(A, "cup", hp)


def test_cap_and_bracket_chapters_render():
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    for kind, obj in (("cap", A.cap_products(2)),
                      ("bracket", A.gerstenhaber_brackets(2))):
        events = products_chapter(A, kind, obj)
        assert events
        html = render_html(events, title=kind, algebra=A)
        assert "HH" in html
        # rendered as Cayley grids (unless the whole family vanishes to one line)
        assert "ql-cayley" in html or "vanish" in html.lower()


def test_bracket_shows_balanced_representatives():
    # GF(7) bracket constants include 5, which displays as the balanced rep -2 (the
    # JSON keeps the raw residue "5"). The balanced-rep legend is stated once.
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    from tests.trace._matrix_grid import cayley_cells
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    obj = A.gerstenhaber_brackets(2)
    # a raw 5 exists in the recorded constants
    assert any("5" in str(c) for t in obj.tables.values()
               for mat in t.constants for row in mat for c in row)
    html = render_html(products_chapter(A, "bracket", obj), title="bracket", algebra=A)
    assert "balanced representatives mod 7" in html
    cells = [c for tbl in cayley_cells(html) for row in tbl for c in row]
    joined = " ".join(cells)
    assert "-2" in joined                       # 5 shown balanced as -2
    assert r"5 \," not in joined                 # never the raw residue in a cell


def test_connes_chapter_carries_the_matrix():
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    from tests.trace._matrix_grid import grids
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    cb = A.connes_differentials(2)
    events = products_chapter(A, "connes_b", cb)
    assert events
    html = render_html(events, title="B", algebra=A)
    # the induced B matrices are rendered as indexed grids
    got = grids(html)
    assert got, "the Connes differential matrices must render as grids"
    # B_0 = [['0', '1']] must be among them
    assert [["0", "1"]] in got


def test_connes_drift_gate_fires_on_dim_mismatch():
    from quiverlab.trace.products import products_chapter
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    cb = A.connes_differentials(2)
    cb.hh_dims[0] = 99            # sabotage the recorded HH_0 dimension
    with pytest.raises(Exception):
        products_chapter(A, "connes_b", cb)


def test_products_back_the_worked_steps_bundle(tmp_path):
    """A products request with artifacts.pdf and no HH/module trace: the first
    product kind backs the worked-steps HTML+JSON bundle (spec.run selection)."""
    from quiverlab.hpc.spec import run as spec_run, parse_request
    req = parse_request({
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": [1],
                    "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                    "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["cup:0..2"],
        "artifacts": {"pdf": True, "tikz": False}})
    res = spec_run(req, tmp_path)
    assert res["meta"]["pdf"] == "worked steps in trace_steps.html"
    html = (tmp_path / "trace_steps.html").read_text(encoding="utf-8")
    assert "cup product" in html.lower()
    assert (tmp_path / "trace.json").exists()
