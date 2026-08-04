"""Report fixes from Marco's 2026-08-03 pass over the desktop-app report
(k[x,y]/(x^2, y^2) over GF(7), every compute kind in one session):

  * a differential that IS the zero matrix printed a full zero grid (his report:
    every Tor differential d_1..d_6) -- a zero map is STATED, ``d_n = 0``, never
    drawn, and never echoed as "same matrix as above";
  * the Ext/Tor sections gave results without saying WHICH module was resolved
    and BY WHICH resolution -- the objects now precede the numbers (the products
    precedent);
  * engine provenance lines named "hanlab engine" without saying what that is --
    every engine note now carries a one-line gloss;
  * the "Worked resolution steps" chapter never said the resolution is of A as
    an A^e-module, nor which resolution was used (bar / Chouhy-Solotar).
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.trace.events import RankStep
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_html import _MatrixEcho, matrix_grid, render_html
from quiverlab.trace.results_html import results_section

pytestmark = [pytest.mark.oracle_selfcert]


def _dual_numbers(p=7):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=["x*x", "y*y", "x*y - y*x"], field=GF(p))


def _loop_x3(p=5):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x*x*x"], field=GF(p))


# --------------------------------------------------------------------------- #
# 1. A zero matrix is stated, never drawn.
# --------------------------------------------------------------------------- #
def test_zero_matrix_grid_prints_zero_not_grid():
    html = matrix_grid([[0, 0, 0], [0, 0, 0]], label="d_{1}")
    assert "ql-matrix" not in html
    assert "d_{1} = 0" in html


def test_zero_matrix_grid_string_entries():
    html = matrix_grid([["0", "0"], ["0", "0"]], label="d_{2}")
    assert "ql-matrix" not in html
    assert "d_{2} = 0" in html


def test_zero_matrix_grid_without_label():
    html = matrix_grid([[0], [0]])
    assert "ql-matrix" not in html


def test_nonzero_matrix_still_renders_grid():
    html = matrix_grid([[0, 1], [0, 0]], label="d_{1}")
    assert "ql-matrix" in html


def test_resolution_differentials_zero_named_never_echoed():
    # d0 nonzero; d1, d3 the SAME zero matrix (the echo trap); d2 zero of another
    # shape; d4 repeats d0 (the echo must still work for NON-zero matrices).
    z12 = [[0, 0]]
    block = {"kind": "projective_resolution", "top": 4,
             "terms": ["(1)"] * 5, "betti": [1] * 5, "summands": ["P_{1}"] * 5,
             "differentials": [
                 {"matrix": [[1, 0]], "rows": 1, "cols": 2},
                 {"matrix": z12, "rows": 1, "cols": 2},
                 {"matrix": [[0], [0]], "rows": 2, "cols": 1},
                 {"matrix": z12, "rows": 1, "cols": 2},
                 {"matrix": [[1, 0]], "rows": 1, "cols": 2},
             ]}
    html = "".join(results_section({"projective_resolution": block}))
    assert "d_{1} = 0" in html
    assert "d_{2} = 0" in html
    assert "d_{3} = 0" in html
    # the zero maps are never cross-referenced...
    assert "d_{3} = d_{1}" not in html
    # ...but a repeated NON-zero matrix still is:
    assert "d_{4} = d_{0}" in html
    # exactly one grid is drawn (d_0; d_4 echoes it)
    assert html.count('<table class="ql-matrix"') == 1


def test_matrix_echo_skips_zero_matrices():
    echo = _MatrixEcho()
    z = RankStep(degree=1, side="chain", nrows=2, ncols=2, rank=0,
                 field="GF(5)", matrix=[["0", "0"], ["0", "0"]])
    assert echo.label_for(z, "b_{1}") is None
    assert echo.label_for(z, "b_{2}") is None          # never recorded, never echoed


def test_tor_reps_zero_differentials_print_no_grid():
    # Over the dual numbers, the minimal resolution of any module tensored with a
    # SIMPLE N has identically-zero differentials: dim Tor_n == dim (P_n (x) N).
    # Marco's report drew every one of them as a zero grid.
    from quiverlab.modules.tor import tor_dims
    A = _dual_numbers()
    M = A.simple(1)
    N = A.simple(1, side="left")
    dims, reps = tor_dims(A, M, N, 3, with_reps=True)
    block = {"kind": "tor", "top": 3, "dims": [int(d) for d in dims]}
    block.update(reps)
    html = "".join(results_section({"tor": block}))
    assert '<table class="ql-matrix"' not in html
    assert "= 0" in html


# --------------------------------------------------------------------------- #
# 2. Ext/Tor name the resolved module + the resolution BEFORE the results.
# --------------------------------------------------------------------------- #
def test_ext_and_tor_blocks_carry_resolved_provenance():
    from quiverlab.hpc.spec import ComputeItem, _dispatch_module
    A = _loop_x3()
    M = A.simple(1)
    N = A.simple(1)
    T = A.simple(1, side="left")
    ext = _dispatch_module(A, ComputeItem("ext", 0, 2), M, N, T)
    tor = _dispatch_module(A, ComputeItem("tor", 0, 2), M, N, T)
    for block in (ext, tor):
        r = block["resolved"]
        assert r["module"] == "M"
        assert r["side"] == "right"
        assert "minimal projective" in r["resolution"]


def test_ext_section_states_resolved_module_before_results():
    block = {"kind": "ext", "top": 1, "dims": [1, 1],
             "target": {"dimvec": {"1": 1}},
             "resolved": {"module": "M", "side": "right",
                          "resolution": "minimal projective resolution"}}
    html = "".join(results_section({"ext": block}))
    stmt = html.find("minimal projective resolution")
    table = html.find("ql-dims")
    assert stmt != -1 and table != -1 and stmt < table
    assert "Hom" in html                       # Ext^n = H^n Hom_A(P_bullet, N)


def test_tor_section_states_resolved_module_before_results():
    block = {"kind": "tor", "top": 1, "dims": [2, 1],
             "target": {"dimvec": {"1": 1}},
             "resolved": {"module": "M", "side": "right",
                          "resolution": "minimal projective resolution"}}
    html = "".join(results_section({"tor": block}))
    stmt = html.find("minimal projective resolution")
    table = html.find("ql-dims")
    assert stmt != -1 and table != -1 and stmt < table
    assert "otimes" in html or "&#8855;" in html   # Tor_n = H_n(P_bullet (x)_A N)


# --------------------------------------------------------------------------- #
# 3. Engine provenance lines explain themselves.
# --------------------------------------------------------------------------- #
def test_engine_note_glosses_hanlab_in_cyclic_homology():
    block = {"kind": "HC_", "top": 1, "dims": [2, 1],
             "engine": "hanlab engine (F_p fast rank)"}
    html = "".join(results_section({"cyclic_homology": block}))
    # public wording only (the bank the engine was ported from is not public):
    assert "HansConjecture" not in html
    assert "Gaussian elimination mod p" in html


def test_engine_note_glosses_hanlab_in_connes_b():
    block = {"kind": "connes_b", "top": 1, "hh_dims": [2, 2],
             "matrices": {"0": [[1, 0], [0, 0]]}, "ranks": {"0": 1},
             "engine": "engine (b,B) GF(7)"}
    html = "".join(results_section({"connes_b": block}))
    assert "(b, B)" in html or "mixed complex" in html


def test_engine_note_glosses_chouhy_solotar_in_hh():
    block = {"kind": "HH^", "top": 1, "dims": [2, 2],
             "engine": "Chouhy-Solotar"}
    html = "".join(results_section({"hh_cohomology": block}))
    # the engine NOTE itself carries the gloss (not just the typing paragraph)
    assert "engine: Chouhy-Solotar &#8212;" in html or "engine: Chouhy-Solotar —" in html


# --------------------------------------------------------------------------- #
# 4. Worked resolution steps: A as an A^e-module + the resolution's name.
# --------------------------------------------------------------------------- #
def test_worked_resolution_steps_name_bimodule_and_chouhy_solotar():
    A = _loop_x3()
    tr = Trace()
    A.hochschild_homology(2, engine="cs", trace=tr)
    html = render_html(list(tr), title="t", algebra=A)
    i = html.find("id='resolution-steps'")
    assert i != -1
    intro = html[i:i + 1600]
    assert "A<sup>e</sup>" in intro
    assert "Chouhy" in intro


def test_worked_resolution_steps_name_bimodule_and_bar():
    A = _loop_x3()
    tr = Trace()
    A.hochschild_homology(2, engine="bar", trace=tr)
    html = render_html(list(tr), title="t", algebra=A)
    i = html.find("id='resolution-steps'")
    assert i != -1
    intro = html[i:i + 1600]
    assert "A<sup>e</sup>" in intro
    assert "bar resolution" in intro
