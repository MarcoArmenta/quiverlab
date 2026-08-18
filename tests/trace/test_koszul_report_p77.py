"""Plan 77 touchpoint 6 -- the `koszul` section of the worked-steps report.

Self-certifying pins: the section prints the ENGINE's generation degrees (never
recomputed here, so the table cannot claim a degree the engine did not produce), a
window-bounded verdict is LABELLED as such rather than typeset as proved, the quadratic
half is named as Plan 27's verdict, a refusal renders as prose instead of numbers, and
the TikZ staircase draws only numbers taken verbatim off the block.
"""
import pytest

from quiverlab import PreprojectiveAlgebra, TruncatedPathAlgebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import koszul_profile_block
from quiverlab.trace.results_html import _HEADINGS, _block_html
from quiverlab.viz.tikz import tikz_koszul

pytestmark = pytest.mark.oracle_selfcert


def _html(A, top=8):
    b = koszul_profile_block(A, top)
    return b, "\n".join(_block_html("koszul", b))


def test_the_heading_names_all_three_generalizations():
    h = _HEADINGS["koszul"]
    assert "N-Koszul" in h and "K" in h and "almost-Koszul" in h


def test_the_table_carries_the_ENGINE_generation_degrees():
    b, html = _html(truncated_polynomial(3, field=QQ))
    assert b["generation_degrees"]["1"] == [0, 1, 3, 4, 6, 7, 9, 10, 12]
    # every degree the engine produced must appear as a cell, in order
    cells = [f"<td>{d}</td>" for d in b["generation_degrees"]["1"]]
    pos = 0
    for cell in cells:
        found = html.find(cell, pos)
        assert found >= 0, f"{cell} missing from the rendered table"
        pos = found + 1
    assert "&#8467;(S<sub>1</sub>)" in html


def test_a_window_bounded_verdict_is_labelled_not_typeset_as_proved():
    # k[x]/x^3 is self-injective, so nothing here is unconditional.
    b, html = _html(truncated_polynomial(3, field=QQ))
    assert b["complete"] is False
    assert "not exact-finite" in html
    assert "certified THROUGH THAT DEGREE" in html
    assert "unconditional" in html


def test_a_complete_window_prints_no_window_caveat():
    b, html = _html(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert b["complete"] is True
    assert "not exact-finite" not in html
    assert "<b>yes</b>" in html                       # N-Koszul decided True


def test_the_quadratic_half_is_named_as_plan_27s_verdict():
    b, html = _html(truncated_polynomial(3, field=QQ))
    assert "Plan-27 verdict" in html
    assert "not Koszul" in html                       # k[x]/x^3 is not quadratic
    assert str(b["quadratic_obstruction"][0]) in html


def test_the_almost_koszul_label_appears_for_a_preprojective():
    b, html = _html(PreprojectiveAlgebra("A3", field=QQ), top=6)
    assert b["almost_koszul"]["verdict"] is True
    assert "(p, q) = (2, 2)" in html
    assert "Brenner" in html


def test_a_refusal_renders_as_prose_not_numbers():
    html = "\n".join(_block_html("koszul", {"kind": "koszul",
                                            "error": "needs the quiver presentation"}))
    assert "needs the quiver presentation" in html
    assert "<table" not in html


def test_the_tikz_staircase_prints_only_engine_numbers():
    b, _ = _html(PreprojectiveAlgebra("A3", field=QQ), top=6)
    tex = tikz_koszul(b)
    assert r"\begin{tikzpicture}" in tex and r"\end{tikzpicture}" in tex
    assert r"\ell = n" in tex                         # the diagonal is drawn
    row = b["generation_degrees"]["1"]
    for n, d in enumerate(row):
        assert "(%d,%d)" % (n, d) in tex
    assert "certified through degree %s" % b["certified_through_degree"] in tex


def test_the_tikz_twin_refuses_an_empty_block_instead_of_drawing_axes():
    tex = tikz_koszul({"generation_degrees": {}})
    assert "no generation-degree table recorded" in tex
    assert r"\draw" not in tex
