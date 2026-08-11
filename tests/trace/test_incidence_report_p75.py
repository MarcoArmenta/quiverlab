"""Plan 75 Task G1 (touchpoint 6) -- the ``incidence_cohomology`` section of the
worked-steps report.

Self-certifying pins: the section STATES the theorem it used before showing numbers, the
face vector it prints is the one the engine returned (never recomputed here, so the
picture can never claim a count the engine did not produce), the contractibility verdict
is rendered three-valued (proved / disproved / merely observed -- never collapsed to a
boolean), the characteristic remark appears ONLY in positive characteristic where it says
something, and a poset-less algebra prints the loud refusal instead of numbers.
"""
import pytest

import quiverlab as ql
from quiverlab.fields import GF, QQ
from quiverlab.hochschild.simplicial import incidence_cohomology_block
from quiverlab.trace.results_html import _HEADINGS, _block_html
from quiverlab.viz.tikz import tikz_order_complex

pytestmark = pytest.mark.oracle_selfcert

# The crown C(2,2): Delta ~ S^1, so H^1 != 0 and there is NO global bound.
_CROWN = [[1, 3], [1, 4], [2, 3], [2, 4]]
# The Boolean lattice B_3: a global bound, so Delta is a cone.
_B3 = [[x, x | (1 << b)] for x in range(8) for b in range(3) if not (x >> b) & 1]
# The minimal triangulation of RP^2 on 6 vertices, as the FACE poset of its simplices:
# H_1(RP^2; Z) = Z/2, so HH^* over GF(2) differs from HH^* over QQ. Built by
# quiverlab's own test helpers is overkill here -- the char remark is exercised on B_3
# over GF(2), where the verdict is "no torsion divisible by 2".


def _html(covers, top=3, field=None, elements=None):
    A = ql.IncidenceAlgebra(covers, elements, field=field)
    b = incidence_cohomology_block(A, top)
    return b, "\n".join(_block_html("incidence_cohomology", b))


def test_the_section_states_the_theorem_before_the_numbers():
    b, html = _html(_B3)
    assert "order complex" in html and "Gerstenhaber" in html and "Cibils" in html
    assert html.index("order complex") < html.index("dim HH^n")
    assert _HEADINGS["incidence_cohomology"].startswith("Hochschild cohomology")


def test_the_dims_table_carries_the_engine_dims():
    b, html = _html(_B3)
    assert b["dims"] == [1, 0, 0, 0]
    for n, d in enumerate(b["dims"]):
        assert "<td>%d</td>" % n in html
    assert "<th>dim HH^n</th>" in html


def test_the_face_vector_printed_is_the_engine_face_vector():
    b, html = _html(_B3)
    assert b["face_vector"] == [8, 19, 18, 6]
    assert "(8, 19, 18, 6)" in html
    # ... and the Euler characteristic is the alternating sum of THAT vector.
    assert b["euler_characteristic"] == sum(
        (-1) ** p * f for p, f in enumerate(b["face_vector"]))
    assert "= %d" % b["euler_characteristic"] in html


def test_contractibility_is_rendered_three_valued():
    # PROVED: B_3 has a global bound, so Delta is a cone.
    b, html = _html(_B3)
    assert b["contractible"] is True and "contractible</b> (proved)" in html
    # DISPROVED: the crown's H^1 is nonzero.
    b, html = _html(_CROWN)
    assert b["contractible"] is False and "not contractible" in html
    # NEITHER: no global bound, and nothing nonzero in the computed range. The section
    # must say "not proved", never claim contractibility from an observation.
    #   1 < 2, 3 < 2, 3 < 4 (a fence): contractible, but not by a global bound.
    b, html = _html([[1, 2], [3, 2], [3, 4]])
    assert b["dims"] == [1, 0, 0, 0]
    assert b["contractible"] is None
    assert "not proved" in html and "contractible</b> (proved)" not in html


def test_the_characteristic_remark_is_silent_in_characteristic_zero():
    """Over QQ, "these dimensions agree with the characteristic-0 ones" is vacuous --
    printing it would read as a computed fact."""
    b, html = _html(_B3, field=QQ)
    assert b["characteristic"] == 0
    assert "agree with the characteristic-0 ones" not in html
    # In positive characteristic the same verdict IS informative, and names the prime.
    b2, html2 = _html(_B3, field=GF(2))
    assert b2["characteristic"] == 2 and b2["char_dependent"] is False
    assert "divisible by 2" in html2


def test_a_poset_less_algebra_prints_the_refusal_not_numbers():
    A = ql.Quiver(vertices=[1, 2], arrows={"a": (1, 2)}).algebra()
    b = incidence_cohomology_block(A, 2)
    html = "\n".join(_block_html("incidence_cohomology", b))
    assert "not computed" in html and "provenance" in html
    assert "dim HH^n" not in html


def test_tikz_order_complex_states_only_what_it_was_given():
    A = ql.IncidenceAlgebra(_CROWN)
    b = incidence_cohomology_block(A, 3)
    tex = tikz_order_complex(A._poset, b["face_vector"])
    assert tex.startswith(r"\begin{tikzpicture}") and tex.rstrip().endswith(
        r"\end{tikzpicture}")
    # one node per element, one edge per cover -- the Hasse diagram, not the complex
    assert tex.count(r"\node[draw, circle]") == len(A._poset.elements)
    assert tex.count(r"\draw (p") == len(A._poset.covers)
    assert "face vector $(%s)$" % ", ".join(str(f) for f in b["face_vector"]) in tex
    # ... and with no face vector supplied it claims no count at all.
    assert "face vector" not in tikz_order_complex(A._poset)
