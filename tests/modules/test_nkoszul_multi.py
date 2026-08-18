# SPDX-License-Identifier: MIT
"""Multi-Koszul (Herscovich): a CONNECTED-graded (A_0 = k) notion; multi-vertex kQ/I
refuses and reports K2 (the settled transfer, Prop. 3.30). Plan 77 / R36."""
import pytest

from quiverlab import PreprojectiveAlgebra, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import multi_koszul_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ),
                               PreprojectiveAlgebra("A3", field=QQ)])
def test_multivertex_refuses_and_points_at_K2(A):
    c = multi_koszul_certificate(A, top=6)
    assert c["applicable"] is False and c["verdict"] is None
    assert "K2" in c["reason"] and "connected" in c["reason"]
    assert c["k2"]["verdict"] is not None or c["k2"]["generator_degrees"] is not None


@lit
def test_local_kxN_is_connected_and_reports_generation_and_K2():
    A = truncated_polynomial(3, field=QQ)         # single vertex -> connected graded
    c = multi_koszul_certificate(A, top=8)
    assert c["applicable"] is True
    assert c["generation_degrees"][1][:4] == [0, 1, 3, 4]     # Berger N=3
    assert c["k2"]["generator_degrees"] == [1, 2]
    # verdict is None: the multi-Koszul DECISION (Herscovich Section 3.2 bimodule
    # Tor/Ext) is scoped out; K2 (Prop. 3.30) + generation degrees are what v1 reports.
    assert c["verdict"] is None
    assert "3.2" in c["status"] and "scope" in c["status"].lower()
    assert "3.30" in c["status"]


@selfcert
def test_two_relation_degrees_local_is_still_connected():
    # The genuinely multi-Koszul-shaped input: ONE vertex, relations in TWO degrees
    # (2 and 3), so Berger's N-homogeneity does not apply but Herscovich's setting
    # does. The recognizer still reports no verdict (the decision is scoped out) and
    # still reports the transferable K2 data -- never a guess.
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y*x"], field=QQ)
    c = multi_koszul_certificate(A, top=5)
    assert c["applicable"] is True and c["verdict"] is None
    assert c["generation_degrees"][1][:2] == [0, 1]
    assert c["k2"]["verdict"] in (True, False, None)


@selfcert
def test_verdict_is_never_fabricated():
    # Across every shipped shape the multi-Koszul verdict is None: the decision engine
    # is not built, and no branch may invent one.
    for A in (truncated_polynomial(2, field=QQ), truncated_polynomial(4, field=QQ),
              linear_path_algebra(2, field=QQ), PreprojectiveAlgebra("A2", field=QQ)):
        assert multi_koszul_certificate(A, top=5)["verdict"] is None
