# SPDX-License-Identifier: MIT
"""(p,q)-almost-Koszul (Brenner-Butler-King): p = top degree, q = e - p (first internal
jump minus p). Reproduces BBK's (h-2, 2) on the Dynkin preprojectives -- exactly the
algebras a Koszul route refuses (Plan 77 / R36)."""
import pytest

from quiverlab import (PreprojectiveAlgebra, TruncatedPathAlgebra,
                       linear_path_algebra, truncated_polynomial)
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import _top_degree, almost_koszul_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("typ,p,h", [("A3", 2, 4), ("A4", 3, 5), ("A5", 4, 6),
                                     ("D4", 4, 6)])
def test_preprojective_is_p_2_almost_koszul(typ, p, h):
    c = almost_koszul_certificate(PreprojectiveAlgebra(typ, field=QQ), top=6)
    assert c["verdict"] is True
    assert c["p"] == p == h - 2                   # BBK: p = h-2 (the top degree)
    assert c["q"] == 2                            # BBK: q = 2
    assert c["break_hom_degree"] == 3             # the first non-linear jump
    assert c["break_internal_degree"] == p + 2    # e = p + q


@lit
@pytest.mark.parametrize("typ,h", [("A2", 3), ("A3", 4), ("A4", 5), ("A5", 6),
                                   ("A6", 7), ("D4", 6)])
def test_top_degree_is_h_minus_2_on_dynkin_preprojectives(typ, h):
    # BBK's p for the Dynkin preprojectives is h-2; the engine reads it off the
    # algebra's own basis (no Groebner completion, so Pi(A5)/Pi(D4) do not raise).
    assert _top_degree(PreprojectiveAlgebra(typ, field=QQ)) == h - 2


@lit
def test_preprojective_A2_is_koszul_not_almost():
    # Pi(A2) = kZ2/rad^2 is genuinely Koszul (linear forever): no break -> None.
    # (BBK's nominal (1,2) label is subsumed by genuine Koszulity.)
    A = PreprojectiveAlgebra("A2", field=QQ)
    c = almost_koszul_certificate(A, top=8)
    assert c["verdict"] is None
    assert c["break_hom_degree"] is None and c["q"] is None
    assert "linear" in c["reason"]
    assert A.ext_algebra(8).koszul is True        # the Koszul field reports it Koszul


@lit
@pytest.mark.parametrize("N", [3, 4, 5])
def test_kxN_is_not_almost_koszul_q_equals_1(N):
    # THE negative pin (the q>=2 gate). k[x]/x^N has top degree p = N-1 and its first
    # jump at e = N, so q = e - p = 1 -> NOT almost-Koszul (it is N-Koszul, caught by
    # n_koszul_certificate). q = 1 is BBK's Koszul-type degenerate boundary.
    c = almost_koszul_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["p"] == N - 1 and c["q"] == 1
    assert c["break_hom_degree"] == 2 and c["break_internal_degree"] == N
    assert c["verdict"] is None


@lit
def test_kx2_is_koszul_so_it_has_NO_jump_and_no_q():
    # PLAN-DOC CORRECTION (measured): the plan parametrized the q==1 negative pin over
    # N = 2,3,4 and asserted q == 1 for EVERY N. For N = 2 that is wrong -- k[x]/x^2 is
    # radical-square-zero, hence Koszul, and its resolution is linear FOREVER
    # (l(n) = n, verified through degree 8), so there is no first jump at all and
    # q = e - p is undefined. It gets the same treatment as Pi(A2): verdict None with
    # q None, never a fabricated 1.
    A = truncated_polynomial(2, field=QQ)
    c = almost_koszul_certificate(A, top=8)
    assert c["p"] == 1
    assert c["break_hom_degree"] is None and c["q"] is None
    assert c["verdict"] is None
    assert A.ext_algebra(8).koszul is True


@selfcert
@pytest.mark.parametrize("typ", ["A3", "A4"])
def test_p75_seam_obstruction_equals_break(typ):
    # THE P75 seam: the ext_algebra.koszul obstruction degree that P75's GHMS route
    # REFUSES on equals the almost-Koszul break degree that P77 CLASSIFIES. This is a
    # STRUCTURAL-CONSISTENCY invariant, not two independent computations -- both read
    # the same first-failure-of-linearity degree, one through the Yoneda engine ("a new
    # generator appears") and one through the resolution ("the first jump above the
    # diagonal").
    c = almost_koszul_certificate(PreprojectiveAlgebra(typ, field=QQ), top=6)
    assert c["seam_obstruction_degree"] == c["break_hom_degree"] == 3


@selfcert
@pytest.mark.parametrize("typ", ["A3", "A4", "A5", "D4"])
def test_the_two_readings_of_q_agree(typ):
    # BBK's definition puts the error of internal degree p+q at homological step q+1,
    # so q = e - p and q = n* - 1 are two readings of the SAME q. verdict True is
    # gated on their agreement; here it is asserted directly.
    c = almost_koszul_certificate(PreprojectiveAlgebra(typ, field=QQ), top=6)
    assert c["q"] == c["linear_steps"] == c["break_hom_degree"] - 1


@selfcert
def test_jump_spacing_is_window_observed_and_is_not_bbk_period():
    # Pi(A3): h = 4, so BBK's periodicity is 2(h-1) = 6, while the OBSERVED spacing of
    # the internal-degree jumps is 3. The engine records the observed spacing and the
    # docstring/reason refuse to call it BBK's period -- this test pins that they are
    # genuinely different numbers, so a future rename cannot quietly conflate them.
    c = almost_koszul_certificate(PreprojectiveAlgebra("A3", field=QQ), top=6)
    assert c["jump_degrees"] == [3, 6]
    assert c["jump_spacing"] == 3
    assert c["jump_spacing"] != 2 * (4 - 1)
    assert "period" not in c.keys()


@selfcert
def test_finite_gldim_monomial_breaks_once_with_q_one():
    # kA4/J^3: only S_1 has a break (the others resolve in <= 1 step and stay linear),
    # so the "e agrees across simples" rule must range over the simples that BREAK --
    # a terminated linear simple imposes no constraint.
    c = almost_koszul_certificate(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert c["p"] == 2 and c["break_hom_degree"] == 2
    assert c["break_internal_degree"] == 3 and c["q"] == 1
    assert c["verdict"] is None                   # q = 1: 3-Koszul, not almost-Koszul


@selfcert
def test_hereditary_has_no_break():
    c = almost_koszul_certificate(linear_path_algebra(3, field=QQ), top=5)
    assert c["verdict"] is None and c["break_hom_degree"] is None
    assert c["seam_obstruction_degree"] is None   # kA3 is Koszul: no obstruction
