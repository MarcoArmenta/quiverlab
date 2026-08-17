# SPDX-License-Identifier: MIT
"""N-Koszul recognizer (Berger) + the Ext-generation cross-check (Plan 77 / R36).
k[x]/(x^N), N>=3 -> N-Koszul (pure resolution, delta(n)); the quadratic case defers
to Plan 27; the preprojectives are quadratic-but-not-Koszul (Plan 27's False)."""
import pytest

from quiverlab import (PreprojectiveAlgebra, Quiver, TruncatedPathAlgebra,
                       linear_path_algebra, truncated_polynomial)
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import berger_degree, n_koszul_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("N", [3, 4, 5])
def test_kxN_is_N_koszul(N):
    c = n_koszul_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["n_homogeneous"] == N
    # k[x]/x^N is self-injective (gl.dim infinite): the Berger pattern is certified
    # THROUGH THE WINDOW, so the three-valued verdict is None, never a bare True.
    assert c["verdict"] is None and c["complete"] is False
    assert c["pure"] is True
    assert c["berger_expected"][:4] == [0, 1, N, N + 1]
    assert c["internal_degrees"][1][:4] == [0, 1, N, N + 1]
    assert "window" in c["reason"]


@lit
def test_cubic_monomial_finite_gldim_is_3koszul():
    c = n_koszul_certificate(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert c["n_homogeneous"] == 3
    assert c["verdict"] is True and c["complete"] is True
    assert c["internal_degrees"][1] == [0, 1, 3]      # Berger delta for N=3


@xeng
def test_quadratic_defers_to_plan27():
    # kA2 (hereditary) is quadratic-trivially -> N=2, verdict == ext_algebra.koszul.
    A = linear_path_algebra(2, field=QQ)
    c = n_koszul_certificate(A, top=6)
    assert c["n_homogeneous"] == 2
    assert c["verdict"] == A.ext_algebra(6).koszul
    assert "Plan 27" in c["reason"]


@lit
@pytest.mark.parametrize("typ", ["A3", "A4"])
def test_preprojective_is_not_N_koszul(typ):
    # Pi(A_{n>=3}) is quadratic (N=2) but NOT Koszul: defers to Plan 27's False.
    A = PreprojectiveAlgebra(typ, field=QQ)
    c = n_koszul_certificate(A, top=6)
    assert c["n_homogeneous"] == 2 and c["verdict"] is False
    assert A.ext_algebra(6).koszul is False


@xeng
@pytest.mark.parametrize("N", [3, 4, 5])
def test_two_certificates_agree(N):
    # oracle_crossengine: Berger's N>=3 theorem says the PURITY of the minimal
    # resolution (internal degrees, read off the Plan-05 resolution) and K2 (the
    # Yoneda generator degrees, read off the Plan-27 Ext-algebra engine) decide the
    # same property. Two independent readings, one verdict.
    c = n_koszul_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["k2_agrees"] is True
    assert c["pure"] is True


@selfcert
def test_berger_closed_form_alternates_one_and_N_minus_one():
    # The 2-N alternation: consecutive differences are 1, N-1, 1, N-1, ...
    for N in (2, 3, 4, 5, 7):
        deltas = [berger_degree(n, N) for n in range(9)]
        jumps = [b - a for a, b in zip(deltas, deltas[1:])]
        assert jumps == [1, N - 1] * 4


@selfcert
def test_non_n_homogeneous_reports_none_and_points_elsewhere():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x - x*x"], field=QQ)
    c = n_koszul_certificate(A, top=6)
    assert c["n_homogeneous"] is None and c["verdict"] is None
    assert "K2" in c["reason"]


@selfcert
def test_mixed_relation_degrees_are_not_n_homogeneous():
    # Homogeneous relations, but in TWO different degrees (2 and 3): N-homogeneity
    # asks for ONE degree, so Berger does not apply and the recognizer says so.
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y*x"], field=QQ)
    c = n_koszul_certificate(A, top=4)
    assert c["n_homogeneous"] is None and c["verdict"] is None
