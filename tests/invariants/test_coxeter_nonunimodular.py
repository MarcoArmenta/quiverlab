"""Characterization + regression tests for coxeter_polynomial's domain (spec section 5
component 8; ledger). The domain follows the COEFFICIENTS (sympy's inference), NOT det C:
a non-unimodular Cartan may be rational (QQ) OR integral (ZZ). These tests pin current
correct behavior so a naive det-based rule cannot creep back in. Fixture D + diag(1,2)."""
import sympy as sp
from quiverlab import Quiver, CC, linear_path_algebra, truncated_polynomial

t = sp.Symbol("t")


def _nonunimodular_rational():
    # C = [[2,1],[0,1]], det = 2; Phi is genuinely rational -> QQ
    return Quiver([1, 2], {"x": (1, 1), "a": (1, 2)}).algebra(
        relations=["x^2", "x*a"], field=CC)


def test_nonunimodular_cartan_rational_is_qq():
    A = _nonunimodular_rational()
    assert A.cartan_matrix() == [[2, 1], [0, 1]]
    p = A.coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t**2 + sp.Rational(3, 2) * t + 1)) == 0
    assert p.domain == sp.QQ                    # rational Coxeter transformation


def test_unimodular_cartan_is_zz():
    p = linear_path_algebra(2).coxeter_polynomial()   # C = [[1,1],[0,1]], det 1
    assert sp.expand(p.as_expr() - (t**2 + t + 1)) == 0
    assert p.domain == sp.ZZ


def test_nonunimodular_but_integral_stays_zz_local():
    # k[x]/(x^2): C = [[2]], det 2 (non-unimodular) but Phi = [-1] integral -> ZZ
    p = truncated_polynomial(2).coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t + 1)) == 0
    assert p.domain == sp.ZZ                    # NON-unimodular yet integral: NOT QQ


def test_nonunimodular_diag_stays_zz():
    # C = diag(1, 2): isolated vertex 1 + a loop y at 2 with y^2 = 0; det 2, Phi = -I,
    # Coxeter polynomial (t+1)^2 -- integral, so ZZ (a |det C| rule would wrongly give QQ)
    A = Quiver([1, 2], {"y": (2, 2)}).algebra(relations=["y^2"], field=CC)
    assert A.cartan_matrix() == [[1, 0], [0, 2]]
    p = A.coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t + 1)**2) == 0
    assert p.domain == sp.ZZ


def test_singular_cartan_still_raises():
    from quiverlab.errors import QuiverlabError
    import pytest
    # kZ_4/rad^2 has det C = 0 -> singular -> loud
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 1)}
               ).algebra(relations=["a*b", "b*c", "c*d", "d*a"], field=CC)
    with pytest.raises(QuiverlabError):
        A.coxeter_polynomial()
