import sympy
import pytest
from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError


def test_cartan_kA2():
    A = linear_path_algebra(2)
    # C[i][j] = dim e_i A e_j = number of basis paths v_i -> v_j
    assert A.cartan_matrix() == [[1, 1], [0, 1]]


def test_cartan_field_independent():
    for field in (CC, GF(2), GF(5)):
        Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
        A = Q.algebra(relations=["a*b"], field=field)
        assert A.cartan_matrix() == [[1, 1, 0], [0, 1, 1], [0, 0, 1]]


def test_cartan_dual_numbers():
    assert truncated_polynomial(2).cartan_matrix() == [[2]]


def test_coxeter_polynomial_A2():
    # Phi = -C^{-T} C for kA_2 has characteristic polynomial t^2 + t + 1 (Coxeter number 3)
    A = linear_path_algebra(2)
    t = sympy.Symbol("t")
    assert sympy.expand(A.coxeter_polynomial().as_expr() - (t**2 + t + 1)) == 0


def test_cartan_without_provenance_fails_loudly():
    # a bare structure-constant algebra has no path basis: the Cartan matrix must
    # refuse loudly with an actionable hint, never guess
    from quiverlab import Algebra
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(QuiverlabError):
        A.cartan_matrix()
