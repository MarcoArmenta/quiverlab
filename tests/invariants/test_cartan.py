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


def test_coxeter_matrix_A2_exact_value():
    # Pin the public API return for kA_2: Phi = -C^{-T} C with C = [[1,1],[0,1]].
    # det C = 1, so entries are exact integers (plain Python ints, not sympy).
    A = linear_path_algebra(2)
    Phi = A.coxeter_matrix()
    assert Phi == [[-1, -1], [1, 0]]
    assert all(type(x) is int for row in Phi for x in row)


def test_nakayama_automorphism_dual_numbers_is_identity():
    # k[x]/(x^2) over GF(5) is symmetric: the Nakayama automorphism is the
    # identity matrix, returned as plain Python ints (no numpy scalar types).
    A = truncated_polynomial(2, field=GF(5))
    nu = A.nakayama_automorphism()
    assert nu == [[1, 0], [0, 1]]
    assert all(type(x) is int for row in nu for x in row)


def test_cartan_without_provenance_fails_loudly():
    # a bare structure-constant algebra has no path basis: the Cartan matrix must
    # refuse loudly with an actionable hint, never guess
    from quiverlab import Algebra
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(QuiverlabError):
        A.cartan_matrix()
