"""Oracle for the Theorem B convention derivation (notes/theorem_B_convention.md).

Verifies the exact relationship: for a self-injective elementary algebra the genuine
Theta-action on HH_0 is the bare Nakayama permutation P_nu with char poly
det(t*I - P_nu), while the classical -C^{-T}C is det(t*I + P_nu); the two differ by
Phi = -P_nu (eigenvalues negated). Settles PLAN.md item E and fixes the orientation
the AR-bimodule action (item A) must use.
"""
import numpy as np
import sympy as sp

from quiverlab.engine.coxeter2 import coxeter_polynomial_from_cartan, cartan_from_raw

# hanlab __init__ aliases, reproduced locally:
coxeter_polynomial = coxeter_polynomial_from_cartan
cartan_matrix = cartan_from_raw

t = sp.symbols('t')


def _perm_charpoly(P):
    return sp.expand((t * sp.eye(P.shape[0]) - sp.Matrix(P.tolist())).det())


def test_kZ3_radsq_relationship():
    # Cartan matrix of kZ_3/rad^2 (vertices + one arrow out of each, rad^2 = 0)
    C = np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]], dtype=np.int64)
    Cs = sp.Matrix(C.tolist())

    # C^{-T} C is the Nakayama permutation P_nu (the 3-cycle) -- signature of self-injectivity
    P3 = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    assert (Cs.inv().T) * Cs == P3

    # classical Coxeter polynomial det(t*I - Phi), Phi = -C^{-T}C  ->  t^3 + 1
    poly, Phi = coxeter_polynomial(C)
    assert sp.expand(poly - (t ** 3 + 1)) == 0
    assert Phi == -P3                                   # Phi = -P_nu

    # genuine Theta on HH_0 = P_nu  ->  t^3 - 1 (the homologically meaningful one)
    assert sp.expand(_perm_charpoly(np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])) - (t ** 3 - 1)) == 0


def test_kZ2_radsq_classical_undefined_but_nakayama_defined():
    # kZ_2/rad^2: C = [[1,1],[1,1]], det C = 0 -> classical -C^{-T}C is undefined,
    # but the Nakayama permutation (transposition) gives a well-defined t^2 - 1.
    C = np.array([[1, 1], [1, 1]], dtype=np.int64)
    poly, Phi = coxeter_polynomial(C)
    assert poly is None and Phi is None
    assert sp.expand(_perm_charpoly(np.array([[0, 1], [1, 0]])) - (t ** 2 - 1)) == 0


def _nakayama_cartan(n, ell):
    """Cartan matrix of kZ_n/rad^ell via its paths (i -> i+s, s = 0..ell-1)."""
    paths = [(i, (i + s) % n) for i in range(n) for s in range(ell)]
    return cartan_matrix(n, None, paths)


def _rot(n, k):
    P = sp.zeros(n, n)
    for j in range(n):
        P[(j + k) % n, j] = 1
    return P


def test_nakayama_cartan_panel_CinvT_C_equals_P_nu():
    """C^{-T} C = P_nu (rotation by ell-1) across a panel of kZ_n/rad^ell with
    det C != 0 -- the ν-symmetry of the Cartan matrix (closes the note's [verify])."""
    checked = 0
    for n in (2, 3, 4, 5):
        for ell in (2, 3):
            C = _nakayama_cartan(n, ell)
            Cs = sp.Matrix(C.tolist())
            if Cs.det() == 0:
                continue
            assert (Cs.inv().T) * Cs == _rot(n, (ell - 1) % n), f"kZ_{n}/rad^{ell}"
            checked += 1
    assert checked >= 5


def test_symmetric_self_injective_nakayama_is_identity():
    """kZ_2/rad^3 is symmetric (Nakayama permutation = rotation by ell-1 = 2 ≡ 0):
    C^{-T} C = I, so Phi = -I (classical (t+1)^n) while genuine Theta = id ((t-1)^n)."""
    n, ell = 2, 3
    C = _nakayama_cartan(n, ell)
    Cs = sp.Matrix(C.tolist())
    assert Cs.det() != 0
    assert (Cs.inv().T) * Cs == sp.eye(n)          # nu = id (symmetric)
    poly, Phi = coxeter_polynomial(C)
    assert Phi == -sp.eye(n)
    assert sp.expand(poly - (t + 1) ** n) == 0
    assert sp.expand(_perm_charpoly(np.eye(n, dtype=np.int64)) - (t - 1) ** n) == 0


def test_symmetric_algebra_sign_relationship():
    # For symmetric algebras C = C^T (invertible) => C^{-T}C = I => Phi = -I, so the
    # classical poly is (t+1)^n while the genuine Theta (nu = id) gives (t-1)^n:
    # they are NOT equal -- they differ by exactly the sign Phi = -P_nu.
    C = np.array([[2, 1, 0], [1, 2, 1], [0, 1, 2]], dtype=np.int64)
    assert np.array_equal(C, C.T)
    poly, Phi = coxeter_polynomial(C)
    assert sp.expand(poly - (t + 1) ** 3) == 0
    assert Phi == -sp.eye(3)
