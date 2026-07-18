"""Non-linear Dynkin path algebras and their Coxeter polynomials (PLAN item C).

`quiver_path_algebra` builds the path algebra of any finite acyclic quiver; the
`dynkin_quiver` presets give the simply-laced A/D/E diagrams.  The oracle is the
classical theory of the Coxeter transformation Phi = -C^{-T}C:

  * the Cartan matrix of an acyclic-tree path algebra has det 1 (hereditary);
  * Phi has finite order equal to the Coxeter number h;
  * the Coxeter polynomial det(tI - Phi) is palindromic, of degree = #vertices, and
    equals the tabulated product of cyclotomic polynomials (orientation-independent);
  * the homology engine sees the hereditary signature HH_0 = #vertices, HH_n = 0
    for n >= 1 (gl.dim 1, tree quiver) on the small cases.
"""

import numpy as np
import sympy as sp
import pytest

from quiverlab.engine.coxeter2 import (
    quiver_path_algebra,
    dynkin_quiver,
    cartan_from_raw,
    coxeter_polynomial_from_cartan,
    coxeter_element_order,
)
from quiverlab.engine.hh_engine import hochschild_homology_dims, check_associative

P = 32003
t = sp.symbols("t")

# the classical (orientation-independent) Coxeter polynomials, as cyclotomic products
EXPECTED = {
    ("A", 4): (t**5 - 1) / (t - 1),
    ("D", 4): (t + 1) ** 2 * (t**2 - t + 1),
    ("D", 5): (t + 1) * (t**4 + 1),
    ("E", 6): (t**2 + t + 1) * (t**4 - t**2 + 1),
    ("E", 7): (t + 1) * (t**6 - t**3 + 1),
    ("E", 8): t**8 + t**7 - t**5 - t**4 - t**3 + t + 1,
}
COX_NUMBER = {("A", 4): 5, ("D", 4): 6, ("D", 5): 8,
              ("E", 6): 12, ("E", 7): 18, ("E", 8): 30}


def _cartan(typ, n):
    nv, arrows, h = dynkin_quiver(typ, n)
    alg, paths = quiver_path_algebra(nv, arrows, name=f"k{typ}{n}")
    se = [(p[0], p[-1]) for p in paths]
    return alg, cartan_from_raw(nv, None, se), h, nv


@pytest.mark.parametrize("typ,n", list(EXPECTED.keys()))
def test_cartan_determinant_is_one(typ, n):
    _, C, _, _ = _cartan(typ, n)
    assert int(sp.Matrix(C.tolist()).det()) == 1


@pytest.mark.parametrize("typ,n", list(EXPECTED.keys()))
def test_coxeter_element_order_is_coxeter_number(typ, n):
    _, C, h, _ = _cartan(typ, n)
    assert coxeter_element_order(C) == h == COX_NUMBER[(typ, n)]


@pytest.mark.parametrize("typ,n", list(EXPECTED.keys()))
def test_coxeter_polynomial_matches_classical(typ, n):
    _, C, _, nv = _cartan(typ, n)
    poly, _ = coxeter_polynomial_from_cartan(C)
    diff = sp.cancel(sp.expand(poly) - EXPECTED[(typ, n)])
    assert diff == 0, f"k{typ}{n}: {sp.factor(poly)} != {sp.factor(EXPECTED[(typ, n)])}"
    assert sp.degree(sp.expand(poly), t) == nv


@pytest.mark.parametrize("typ,n", list(EXPECTED.keys()))
def test_coxeter_polynomial_is_palindromic(typ, n):
    _, C, _, _ = _cartan(typ, n)
    poly, _ = coxeter_polynomial_from_cartan(C)
    coeffs = sp.Poly(sp.expand(poly), t).all_coeffs()
    assert coeffs == coeffs[::-1], f"k{typ}{n} Coxeter poly not self-reciprocal"


@pytest.mark.parametrize("typ,n", [("A", 4), ("D", 4)])
def test_hereditary_homology_signature(typ, n):
    """Small Dynkin path algebras: associative, HH_0 = #vertices, HH_{>=1} = 0."""
    alg, _, _, nv = _cartan(typ, n)
    assert check_associative(alg)[0]
    hh = hochschild_homology_dims(alg, 2, primes=(P,))[P]
    assert hh == [nv, 0, 0]


def test_quiver_path_algebra_generic_acyclic():
    """A non-Dynkin acyclic quiver (a 3-vertex 'V'): builds and is associative."""
    alg, paths = quiver_path_algebra(3, [(0, 1), (0, 2)], name="V")
    # paths: e_0, 0->1, 0->2, e_1, e_2  => dim 5
    assert alg.m == 5
    assert check_associative(alg)[0]
    hh0 = hochschild_homology_dims(alg, 0, primes=(P,))[P][0]
    assert hh0 == 3  # three vertices
