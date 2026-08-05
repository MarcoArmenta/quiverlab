"""Worked examples of Armenta, 'The Coxeter transformation as an
automorphism of the Tamarkin-Tsygan calculus' (arXiv:2606.15595), pinned
on the existing Coxeter surface. D4 vs A4: isomorphic TT calculi,
distinguished by the Coxeter polynomial; the 8-vertex cospectral trees:
equal Coxeter polynomials, not derived equivalent (honest-scope demo:
these invariants are necessary conditions only)."""
import sympy as sp
import pytest

from quiverlab import Quiver

pytestmark = pytest.mark.oracle_literature

t = sp.Symbol("t")


def test_d4_coxeter_polynomial():
    # D4: three arrows into the central vertex 1 (paper Ex 5.7)
    A = Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra()
    p = A.coxeter_polynomial().as_expr()
    assert sp.expand(p - (t**4 + t**3 + t + 1)) == 0
    assert sp.factor(p) == sp.factor((t + 1) ** 2 * (t**2 - t + 1))


def test_a4_coxeter_polynomial():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra()
    p = A.coxeter_polynomial().as_expr()
    assert sp.expand(p - (t**4 + t**3 + t**2 + t + 1)) == 0   # Phi_5


def test_d4_a4_differ_by_t_squared():
    D4 = Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra()
    A4 = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra()
    diff = sp.expand(A4.coxeter_polynomial().as_expr()
                     - D4.coxeter_polynomial().as_expr())
    assert diff == t**2                                        # paper Ex 5.7


def _spider():
    # 8 vertices: center 1 with arms of lengths 3,1,1,1,1 (degree seq 5,2,2,1x5)
    return Quiver(list(range(1, 9)),
                  {"a1": (1, 2), "a2": (2, 3), "a3": (3, 4),
                   "b": (1, 5), "c": (1, 6), "d": (1, 7), "e": (1, 8)}).algebra()


def _double_star():
    # 8 vertices: centers 1,2 joined; three leaves on each (degree seq 4,4,1x6)
    return Quiver(list(range(1, 9)),
                  {"m": (1, 2), "a": (1, 3), "b": (1, 4), "c": (1, 5),
                   "d": (2, 6), "e": (2, 7), "f": (2, 8)}).algebra()


def test_cospectral_trees_share_coxeter_polynomial():
    ps = _spider().coxeter_polynomial().as_expr()
    pd = _double_star().coxeter_polynomial().as_expr()
    target = sp.expand((t + 1) ** 4 * (t**4 - 3 * t**3 + t**2 - 3 * t + 1))
    assert sp.expand(ps - target) == 0 and sp.expand(pd - target) == 0
