"""The Gerstenhaber bracket on HH^* (the Lie side of the TT calculus, PLAN item D).

Oracles (all asserted directly on basis / representative cochains, the project's
rigorous descent style):

  * [m, m] = 0 for the multiplication 2-cochain m  (<=> associativity);
  * graded antisymmetry  [f, g] = -(-1)^{(p-1)(q-1)} [g, f]  (chain level);
  * descent to HH^*: the bracket of two cocycles is a cocycle, and the bracket of
    a coboundary with a cocycle is a coboundary -- so [,] is well-defined on classes;
  * on HH^1 the bracket is the commutator of (outer) derivations  D_f D_g - D_g D_f;
  * graded Jacobi on HH^1 classes.
"""

import itertools

import numpy as np
import pytest

# hanlab sys.path shim dropped: quiverlab uses absolute package imports.
from quiverlab.engine.hh_engine import truncated_polynomial, two_gen_local
from quiverlab.engine.scan3 import cochain_basis, coboundary_matrix, quantum_ci
from quiverlab.engine.coxeter import nullspace_mod_p, colspace_basis_mod_p, solve_mod_p
from quiverlab.engine.tt_calculus import (
    circle_cochain,
    gerstenhaber_bracket_cochain,
    multiplication_cochain,
    gerstenhaber_bracket_matrix,
    cochain_to_dict,
    cup_cochain,
)

P = 32003


def _zoo():
    return [
        truncated_polynomial(2),
        truncated_polynomial(3),
        two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x2,y2)"),
        quantum_ci(3),
    ]


def _delta(alg, p, vec):
    bp = cochain_basis(alg, p)
    idx = {g: i for i, g in enumerate(cochain_basis(alg, p + 1))}
    return coboundary_matrix(alg, p, bp, idx) @ np.asarray(vec, dtype=np.int64)


def _cocycles(alg, p, prime):
    bp = cochain_basis(alg, p)
    idx = {g: i for i, g in enumerate(cochain_basis(alg, p + 1))}
    return nullspace_mod_p(coboundary_matrix(alg, p, bp, idx), prime)


def _coboundaries(alg, p, prime):
    if p == 0:
        return np.zeros((len(cochain_basis(alg, 0)), 0), dtype=np.int64)
    bpm = cochain_basis(alg, p - 1)
    idx = {g: i for i, g in enumerate(cochain_basis(alg, p))}
    return colspace_basis_mod_p(coboundary_matrix(alg, p - 1, bpm, idx), prime)


def _rand(alg, p, rng):
    return rng.integers(-2, 3, size=len(cochain_basis(alg, p)))


# ---------------------------------------------------------------------------
# [m, m] = 0  <=>  associativity
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", _zoo(), ids=lambda a: a.name)
def test_bracket_of_multiplication_vanishes(alg):
    m = multiplication_cochain(alg)
    mm = gerstenhaber_bracket_cochain(alg, 2, 2, m, m)
    assert np.all(mm % P == 0), f"{alg.name}: [m, m] != 0"


# ---------------------------------------------------------------------------
# graded antisymmetry (chain level)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", _zoo(), ids=lambda a: a.name)
def test_graded_antisymmetry(alg):
    rng = np.random.default_rng(0)
    for (p, q) in [(1, 1), (1, 2), (2, 1), (2, 2), (1, 3), (3, 1)]:
        f = _rand(alg, p, rng)
        g = _rand(alg, q, rng)
        b1 = gerstenhaber_bracket_cochain(alg, p, q, f, g)
        b2 = gerstenhaber_bracket_cochain(alg, q, p, g, f)
        sign = 1 if ((p - 1) * (q - 1)) % 2 == 0 else -1
        # [f,g] = -(-1)^{(p-1)(q-1)} [g,f]  <=>  [f,g] + sign*[g,f] = 0
        assert np.all((b1 + sign * b2) % P == 0), f"{alg.name} (p,q)=({p},{q})"


# ---------------------------------------------------------------------------
# descent to cohomology
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", [truncated_polynomial(3), quantum_ci(3)],
                         ids=lambda a: a.name)
def test_bracket_descends_to_cohomology(alg):
    for (p, q) in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        Zp = _cocycles(alg, p, P)
        Zq = _cocycles(alg, q, P)
        Bp = _coboundaries(alg, p, P)
        Bpq = _coboundaries(alg, p + q - 1, P)
        # bracket of cocycles is a cocycle
        for i in range(Zp.shape[1]):
            for j in range(Zq.shape[1]):
                br = gerstenhaber_bracket_cochain(alg, p, q, Zp[:, i], Zq[:, j])
                assert np.all(_delta(alg, p + q - 1, br) % P == 0)
        # bracket of a coboundary with a cocycle is a coboundary
        for i in range(Bp.shape[1]):
            for j in range(Zq.shape[1]):
                br = gerstenhaber_bracket_cochain(alg, p, q, Bp[:, i], Zq[:, j]) % P
                if Bpq.shape[1] == 0:
                    assert np.all(br == 0)
                else:
                    assert solve_mod_p(Bpq, br, P) is not None


# ---------------------------------------------------------------------------
# on HH^1: bracket = commutator of derivations
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", _zoo(), ids=lambda a: a.name)
def test_bracket_on_C1_is_derivation_commutator(alg):
    def commutator(f, g):
        fd = cochain_to_dict(alg, 1, f)
        gd = cochain_to_dict(alg, 1, g)

        def apply(d, av):
            out = np.zeros(alg.m, dtype=np.int64)
            for c in np.nonzero(av)[0]:
                if c == alg.t:
                    continue
                val = d.get((int(c),))
                if val is not None:
                    out = out + int(av[c]) * val
            return out

        out = {}
        zero = np.zeros(alg.m, dtype=np.int64)
        for r in alg.R:
            c = apply(fd, gd.get((r,), zero)) - apply(gd, fd.get((r,), zero))
            if np.any(c):
                out[(r,)] = c
        from quiverlab.engine.tt_calculus import dict_to_cochain
        return dict_to_cochain(alg, 1, out)

    rng = np.random.default_rng(2)
    for _ in range(5):
        f = _rand(alg, 1, rng)
        g = _rand(alg, 1, rng)
        assert np.all((gerstenhaber_bracket_cochain(alg, 1, 1, f, g)
                       - commutator(f, g)) % P == 0)


# ---------------------------------------------------------------------------
# graded Jacobi on HH^1 classes
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", [truncated_polynomial(3), quantum_ci(3)],
                         ids=lambda a: a.name)
def test_graded_jacobi_on_HH1(alg):
    Z = _cocycles(alg, 1, P)
    B1 = _coboundaries(alg, 1, P)
    br = lambda x, y: gerstenhaber_bracket_cochain(alg, 1, 1, x, y)
    n = Z.shape[1]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                f, g, h = Z[:, i], Z[:, j], Z[:, k]
                s = (br(br(f, g), h) + br(br(g, h), f) + br(br(h, f), g)) % P
                if B1.shape[1] == 0:
                    assert np.all(s == 0)
                else:
                    assert solve_mod_p(B1, s, P) is not None


# ---------------------------------------------------------------------------
# graded (Koszul-signed) Jacobi in a mixed / higher degree
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", [truncated_polynomial(3), quantum_ci(3)],
                         ids=lambda a: a.name)
def test_graded_jacobi_mixed_degree(alg):
    """Graded Jacobi identity for the Gerstenhaber bracket in a MIXED degree (p,q,r)=(2,2,1)
    ≠ (1,1,1) (test_graded_jacobi_on_HH1 only covers (1,1,1), where every Koszul sign is +1).
    With |f| = p−1 the shifted degree, the identity is

        (−1)^{|f||h|}[[f,g],h] + (−1)^{|g||f|}[[g,h],f] + (−1)^{|h||g|}[[h,f],g] = 0,

    which for (2,2,1) reads  [[f,g],h] − [[g,h],f] + [[h,f],g] = 0  — one genuine −1 Koszul
    sign. The Gerstenhaber bracket is the graded commutator of a pre-Lie (comp) product, so
    this holds EXACTLY at the cochain level (stronger than mod-coboundaries). Non-vacuous:
    HH^1 and HH^2 are both nonzero here (so the cocycle reps are genuine classes), and the
    WRONG all-+ sign is asserted non-identically-zero, pinning the Koszul sign load-bearing."""
    p, q, r = 2, 2, 1
    a, b, c = p - 1, q - 1, r - 1
    s_fgh = 1 if (a * c) % 2 == 0 else -1        # sign on [[f,g],h]
    s_ghf = 1 if (b * a) % 2 == 0 else -1        # sign on [[g,h],f]  -> −1 for (2,2,1)
    s_hfg = 1 if (c * b) % 2 == 0 else -1        # sign on [[h,f],g]
    assert (s_fgh, s_ghf, s_hfg) == (1, -1, 1)   # the mixed-degree Koszul signs, spelled out
    Zp, Zq, Zr = _cocycles(alg, p, P), _cocycles(alg, q, P), _cocycles(alg, r, P)
    # non-vacuous: the input degrees carry nonzero cohomology classes (dim Z − dim B > 0)
    assert Zp.shape[1] - _coboundaries(alg, p, P).shape[1] > 0   # HH^2 ≠ 0
    assert Zr.shape[1] - _coboundaries(alg, r, P).shape[1] > 0   # HH^1 ≠ 0
    br = lambda pp, qq, x, y: gerstenhaber_bracket_cochain(alg, pp, qq, x, y)
    wrong_sign_seen_nonzero = False
    for i in range(Zp.shape[1]):
        for j in range(Zq.shape[1]):
            for k in range(Zr.shape[1]):
                f, g, h = Zp[:, i], Zq[:, j], Zr[:, k]
                t1 = br(p + q - 1, r, br(p, q, f, g), h)   # [[f,g],h]
                t2 = br(q + r - 1, p, br(q, r, g, h), f)   # [[g,h],f]
                t3 = br(r + p - 1, q, br(r, p, h, f), g)   # [[h,f],g]
                jac = (s_fgh * t1 + s_ghf * t2 + s_hfg * t3) % P
                assert np.all(jac == 0), f"{alg.name}: graded Jacobi (2,2,1) fails"
                if np.any((t1 + t2 + t3) % P):             # wrong all-+ sign combination
                    wrong_sign_seen_nonzero = True
    assert wrong_sign_seen_nonzero, "Koszul sign must be load-bearing (all-+ should not vanish)"


# ---------------------------------------------------------------------------
# Gerstenhaber algebra: the bracket is a graded derivation of the cup product
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", [truncated_polynomial(3), quantum_ci(3)],
                         ids=lambda a: a.name)
def test_poisson_identity_on_classes(alg):
    """[a, b ⌣ c] = [a,b] ⌣ c + (-1)^{(p-1)q} b ⌣ [a,c] on HH^* (the Gerstenhaber
    algebra / graded Poisson identity), checked on cocycle representatives mod
    coboundaries."""
    def in_span(B, v):
        v = v % P
        if B.shape[1] == 0:
            return not np.any(v)
        return solve_mod_p(B, v, P) is not None

    for (p, q, r) in [(1, 1, 1), (1, 1, 2), (2, 1, 1), (1, 2, 1)]:
        Za, Zb, Zc = _cocycles(alg, p, P), _cocycles(alg, q, P), _cocycles(alg, r, P)
        Bres = _coboundaries(alg, p + q + r - 1, P)
        sign = 1 if ((p - 1) * q) % 2 == 0 else -1
        for i in range(Za.shape[1]):
            for j in range(Zb.shape[1]):
                for k in range(Zc.shape[1]):
                    a, b, c = Za[:, i], Zb[:, j], Zc[:, k]
                    lhs = gerstenhaber_bracket_cochain(
                        alg, p, q + r, a, cup_cochain(alg, q, r, b, c))
                    t1 = cup_cochain(alg, p + q - 1, r,
                                     gerstenhaber_bracket_cochain(alg, p, q, a, b), c)
                    t2 = cup_cochain(alg, q, p + r - 1, b,
                                     gerstenhaber_bracket_cochain(alg, p, r, a, c))
                    assert in_span(Bres, lhs - (t1 + sign * t2))


# ---------------------------------------------------------------------------
# the induced bracket matrix on classes
# ---------------------------------------------------------------------------
def test_induced_bracket_matrix_antisymmetric_on_HH1():
    """[.,.] : HH^1 (x) HH^1 -> HH^1 is antisymmetric (a Lie algebra)."""
    A = truncated_polynomial(3)
    C, d1, d1b, dout = gerstenhaber_bracket_matrix(A, 1, 1, P)
    assert d1 == d1b
    # antisymmetry of the induced structure constants: C[:,i,j] = -C[:,j,i]
    for i in range(d1):
        for j in range(d1):
            assert np.all((C[:, i, j] + C[:, j, i]) % P == 0)
    # self-bracket is zero
    for i in range(d1):
        assert np.all(C[:, i, i] % P == 0)
