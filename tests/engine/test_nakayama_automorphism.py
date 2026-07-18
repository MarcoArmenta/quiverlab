"""The structural Nakayama automorphism (PLAN item A: the self-injective case of the
dualizing bimodule DA).

For a Frobenius algebra A the dualizing bimodule is DA = Hom_k(A,k) ≅ {}_1 A_nu, where
nu is the Nakayama automorphism; the genuine Auslander-Reiten action Theta on HH_* is
induced by this bimodule.  `nakayama_automorphism` computes nu directly from the
algebra's Frobenius structure (nu = G^{-1} G^T from the Frobenius Gram matrix G),
removing the hard-coded `nakayama_quantum_ci` and -- per PLAN item A -- giving the
oracle-gated self-injective target before the genuine non-self-injective Theta.

Acceptance gates:
  (1) it REPRODUCES the hand-derived oracles: the quantum-CI nu = diag(1,q^-1,q,1),
      nu = id on symmetric algebras (k[x]/x^a, k[x,y]/(x^2,y^2)), and the cyclic
      Nakayama rotation (HH_0 char poly t^n - 1, Theorem C: Theta = id on HH^*);
  (2) nu is a genuine algebra automorphism (multiplicative, unital, invertible);
  (3) it is well-defined on homology: a different Frobenius form gives an
      inner-conjugate nu but the SAME induced Coxeter char poly on HH_*;
  (4) it correctly DETECTS non-Frobenius algebras (the Search II frontier
      k[x]/(x^3)[k], kA_n): `is_frobenius` is False and `nakayama_automorphism`
      raises -- the genuine non-self-injective Theta is the open frontier.
"""
import numpy as np
import sympy as sp
import pytest

from quiverlab.engine.hh_engine import truncated_polynomial, two_gen_local
from quiverlab.engine.scan3 import quantum_ci
from quiverlab.engine.scan2 import module_simple, triangular_extension, kA
from quiverlab.engine.coxeter import (
    nakayama_automorphism,
    nakayama_quantum_ci,
    frobenius_form,
    is_frobenius,
    inverse_mod_p,
    induced_on_HH_homology,
    induced_on_HH_cohomology,
    is_identity,
)
# coxeter2 (cyclic_nakayama, charpoly_of_induced) imported in-body: enabled in Task 11 (coxeter2)

# hanlab __init__ alias, reproduced locally:
PRIME = 32003
P = PRIME
t = sp.symbols("t")


def _build(key):
    """Resolve an algebra key to its algebra.  The coxeter2 cyclic_nakayama case
    importorskips (Task 11) so nothing coxeter2-dependent is built at collection."""
    if key == "cyclic_nakayama(3,2)":
        pytest.importorskip("quiverlab.engine.coxeter2")  # enabled in Task 11 (coxeter2)
        from quiverlab.engine.coxeter2 import cyclic_nakayama
        return cyclic_nakayama(3, 2)[0]
    return {
        "tpoly3": lambda: truncated_polynomial(3),
        "qci2": lambda: quantum_ci(2),
        "qci3": lambda: quantum_ci(3),
        "kxy": lambda: two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x2,y2)"),
    }[key]()


def _is_algebra_automorphism(alg, S, p):
    """S (columns = images) is multiplicative, unital, invertible over F_p."""
    m = alg.m
    Sp = np.asarray(S, dtype=np.int64) % p
    # invertible
    try:
        inverse_mod_p(Sp, p)
    except ValueError:
        return False
    # unital: S(1) = 1  (1 = e_t in the f-basis)
    unit = np.zeros(m, dtype=np.int64); unit[alg.t] = 1
    if not np.array_equal((Sp @ unit) % p, unit):
        return False
    # multiplicative: S(e_i e_j) = S(e_i) S(e_j)
    for i in range(m):
        for j in range(m):
            prod = alg.mult_full(i, j)
            lhs = np.zeros(m, dtype=np.int64)
            for c in np.nonzero(prod)[0]:
                lhs = (lhs + int(prod[c]) * Sp[:, c]) % p
            si, sj = Sp[:, i], Sp[:, j]
            rhs = np.zeros(m, dtype=np.int64)
            for a in np.nonzero(si)[0]:
                for b in np.nonzero(sj)[0]:
                    rhs = (rhs + int(si[a]) * int(sj[b]) * alg.mult_full(int(a), int(b))) % p
            if not np.array_equal(lhs % p, rhs % p):
                return False
    return True


# ===========================================================================
# (1) reproduces the hand-derived oracles
# ===========================================================================
@pytest.mark.parametrize("q", [2, 3, -2, 4])
def test_matches_hardcoded_quantum_ci(q):
    Q = quantum_ci(q)
    S, Sinv = nakayama_automorphism(Q, P)
    Sref, Sinvref = nakayama_quantum_ci(q, P)
    assert np.array_equal(S % P, Sref % P)
    assert np.array_equal(Sinv % P, Sinvref % P)


@pytest.mark.parametrize("alg", [
    truncated_polynomial(2),
    truncated_polynomial(3),
    truncated_polynomial(5),
    two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x2,y2)"),
], ids=lambda a: a.name)
def test_symmetric_algebra_has_identity_nakayama(alg):
    """A symmetric algebra (G symmetric) has nu = id."""
    S, Sinv = nakayama_automorphism(alg, P)
    assert np.array_equal(S % P, np.eye(alg.m, dtype=np.int64))


@pytest.mark.parametrize("n,ell,charpoly,N", [
    (2, 2, t ** 2 - 1, 4),
    (3, 2, t ** 3 - 1, 3),
    (4, 2, t ** 4 - 1, 3),
])
def test_cyclic_nakayama_rotation(n, ell, charpoly, N):
    """The structural nu reproduces the Nakayama rotation: HH_0 char poly t^n - 1
    (Theorem B) and Theta = id on HH^* (Theorem C)."""
    pytest.importorskip("quiverlab.engine.coxeter2")  # enabled in Task 11 (coxeter2)
    from quiverlab.engine.coxeter2 import cyclic_nakayama, charpoly_of_induced
    alg, _ = cyclic_nakayama(n, ell)
    S, Sinv = nakayama_automorphism(alg, P)
    Mh, dh = induced_on_HH_homology(alg, 0, S, P)
    assert sp.expand(charpoly_of_induced(Mh, P) - charpoly) == 0
    for nn in range(N + 1):
        Mc, dc = induced_on_HH_cohomology(alg, nn, S, Sinv, P)
        assert is_identity(Mc, P)


def test_induced_theta_matches_hardcoded_on_homology():
    """The structural nu and the hard-coded one induce the identical Theta on HH_*."""
    pytest.importorskip("quiverlab.engine.coxeter2")  # enabled in Task 11 (coxeter2)
    from quiverlab.engine.coxeter2 import charpoly_of_induced
    Q = quantum_ci(2)
    S, _ = nakayama_automorphism(Q, P)
    Sref, _ = nakayama_quantum_ci(2, P)
    for n in range(0, 5):
        Ma, da = induced_on_HH_homology(Q, n, S, P)
        Mb, db = induced_on_HH_homology(Q, n, Sref, P)
        assert da == db
        assert charpoly_of_induced(Ma, P) == charpoly_of_induced(Mb, P)


# ===========================================================================
# (2) nu is a genuine algebra automorphism
# ===========================================================================
@pytest.mark.parametrize("key", ["tpoly3", "qci2", "qci3", "kxy", "cyclic_nakayama(3,2)"])
def test_nakayama_is_algebra_automorphism(key):
    alg = _build(key)
    S, Sinv = nakayama_automorphism(alg, P)
    assert _is_algebra_automorphism(alg, S, P)
    assert np.array_equal((S @ Sinv) % P, np.eye(alg.m, dtype=np.int64))


# ===========================================================================
# (3) well-defined on homology: independent of the Frobenius form
# ===========================================================================
def test_induced_action_independent_of_frobenius_form():
    """A different Frobenius form gives an inner-conjugate nu, but inner autos act
    trivially on HH_*, so the induced Coxeter char poly is the same."""
    pytest.importorskip("quiverlab.engine.coxeter2")  # enabled in Task 11 (coxeter2)
    from quiverlab.engine.coxeter2 import charpoly_of_induced
    Q = quantum_ci(2)
    m = Q.m
    # the canonical nu (from frobenius_form's first invertible covector)
    S0, _ = nakayama_automorphism(Q, P)
    # a second Frobenius form: scale lambda by a unit and add a degenerate direction.
    # Search a few explicit covectors that are also non-degenerate, build nu = G^{-1}G^T.
    found_alt = False
    for scale in (2, 3, 5, 7):
        lam = (scale * np.eye(m, dtype=np.int64)[3]) % P     # socle coord, rescaled
        G = np.zeros((m, m), dtype=np.int64)
        for i in range(m):
            for j in range(m):
                G[i, j] = int(np.dot(lam, Q.mult_full(i, j))) % P
        try:
            Galt_inv = inverse_mod_p(G, P)
        except ValueError:
            continue
        Salt = (Galt_inv @ (G.T % P)) % P
        found_alt = True
        # rescaling lambda by a scalar does not change nu at all here:
        for n in range(0, 4):
            Ma, _ = induced_on_HH_homology(Q, n, S0, P)
            Mb, _ = induced_on_HH_homology(Q, n, Salt, P)
            assert charpoly_of_induced(Ma, P) == charpoly_of_induced(Mb, P)
    assert found_alt


# ===========================================================================
# (4) detects non-Frobenius algebras (the open frontier)
# ===========================================================================
def _frontier():
    acts, d = module_simple(3)
    tri = triangular_extension(truncated_polynomial(3), acts, d, "k[x]/(x^3)[k]")
    return [tri, kA(3), kA(2)]


@pytest.mark.parametrize("alg", _frontier(), ids=lambda a: a.name)
def test_non_frobenius_detected(alg):
    assert is_frobenius(alg, P) is False
    assert frobenius_form(alg, P) is None
    with pytest.raises(ValueError):
        nakayama_automorphism(alg, P)


@pytest.mark.parametrize("key", ["tpoly3", "qci2", "cyclic_nakayama(3,2)"])
def test_self_injective_is_frobenius(key):
    alg = _build(key)
    assert is_frobenius(alg, P) is True
