"""Item 4 (TODO P2): regression-pin Theorem B on multi-vertex algebras.

Theorem B (the project's headline math): the AR/Coxeter-induced automorphism Theta acts on
HH_* with characteristic polynomial equal to the Coxeter polynomial -- and on a self-injective
Nakayama algebra Theta is the Nakayama rotation, so on HH_0 the genuine char poly is that of the
bare Nakayama permutation P_nu, namely t^n - 1 (det(t I - P_nu)). Theorem C says Theta = id on
HH^*. The `coxeter2 nakayama` driver computes both exactly but nothing asserted them; the pure
Cartan/permutation side is pinned in test_theoremB_convention.py, the INDUCED-ACTION side (through
the actual HH engine) is pinned here.

See notes/theorem_B_convention.md for why the homologically meaningful poly is det(t I - P_nu)
(= t^n - 1) while the classical -C^{-T}C gives det(t I + P_nu) (= t^n + 1): they differ by the
sign Phi = -P_nu.
"""
import numpy as np
import sympy as sp
import pytest
from quiverlab.engine.hh_engine import hochschild_homology_dims
from quiverlab.engine.coxeter2 import cyclic_nakayama, rotation_matrix_orig, auto_to_f_basis, charpoly_of_induced
from quiverlab.engine.coxeter import induced_on_HH_homology, induced_on_HH_cohomology, is_identity

# hanlab __init__ aliases, reproduced locally:
homology_dims = hochschild_homology_dims
PRIME = 32003

P = PRIME
t = sp.symbols('t')


def _nakayama_setup(n, ell):
    """Build kZ_n/rad^ell and the Nakayama rotation Theta (and its inverse) in the f-basis."""
    alg, idx = cyclic_nakayama(n, ell)
    unit = np.zeros(alg.m, dtype=np.int64)
    for i in range(n):
        unit[idx(i, 0)] = 1
    k = (ell - 1) % n
    S = auto_to_f_basis(alg, unit, rotation_matrix_orig(n, ell, k))
    Sinv = auto_to_f_basis(alg, unit, rotation_matrix_orig(n, ell, (-k) % n))
    return alg, S, Sinv


@pytest.mark.parametrize("n,ell,hh_dims,hh0_charpoly", [
    (2, 2, [2, 1, 1, 1, 1], t ** 2 - 1),   # kZ_2/rad^2: Nakayama transposition  -> t^2 - 1
    (3, 2, [3, 0, 1, 1],     t ** 3 - 1),   # kZ_3/rad^2: Nakayama 3-cycle        -> t^3 - 1
])
def test_theoremB_charpoly_on_HH0(n, ell, hh_dims, hh0_charpoly):
    alg, S, Sinv = _nakayama_setup(n, ell)
    N = len(hh_dims) - 1
    # homology dimensions are the frozen oracle
    assert homology_dims(alg, N)[P] == hh_dims
    # char poly of the induced Theta on HH_0 is t^n - 1 (the genuine Nakayama permutation)
    Mh, dh = induced_on_HH_homology(alg, 0, S, P)
    assert sp.expand(charpoly_of_induced(Mh, P) - hh0_charpoly) == 0


@pytest.mark.parametrize("n,ell,N", [(2, 2, 4), (3, 2, 3)])
def test_theoremC_theta_is_identity_on_cohomology(n, ell, N):
    # Theorem C: Theta acts as the identity on HH^n in every degree (invisible to cohomology).
    alg, S, Sinv = _nakayama_setup(n, ell)
    for nn in range(N + 1):
        Mc, dc = induced_on_HH_cohomology(alg, nn, S, Sinv, P)
        assert is_identity(Mc, P), f"Theta != id on HH^{nn} of kZ_{n}/rad^{ell}"


def test_theoremB_nontrivial_certificate_kZ3():
    # Asymmetric non-vanishing witness: on kZ_3/rad^2, HH_2 != 0 carries a NONTRIVIAL Theta
    # action (char poly t - 1 on a 1-dim space is trivial, so check a degree where it bites:
    # HH_0 has the nontrivial 3-cycle action while HH^0 = center is fixed).
    alg, S, Sinv = _nakayama_setup(3, 2)
    Mh, dh = induced_on_HH_homology(alg, 0, S, P)
    Mc, dc = induced_on_HH_cohomology(alg, 0, S, Sinv, P)
    assert dh == 3 and not is_identity(Mh, P)   # nontrivial Coxeter/Nakayama action on HH_0
    assert is_identity(Mc, P)                    # but identity on HH^0 (Thm C)
