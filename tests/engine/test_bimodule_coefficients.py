"""Hochschild (co)homology with bimodule coefficients (PLAN item A step 2).

The genuine non-self-injective Auslander-Reiten action Theta is built from the
dualizing bimodule DA = Hom_k(A,k).  This is the substrate: HH_*(A, M) / HH^*(A, M)
for an arbitrary A-bimodule M, on the normalized bar complex.

Acceptance gates:
  (1) reproduction -- for the regular bimodule M = A, HH_*(A,A) and HH^*(A,A) equal
      the standard `hh_engine` / `scan3` engines exactly (all primes);
  (2) DA is a genuine bimodule (`check_bimodule`), and its unit/duality structure
      holds;
  (3) the two classical dualities:
          dim HH_n(A, DA) = dim HH^n(A)   and   dim HH^n(A, DA) = dim HH_n(A),
      checked across the zoo incl. the asymmetric quantum CI (where HH_*(A,A) and
      HH^*(A,A) differ, so the swap is a real test, not a tautology);
  (4) HH_0(A, DA) = HH^0(A) = Z(A) in dimension.
"""
import numpy as np
import pytest

from quiverlab.engine.hh_engine import (
    truncated_polynomial,
    two_gen_local,
    hochschild_homology_dims,
)
from quiverlab.engine.scan3 import quantum_ci, hochschild_cohomology_dims
from quiverlab.engine.scan2 import module_simple, triangular_extension
from quiverlab.engine.bimodule import (
    regular_bimodule,
    dual_bimodule,
    check_bimodule,
    hochschild_homology_with_coefficients,
    hochschild_cohomology_with_coefficients,
)

# hanlab __init__ aliases, reproduced locally:
homology_dims = hochschild_homology_dims        # dim HH_n(A)
cohomology_dims = hochschild_cohomology_dims    # dim HH^n(A)
PRIME = 32003
P = PRIME


# Panel kept to dim <= 4 so the degree-4 bar complexes (dim C^n = m*(m-1)^n) stay cheap;
# the multi-vertex case is cyclic_nakayama(2,2) (dim 4, two idempotents 1_A = e_0 + e_1).
# String keys (not built algebras) so the coxeter2-dependent multi-vertex case is NOT
# constructed at collection time; `_build` importorskips coxeter2 (Task 11) for it, so the
# single-idempotent cases still run now and the multi-vertex case lights up when it lands.
_PANEL_KEYS = [
    "k[x]/(x^2)", "k[x]/(x^3)", "k[x,y]/(x2,y2)", "qci2", "qci3",
    "cyclic_nakayama(2,2)",
]


def _build(key):
    """Resolve a panel key to its algebra."""
    if key == "cyclic_nakayama(2,2)":
        pytest.importorskip("quiverlab.engine.coxeter2")  # enabled in Task 11 (coxeter2)
        from quiverlab.engine.coxeter2 import cyclic_nakayama
        return cyclic_nakayama(2, 2)[0]
    return {
        "k[x]/(x^2)": lambda: truncated_polynomial(2),
        "k[x]/(x^3)": lambda: truncated_polynomial(3),
        "k[x,y]/(x2,y2)": lambda: two_gen_local(
            [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x2,y2)"),
        "qci2": lambda: quantum_ci(2),
        "qci3": lambda: quantum_ci(3),
    }[key]()


DEG = 3   # depth for the (co)homology cross-checks (bar complex blows up with n)


# ===========================================================================
# (1) reproduction: M = A gives the standard engines
# ===========================================================================
@pytest.mark.parametrize("key", _PANEL_KEYS)
def test_regular_bimodule_reproduces_homology(key):
    alg = _build(key)
    M = regular_bimodule(alg)
    got = hochschild_homology_with_coefficients(alg, M, DEG, primes=(P, 2, 3))
    ref = homology_dims(alg, DEG, primes=(P, 2, 3))
    assert got == ref


@pytest.mark.parametrize("key", _PANEL_KEYS)
def test_regular_bimodule_reproduces_cohomology(key):
    alg = _build(key)
    M = regular_bimodule(alg)
    got = hochschild_cohomology_with_coefficients(alg, M, DEG, primes=(P, 2, 3))
    ref = cohomology_dims(alg, DEG, primes=(P, 2, 3))
    assert got == ref


# ===========================================================================
# (2) DA is a genuine bimodule
# ===========================================================================
@pytest.mark.parametrize("key", _PANEL_KEYS)
def test_dual_bimodule_axioms(key):
    alg = _build(key)
    assert check_bimodule(dual_bimodule(alg), P)
    assert check_bimodule(regular_bimodule(alg), P)


def test_dual_bimodule_dimension():
    """dim DA = dim A."""
    for key in _PANEL_KEYS:
        alg = _build(key)
        assert dual_bimodule(alg).mu == alg.m


# ===========================================================================
# (3) the classical dualities  HH_n(A,DA) = HH^n(A),  HH^n(A,DA) = HH_n(A)
# ===========================================================================
@pytest.mark.parametrize("key", _PANEL_KEYS)
def test_homology_with_DA_equals_cohomology(key):
    alg = _build(key)
    DA = dual_bimodule(alg)
    hom_DA = hochschild_homology_with_coefficients(alg, DA, DEG, primes=(P,))[P]
    coh = cohomology_dims(alg, DEG, primes=(P,))[P]
    assert hom_DA == coh


@pytest.mark.parametrize("key", _PANEL_KEYS)
def test_cohomology_with_DA_equals_homology(key):
    alg = _build(key)
    DA = dual_bimodule(alg)
    coh_DA = hochschild_cohomology_with_coefficients(alg, DA, DEG, primes=(P,))[P]
    hom = homology_dims(alg, DEG, primes=(P,))[P]
    assert coh_DA == hom


def test_duality_is_nontrivial_on_quantum_ci():
    """On the quantum CI HH_*(A,A) != HH^*(A,A), so DA genuinely swaps them:
    HH_*(A,DA) tracks cohomology, HH^*(A,DA) tracks homology."""
    Q = quantum_ci(2)
    DA = dual_bimodule(Q)
    hom_A = homology_dims(Q, DEG, primes=(P,))[P]
    coh_A = cohomology_dims(Q, DEG, primes=(P,))[P]
    assert hom_A != coh_A                              # the asymmetry is real
    assert hochschild_homology_with_coefficients(Q, DA, DEG, primes=(P,))[P] == coh_A
    assert hochschild_cohomology_with_coefficients(Q, DA, DEG, primes=(P,))[P] == hom_A


# ===========================================================================
# (4) degree 0
# ===========================================================================
@pytest.mark.parametrize("key", _PANEL_KEYS)
def test_degree0_DA_homology_is_center_dimension(key):
    """dim HH_0(A, DA) = dim HH^0(A) = dim Z(A)."""
    alg = _build(key)
    DA = dual_bimodule(alg)
    h0 = hochschild_homology_with_coefficients(alg, DA, 0, primes=(P,))[P][0]
    z = cohomology_dims(alg, 0, primes=(P,))[P][0]
    assert h0 == z


# ===========================================================================
# non-self-injective frontier: the dualities still hold (DA always exists)
# ===========================================================================
def test_dualities_hold_on_non_frobenius_frontier():
    """DA = Hom_k(A,k) exists for every A; the dualities hold even where there is
    no Nakayama automorphism (the open frontier k[x]/(x^3)[k])."""
    pytest.importorskip("quiverlab.engine.coxeter")  # enabled in Task 6 (coxeter)
    from quiverlab.engine.coxeter import is_frobenius
    acts, d = module_simple(3)
    A = triangular_extension(truncated_polynomial(3), acts, d, "k[x]/(x^3)[k]")
    assert is_frobenius(A, P) is False
    DA = dual_bimodule(A)
    assert check_bimodule(DA, P)
    hom_DA = hochschild_homology_with_coefficients(A, DA, 3, primes=(P,))[P]
    coh = cohomology_dims(A, 3, primes=(P,))[P]
    assert hom_DA == coh
