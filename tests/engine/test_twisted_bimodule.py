"""The twisted bimodule {}_phi A_psi and the self-injective dualizing bimodule
DA ≅ {}_1 A_ν  (PLAN item A step 2, the oracle-gated first target).

For a Frobenius algebra A with Nakayama automorphism ν, the dualizing bimodule is
`DA ≅ {}_1 A_ν`.  This is the self-injective case of the genuine non-self-injective
DA-twist Θ: here we verify, computationally, that the twisted bimodule built from ν
reproduces DA -- `HH_*(A, {}_1 A_ν) = HH_*(A, DA)` across the self-injective zoo --
and that the *wrong* twist conventions disagree (so the test bites), and that
`{}_phi A_psi` is a genuine bimodule for any pair of automorphisms.

The genuine non-self-injective Θ (where no Nakayama ν exists, e.g. k[x]/(x^3)[k])
remains the open headline of PLAN item A.

Nakayama/identity automorphisms live in coxeter (Task 6) and the multi-vertex
cyclic_nakayama zoo entries in coxeter2 (Task 11); the tests that need them
importorskip in-body so the coxeter-free identity-twist reproduction runs now.
"""

import numpy as np
import pytest

from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.scan3 import quantum_ci
from quiverlab.engine.bimodule import (
    twisted_bimodule,
    dual_bimodule,
    check_bimodule,
    hochschild_homology_with_coefficients,
)

P = 32003


# The self-injective zoo, as string keys so the coxeter2-dependent cyclic_nakayama
# cases are not built at collection time; `_zoo_build` importorskips coxeter2 (Task 11)
# for them.  Depths N match the bank panel.
_ZOO_KEYS = ["tpoly2", "qci3", "cyclic_nakayama(3,2)", "cyclic_nakayama(2,2)"]
_ZOO_DEPTH = {"tpoly2": 5, "qci3": 4, "cyclic_nakayama(3,2)": 3, "cyclic_nakayama(2,2)": 4}


def _zoo_build(key):
    """Resolve a self-injective-zoo key to its algebra."""
    if key.startswith("cyclic_nakayama"):
        pytest.importorskip("quiverlab.engine.coxeter2")  # enabled in Task 11 (coxeter2)
        from quiverlab.engine.coxeter2 import cyclic_nakayama
        n = 3 if key == "cyclic_nakayama(3,2)" else 2
        return cyclic_nakayama(n, 2)[0]
    return {"tpoly2": lambda: truncated_polynomial(2),
            "qci3": lambda: quantum_ci(3)}[key]()


@pytest.mark.parametrize("key", _ZOO_KEYS)
def test_DA_iso_twisted_by_nakayama(key):
    """DA ≅ {}_1 A_ν : HH_*(A, {}_1 A_ν) = HH_*(A, DA) for self-injective A."""
    pytest.importorskip("quiverlab.engine.coxeter")  # enabled in Task 6 (coxeter)
    from quiverlab.engine.coxeter import nakayama_automorphism
    alg = _zoo_build(key)
    N = _ZOO_DEPTH[key]
    S, _ = nakayama_automorphism(alg, P)
    A_nu = twisted_bimodule(alg, psi=S, name="_1A_nu")
    assert check_bimodule(A_nu, P)
    h_nu = hochschild_homology_with_coefficients(alg, A_nu, N, primes=(P,))[P]
    h_DA = hochschild_homology_with_coefficients(alg, dual_bimodule(alg), N, primes=(P,))[P]
    assert h_nu == h_DA


def test_wrong_twist_disagrees_for_nonsymmetric():
    """The convention bites: for a genuinely non-symmetric self-injective algebra
    (quantum CI), the wrong-side twist {}_ν A_1 does NOT reproduce DA."""
    pytest.importorskip("quiverlab.engine.coxeter")  # enabled in Task 6 (coxeter)
    from quiverlab.engine.coxeter import nakayama_automorphism
    alg = quantum_ci(3)
    S, _ = nakayama_automorphism(alg, P)
    h_DA = hochschild_homology_with_coefficients(alg, dual_bimodule(alg), 4, primes=(P,))[P]
    wrong = twisted_bimodule(alg, phi=S, name="_nu A_1")
    h_wrong = hochschild_homology_with_coefficients(alg, wrong, 4, primes=(P,))[P]
    assert h_wrong != h_DA


def test_symmetric_algebra_DA_is_regular():
    """A symmetric algebra has ν = id, so {}_1 A_ν = A and HH_*(A, DA) = HH_*(A)."""
    pytest.importorskip("quiverlab.engine.coxeter")  # enabled in Task 6 (coxeter)
    from quiverlab.engine.coxeter import nakayama_automorphism
    alg = truncated_polynomial(2)  # symmetric
    S, _ = nakayama_automorphism(alg, P)
    assert np.array_equal(S % P, np.eye(alg.m, dtype=np.int64))  # ν = id
    A_nu = twisted_bimodule(alg, psi=S)
    from quiverlab.engine.bimodule import regular_bimodule
    h_nu = hochschild_homology_with_coefficients(alg, A_nu, 5, primes=(P,))[P]
    h_reg = hochschild_homology_with_coefficients(alg, regular_bimodule(alg), 5, primes=(P,))[P]
    assert h_nu == h_reg


def test_twisted_bimodule_is_a_bimodule_for_any_automorphism():
    """{}_phi A_psi satisfies the bimodule axioms for any pair of automorphisms."""
    pytest.importorskip("quiverlab.engine.coxeter")  # enabled in Task 6 (coxeter)
    from quiverlab.engine.coxeter import nakayama_automorphism, identity_auto
    alg = quantum_ci(2)
    S, Sinv = nakayama_automorphism(alg, P)
    Id, _ = identity_auto(alg.m, P)
    for phi, psi in [(None, None), (S, None), (None, S), (S, Sinv), (Sinv, S)]:
        M = twisted_bimodule(alg, phi=phi, psi=psi)
        assert check_bimodule(M, P)


def test_identity_twist_is_regular_bimodule():
    """{}_1 A_1 is the regular bimodule A (reproduces HH_*(A,A) = HH_*(A))."""
    alg = quantum_ci(3)
    A1 = twisted_bimodule(alg)
    h = hochschild_homology_with_coefficients(alg, A1, 4, primes=(P,))[P]
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    assert h == hochschild_homology_dims(alg, 4, primes=(P,))[P]
