"""Regression oracle for the automorphism / Coxeter route.

Pins Theorem C (the calculus automorphism is the identity on HH^*), the asymmetric
non-vanishing certificate (nontrivial on HH_*), the chain-map property that any future
AR-bimodule action must also satisfy, and the classical Coxeter polynomial of kA_3.
"""
import numpy as np
import pytest
from sympy import symbols, expand

from quiverlab.engine.scan3 import quantum_ci
from quiverlab.engine.coxeter import (
    identity_auto,
    is_identity,
    nakayama_quantum_ci,
    induced_on_HH_homology,
    induced_on_HH_cohomology,
)

# hanlab __init__ aliases, reproduced locally:
induced_homology = induced_on_HH_homology
induced_cohomology = induced_on_HH_cohomology
PRIME = 32003
P = PRIME


# ---- identity automorphism induces the identity on every HH_n and HH^n ----
@pytest.mark.parametrize("n", [0, 1, 2, 3])
def test_identity_automorphism_induces_identity(n):
    A = quantum_ci(2)
    S, Sinv = identity_auto(A.m, P)
    Mh, dh = induced_homology(A, n, S, P)
    Mc, dc = induced_cohomology(A, n, S, Sinv, P)
    assert is_identity(Mh, P)
    assert is_identity(Mc, P)


# ---- Theorem C certificate: Nakayama auto is NONtrivial on HH_2, identity on HH^2 ----
def test_theoremC_certificate_quantum_ci():
    Q = quantum_ci(2)
    S, Sinv = nakayama_quantum_ci(2, P)
    Mh, dh = induced_homology(Q, 2, S, P)
    Mc, dc = induced_cohomology(Q, 2, S, Sinv, P)
    assert not is_identity(Mh, P)   # non-vanishing witness, visible on homology
    assert is_identity(Mc, P)       # invisible to cohomology


# ---- chain-map commutation b_n . Sigma_n == Sigma_{n-1} . b_n  (guards B1 / item A) ----
@pytest.mark.parametrize("n", [1, 2])
def test_sigma_is_chain_map(n):
    from quiverlab.engine.hh_engine import differential_matrix
    from quiverlab.engine.coxeter import sigma_chain_matrix
    Q = quantum_ci(2)
    S, Sinv = nakayama_quantum_ci(2, P)
    Sig_n, basis_n, idx_n = sigma_chain_matrix(Q, n, S, P)
    Sig_n1, basis_n1, idx_n1 = sigma_chain_matrix(Q, n - 1, S, P)
    Bn = differential_matrix(Q, n, basis_n, idx_n1)
    lhs = (Bn @ Sig_n) % P
    rhs = (Sig_n1 @ Bn) % P
    assert np.array_equal(lhs, rhs)


# ---- classical Coxeter polynomial of the hereditary algebra kA_3 ----
def test_coxeter_polynomial_kA3():
    pytest.importorskip("quiverlab.engine.coxeter2")  # enabled in Task 11 (coxeter2)
    from quiverlab.engine.coxeter2 import cartan_from_raw, coxeter_polynomial_from_cartan
    paths = [(i, j) for i in range(3) for j in range(i, 3)]
    C = cartan_from_raw(3, None, paths)
    poly, Phi = coxeter_polynomial_from_cartan(C)
    t = symbols('t')
    assert expand(poly - (t + 1) * (t ** 2 + 1)) == 0
