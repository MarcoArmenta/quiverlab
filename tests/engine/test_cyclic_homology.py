"""Connes' B operator and cyclic homology HC_* (PLAN item D).

The rigorous oracle is the **mixed-complex identities** b^2 = 0, B^2 = 0,
bB + Bb = 0, asserted directly on the matrices across the zoo and several primes:
they prove (C_*, b, B) is a mixed complex, so its (b, B) total complex computes a
genuine cyclic homology.  The *specific* Connes convention is then pinned by closed
forms that the construction must reproduce:

  * HC_*(k) = [1, 0, 1, 0, ...]                       (the ground field);
  * HC_*(k x k) = [2, 0, 2, 0, ...]                   (separable: periodicity S
    pushes HH_0 into every even degree, HH_n = 0 for n > 0);
  * HC_1(k[x]/x^2) = 0                                (= Omega^1 / dA, char 0);
  * HC_*(k[x]/x^a) = [a, 0, a, 0, ...] (char-0 proxy) (truncated polynomials;
    consistent with Goodwillie: HP_*(A) = HP_*(A_red) = HP_*(k));
  * HC_0(A) = dim HH_0(A) across the zoo               (always).
"""

import numpy as np
import pytest

# hanlab sys.path shim dropped: quiverlab uses absolute package imports.
from quiverlab.engine.hh_engine import (
    Algebra,
    truncated_polynomial,
    two_gen_local,
    hochschild_homology_dims,
)
from quiverlab.engine.cyclic import (
    connes_B_matrix,
    cyclic_homology_dims,
    check_mixed_complex,
)
from quiverlab.engine.scan3 import quantum_ci

P = 32003


def kk_semisimple():
    """A = k x k (two orthogonal idempotents); separable, unit = e_0 + e_1."""
    T = np.zeros((2, 2, 2), dtype=np.int64)
    T[0, 0, 0] = 1
    T[1, 1, 1] = 1
    return Algebra(2, T, np.array([1, 1], dtype=np.int64), name="kxk")


def _zoo():
    return [
        truncated_polynomial(1),
        truncated_polynomial(2),
        truncated_polynomial(3),
        kk_semisimple(),
        two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x2,y2)"),
        two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, -3], "qci3"),
    ]


# ---------------------------------------------------------------------------
# The mixed-complex identities -- the rigorous descent oracle
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("alg", _zoo(), ids=lambda a: a.name)
@pytest.mark.parametrize("prime", [32003, 2, 3, 5])
def test_mixed_complex_identities(alg, prime):
    """b^2 = 0, B^2 = 0, bB + Bb = 0 on the bar complex over F_p."""
    fail = check_mixed_complex(alg, 4, prime=prime)
    assert fail is None, f"{alg.name}: mixed-complex identity {fail} failed at p={prime}"


def test_B_squared_zero_explicit():
    """B_{n+1} . B_n = 0 directly on the matrices (k[x]/x^2)."""
    A = truncated_polynomial(2)
    for n in range(0, 4):
        Bn = connes_B_matrix(A, n)
        Bnp1 = connes_B_matrix(A, n + 1)
        assert np.all((Bnp1 @ Bn) % P == 0)


def test_B_kills_unit_in_A_slot():
    """B(1 (x) r_1 (x) ... ) = 0: the extra-degeneracy image is annihilated."""
    A = truncated_polynomial(3)
    from quiverlab.engine.hh_engine import cn_basis
    basis2 = cn_basis(A, 2)
    B = connes_B_matrix(A, 2)
    for j, gen in enumerate(basis2):
        if gen[0] == A.t:            # A-slot is the unit
            assert np.all(B[:, j] == 0)


# ---------------------------------------------------------------------------
# Closed-form cyclic homology
# ---------------------------------------------------------------------------
def test_HC_of_ground_field():
    """HC_*(k) = [1, 0, 1, 0, ...]."""
    k = truncated_polynomial(1)
    for prime in (P, 2, 3, 5):
        hc = cyclic_homology_dims(k, 7, primes=(prime,))[prime]
        assert hc == [1, 0, 1, 0, 1, 0, 1, 0]


def test_HC_of_kxk_separable():
    """HC_*(k x k) = [2, 0, 2, 0, ...]: separable, periodicity lifts HH_0."""
    kk = kk_semisimple()
    # separable: HH_n = 0 for n > 0
    hh = hochschild_homology_dims(kk, 5, primes=(P,))[P]
    assert hh == [2, 0, 0, 0, 0, 0]
    for prime in (P, 2, 3, 5):
        hc = cyclic_homology_dims(kk, 7, primes=(prime,))[prime]
        assert hc == [2, 0, 2, 0, 2, 0, 2, 0]


def test_HC1_truncated_is_zero_char0():
    """HC_1(k[x]/x^a) = Omega^1/dA = 0 in char 0 (independent de Rham check)."""
    for a in (2, 3, 4):
        hc = cyclic_homology_dims(truncated_polynomial(a), 1, primes=(P,))[P]
        assert hc[1] == 0


def test_HC_truncated_polynomial_char0():
    """HC_*(k[x]/x^a) = [a, 0, a, 0, ...] over the char-0 proxy (regression pin)."""
    for a in (2, 3, 4):
        hc = cyclic_homology_dims(truncated_polynomial(a), 6, primes=(P,))[P]
        assert hc == [a, 0, a, 0, a, 0, a]


def test_HC0_equals_HH0_across_zoo():
    """HC_0(A) = dim HH_0(A) for every algebra (always true)."""
    for alg in _zoo():
        hh0 = hochschild_homology_dims(alg, 0, primes=(P,))[P][0]
        hc0 = cyclic_homology_dims(alg, 0, primes=(P,))[P][0]
        assert hc0 == hh0, alg.name


def test_HC_quantum_ci_runs_and_is_consistent():
    """The non-monomial quantum CI: HC_0 = HH_0 and the bicomplex is a mixed complex."""
    A = quantum_ci(3)
    assert check_mixed_complex(A, 3) is None
    hh0 = hochschild_homology_dims(A, 0, primes=(P,))[P][0]
    hc = cyclic_homology_dims(A, 3, primes=(P,))[P]
    assert hc[0] == hh0
