"""Regression oracle for the Hochschild (co)homology engine.

Every value frozen here matches a closed-form result from the literature or the
engine's own validation block; together they are the ground-truth oracle that the
bug fixes and the B1 Chouhy-Solotar resolution backend must reproduce exactly.
These tests are expected to PASS against the current bar-complex engine.
"""
import pytest

# cohomology_dims and quantum_ci live in scan3 (Task 4); every test in this oracle
# computes HH_* or HH^*, so the whole module self-heals once scan3 lands. (The
# homology backend quiverlab.engine.resolutions was pulled forward into Task 3.)
pytest.importorskip("quiverlab.engine.scan3")        # cohomology_dims, quantum_ci (Task 4)

from quiverlab.engine.hh_engine import (
    truncated_polynomial,
    two_gen_local,
    hochschild_homology_dims,
)
from quiverlab.engine.scan3 import hochschild_cohomology_dims, quantum_ci

# hanlab __init__ aliases, reproduced locally (see hanlab/__init__.py):
homology_dims = hochschild_homology_dims
cohomology_dims = hochschild_cohomology_dims
PRIME = 32003

P = PRIME


def kxy():
    """Commutative complete intersection k[x,y]/(x^2,y^2): x^2=0, y^2=0, yx=xy."""
    return two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x^2,y^2)")


# ---- Hochschild homology, characteristic-0 proxy (p = 32003) ----
@pytest.mark.parametrize("alg_fn,N,expected", [
    (lambda: truncated_polynomial(2), 9, [2, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
    (lambda: truncated_polynomial(3), 8, [3, 2, 2, 2, 2, 2, 2, 2, 2]),
    (lambda: truncated_polynomial(4), 5, [4, 3, 3, 3, 3, 3]),
    (kxy, 6, [4, 4, 5, 6, 7, 8, 9]),
    (lambda: quantum_ci(2), 6, [3, 2, 2, 2, 2, 2, 2]),
    (lambda: quantum_ci(3), 6, [3, 2, 2, 2, 2, 2, 2]),
])
def test_homology_char0(alg_fn, N, expected):
    assert homology_dims(alg_fn(), N)[P] == expected


# ---- characteristic-specific behaviour (small primes carry real signal) ----
def test_homology_char2_truncated2_doubles():
    # k[x]/(x^2): char-2 divided-power jump, HH_n = 2 for all n
    assert homology_dims(truncated_polynomial(2), 9)[2] == [2] * 10


def test_homology_char3_truncated3():
    assert homology_dims(truncated_polynomial(3), 8)[3] == [3] * 9


def test_homology_char2_kxy():
    assert homology_dims(kxy(), 6)[2] == [4, 8, 12, 16, 20, 24, 28]


# ---- the homology/cohomology ASYMMETRY (the program's scientific core) ----
@pytest.mark.parametrize("c", [2, 3])
def test_quantum_ci_asymmetry(c):
    Q = quantum_ci(c)
    assert homology_dims(Q, 6)[P] == [3, 2, 2, 2, 2, 2, 2]    # homology persists
    assert cohomology_dims(Q, 6)[P] == [2, 2, 1, 0, 0, 0, 0]  # cohomology dies


# ---- symmetric algebras: HH^n == HH_n in every degree ----
@pytest.mark.parametrize("alg_fn,N", [
    (lambda: truncated_polynomial(2), 7),
    (lambda: truncated_polynomial(3), 5),
    (kxy, 5),
])
def test_symmetric_homology_equals_cohomology(alg_fn, N):
    A = alg_fn()
    assert homology_dims(A, N)[P] == cohomology_dims(A, N)[P]


# ---- Kunneth: HH_*(A (x) B) = HH_*(A) (x) HH_*(B) ----
def test_kunneth_tensor_square_equals_commutative_ci():
    from quiverlab.engine.scan2 import tensor_product
    t = tensor_product(truncated_polynomial(2), truncated_polynomial(2))
    assert homology_dims(t, 6)[P] == homology_dims(kxy(), 6)[P] == [4, 4, 5, 6, 7, 8, 9]


# ---- characteristic-p COHOMOLOGY (item 9): the dual side also carries torsion signal ----
# Homology char-p is frozen above; cohomology was only pinned at the char-0 proxy. Freeze a
# small-prime cohomology sequence too. For a SYMMETRIC algebra HH^n == HH_n gives an
# independent cross-check of the value.
def test_cohomology_char2_truncated2():
    A = truncated_polynomial(2)
    assert cohomology_dims(A, 7)[2] == [2] * 8
    assert cohomology_dims(A, 7)[2] == homology_dims(A, 7)[2]   # symmetric: HH^ == HH_


def test_cohomology_char3_truncated3():
    A = truncated_polynomial(3)
    assert cohomology_dims(A, 5)[3] == [3] * 6
    assert cohomology_dims(A, 5)[3] == homology_dims(A, 5)[3]


def test_cohomology_char2_kxy():
    # symmetric, so HH^n == HH_n; char-2 sequence is 4,8,12,... (divided-power jump)
    assert cohomology_dims(kxy(), 5)[2] == [4, 8, 12, 16, 20, 24]


def test_cohomology_quantum_ci_asymmetry_is_char0_specific():
    # The HH^n -> 0 collapse of the quantum CI is a CHARACTERISTIC-0 (large-prime) phenomenon.
    # At p=2 the parameter c=2 == 0 degenerates the algebra, so cohomology does NOT die:
    Q2 = quantum_ci(2)
    assert Q2 and cohomology_dims(Q2, 6)[PRIME] == [2, 2, 1, 0, 0, 0, 0]   # char 0: dies
    assert cohomology_dims(Q2, 6)[2] == [2, 2, 3, 5, 7, 9, 11]             # p=2: does not


def test_cohomology_quantum_ci3_commutative_at_char2():
    # c=3 == 1 (mod 2), so qCI(3) becomes the commutative k[x,y]/(x^2,y^2) at p=2: its char-2
    # cohomology equals that of kxy (a cross-check the small prime makes visible).
    assert cohomology_dims(quantum_ci(3), 6)[2] == [4, 8, 12, 16, 20, 24, 28]
    assert cohomology_dims(kxy(), 6)[2] == [4, 8, 12, 16, 20, 24, 28]
