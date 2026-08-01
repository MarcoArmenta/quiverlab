"""Plan 35 literature pins for the HH product surface.

k[x]/(x^2), the dual numbers. Over GF(32003) (a char-0-shaped prime, big enough
that 2 is a unit) the Hochschild (co)homology is the classical
    HH^0 = Z(A) = A            (commutative => the whole algebra),   dim 2
    HH^n = k  for every n >= 1,                                      dim 1
so the cup ring is k[u] (deg-2 generator u) tensor an exterior deg-1 class, with
the char-0 relation that the odd generator squares to zero. In char 2 that odd
square SURVIVES (the classical char-2 phenomenon). QuantumCI q=1 over GF(2) is
the (2,2) commutative complete intersection whose HH dims are the BGMS pin
[4, 8, 12, 16].

ENGINE-VERIFIED (2026-08-01): the brief's stated degree-0 dims (1 for HH^0 and
HH_0) DISAGREE with the engine, which returns 2 -- and the engine wins (the
CRS-2004 precedent). The reason is textbook: for a COMMUTATIVE algebra
HH^0 = Z(A) = A and HH_0 = A/[A,A] = A, both of dimension dim A = 2 for the dual
numbers, not 1. The brief's "dim HH^n = 1 for all n" holds only for n >= 1. The
engine dims are oracle-gated upstream (tests/test_tt_calculus.py,
tests/engine/test_gerstenhaber.py, and the Task-6 identity batteries):
    A.hochschild_cohomology(4) -> HH^0=2 HH^1=1 HH^2=1 HH^3=1 HH^4=1
    A.hochschild_homology(4)   -> HH_0=2 HH_1=1 HH_2=1 HH_3=1 HH_4=1
The two degree-0 expectations below are the corrected (verified) values.
"""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_literature]


def _entry(t, k, i, j):
    return int(t.constants[k][i][j])


def test_dual_numbers_cup_ring_char0_shape():
    A = ql.truncated_polynomial(2, field=ql.GF(32003))
    hp = A.cup_products(4)
    dims = {n: hp.tables[(0, n)].dims[1] for n in range(5)}
    # CORRECTION (engine wins over the brief): HH^0 = Z(A) = A has dim 2 for the
    # commutative dual numbers; the brief said 1. HH^n = k (dim 1) for n >= 1.
    assert dims == {0: 2, 1: 1, 2: 1, 3: 1, 4: 1}
    # even generators compose to a nonzero multiple of the even basis class
    assert _entry(hp.tables[(2, 2)], 0, 0, 0) != 0
    # odd squares vanish in char != 2 (graded commutativity forces 2x = 0)
    assert _entry(hp.tables[(1, 1)], 0, 0, 0) == 0
    assert _entry(hp.tables[(3, 1)], 0, 0, 0) == 0


def test_dual_numbers_char2_odd_squares_survive():
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    hp = A.cup_products(2)
    # char 2: HH^* of k[x]/(x^2) is the polynomial-like ring where the odd
    # generator squares NONZERO (the classical char-2 phenomenon)
    assert _entry(hp.tables[(1, 1)], 0, 0, 0) != 0


def test_qci_dims_line_up_with_bgms():
    A = ql.QuantumCI(field=ql.GF(2), q=1)
    hp = A.cup_products(2)
    assert [hp.tables[(0, n)].dims[1] for n in range(3)] == [4, 8, 12]


def test_connes_b_rank_dual_numbers():
    A = ql.truncated_polynomial(2, field=ql.GF(32003))
    cb = A.connes_differentials(4)
    # CORRECTION (engine wins over the brief): HH_0 = A/[A,A] = A has dim 2 for
    # the commutative dual numbers; the brief said 1. HH_n = k (dim 1) for n >= 1.
    # Connes B alternates iso/zero along the SBI pattern: rank B_0 = 1 here.
    assert cb.hh_dims == [2, 1, 1, 1, 1]
    assert cb.ranks[0] == 1
