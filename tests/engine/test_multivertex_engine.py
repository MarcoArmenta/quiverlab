"""Item 5 (TODO P2): the multi-idempotent basis transport (1_A = sum e_i).

CLAUDE.md calls the unit-adapted basis "what makes multi-vertex algebras work", but
tests/test_basis_transport.py only ever uses a single-idempotent unit (t=0, B=I). Here
we exercise the genuine `1_A = sum_i e_i` change-of-basis (B with column t = the full
unit vector, so B != I) end-to-end: the constructor must succeed, e_t must act as a
two-sided identity in the f-basis, the reduced index set must drop exactly t, and the
resulting Hochschild homology must agree with an independent oracle (Bardzell) and with
closed-form values.
"""
import numpy as np
import pytest

# This oracle exercises the multi-vertex unit-adaptation end-to-end against the
# coxeter2 cyclic-Nakayama builder (Task 11) and the Bardzell resolution oracle
# (Task 10) -- and the parametrize builds cyclic_nakayama algebras -- so it
# self-heals as a whole once both land. (The bar-homology backend
# quiverlab.engine.resolutions was pulled forward into Task 3.)
pytest.importorskip("quiverlab.engine.resolutions_bardzell")  # Bardzell oracle (Task 10)
pytest.importorskip("quiverlab.engine.coxeter2")              # cyclic_nakayama (Task 11)

from quiverlab.engine.hh_engine import Algebra, check_associative, hochschild_homology_dims
from quiverlab.engine.coxeter2 import cyclic_nakayama
from quiverlab.engine.scan2 import kA
from quiverlab.engine.resolutions_bardzell import BardzellResolution, MonomialPresentation

# hanlab __init__ alias, reproduced locally:
homology_dims = hochschild_homology_dims
PRIME = 32003
P = PRIME


def _is_two_sided_identity(A, idx):
    I = np.eye(A.m, dtype=np.int64)
    return all(np.array_equal(A.T[idx, j, :], I[j]) and
               np.array_equal(A.T[j, idx, :], I[j]) for j in range(A.m))


@pytest.mark.parametrize("builder,label,dim", [
    (lambda: cyclic_nakayama(2, 2)[0], "kZ_2/rad^2", 4),
    (lambda: cyclic_nakayama(3, 2)[0], "kZ_3/rad^2", 6),
    (lambda: kA(3), "kA_3", 6),
])
def test_multivertex_unit_adaptation_invariants(builder, label, dim):
    A = builder()
    assert A.m == dim
    ok, _ = check_associative(A)
    assert ok, f"{label} not associative"
    # 1_A is a single basis vector e_t that acts as a two-sided identity
    assert _is_two_sided_identity(A, A.t)
    # reduced space A/k.1 = "drop coordinate t": R is all indices except t
    assert A.t not in A.R
    assert sorted(A.R + [A.t]) == list(range(A.m))
    assert len(A.R) == A.m - 1


def test_genuine_multi_idempotent_change_of_basis():
    # k x k with 1 = e_0 + e_1 (TWO idempotents): the change-of-basis B has column t equal
    # to [1, 1], so B != I -- this is the transport single-idempotent tests never reach.
    m = 2
    T = np.zeros((m, m, m), dtype=np.int64)
    T[0, 0, 0] = 1
    T[1, 1, 1] = 1                          # e_i e_i = e_i (no cross terms)
    A = Algebra(m, T, np.array([1, 1], dtype=np.int64), name="k x k")
    ok, _ = check_associative(A)
    assert ok
    assert _is_two_sided_identity(A, A.t)
    # k x k is semisimple: HH_0 = k^2, HH_n = 0 for n >= 1
    assert homology_dims(A, 3)[P] == [2, 0, 0, 0]


def test_multivertex_homology_matches_bardzell_oracle():
    # The transported multi-vertex algebra must give the SAME homology as the independent
    # Bardzell backend on the same monomial algebra.
    A, _ = cyclic_nakayama(2, 2)
    pres = MonomialPresentation.cyclic_nakayama(2, 2)
    bard = BardzellResolution(pres)
    assert homology_dims(A, 4)[P] == homology_dims(A, 4, resolution=bard)[P] == [2, 1, 1, 1, 1]
