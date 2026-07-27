"""Oracle + guard for cartan_from_raw and auto_to_f_basis (PLAN.md hardening).

auto_to_f_basis inverts the change-of-basis exactly (Sherman-Morrison integer
inverse, valid only for unit[t] == 1); a non-unimodular B must be rejected with a
loud QuiverlabError -- not a bare assert, which python -O would strip, silently
returning a wrong matrix in an otherwise exact pipeline.
"""
import numpy as np
import pytest
from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.coxeter2 import cartan_from_raw, auto_to_f_basis
from quiverlab.errors import QuiverlabError


@pytest.mark.oracle_selfcert
def test_cartan_from_raw_kZ3_radsq():
    # paths as (start, end): vertices (i,i) and arrows i -> i+1
    paths = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)]
    C = cartan_from_raw(3, None, paths)
    assert np.array_equal(C, np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]]))


@pytest.mark.oracle_selfcert
def test_auto_to_f_basis_identity_is_exact():
    A = truncated_polynomial(2)                 # unit = e_0, t = 0, so B = I
    unit = np.eye(A.m, dtype=np.int64)[A.t]
    out = auto_to_f_basis(A, unit, np.eye(A.m, dtype=np.int64))
    assert np.array_equal(out, np.eye(A.m, dtype=np.int64))


def test_auto_to_f_basis_rejects_non_unimodular():
    # unit[t] = 2 makes the change-of-basis B non-unimodular; the Sherman-Morrison
    # inverse only holds for unit[t] == 1, so this must refuse loudly (-O safe).
    A = truncated_polynomial(2)                 # t = 0
    bad_unit = np.array([2, 0], dtype=np.int64)
    with pytest.raises(QuiverlabError, match="not unimodular"):
        auto_to_f_basis(A, bad_unit, np.eye(A.m, dtype=np.int64))
