"""Oracle + guard for cartan_from_raw and auto_to_f_basis (PLAN.md hardening).

test_auto_to_f_basis_rejects_non_unimodular is RED until the unimodularity guard is
added: auto_to_f_basis rounds a floating-point matrix inverse to integers, and for a
non-unimodular change-of-basis that silently returns a WRONG integer matrix in an
otherwise exact pipeline. It must fail loudly instead.
"""
import numpy as np
import pytest
from quiverlab.engine.hh_engine import truncated_polynomial
from quiverlab.engine.coxeter2 import cartan_from_raw, auto_to_f_basis


def test_cartan_from_raw_kZ3_radsq():
    # paths as (start, end): vertices (i,i) and arrows i -> i+1
    paths = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)]
    C = cartan_from_raw(3, None, paths)
    assert np.array_equal(C, np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]]))


def test_auto_to_f_basis_identity_is_exact():
    A = truncated_polynomial(2)                 # unit = e_0, t = 0, so B = I
    unit = np.eye(A.m, dtype=np.int64)[A.t]
    out = auto_to_f_basis(A, unit, np.eye(A.m, dtype=np.int64))
    assert np.array_equal(out, np.eye(A.m, dtype=np.int64))


def test_auto_to_f_basis_rejects_non_unimodular():
    # unit[t] = 2 makes the change-of-basis B non-unimodular; the float-rounded
    # inverse is then wrong and must be rejected, not returned silently.
    A = truncated_polynomial(2)                 # t = 0
    bad_unit = np.array([2, 0], dtype=np.int64)
    with pytest.raises(AssertionError):
        auto_to_f_basis(A, bad_unit, np.eye(A.m, dtype=np.int64))
