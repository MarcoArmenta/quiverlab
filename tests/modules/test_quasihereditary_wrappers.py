"""Algebra-level wrappers for the quasi-hereditary surface (Plan 47)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import is_isomorphic

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def test_wrappers_delegate():
    A = _a3()
    assert A.is_quasi_hereditary().is_quasi_hereditary is True
    assert all(A.standard_modules()[v].dim == 1 for v in (1, 2, 3))
    R = A.recollement([2])
    assert A.corner_algebra([1, 3]).dim == 3               # the eAe trap via the wrapper
    assert A.quotient_by_idempotent([2]).dim == 2
    assert R.eAe.dim == 1
    T = A.characteristic_tilting()
    from quiverlab.modules.morphism import direct_sum
    DA, _, _ = direct_sum(*[A.injective(v) for v in (1, 2, 3)])
    assert is_isomorphic(T, DA)


def test_quasi_hereditary_block_shape():
    from quiverlab.modules.quasihereditary import quasi_hereditary_block
    b = quasi_hereditary_block(_a3())
    assert b["kind"] == "quasi_hereditary"
    assert b["is_quasi_hereditary"] is True
    assert b["gl_dim"] == {"value": 1, "exact": True}
    assert set(b["standard_dims"]) == {"1", "2", "3"}
    assert "order-dependent" in b["order_note"]
    assert b["references"] == ["dlab_ringel", "assem_book"]
