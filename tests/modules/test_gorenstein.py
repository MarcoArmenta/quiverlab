"""Gorenstein dimension via the injective dimension of the regular module on both
sides. Self-injective => 0; hereditary => gl.dim; both cases Gorenstein."""
import pytest

from quiverlab import Quiver, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import gorenstein_dimension, is_gorenstein

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_self_injective_gorenstein_zero():
    A = truncated_polynomial(3, field=QQ)              # symmetric => id(A) = 0
    gd = gorenstein_dimension(A)
    assert gd.right_id == 0 and gd.left_id == 0
    assert gd.is_gorenstein is True
    assert is_gorenstein(A) is True


def test_hereditary_gorenstein_equals_gldim():
    A = _kA2()                                          # gl.dim 1
    gd = gorenstein_dimension(A)
    assert gd.right_id == 1 and gd.left_id == 1
    assert gd.is_gorenstein is True
    assert int(A.global_dimension()) == 1


def test_wrapper_matches():
    A = _kA2()
    assert A.is_gorenstein() is True
