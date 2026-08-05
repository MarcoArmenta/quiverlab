"""Finitistic-dimension bounds. Exact when gl.dim is finite (findim = gl.dim);
honest [lower, upper] otherwise. Barrios-Mata truncated closed forms."""
import pytest

from quiverlab import Quiver, TruncatedPathAlgebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import finitistic_dimension_bounds

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_hereditary_exact():
    A = _kA2()                                          # gl.dim 1, findim 1
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 1 and fb.upper == 1 and fb.exact is True


def test_nakayama_exact_equals_gldim():
    A = TruncatedPathAlgebra("A3", 2, field=QQ)         # gl.dim 2, findim 2
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 2 and fb.upper == 2 and fb.exact is True


def test_self_injective_lower_zero_upper_valid():
    # k[x]/(x^3): findim = 0 (only projectives have finite pd), gl.dim infinite.
    A = truncated_polynomial(3, field=QQ)
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 0
    assert fb.exact is False
    # upper is a valid bound (>= lower) or an honest None
    assert fb.upper is None or fb.upper >= fb.lower


def test_wrapper_matches():
    A = _kA2()
    fb = A.finitistic_dimension_bounds()
    assert fb.exact and fb.lower == 1
