"""loewy_length / complexity / center (spec section 3.5). Fixtures A, B; k[x]/(x^n)."""
import pytest
from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial
from quiverlab.errors import FieldError


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_loewy_length():
    assert linear_path_algebra(2).loewy_length() == 2        # rad^2 = 0
    assert _square().loewy_length() == 3                     # rad^3 = 0, rad^2 = <a*b>
    assert truncated_polynomial(4).loewy_length() == 4       # k[x]/(x^4): rad^4 = 0
    assert truncated_polynomial(2, field=GF(3)).loewy_length() == 2


def test_center_dimension_and_basis():
    d, basis = linear_path_algebra(2).center()
    assert d == 1                                            # connected -> center = k*1
    assert len(basis) == 1
    assert _square().center()[0] == 1
    # k[x]/(x^2) is commutative -> center is the whole algebra
    assert truncated_polynomial(2).center()[0] == 2


def test_center_over_gfp():
    assert linear_path_algebra(2, field=GF(7)).center()[0] == 1


def test_complexity_gfp():
    # finite gl.dim (kA_2) -> minimal A^e resolution terminates -> complexity 0
    assert linear_path_algebra(2, field=GF(32003)).complexity(6) == 0
    # k[x]/(x^2) self-injective -> constant-rank periodic resolution -> complexity 1
    assert truncated_polynomial(2, field=GF(32003)).complexity(6) == 1


def test_complexity_cc_loud():
    with pytest.raises(FieldError):
        truncated_polynomial(2, field=CC).complexity(4)
