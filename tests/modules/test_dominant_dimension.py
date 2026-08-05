"""Dominant dimension: leading projective terms of the injective coresolution of
the regular module. Self-injective => infinity; hereditary kA2 => 1."""
import pytest

from quiverlab import Quiver, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import dominant_dimension

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_self_injective_is_infinite():
    A = truncated_polynomial(3, field=QQ)              # k[x]/(x^3): self-injective
    dd = dominant_dimension(A)
    assert dd.infinite is True
    assert dd.value is None


def test_hereditary_kA2_domdim_one():
    # A = P1 (+) P2; E(A) = P1 (+) P1 (projective, counts), E^1 = S1 (not
    # projective, stops) => domdim = 1.
    A = _kA2()
    dd = dominant_dimension(A)
    assert dd.value == 1 and dd.exact is True and dd.infinite is False
    assert int(dd) == 1 and dd == 1


def test_wrapper_matches():
    A = _kA2()
    assert A.dominant_dimension() == dominant_dimension(A)
