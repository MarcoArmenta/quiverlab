"""Finite repetitive-algebra slices (Plan 44 / C7). copies=1 is A itself (byte-Cartan);
the general slice is certified by dim == (2*copies - 1)*dim A. The full repetitive algebra
is infinite-dimensional -- out of scope, stated on the verification page."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.repetitive import repetitive_slice

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_one_copy_is_A():
    A = _kA2()
    R1 = repetitive_slice(A, 1)
    assert R1.dim == A.dim
    assert R1.cartan_matrix() == A.cartan_matrix()          # copies=1 == A


def test_two_slice_dimension():
    A = _kA2()                                              # dim 3
    R2 = repetitive_slice(A, 2)
    assert R2.dim == (2 * 2 - 1) * A.dim == 9               # 2 A's + 1 D(A)


@pytest.mark.parametrize("copies", [1, 2, 3])
def test_slice_dimension_certificate(copies):
    A = _kA2()
    R = repetitive_slice(A, copies)
    assert R.dim == (2 * copies - 1) * A.dim


def test_bad_copies_refused():
    A = _kA2()
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError):
        repetitive_slice(A, 0)
