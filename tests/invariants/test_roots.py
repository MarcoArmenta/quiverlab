"""Root counts pinned to the classical tables (Bourbaki; ASS VII):
A_n: n(n+1)/2, D_n: n(n-1), E6: 36, E7: 63, E8: 120. Gabriel: for a
Dynkin quiver these are exactly the dimension vectors of the
indecomposables."""
import pytest

from quiverlab.families.dynkin import dynkin_quiver
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.roots import positive_roots

pytestmark = pytest.mark.oracle_literature


@pytest.mark.parametrize("typ,count", [
    ("A2", 3), ("A3", 6), ("A5", 15), ("D4", 12), ("D5", 20), ("E6", 36),
])
def test_root_counts(typ, count):
    A = dynkin_quiver(typ).algebra()
    roots = positive_roots(A)
    assert len(roots) == count
    assert len(set(roots)) == count            # no duplicates
    from quiverlab.invariants.forms import tits_form
    assert all(tits_form(A, list(r)) == 1 for r in roots)


def test_a2_roots_explicit():
    A = dynkin_quiver("A2").algebra()
    assert sorted(positive_roots(A)) == [(0, 1), (1, 0), (1, 1)]


def test_highest_root_d4():
    A = dynkin_quiver("D4").algebra()
    # generator center is vertex 2 (index 1); highest root (1, 2, 1, 1)
    assert sorted(max(positive_roots(A), key=sum)) == [1, 1, 1, 2]


def test_non_dynkin_refused():
    from quiverlab import Quiver
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra()
    with pytest.raises(QuiverlabError, match="Dynkin"):
        positive_roots(K)
