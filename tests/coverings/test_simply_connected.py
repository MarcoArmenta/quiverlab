"""Three-valued is_simply_connected. The CORRECTED square oracle drives the
True/False split; trees True; disconnected/cyclic/nontrivial-pi1ab False;
bypass-bearing perfect corner => None (Adian-Rabin honesty)."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.invariants.coverings import (
    bypasses, fundamental_group, is_simply_connected)

pytestmark = pytest.mark.oracle_literature


def test_square_with_relation_true():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=["a*b - c*d"], field=GF(7))
    assert is_simply_connected(A).verdict is True          # no bypasses + pi1 trivial (R3)


def test_square_without_relation_false():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=[], field=GF(7))
    r = is_simply_connected(A)
    assert r.verdict is False and r.witness["kind"] == "nontrivial_pi1ab"


def test_tree_true():
    assert is_simply_connected(linear_path_algebra(5, field=GF(7))).verdict is True


def test_loop_not_triangular_false():
    A = truncated_polynomial(3, field=GF(7))               # loop: not triangular
    r = is_simply_connected(A)
    assert r.verdict is False and r.witness["kind"] == "oriented_cycle"


def test_disconnected_false():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)}).algebra(field=GF(5))
    r = is_simply_connected(A)
    assert r.verdict is False and r.witness["kind"] == "disconnected"


def test_never_true_from_failed_search():
    # ADMISSIBLE bypass example (B3): alpha:1->3 bypasses the path a*b:1->2->3; the
    # suffix arrow s:3->4 lets the ADMISSIBLE (rad^2) relation alpha*s - a*b*s kill the
    # triangle cycle => pi1^ab = 0 (NOT caught by the nontrivial-pi1ab False guard).
    # With only the cheap routes (strong=False): NOT a tree (R1 off), a bypass is
    # present (R3 off) => the honest verdict is None -- the cheap routes can NEVER
    # fabricate True from a failed search.
    Q = Quiver([1, 2, 3, 4], {"al": (1, 3), "a": (1, 2), "b": (2, 3), "s": (3, 4)})
    A = Q.algebra(relations=["al*s - a*b*s"], field=GF(7))
    assert fundamental_group(A).free_rank == 0             # pi1^ab = 0 (cycle killed)
    assert bypasses(A)                                     # (al, a*b) is a bypass
    assert is_simply_connected(A, strong=False).verdict is None   # never a search-True
