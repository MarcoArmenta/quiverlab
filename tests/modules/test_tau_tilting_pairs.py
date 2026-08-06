"""Support tau-tilting pairs (Plan 45 / C4). Self-cert: the four axioms are enforced on
construction; (A,0) and (0,A) are the extreme pairs; a non-tau-rigid or wrong-rank input
refuses loudly."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.tautilting.pairs import (initial_pair, make_pair, terminal_pair)

pytestmark = pytest.mark.oracle_selfcert


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def test_initial_pair_is_the_free_module():
    A = _kA2()
    p = initial_pair(A)
    assert len(p.summands) == 2 and p.support == frozenset()      # (P1 (+) P2, {})
    assert p.g_matrix() == [[1, 0], [0, 1]]


def test_terminal_pair_kills_everything():
    A = _kA2()
    p = terminal_pair(A)
    assert p.summands == () and p.support == frozenset({1, 2})
    assert p.g_matrix() == [[-1, 0], [0, -1]]                     # columns -e_1, -e_2


def test_apr_style_pair_validates():
    # (P1 (+) S1, {}) : over kA2 this is a genuine tau-tilting pair (the APR tilt).
    A = _kA2()
    p = make_pair(A, [A.projective(1), A.simple(1)], support=[])
    assert len(p.summands) == 2 and p.support == frozenset()


def test_support_pair_validates():
    # (S1, {2}) : S1 tau-rigid, P_2 killed with Hom(P_2, S1) = 0, |M|+|supp| = 1+1 = 2.
    A = _kA2()
    p = make_pair(A, [A.simple(1)], support=[2])
    assert p.support == frozenset({2})


def test_non_tau_rigid_refused():
    A = _kA2()
    # wrong rank is the cleanest loud case: two summands claimed with a nonempty support
    # => rank 3 != 2.
    with pytest.raises(QuiverlabError):
        make_pair(A, [A.projective(1), A.simple(1)], support=[2])   # rank 3 != 2


def test_support_axiom_violation_refused():
    # (P1, {1}) : Hom(P_1, P_1) != 0, so the support axiom fails -- loud.
    A = _kA2()
    with pytest.raises(QuiverlabError):
        make_pair(A, [A.projective(1)], support=[1])
