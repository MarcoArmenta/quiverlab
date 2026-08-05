"""Crosscheck plumbing with a canned M2 transcript -- no M2 binary."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.fast


def _dual_numbers():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(7))


def test_graded_dims_agree_path(monkeypatch):
    # k[x]/(x^2): graded dims 1,1,0,0,0
    monkeypatch.setattr(session, "run_script",
                        lambda s, timeout=120: "<<QL>> 0 1\n<<QL>> 1 1\n"
                                               "<<QL>> 2 0\n<<QL>> 3 0\n<<QL>> 4 0\n")
    rep = cc.crosscheck_graded_dims(_dual_numbers(), top=4)
    assert rep.agree and rep.ours == rep.qpa == [1, 1, 0, 0, 0]
    rep.assert_agree()


def test_disagreement_is_loud(monkeypatch):
    monkeypatch.setattr(session, "run_script",
                        lambda s, timeout=120: "<<QL>> 0 1\n<<QL>> 1 2\n"
                                               "<<QL>> 2 0\n<<QL>> 3 0\n<<QL>> 4 0\n")
    rep = cc.crosscheck_graded_dims(_dual_numbers(), top=4)
    assert not rep.agree
    with pytest.raises(AssertionError):
        rep.assert_agree()


def test_dispatcher_unknown_subject():
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="graded_dims"):
        cc.crosscheck(_dual_numbers(), "hochschild")   # honest: no M2 HH


def test_commutative_guard():
    from quiverlab.errors import QuiverlabError
    # A finite-dimensional, genuinely NON-commutative single-vertex algebra
    # whose presentation carries NO commutator relation (x*y = 0 but y*x != 0,
    # basis 1,x,y,yx). NOTE: the plan's k<x,y>/(x^2,y^2) is INFINITE-dimensional
    # (the free product of two dual numbers -- alternating words never die) and
    # would raise NotFiniteDimensionalError at CONSTRUCTION, never reaching the
    # commutativity guard this test exercises.
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y"], field=GF(7))
    with pytest.raises(QuiverlabError, match="commutative"):
        cc.crosscheck_commutative_ext(A, ["x", "y"], ["x^2", "y^2"], top=4)
