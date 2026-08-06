"""Live M2 battery: Ext_A(k,k) dims of commutative examples vs the minimal
free resolution M2 computes -- a fully independent homological route."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.skipif(session.should_skip_m2(),
                                reason="Macaulay2 not installed")


def _comm(relations, extra, p):
    """k[x,y]/(relations) as a quiver algebra: loops + commutator."""
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=relations + ["x*y-y*x"] + extra, field=GF(p))


CASES = [
    ("kxy_x2y2_gf7",   ["x*x", "y*y"],        [], 7, ["x^2", "y^2"], 6),
    ("kxy_x2y2_gf3",   ["x*x", "y*y"],        [], 3, ["x^2", "y^2"], 6),
    ("kxy_x3y2",       ["x*x*x", "y*y"],      [], 5, ["x^3", "y^2"], 6),
    ("kxy_x2y3",       ["x*x", "y*y*y"],      [], 32003, ["x^2", "y^3"], 6),
]


@pytest.mark.parametrize("name,rels,extra,p,m2rels,top", CASES,
                         ids=[c[0] for c in CASES])
def test_ext_k_matches_m2(name, rels, extra, p, m2rels, top):
    A = _comm(rels, extra, p)
    cc.crosscheck_commutative_ext(A, ["x", "y"], m2rels, top=top).assert_agree()
