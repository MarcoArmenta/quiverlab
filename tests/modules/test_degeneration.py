"""Degeneration = hom order for representation-finite algebras (Plan 49 / C8).
Literature (hand-derived): kA2 d=(1,1) is the 2-chain S1(+)S2 <_deg P1; kA3
d=(1,1,1) is the diamond (semisimple bottom, two incomparable middles, uniserial
top). Self-cert: the relation is a partial order; orbit dim strictly increases up
each cover; the unique maximum is the generic (rigid) module."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.geometry import is_rigid, orbit_dimension
from quiverlab.modules.degeneration import degeneration_order

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def _reachability(P):
    """Reflexive-transitive closure of the covers as a bool matrix leq[a][b]
    ('a degenerates from b': a is at or below b)."""
    n = len(P.vertices)
    leq = [[i == j for j in range(n)] for i in range(n)]
    for lo, hi in P.covers:
        leq[lo][hi] = True
    changed = True
    while changed:
        changed = False
        for a in range(n):
            for b in range(n):
                if leq[a][b]:
                    for c in range(n):
                        if leq[b][c] and not leq[a][c]:
                            leq[a][c] = True
                            changed = True
    return leq


@lit
def test_kA2_11_is_a_two_chain():
    A = _kA2()
    P = degeneration_order(A, {1: 1, 2: 1})
    assert P.is_complete and len(P.vertices) == 2      # {P1}, {S1 (+) S2}
    assert len(P.covers) == 1                          # a single cover
    lo, hi = P.covers[0]
    # the lower (degenerate) class is S1 (+) S2 (two summands); the upper is P1
    assert len(P.vertices[lo]["summands"]) == 2
    assert len(P.vertices[hi]["summands"]) == 1
    assert P.vertices[hi]["orbit_dim"] > P.vertices[lo]["orbit_dim"]
    assert P.vertices[hi]["is_generic"] and is_rigid(P.vertices[hi]["module"])


@lit
def test_kA3_111_is_the_diamond():
    A = linear_path_algebra(3, field=QQ)
    P = degeneration_order(A, {1: 1, 2: 1, 3: 1})
    assert P.is_complete and len(P.vertices) == 4       # [1,3]; [1,2](+)S3; S1(+)[2,3]; S1(+)S2(+)S3
    # diamond: 4 covers (bottom->2 middles, 2 middles->top), 2 incomparable middles
    assert len(P.covers) == 4
    orbit = sorted(v["orbit_dim"] for v in P.vertices)
    assert orbit == [0, 1, 1, 2]                        # Z=0, A=B=1, G=2 (open)
    tops = [v for v in P.vertices if v["is_generic"]]
    assert len(tops) == 1 and tops[0]["orbit_dim"] == 2
    # the two middles are incomparable (a genuine diamond, not a chain)
    leq = _reachability(P)
    middles = [v["index"] for v in P.vertices if v["orbit_dim"] == 1]
    a, b = middles
    assert not leq[a][b] and not leq[b][a]


@selfcert
def test_orbit_dim_matches_geometry_invariant():
    # independent oracle: the per-class orbit dim (via the Hom matrix) equals
    # invariants.geometry.orbit_dimension on the assembled module.
    A = linear_path_algebra(3, field=QQ)
    P = degeneration_order(A, {1: 1, 2: 1, 3: 1})
    for v in P.vertices:
        assert v["orbit_dim"] == orbit_dimension(v["module"])


@selfcert
def test_hom_order_is_a_partial_order():
    A = linear_path_algebra(3, field=QQ)
    P = degeneration_order(A, {1: 1, 2: 1, 3: 1})
    leq = _reachability(P)                              # reflexive-transitive closure
    n = len(P.vertices)
    for a in range(n):
        assert leq[a][a]                               # reflexive
    for a in range(n):
        for b in range(n):
            if a != b and leq[a][b]:
                assert not leq[b][a]                   # antisymmetric


@selfcert
def test_orbit_dim_monotone_up_covers():
    A = linear_path_algebra(3, field=QQ)
    P = degeneration_order(A, {1: 1, 2: 1, 3: 1})
    for lo, hi in P.covers:
        assert P.vertices[lo]["orbit_dim"] < P.vertices[hi]["orbit_dim"]


@selfcert
def test_representation_infinite_refused():
    # 2-Kronecker is representation-infinite: knit cannot close -> loud, not a poset.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    P = degeneration_order(A, {1: 2, 2: 2}, budget=40)
    assert P.is_complete is False and P.status in ("budget", "error", "unsupported")
