# SPDX-License-Identifier: MIT
"""Cluster-tilting objects ARE support tau-tilting pairs (Adachi-Iyama-Reiten), and the
certification of the COUNT is a per-input RUNTIME property (Plan 79 / R31, Task 2).

MEASURED RUNTIMES on this tree (single-threaded, no contention noted otherwise):
A2 ~0s, A3 ~0.4s, A4 ~9s, D4 ~20s, A5 ~160s, D5 ~842s (under contention; the plan doc
measured ~284s). The A5/D5 cases are DEEP-bucket by directory and marked `slow` so the
routine deep run is not held hostage to a 14-minute BFS.
"""
import pytest

from quiverlab import Quiver
from quiverlab.cluster import ClusterCategory
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


def _catalan(m):
    num = den = 1
    for i in range(m):
        num *= (2 * m - i)
        den *= (i + 1)
    return num // den // (m + 1)


@lit
@xeng
@pytest.mark.parametrize("typ,n,expected", [("A2", 2, 5), ("A3", 3, 14), ("A4", 4, 42)])
def test_type_A_cluster_number_is_catalan(typ, n, expected):
    # THE cross-engine backbone: the P45 exchange graph is built by tau-tilting mutation
    # and knows nothing about cluster combinatorics, yet by AIR its vertex count IS the
    # cluster number -- Catalan(n+1) in type A.
    info = ClusterCategory(typ).num_cluster_tilting()
    assert info["count"] == expected == _catalan(n + 1)
    assert info["certified"] is True and info["status"] == "complete"


@lit
@selfcert
def test_D4_is_50_and_is_CERTIFIED_BY_THE_N_REGULARITY_RECOVERY():
    # Buan-Marsh: the D4 cluster number is 50. The BFS still reports status='error' (a
    # mutate boundary), so the count is certified ONLY through the Plan-63 n-regularity
    # recovery -- the recovery is load-bearing here, not decoration.
    info = ClusterCategory("D4").num_cluster_tilting()
    assert info["count"] == 50 and info["certified"] is True
    assert info["status"] == "error"
    assert "n-regular" in info["note"] and "Plan-63" in info["note"]


@xeng
@pytest.mark.parametrize("typ", ["A3", "D4"])
def test_the_count_equals_the_exchange_graph_vertex_count(typ):
    from quiverlab.tautilting.mutation import exchange_graph
    C = ClusterCategory(typ)
    info = C.num_cluster_tilting()
    eg = exchange_graph(C.algebra)
    assert info["count"] == len(eg.vertices)


@selfcert
@pytest.mark.parametrize("typ", ["A3", "D4"])
def test_every_cluster_tilting_object_has_exactly_n_summands(typ):
    # T = M (+) P[1] is basic with n indecomposable summands: |M| + |support| = n, the
    # AIR full-rank axiom read through the cluster category.
    C = ClusterCategory(typ)
    objs = C.cluster_tilting_objects()
    assert len(objs) == C.num_cluster_tilting()["count"]
    for o in objs:
        assert o["rank"] == C.n
        assert len(o["summands"]) + len(o["support"]) == C.n


@selfcert
def test_G3a_over_budget_is_NOT_an_infiniteness_claim():
    # THE split that protects a rep-finite algebra from being slandered as infinite.
    # A5 is Dynkin (rep-FINITE) with cluster number 132; with a tiny budget the BFS stops
    # at 'budget', and the refusal must say OVER-BUDGET and quote 132 -- never "infinite".
    info = ClusterCategory("A5").num_cluster_tilting(budget_pairs=8)
    assert info["status"] == "budget"
    assert info["certified"] is False and info["count"] is None
    assert info["dynkin_type"] == "A5"
    assert info["expected_count"] == 132
    assert "132" in info["note"] and "budget" in info["note"].lower()
    # The note may DENY infiniteness ("NOT an infiniteness claim") -- what it must never
    # do is ASSERT it. G3b's affirmative wording is "infinitely many", so that exact
    # claim is what this forbids; a blunt "infinit" substring check would fail on the
    # denial itself and would have to be weakened, losing the point.
    assert "infinitely many" not in info["note"]
    assert "NOT an infiniteness claim" in info["note"]


@selfcert
def test_G3b_an_affine_quiver_DOES_claim_infiniteness():
    # The contrast case: the Kronecker ~A_1 is representation-INFINITE, and there the
    # claim is made -- off the Gabriel/Dynkin type certificate, not off the budget stop.
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    info = ClusterCategory(Q, field=QQ).num_cluster_tilting(budget_pairs=8)
    assert info["certified"] is False and info["count"] is None
    assert "infinitely many" in info["note"] and "Gabriel" in info["note"]


@selfcert
def test_an_uncertified_count_never_leaks_the_discovered_number():
    # The whole point of the gate: on ANY refusal `count` is None. The literature value
    # lives in `expected_count` and is explicitly NOT what the engine certified.
    info = ClusterCategory("A5").num_cluster_tilting(budget_pairs=8)
    assert info["count"] is None
    assert info["expected_count"] == 132          # quoted to size the budget, not a result
    with pytest.raises(QuiverlabError, match="not certified"):
        ClusterCategory("A5").cluster_tilting_objects(budget_pairs=8)


@selfcert
def test_the_refusal_never_hands_out_a_partial_enumeration():
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    C = ClusterCategory(Q, field=QQ)
    with pytest.raises(QuiverlabError):
        C.cluster_tilting_objects(budget_pairs=8)


@lit
@xeng
@pytest.mark.slow
def test_A5_cluster_number_is_132():
    # ~160 s: deep + slow (opt-in), so the routine deep run is not held up.
    info = ClusterCategory("A5").num_cluster_tilting()
    assert info["count"] == 132 == _catalan(6)
    assert info["certified"] is True


@lit
@selfcert
@pytest.mark.slow
def test_D5_is_CERTIFIED_OR_REFUSED_and_never_in_between():
    """The plan's both-regimes gate -- and a PLAN-DOC CORRECTION.

    The plan doc predicted TWO possible outcomes for D5: pre-P65 a refusal
    (`certified False`, count None, because the discovered graph was NOT 5-regular), or
    post-P65 a native `status == "complete"` with count 182. MEASURED on this tree the
    outcome is a THIRD regime the doc did not anticipate: `status == "error"` (the mutate
    boundary still bites) but the discovered graph IS 5-regular, so the Plan-63 recovery
    certifies count == 182 -- the Buan-Marsh cluster number for D5.

    So the assertion is written on the INVARIANT rather than on a predicted regime: the
    count is either certified and equal to the literature value, or it is refused with
    count None -- and a refusal must never be an infiniteness claim (D5 is Dynkin, hence
    representation-finite). ~842 s under contention.
    """
    info = ClusterCategory("D5").num_cluster_tilting()
    if info["certified"]:
        assert info["count"] == 182               # Buan-Marsh cluster number, type D5
        assert info["status"] in ("complete", "error")
    else:
        assert info["count"] is None
        assert "infinitely many" not in info["note"]
