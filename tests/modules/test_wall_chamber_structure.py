"""The wall-and-chamber structure via bricks (Plan 63 / R25). Literature: kA2 = 5 chambers
/ 3 walls (the P1 wall a ray); DIJ brick-finite <=> tau-tilting-finite gate. Cross-engine:
#chambers = #support tau-tilting. Self-cert: each exchange-edge facet is shared by exactly
two chambers; every wall is coplanar (theta.dim B = 0 on its facets)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_kA2_five_chambers_three_walls():
    A = _kA2()
    wc = A.wall_chamber_structure()
    assert wc["complete"] and wc["status"] == "complete"
    assert wc["num_chambers"] == 5 and wc["num_walls"] == 3     # 3 bricks -> 3 walls
    assert wc["render"] == "fan2d"
    # exactly one wall is a proper ray (D(P1)); the other two are full lines (D(S1),D(S2))
    rays = [w for w in wc["walls"] if not w["is_full_hyperplane"]]
    lines = [w for w in wc["walls"] if w["is_full_hyperplane"]]
    assert len(rays) == 1 and len(lines) == 2
    assert sorted(rays[0]["brick_dimvec"].items()) == [(1, 1), (2, 1)]   # the P1 ray


@xeng
@pytest.mark.parametrize("n, chambers", [(2, 5), (3, 14)])   # HEREDITARY kA_n: Catalan(n+1)
def test_chamber_count_equals_support_tau_tilting(n, chambers):
    # linear_path_algebra(n) is HEREDITARY kA_n (no relations) -> Catalan(n+1) chambers. This
    # is a DIFFERENT algebra/fact from the rad^2-Nakayamas (linear = 12, cyclic = 14 chambers).
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(n, field=QQ)
    wc = A.wall_chamber_structure()
    eg = exchange_graph(A)
    assert wc["complete"]
    assert wc["num_chambers"] == len(eg.vertices) == chambers
    assert wc["counts"]["chambers"] == wc["counts"]["s_tau_tilt"] == chambers


@selfcert
def test_each_facet_bounds_exactly_two_chambers():
    # the fan property, FALSIFIABLY (W2): the grouped facets are in BIJECTION with the
    # exchange graph's edges -- every exchange edge appears as exactly one facet under exactly
    # one wall, every facet is a real edge between two DISTINCT chambers. A grouping that
    # dropped/duplicated an edge, or invented a non-edge facet, fails here (the old version
    # only checked i != j, which is tautological for an exchange edge).
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(3, field=QQ)
    wc = A.wall_chamber_structure()
    eg = exchange_graph(A)
    edges = {tuple(sorted(e)) for e in eg.arrows}           # eg.arrows keyed by (i<j)
    facets = [tuple(sorted(f)) for w in wc["walls"] for f in w["facets"]]
    for (i, j) in facets:
        assert i != j and 0 <= i < wc["num_chambers"] and 0 <= j < wc["num_chambers"]
    assert len(facets) == len(set(facets)) == len(edges)    # bijection: no drop, no dup
    assert set(facets) == edges                             # every facet is a real edge


@selfcert
def test_walls_are_coplanar_theta_dot_dimB_zero_on_facets():
    # every facet grouped under a wall lies in the hyperplane theta . dim B = 0: the facet
    # vector (a shared g-vector) is orthogonal to the brick dim-vector.
    A = _kA2()
    wc = A.wall_chamber_structure()
    verts = [1, 2]
    for w in wc["walls"]:
        nrm = [w["brick_dimvec"][v] for v in verts]
        # the wall's own extreme rays (n<=2) satisfy the equality + all inequalities
        for r in (w["rays"] or []):
            from fractions import Fraction
            theta = [Fraction(x) for x in r]
            assert sum(theta[k] * nrm[k] for k in range(2)) == 0


@lit
def test_brick_infinite_is_bounded_region_not_complete():
    # the 2-Kronecker is tau-tilting-infinite / brick-infinite (DIJ): the structure cannot
    # close -> honest truncation, a bounded region, NO counts. Budget 10 (not the plan's 40):
    # the 2-Kronecker's preprojectives grow without bound, so each mutation on the larger
    # modules is costlier -- exchange_graph is ~6s at budget 10 but ~47s at budget 15 (measured),
    # super-linear; budget 40 would run for hours. Budget 10 trips the cap identically and is
    # the honest bounded-region witness. See the plan's "adjust to reality" (Task 2).
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    wc = K.wall_chamber_structure(budget_pairs=10)
    assert wc["complete"] is False and wc["status"] == "budget"
    assert wc["counts"] is None and wc["green_count"] is None
    assert wc["truncation"]                          # a non-empty honest note
    assert wc["num_chambers"] > 0                    # but the bounded region IS shown


def _kZ2_radsq():
    # self-injective Nakayama N2^2: 1 <-> 2, BOTH length-2 paths killed (rad^2 = 0).
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


@lit
def test_kZ2_radsq_four_walls_two_share_a_hyperplane():
    # THE primary discriminating oracle (BST Remark 3.19). kZ2/rad^2 is self-injective with 6
    # support tau-tilting pairs and 4 bricks: S1=(1,0), S2=(0,1) simple, and TWO
    # non-isomorphic bricks P1, P2 BOTH of dim-vector (1,1) with submodule-dimvec sets
    # {(0,0),(1,0),(1,1)} and {(0,0),(0,1),(1,1)} -> D(P1), D(P2) are OPPOSITE half-rays of
    # the ONE line theta_1 + theta_2 = 0. A grouping keyed on dim-vector (up to scaling) would
    # MERGE them into 3 walls; the correct ISO-CLASS grouping keeps 4. (Engine-verified in the
    # fix round: 6 pairs, 4 bricks, the two submodule sets above.)
    from fractions import Fraction
    A = _kZ2_radsq()
    wc = A.wall_chamber_structure()
    assert wc["complete"] and wc["num_chambers"] == 6
    assert wc["num_walls"] == 4                       # NOT 3 -- P1, P2 stay separate
    # exactly two walls carry the dim-vector (1,1); both lie in the ONE line theta_1+theta_2=0
    dim11 = [w for w in wc["walls"] if sorted(w["brick_dimvec"].items()) == [(1, 1), (2, 1)]]
    assert len(dim11) == 2
    for w in dim11:
        assert list(w["equality"]) == [1, 1]         # the shared hyperplane normal
        assert w["is_full_hyperplane"] is False      # each is a proper RAY (half of the line)
        assert len(w["rays"]) == 1                    # one extreme ray each (a half-ray)
    # their extreme rays are the two OPPOSITE directions (1,-1) and (-1,1), NOT collapsed
    rays = {tuple(int(Fraction(x)) for x in w["rays"][0]) for w in dim11}
    assert rays == {(1, -1), (-1, 1)}                 # distinct opposite half-rays


@selfcert
@pytest.mark.parametrize("build", [
    _kA2,                                             # thin: kA2, 3 bricks (distinct dimvecs)
    _kZ2_radsq,                                       # non-thin: kZ2/rad^2, 4 bricks (P1,P2)
])
def test_counts_consistency_walls_equal_bricks(build):
    # M2 counts arbiter. ONE canonical brick set: torsion.bricks(A) -- the shipped
    # ISO-CLASS-deduped enumerator (chosen over an inline dim-vector dedup precisely because
    # it disambiguates same-dim-vector bricks by torsion-class membership, torsion.py:142, so
    # kZ2/rad^2's P1,P2 count as 2). In the complete case #walls == #bricks, and the walls'
    # brick dim-vector MULTISET matches torsion.bricks(A)'s (on kZ2 the (1,1) key has
    # multiplicity 2 -- a dim-vector-merged grouping would show 1 here and fail).
    from collections import Counter
    from quiverlab.tautilting.torsion import bricks as torsion_bricks
    A = build()
    wc = A.wall_chamber_structure()
    assert wc["complete"]
    ref = torsion_bricks(A)
    assert wc["counts"]["walls"] == wc["counts"]["bricks"] == len(ref) == wc["num_walls"]
    wall_dvs = Counter(tuple(sorted(w["brick_dimvec"].items())) for w in wc["walls"])
    ref_dvs = Counter(tuple(sorted(B.dimension_vector().items())) for B in ref)
    assert wall_dvs == ref_dvs


# --------------------------------------------------------------------------- #
# Plan 63 fix (adjudicated 2026-08-07): status="error" is NEVER read as tau-tilting-infinite.
# --------------------------------------------------------------------------- #
@lit
def test_d4_star_recovered_complete_despite_mutate_error():
    # The D4 star {1->2, 3->1, 4->1} over QQ is HEREDITARY of underlying type D4, hence
    # tau-tilting FINITE: Cat(D4) = 50 support tau-tilting modules (50 chambers) and 12 bricks
    # = 12 positive roots = 12 walls. The shipped exchange-graph BFS trips a spurious P45
    # mutate error (status="error"), but every one of the 50 discovered chambers has exactly
    # n = 4 neighbours, so the graph is CLOSED (each edge is rediscovered from its other
    # endpoint) and the structure is RECOVERED as complete -- NEVER the false "tau-tilting-
    # infinite" truncation. This pins the recovered OUTCOME, so it survives Plan 65's root fix
    # of mutate (which will simply report status="complete" and skip the certificate).
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 1), "c": (4, 1)}).algebra(
        relations=[], field=QQ)
    wc = A.wall_chamber_structure()
    # OUTCOME (survives P65): complete, the Catalan count, walls == bricks, no false claim.
    assert wc["complete"] is True
    assert wc["num_chambers"] == 50                       # Cat(D4)
    assert wc["counts"]["chambers"] == wc["counts"]["s_tau_tilt"] == 50
    assert wc["num_walls"] == wc["counts"]["walls"] == wc["counts"]["bricks"] == 12
    assert wc["truncation"] is None                       # NOT the false infinite claim
    assert "tau-tilting-infinite" not in (wc.get("truncation") or "")
    # MECHANISM (today only): the recovery fired, so the honest n-regularity note is present and
    # the green count is omitted rather than a false 0. After Plan 65 fixes mutate, status is
    # "complete", there is no note, and green_count is the real integer -- the outcome holds.
    if wc["status"] == "error":
        assert wc.get("note") and "n-regularity" in wc["note"]
        assert wc["green_count"] is None


@selfcert
def test_error_graph_not_n_regular_refuses_loudly(monkeypatch):
    # ruling 1b: a status="error" exchange graph that is NOT n-regular cannot be certified
    # complete -- the structure is UNRELIABLE and must be refused LOUDLY, never dressed up as
    # the tau-tilting-infinite truncation or a bounded sub-fan. Injected via a fake exchange
    # graph (both chambers at degree 1 != n = 2) so no ~20s BFS is spent.
    import quiverlab.tautilting.mutation as mut
    from quiverlab.errors import QuiverlabError

    class _FakeEG:
        is_complete = False
        status = "error"
        vertices = [
            {"g_matrix": [[1, 0], [0, 1]], "label": "(A,0)", "support": [1, 2],
             "is_initial": True},
            {"g_matrix": [[-1, 0], [0, 1]], "label": "x", "support": [2], "is_initial": False},
        ]
        adj = {0: {1}, 1: {0}}                            # both degree 1 != n = 2: NOT n-regular
        arrows = {(0, 1): {"brick": {1: 1, 2: 0}, "brick_name": None}}

    monkeypatch.setattr(mut, "exchange_graph", lambda A, budget_pairs=512: _FakeEG())
    with pytest.raises(QuiverlabError) as ei:
        _kA2().wall_chamber_structure()                   # n = 2
    msg = str(ei.value)
    assert "unreliable" in msg
    assert "tau-tilting-infinite" not in msg and "brick-infinite" not in msg


@selfcert
def test_closed_by_n_regularity_certificate():
    # the certificate itself, FALSIFIABLY: n-regular closed graph -> True; any vertex short of
    # n neighbours -> False. (Recomputed from adj because eg.n_regular is gated on is_complete
    # in mutation.py and is False on an error graph by construction.)
    from quiverlab.tautilting.wallchamber import _closed_by_n_regularity

    class _EG:
        def __init__(self, adj):
            self.adj = adj
            self.vertices = list(range(len(adj)))

    assert _closed_by_n_regularity(_EG({0: {1, 2}, 1: {0, 2}, 2: {0, 1}}), 2) is True
    assert _closed_by_n_regularity(_EG({0: {1}, 1: {0}}), 2) is False   # degree 1 < n = 2


@selfcert
def test_budget_partial_num_walls_is_none_with_group_count():
    # ruling 3: a GENUINE budget cap must NOT present its group-by-normal count as a definitive
    # num_walls (the true wall count is unknown / unbounded). num_walls is None; the discovered
    # wall-normal groups are counted honestly under partial_wall_groups. 2-Kronecker, budget 8
    # (super-linear BFS -> a small budget keeps this fast).
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    wc = K.wall_chamber_structure(budget_pairs=8)
    assert wc["complete"] is False and wc["status"] == "budget"
    assert wc["num_walls"] is None
    assert isinstance(wc["partial_wall_groups"], int) and wc["partial_wall_groups"] > 0
    assert wc["counts"] is None and wc["green_count"] is None and wc["truncation"]


@selfcert
def test_n1_degenerate_wall_and_zero_walls():
    # ruling 4 / commit 7e5b59e: for a single-vertex algebra the stability space is R and
    # D(B) = {theta : theta . dim B = 0} = {0} -- a DEGENERATE point, NOT the full hyperplane.
    # wall_of_brick reports is_full_hyperplane False + no rays; k[x]/(x^2) has no exchange-edge
    # brick (End of the regular module is not k), so the structure has num_walls 0.
    from quiverlab.tautilting.wallchamber import wall_of_brick
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=QQ)
    w = wall_of_brick(A.simple(1))
    assert w.is_full_hyperplane is False
    assert w.rays == ()
    wc = A.wall_chamber_structure()
    assert wc["complete"] and wc["num_walls"] == 0
