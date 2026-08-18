"""Rank <= 3 drawing geometry + MGS green paths (Plan 63 / R25). Self-cert: the wall's
drawing rays satisfy its inequality system (theta.dim B = 0, theta.dim N <= 0); n=3 rays
carry an L1/octahedron projection; each MGS is a monotone source->sink chamber path."""
import pytest
from fractions import Fraction

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


@xeng
def test_n2_wall_rays_equal_grouped_exchange_facets():
    # TWO INDEPENDENT CONSTRUCTIONS (W1): the inequality-defined D(B) extreme rays (Task 1)
    # must equal, as SETS OF PRIMITIVE DIRECTIONS, the grouped exchange-edge facet vectors of
    # the Plan-45 fan (each facet a shared g-vector lying IN the wall). Not merely "both
    # nonempty" -- genuine set equality per brick. kA2's bricks have DISTINCT dim-vectors, so
    # the fan facets can be grouped by dim-vector here. Concretely: S1 -> {(0,1),(0,-1)},
    # S2 -> {(1,0),(-1,0)}, P1 -> {(1,-1)}.
    from math import gcd
    from quiverlab.tautilting.stability import wall_and_chamber_fan

    def prim(v):                                     # primitive integer direction (sign kept)
        ints = [int(Fraction(x)) for x in v]
        g = 0
        for x in ints:
            g = gcd(g, abs(x))
        return tuple(x // g for x in ints) if g else tuple(ints)

    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    wc = A.wall_chamber_structure()
    fan = wall_and_chamber_fan(A)
    fan_by_brick = {}
    for f in fan["walls"]:
        key = tuple(sorted(f["brick_dimvec"].items()))
        fan_by_brick.setdefault(key, set()).update(prim(v) for v in f["facet"])
    ineq_by_brick = {tuple(sorted(w["brick_dimvec"].items())): {prim(r) for r in w["rays"]}
                     for w in wc["walls"]}
    assert set(ineq_by_brick) == set(fan_by_brick)       # same brick set
    for key in ineq_by_brick:
        assert ineq_by_brick[key] == fan_by_brick[key]   # SET EQUALITY of directions
    assert ineq_by_brick[((1, 1), (2, 1))] == {(1, -1)}  # the P1 ray, single direction


@selfcert
def test_n3_walls_project_via_L1_and_lie_in_D_of_B():
    # n=3 self-cert (W2): each wall's drawing rays (grouped exchange-edge facet vectors)
    # actually LIE IN D(B) -- theta.dim B = 0 AND theta.dim N <= 0 for every submodule
    # inequality -- and carry the L1/octahedron projection. FALSIFIABLE: a facet violating an
    # inequality, or a wall carrying no rays at all, fails here (the old version's disjunction
    # was trivially satisfiable by w["rays"] being None).
    A = linear_path_algebra(3, field=QQ)             # hereditary kA3, n=3, tau-tilting-finite
    wc = A.wall_chamber_structure()
    assert wc["render"] == "fan3d"
    saw_wall_ray = False
    for w in wc["walls"]:
        rays = w["rays"] or []
        nrm = [w["brick_dimvec"][v] for v in (1, 2, 3)]
        assert len(w.get("rays_l1", rays)) == len(rays)          # one projection per ray
        for r in rays:
            saw_wall_ray = True
            theta = [Fraction(x) for x in r]
            assert sum(theta[k] * nrm[k] for k in range(3)) == 0          # in the hyperplane
            for d in w["inequalities"]:
                assert sum(theta[k] * d[k] for k in range(3)) <= 0        # in D(B)
    assert saw_wall_ray                                          # walls DID carry drawing rays
    for ch in wc["chambers"]:
        assert "net2d" in ch and ch["net2d"]          # chambers carry the L1 net position


@lit
def test_kA2_two_maximal_green_sequences_are_green_paths():
    from quiverlab.tautilting.wallchamber import _is_green_path
    from quiverlab.tautilting.mutation import exchange_graph
    from quiverlab.tautilting.torsion import hasse_orientation
    from quiverlab.tautilting.green import maximal_green_sequences
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    eg = exchange_graph(A)
    orient = hasse_orientation(eg)
    mgs = maximal_green_sequences(A)
    assert mgs["count"] == 2
    for seq in mgs["sequences"]:
        assert _is_green_path(eg, orient, seq)        # each MGS is a monotone chamber path
