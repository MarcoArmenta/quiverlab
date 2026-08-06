"""Torsion lattice + bricks/semibricks (Plan 45 / C4). Literature: the AIR four-way count
identity #stau-tilt = #f.f. torsion classes = #semibricks = Catalan(n+1) on kA2/kA3 (the
silting leg is added in Task G). Self-cert: (A,0) is the unique Hasse source, (0,A) the
unique sink; every brick label has end_dim 1; semibricks are pairwise Hom-orthogonal."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import end_dim, hom_dim, is_isomorphic
from quiverlab.tautilting.mutation import exchange_graph
from quiverlab.tautilting.torsion import (bricks, hasse_orientation, semibricks,
                                          torsion_class_data)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kz2_radsq():
    """The symmetric Nakayama algebra kZ2/rad^2: quiver 1<->2 (a:1->2, b:2->1), relations
    a*b, b*a (left-to-right composition, so a*b is the 2-cycle 1->2->1). dim 4,
    tau-tilting-finite, self-injective. Its NON-THIN structure -- P1, P2 are two
    non-isomorphic bricks of the SAME dim-vector (1,1) -- is the regression the first-match
    brick labelling silently collapsed."""
    return Quiver(vertices=[1, 2],
                  arrows={"a": (1, 2), "b": (2, 1)}).algebra(relations=["a*b", "b*a"],
                                                             field=QQ)


@selfcert
def test_hasse_has_unique_source_and_sink():
    A = linear_path_algebra(3, field=QQ)
    eg = exchange_graph(A)
    o = hasse_orientation(eg)
    indeg = {i: 0 for i in range(len(eg.vertices))}
    outdeg = dict(indeg)
    for (i, j), d in o.items():
        a, b = (i, j) if d == "down" else (j, i)     # a -> b downward
        outdeg[a] += 1
        indeg[b] += 1
    sources = [i for i in indeg if indeg[i] == 0]
    sinks = [i for i in outdeg if outdeg[i] == 0]
    assert len(sources) == 1 and len(sinks) == 1
    assert eg.vertices[sources[0]]["is_initial"]      # (A,0) on top


@selfcert
def test_brick_labels_are_bricks():
    A = linear_path_algebra(3, field=QQ)
    for B in bricks(A):
        assert end_dim(B) == 1                         # End(B) = k over QQ (alg closed base)


@selfcert
def test_semibricks_are_hom_orthogonal():
    A = linear_path_algebra(2, field=QQ)
    for sb in semibricks(A):
        mods = list(sb)
        for i in range(len(mods)):
            for j in range(len(mods)):
                if i != j:
                    assert hom_dim(mods[i], mods[j]) == 0    # pairwise Hom-orthogonal


@lit
@pytest.mark.parametrize("n, count", [(2, 5), (3, 14)])
def test_air_four_way_identity_torsion_and_semibricks(n, count):
    # THE AIR four-way identity on ONE run: #support tau-tilting pairs == #f.f. torsion
    # classes == #2-term silting == #semibricks == Catalan(n+1). (Task G adds the silting
    # leg so the full four-way count is pinned here.)
    from quiverlab.tautilting.silting import silting_count
    A = linear_path_algebra(n, field=QQ)
    eg = exchange_graph(A)
    n_pairs = len(eg.vertices)
    n_torsion = len({tuple(torsion_class_data(v["pair"])["gen_dimvecs"])
                     for v in eg.vertices})
    n_silting = silting_count(A)["count"]
    n_semibricks = len(semibricks(A))
    assert n_pairs == n_torsion == n_silting == n_semibricks == count


@xeng
def test_torsion_classes_are_distinct_per_pair():
    # the pair -> Gen(M) map is injective (AIR bijection): fingerprints all distinct.
    A = linear_path_algebra(3, field=QQ)
    eg = exchange_graph(A)
    fps = [tuple(torsion_class_data(v["pair"])["gen_dimvecs"]) for v in eg.vertices]
    assert len(set(fps)) == len(fps)


@lit
def test_kz2_radsq_non_thin_bricks_and_four_way_identity():
    # THE non-thin gate: kZ2/rad^2 has FOUR bricks -- S1, S2 (dim (1,0),(0,1)) and the two
    # non-isomorphic projective-injectives P1, P2 (BOTH dim (1,1)). A dim-vector-keyed brick
    # count silently collapses P1, P2 (returns 3) and merges the semibricks {P1}, {P2}
    # (returns 5), breaking the AIR four-way identity in the user-facing payload. With
    # iso-class identification: 4 bricks (pairwise non-isomorphic, each End = k) and the
    # four-way count identity #pairs = #torsion = #silting = #semibricks = 6 holds.
    from quiverlab.tautilting.block import tau_tilting_block
    A = _kz2_radsq()
    assert A.dim == 4
    br = bricks(A)
    assert len(br) == 4
    assert all(end_dim(B) == 1 for B in br)                       # every brick End_A(B)=k
    for i in range(len(br)):
        for j in range(i + 1, len(br)):
            assert not is_isomorphic(br[i], br[j])               # pairwise non-isomorphic
    assert len(semibricks(A)) == 6
    counts = tau_tilting_block(A)["counts"]
    assert counts == {"pairs": 6, "torsion": 6, "silting": 6, "semibricks": 6}


@selfcert
def test_kz2_radsq_same_dimvec_walls_carry_distinct_bricks():
    # Regression for the wall mislabelling: the two exchange walls whose King normal is the
    # dim-vector (1,1) must carry NON-isomorphic brick labels (one P1, one P2) -- a first
    # dim-vector match gives both the same module. Across ALL edges the labels realise
    # exactly the four distinct brick iso-classes.
    from quiverlab.tautilting.torsion import _edge_brick, _torsion_universe
    A = _kz2_radsq()
    eg = exchange_graph(A)
    uni = _torsion_universe(A)

    def label(e):
        return _edge_brick(A, eg.vertices[e[0]]["pair"], eg.vertices[e[1]]["pair"],
                           eg.arrows[e]["brick"], uni)

    walls_11 = [e for e in eg.arrows if eg.arrows[e]["brick"] == {1: 1, 2: 1}]
    assert len(walls_11) == 2
    b0, b1 = label(walls_11[0]), label(walls_11[1])
    assert not is_isomorphic(b0, b1)                              # distinct wall labels
    distinct = []
    for e in eg.arrows:
        B = label(e)
        if not any(is_isomorphic(B, C) for C in distinct):
            distinct.append(B)
    assert len(distinct) == 4                                     # 4 iso-classes on the walls


@xeng
def test_kz2_radsq_matches_over_gf2():
    # MEDIUM-1 honest contract: the engine does NOT proactively refuse over char <= dim. On
    # kZ2/rad^2 every module involved is a brick or splits, so the GF(2) run (char 2 <= dim 4)
    # computes and AGREES with the certified QQ result -- it refuses only where a
    # decompose/is_isomorphic certificate is actually needed and unavailable, never up front.
    from quiverlab.fields import GF
    from quiverlab.tautilting.block import tau_tilting_block
    A_gf2 = Quiver(vertices=[1, 2],
                   arrows={"a": (1, 2), "b": (2, 1)}).algebra(relations=["a*b", "b*a"],
                                                              field=GF(2))
    counts = tau_tilting_block(A_gf2)["counts"]
    assert counts == {"pairs": 6, "torsion": 6, "silting": 6, "semibricks": 6}
