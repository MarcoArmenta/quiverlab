"""Consistency + literature oracles (Plan 63 / R25). #chambers = #support tau-tilting on
the tau-tilting-finite zoo (cross-engine); the Kaipel-Treffinger Ex. 13 (kA2) + Ex. 15
(CYCLIC rad^2-Nakayama N3^2, 14 chambers) worked examples (lit); the linear kA3/rad^2 = 12
cross-engine value (NOT KT); the DIJ brick-finite gate."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ

xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


@xeng
@pytest.mark.parametrize("A", [
    linear_path_algebra(2, field=QQ),                                    # hereditary kA2 (5)
    linear_path_algebra(3, field=QQ),                                    # hereditary kA3 (14)
    Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(              # LINEAR rad^2-Nakayama
        relations=["a*b"], field=QQ),                                    # 12 pairs (NOT KT Ex.15)
])
def test_chamber_count_equals_exchange_graph(A):
    from quiverlab.tautilting.mutation import exchange_graph
    wc = A.wall_chamber_structure()
    if wc["complete"]:
        assert wc["num_chambers"] == len(exchange_graph(A).vertices)


@xeng
def test_linear_kA3_radsq_is_12_not_14():
    # H1 pin: the LINEAR 1->2->3 mod rad^2 (radical-square-zero Nakayama) has 12 support
    # tau-tilting pairs -- NOT 14. Engine-verified in the fix round (12 pairs, 5 bricks). This
    # is an honest cross-engine consistency value, NOT attributed to KT (KT Example 15 is the
    # CYCLIC N3^2, tested separately). Guards against the earlier mis-transcription that
    # pinned this algebra at 14.
    from quiverlab.tautilting.mutation import exchange_graph
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ)
    wc = A.wall_chamber_structure()
    assert wc["complete"]
    assert wc["num_chambers"] == len(exchange_graph(A).vertices) == 12   # NOT 14


@lit
def test_kaipel_treffinger_example_13_kA2_walls():
    # KT arXiv:2302.12699 Example 13 (VERBATIM): D(S1)={(0,y)} full line, D(S2)={(x,0)}
    # full line, D(P1)={(x,-x): x>=0} the RAY. 3 walls, 5 chambers.
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    wc = A.wall_chamber_structure()
    assert wc["num_chambers"] == 5 and wc["num_walls"] == 3
    rays = [w for w in wc["walls"] if not w["is_full_hyperplane"]]
    assert len(rays) == 1 and sorted(rays[0]["brick_dimvec"].items()) == [(1, 1), (2, 1)]


@lit
def test_kaipel_treffinger_example_15_cyclic_radsq_nakayama_14_chambers():
    # KT Example 15 is the CYCLIC radical-square-zero Nakayama N3^2 (Q: 1->2->3->1 mod rad^2,
    # self-injective): its module list contains 3/1 (top S3, socle S1), which needs the arrow
    # 3->1, so KT's algebra is the 3-CYCLE -- NOT the linear 1->2->3 (that one has 12 chambers,
    # see test_linear_kA3_radsq_is_12_not_14). Engine-verified in the fix round: 14 pairs
    # (complete, n_regular), 6 bricks. This test ALSO covers self-injective input (M3). A
    # rank-3 oracle (render fan3d). (NOT "= Catalan(4)": a numerical coincidence, not the
    # hereditary-kA3 Catalan fact.)
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}).algebra(
        relations=["a*b", "b*c", "c*a"], field=QQ)   # 3-cycle, all length-2 paths = 0 (rad^2)
    wc = A.wall_chamber_structure()
    assert wc["complete"] and wc["num_chambers"] == 14 and wc["render"] == "fan3d"
