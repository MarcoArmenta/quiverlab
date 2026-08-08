"""Left/right parts L_A, R_A via the closed-under-predecessors pd/id sweep (Plan 55 / R15).
Literature: over a hereditary algebra (pd/id <= 1 everywhere) L_A = R_A = ind A, empty
complement; the ACLV Example 2.2(b) rad^2=0 linear Nakayama splits as L_A={S1,S2,P2,P3},
R_A={S4,S5,P4,P5}, complement={S3} (pd 2, id 2) -- ada with NON-empty complement. Self-cert:
L_A closed under predecessors, intersection/complement consistent. Cross-engine: the Hom
closure and AR-quiver reachability predecessor relations coincide (rad^infty = 0)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _radsq_nakayama_a5():
    # ACLV Example 2.2(b): 1 <- 2 <- 3 <- 4 <- 5 bound by rad^2 = 0.
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


@lit
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_hereditary_both_parts_are_total(n, count):
    A = linear_path_algebra(n, field=QQ)          # kA_n hereditary, pd/id <= 1 everywhere
    atlas = left_right_parts(A)
    assert atlas.is_complete and atlas.status == "complete"
    assert atlas.universe_size == count == n * (n + 1) // 2
    assert len(atlas.left) == len(atlas.right) == count      # L_A = R_A = ind A
    assert atlas.complement == []                            # hereditary => empty complement


@lit
def test_aclv_example_22b_complement_is_S3():
    A = _radsq_nakayama_a5()
    atlas = left_right_parts(A)
    assert atlas.is_complete and atlas.universe_size == 9
    assert A.global_dimension().value == 4                   # sharp (ACLV Rem 2.7a)
    # S3 (dim-vector {3:1}) is in NEITHER part -- the ada-with-non-empty-complement datum.
    comp = [r["dimvec"] for r in atlas.complement]
    assert comp == [{1: 0, 2: 0, 3: 1, 4: 0, 5: 0}]          # exactly S3
    assert {r["name"] for r in atlas.complement} == {"S_3"}
    left_names = {r["name"] for r in atlas.left}
    right_names = {r["name"] for r in atlas.right}
    assert "S_3" not in left_names and "S_3" not in right_names   # positive: S3 in neither
    assert len(atlas.left) == 4 and len(atlas.right) == 4        # proper subsets (of 9)
    # placement (up to the projective/injective side convention -- see note):
    assert atlas.projective_placement[1] in ("L", "both")
    assert atlas.projective_placement[3] in ("L", "both")
    assert atlas.projective_placement[4] in ("R", "both")
    assert atlas.projective_placement[5] in ("R", "both")


@selfcert
def test_intersection_and_complement_consistent():
    A = _radsq_nakayama_a5()
    atlas = left_right_parts(A)
    L = {r["index"] for r in atlas.left}
    R = {r["index"] for r in atlas.right}
    assert {r["index"] for r in atlas.intersection} == (L & R)   # here == empty set
    assert {r["index"] for r in atlas.complement} == set(range(atlas.universe_size)) - (L | R)


@selfcert
def test_left_part_closed_under_predecessors():
    A = linear_path_algebra(4, field=QQ)
    atlas = left_right_parts(A)
    leq = atlas._leq
    left = {r["index"] for r in atlas.left}
    for j in left:
        for i in range(atlas.universe_size):
            if leq[i][j]:
                assert i in left


@xeng
def test_predecessor_relation_two_routes_agree():
    # Hom-nonzero transitive closure == AR-quiver reachability (rad^infty = 0, rep-finite).
    from quiverlab.modules.left_right import _ar_reachability, _leq_matrix, _universe
    A = _radsq_nakayama_a5()
    ar, U, _names, _recs = _universe(A, 256, 4096)
    assert ar.is_complete
    assert _leq_matrix(U) == _ar_reachability(ar, len(U))


@selfcert
def test_selfinjective_refused_loudly():
    from quiverlab import truncated_polynomial
    A = truncated_polynomial(3, field=QQ)                    # k[x]/(x^3), self-injective
    atlas = left_right_parts(A)
    assert atlas.is_complete is False and atlas.status == "unsupported"
    assert atlas.left == [] and atlas.right == []            # never a partial atlas


@selfcert
def test_rep_infinite_refused_loudly():
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    atlas = left_right_parts(A, budget=40)                   # 2-Kronecker, rep-infinite
    assert atlas.is_complete is False and atlas.status in ("budget", "error", "unsupported")


@selfcert
def test_aclv_22c_nonhereditary_rep_infinite_refused_fast():
    # ACLV Example 2.2(c): 1 => 2 => 3 => 4 (double arrows) bound by rad^2 = 0 -- a
    # MATHEMATICALLY ada but REPRESENTATION-INFINITE algebra (the plan's named refusal
    # oracle, acceptance #6). It is NON-hereditary, so the hereditary fast guard misses it;
    # the knit's per-module almost-split cost makes even budget_dim=16 take ~120s (measured),
    # so a fast SUFFICIENT rep-infinite certificate is required -- Gabriel's separated-quiver
    # criterion (rad^2=0, separated quiver not a disjoint union of Dynkin) fires INSTANTLY.
    import time
    Q = Quiver([1, 2, 3, 4],
               {"a1": (1, 2), "b1": (1, 2), "a2": (2, 3), "b2": (2, 3),
                "a3": (3, 4), "b3": (3, 4)})
    A = RadicalSquareZero(Q, field=QQ)
    assert A.is_hereditary() is False               # bypasses the 2-Kronecker (hereditary) guard
    t = time.time()
    atlas = left_right_parts(A)
    elapsed = time.time() - t
    assert atlas.is_complete is False               # a LOUD bounded refusal, never a partial atlas
    assert atlas.status in ("budget", "unsupported")
    assert atlas.left == [] and atlas.right == [] and atlas.complement == []
    assert elapsed < 15.0, f"2.2(c) refusal took {elapsed:.1f}s (must be fast, not the ~120s knit)"


@selfcert
def test_delegate_matches_free_function():
    A = linear_path_algebra(2, field=QQ)
    assert [r["index"] for r in A.left_right_parts().left] == \
           [r["index"] for r in left_right_parts(A).left]


@selfcert
def test_addendum_pd_id_vectors_index_aligned():
    # P61 contract: pd_le_1[i] / id_le_1[i] index-aligned with _modules; assigned from the
    # same pd/id sweep. On the rad^2=0 A5 exactly S3 has neither pd<=1 nor id<=1.
    A = _radsq_nakayama_a5()
    atlas = left_right_parts(A)
    assert len(atlas.pd_le_1) == len(atlas.id_le_1) == atlas.universe_size == len(atlas._modules)
    # every X in L_A has pd<=1 (closed-under-predecessors => X itself has pd<=1)
    for r in atlas.left:
        assert atlas.pd_le_1[r["index"]] is True
    for r in atlas.right:
        assert atlas.id_le_1[r["index"]] is True
    # S3 is the unique indecomposable with pd>=2 AND id>=2
    (s3,) = [r["index"] for r in atlas.complement]
    assert atlas.pd_le_1[s3] is False and atlas.id_le_1[s3] is False
