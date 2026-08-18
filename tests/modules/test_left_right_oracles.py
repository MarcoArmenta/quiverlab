"""Cross-engine + literature oracles for the left/right parts (Plan 55 / R15). Cross-engine:
D R_A = L_{A^op} (ACLV, verbatim) via the shipped opposite + duality, run on the ACLV 2.2(b)
example where R_A is a PROPER subset (4 of 9) so the duality is NON-vacuous; the hereditary
case is a smoke test only (both parts total, match trivial). Literature: gl.dim <= 1 => both
parts total. Honest scope: self-injective and rep-infinite refuse loudly. QQ-scope (M1)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


def _radsq_nakayama_a5():
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


def _dvms(records):
    return sorted(tuple(sorted(r["dimvec"].items())) for r in records)


@xeng
def test_D_of_right_equals_left_of_opposite_on_proper_subset():
    # ACLV: D R_A = L_{A^op}. On the rad^2=0 A5, R_A is a PROPER subset (4 of 9). The check is
    # dim-vector MULTISET equality; it is membership-exact HERE because the compared parts have
    # PAIRWISE-DISTINCT dim-vectors (asserted below), so the multiset pins the exact module set.
    A = _radsq_nakayama_a5()
    atlas_A = left_right_parts(A)
    atlas_op = left_right_parts(A.opposite())
    assert len(atlas_A.right) == 4 and len(atlas_A.right) < atlas_A.universe_size  # non-vacuous
    # membership-exactness precondition: no two modules in the compared parts share a dim-vector
    for part in (atlas_A.right, atlas_A.left, atlas_op.left, atlas_op.right):
        dvs = [tuple(sorted(r["dimvec"].items())) for r in part]
        assert len(dvs) == len(set(dvs)), "dim-vectors not pairwise distinct -- multiset " \
                                          "equality would not be membership-exact"
    assert _dvms(atlas_A.right) == _dvms(atlas_op.left)         # D R_A = L_{A^op}
    assert _dvms(atlas_A.left) == _dvms(atlas_op.right)         # D L_A = R_{A^op}


@xeng
def test_D_duality_smoke_hereditary():
    # hereditary: both parts total; the duality holds trivially -- a smoke case only.
    A = linear_path_algebra(3, field=QQ)
    assert _dvms(left_right_parts(A).right) == _dvms(left_right_parts(A.opposite()).left)


@lit
def test_gldim_le_1_forces_total_parts():
    A = linear_path_algebra(3, field=QQ)            # gl.dim = 1
    assert A.global_dimension().value <= 1
    atlas = left_right_parts(A)
    assert len(atlas.left) == len(atlas.right) == atlas.universe_size
    assert atlas.complement == []


@selfcert
def test_selfinjective_and_rep_infinite_refuse():
    from quiverlab import truncated_polynomial
    assert left_right_parts(truncated_polynomial(4, field=QQ)).status == "unsupported"
    K2 = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    assert left_right_parts(K2, budget=40).is_complete is False
