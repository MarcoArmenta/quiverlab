"""Plan 61 literature + honest-scope oracle battery (beyond the Task A/B pins).

  * strict weakly shod => HH^i(A) = 0 for i >= 2 (survey S5.1 [43]): asserted on the rad^2=0
    A5 fixture, which IS strict (weakly shod True, quasi-tilted False). The survey warns
    quasi-tilted algebras CAN have HH^2 != 0, so this is pinned only on the strict member.
  * ada => gl.dim <= 4 (ACLV Corollary 2.6(b)) and ada => (pd M <= 2 or id M <= 1) for every
    indecomposable M (Corollary 2.6(a)): asserted over the in-scope ada zoo members.

QQ-scope (M1). Values recomputed live before pinning (A5: HH dims [1,0,0,0])."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.injective import injective_resolution
from quiverlab.modules.left_right import left_right_parts
from quiverlab.modules.recognizers_ladder import recognizer_ladder

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _radsq_nakayama_a5():
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


def _comm_square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                  ).algebra(relations=["a*b - c*d"], field=QQ)


def _ada_members():
    return [_radsq_nakayama_a5(), _comm_square(),
            linear_path_algebra(3, field=QQ), linear_path_algebra(4, field=QQ)]


@lit
def test_strict_weakly_shod_hh_geq2_vanishes():
    # Survey S5.1: a strict weakly shod algebra (weakly shod but NOT quasi-tilted) has
    # H^i(A) = 0 for i >= 2. The rad^2=0 A5 fixture qualifies (weakly_shod True, quasi_tilted
    # False). Recomputed live: dim HH = [1, 0, 0, 0].
    A = _radsq_nakayama_a5()
    L = recognizer_ladder(A)
    assert L.verdict("weakly_shod") is True and L.verdict("quasi_tilted") is False  # strict
    hh = A.hochschild_cohomology(3, verbose=False)
    assert hh[2] == 0 and hh[3] == 0


@selfcert
def test_ada_implies_gldim_le_4():
    # ACLV Corollary 2.6(b): ada => gl.dim <= 4 (sharp; A5 attains 4).
    for A in _ada_members():
        L = recognizer_ladder(A)
        if L.verdict("ada"):
            assert L.gldim is not None and L.gldim <= 4


@selfcert
def test_ada_implies_pd_le_2_or_id_le_1():
    # ACLV Corollary 2.6(a): ada => every indecomposable has pd <= 2 or id <= 1.
    for A in _ada_members():
        L = recognizer_ladder(A)
        if not L.verdict("ada"):
            continue
        atlas = left_right_parts(A)
        for M in atlas._modules:
            pd_le_2 = M.projective_resolution(3).betti(3) == 0
            id_le_1 = injective_resolution(M, 2).betti(2) == 0
            assert pd_le_2 or id_le_1, M.dimension_vector()
