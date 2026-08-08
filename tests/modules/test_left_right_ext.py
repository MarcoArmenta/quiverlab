"""Ext-injectives of add L_A via the tau^{-1} criterion (Plan 55 / R15, ACT [13](3.4)).
Self-cert: an injective module in L_A is always Ext-injective (tau^{-1} I = 0); every
Ext-injective lies in L_A. Literature: over kA_n, L_A = ind A, so X is Ext-injective in
add L_A iff tau^{-1}X = 0 iff X is injective (there are exactly n such). QQ-scope (M1)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_ext_injectives_are_in_left():
    A = linear_path_algebra(3, field=QQ)
    atlas = left_right_parts(A)
    left_ix = {r["index"] for r in atlas.left}
    assert all(r["index"] in left_ix for r in atlas.ext_injectives_left)    # subset of L_A


@lit
def test_kA3_ext_injectives_count_is_n():
    # kA_n hereditary => L_A = ind A, so Ext-injective in add L_A iff tau^{-1}X = 0 iff X
    # injective. kA_3 has exactly 3 indecomposable injectives.
    A = linear_path_algebra(3, field=QQ)
    assert len(left_right_parts(A).ext_injectives_left) == 3


@selfcert
def test_dual_ext_projectives_are_in_right():
    A = linear_path_algebra(3, field=QQ)
    atlas = left_right_parts(A)
    right_ix = {r["index"] for r in atlas.right}
    assert all(r["index"] in right_ix for r in atlas.ext_projectives_right)


@lit
def test_kA3_ext_projectives_count_is_n():
    # dually: over kA_n, X is Ext-projective in add R_A iff tau X = 0 iff X projective;
    # kA_3 has exactly 3 indecomposable projectives.
    A = linear_path_algebra(3, field=QQ)
    assert len(left_right_parts(A).ext_projectives_right) == 3
