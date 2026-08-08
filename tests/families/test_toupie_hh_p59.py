"""Toupie Hochschild cohomology (Plan 59 / R35). Literature (SAFE pin): the a-Kronecker
HH^* = [1, a^2-1, 0, ...] (ALS 2020; the SAME VALUE the m-Kronecker engine test pins,
recomputed here independently over QQ through Algebra.hochschild_cohomology at top=4).
Cross-engine: bar == CS degreewise where bar survives (non-hereditary toupies blow the
bar window past low degree). sl_a (char 0 only): a = # DIRECT source->sink arrows;
dim HH^1 >= a^2-1, equality on the a-Kronecker. ALL VALUES BELOW VERIFIED LIVE on dev.

Deviation from the plan sketch (verified 2026-08-07): engine="bar" at top=4 on the
3-Kronecker (dim 5) blows past max_cells at d^4 (the bar oracle is exponential
regardless of the vanishing homology of a hereditary algebra), so the bar == CS
degreewise comparison is pinned WHERE BAR SURVIVES -- top=4 on the 2-Kronecker, top=3 on
the 3-Kronecker -- and the full CS profile is pinned to top=4. This is exactly the
"bar == CS where bar survives" honest-scope the plan mandates."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.families.toupie import (ToupieAlgebra, toupie_direct_arrow_count,
                                       toupie_sl_a_lower_bound)

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("a", [2, 3, 4])
def test_a_kronecker_hh_is_closed_form(a):
    A = ToupieAlgebra([1] * a, field=QQ)                   # a-Kronecker, hereditary
    hh = A.hochschild_cohomology(4, verbose=False).dims    # VERIFIED: a=2->[1,3,..], a=3->[1,8,..]
    assert hh == [1, a * a - 1, 0, 0, 0]                   # ALS 2020: [1, a^2-1, 0, ..]


@lit
@pytest.mark.parametrize("a", [2, 3, 4])
def test_a_kronecker_hh1_is_dim_sl_a(a):
    # char 0: HH^1 = sl_a exactly on the a-Kronecker (all a arrows are direct).
    A = ToupieAlgebra([1] * a, field=QQ)
    hh1 = A.hochschild_cohomology(1, verbose=False).dims[1]
    assert toupie_direct_arrow_count(A) == a               # all branches length 1
    assert hh1 == a * a - 1 == toupie_sl_a_lower_bound(A)  # dim sl_a = a^2 - 1


@xeng
def test_bar_equals_cs_on_2_kronecker_empty_ideal():
    # CS handles the EMPTY ideal (VERIFIED): 2-Kronecker bar survives to top=4.
    A = ToupieAlgebra([1, 1], field=QQ)                    # dim 4, hereditary
    bar = A.hochschild_cohomology(4, engine="bar", verbose=False).dims
    cs = A.hochschild_cohomology(4, engine="cs", auto_cs=True, verbose=False).dims
    cs_plain = A.hochschild_cohomology(4, engine="cs", verbose=False).dims
    assert bar == cs == cs_plain == [1, 3, 0, 0, 0]        # VERIFIED live


@xeng
def test_bar_equals_cs_on_3_kronecker_where_bar_survives():
    # 3-Kronecker (dim 5): bar survives only to top=3 (d^4 blows past max_cells). Compare
    # bar == CS degreewise where bar SURVIVES (top=3), and pin the full CS profile.
    A = ToupieAlgebra([1, 1, 1], field=QQ)                 # dim 5, hereditary
    bar3 = A.hochschild_cohomology(3, engine="bar", verbose=False).dims
    cs3 = A.hochschild_cohomology(3, engine="cs", auto_cs=True, verbose=False).dims
    cs_plain3 = A.hochschild_cohomology(3, engine="cs", verbose=False).dims
    assert bar3 == cs3 == cs_plain3 == [1, 8, 0, 0]        # VERIFIED live (empty ideal, CS ok)
    full = A.hochschild_cohomology(4, engine="cs", auto_cs=True, verbose=False).dims
    assert full == [1, 8, 0, 0, 0]                         # CS past the bar window


@xeng
def test_bar_equals_cs_on_a_NON_hereditary_toupie():
    # A genuine NON-monomial cross-branch toupie. Bar blows up past deg 2, so compare
    # where bar SURVIVES: top=2. VERIFIED: both [1,0,0].
    A = ToupieAlgebra([2, 2], relations=["p0 - p1"], field=QQ)   # commutative square, dim 9
    bar = A.hochschild_cohomology(2, engine="bar", verbose=False).dims
    cs = A.hochschild_cohomology(2, engine="cs", auto_cs=True, verbose=False).dims
    assert bar == cs == [1, 0, 0]                          # VERIFIED live (a_direct=0 -> no sl_a)
    full = A.hochschild_cohomology(4, engine="cs", auto_cs=True, verbose=False).dims
    assert full == [1, 0, 0, 0, 0]                         # CS past the bar window (VERIFIED)


@xeng
def test_sl_a_lower_bound_with_direct_arrows():
    # char 0: dim HH^1 >= dim sl_a = a_direct^2 - 1, a_direct = # direct source->sink arrows.
    # mixed toupie: 2 DIRECT arrows + one length-2 branch (relation-free). VERIFIED [1,6,..].
    A = Quiver([0, "v", 1],
               {"d0": (0, 1), "d1": (0, 1), "b1": (0, "v"), "b2": ("v", 1)}
               ).algebra(relations=[], field=QQ)
    assert toupie_direct_arrow_count(A) == 2
    hh1 = A.hochschild_cohomology(3, engine="cs", auto_cs=True, verbose=False).dims[1]
    assert hh1 == 6 >= toupie_sl_a_lower_bound(A) == 2 * 2 - 1   # 6 >= 3 (strict inclusion)


def test_sl_a_refused_off_char_zero():
    A = ToupieAlgebra([1, 1], field=GF(7))
    with pytest.raises(QuiverlabError):
        toupie_sl_a_lower_bound(A)                         # sl_a inclusion is char-0 only


@selfcert
def test_sl_a_lower_bound_clamped_at_zero_direct_arrows():
    # a toupie with NO direct source->sink arrows (two length-2 branches): a_direct = 0,
    # so dim sl_0 = 0 -- the public function clamps max(a^2-1, 0), never the bogus -1.
    A = ToupieAlgebra([2, 2], field=QQ)
    assert toupie_direct_arrow_count(A) == 0
    assert toupie_sl_a_lower_bound(A) == 0
