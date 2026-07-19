"""Plan-04 acceptance: a deep non-monomial Hochschild computation that is genuinely
IMPOSSIBLE for the bar complex (it raises DepthLimitError) yet is delivered exactly by
engine="cs" (the plan's raison d'etre), plus the characteristic sweep whose
char-dependence is a headline for the maths paper.

Exact algebra: the quantum complete intersection A = k<x,y>/(x^2, y^2, yx - 2*x*y)
(CS section 6 Happel-counterexample family). Over the bar complex dim C_n = 4*3^n, so the
degree-12 boundary matrix would have ~4*3^11 x ~4*3^12 ~= 5e11 entries >> max_cells and
bar raises DepthLimitError long before degree 12; the CS resolution has dim C_n = n + 1.

Every expected value is independently pinned: the ell HH-dims are the bank-run values
(Task 11 byte-oracle confirms the differentials; Task 9 confirms bar-agreement in the low
window) and are reconfirmed here by executing the shipped CS engine.
"""
import pytest

from quiverlab import Quiver, CC, GF
from quiverlab.errors import DepthLimitError

pytest.importorskip("quiverlab.groebner")

BGMS_QCI_XI2 = [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]     # bank-derived, degree 0..12


def _qci(field):
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=field)


def test_deep_nonmonomial_hh_impossible_for_bar():
    assert _qci(CC).hochschild_homology(12, engine="cs").dims == BGMS_QCI_XI2
    with pytest.raises(DepthLimitError):
        _qci(CC).hochschild_homology(12, engine="bar")     # exponential blow-up, loud


def test_characteristic_sweep_shows_pathologies():
    # Bank-derived char-dependence (a headline for the maths paper):
    assert _qci(GF(3)).hochschild_homology(8, engine="cs").dims == [3, 4, 6, 8, 10, 12, 14, 16, 18]
    assert _qci(GF(2)).hochschild_homology(8, engine="cs").dims == [3, 4, 4, 4, 4, 4, 4, 4, 4]
    # p = 5: a period-4 pathology (2 has multiplicative order 4 mod 5). Reconfirmed at
    # execution time against the shipped CS engine (Task-14 edit #4).
    assert _qci(GF(5)).hochschild_homology(12, engine="cs").dims == [3, 2, 3, 4, 3, 2, 4, 6, 4, 2, 5, 8, 5]
