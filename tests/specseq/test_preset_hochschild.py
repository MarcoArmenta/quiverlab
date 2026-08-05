"""Hochschild (b,B) spectral sequence. Cross-engine: E_inf total == HC dims
(the SBI/free-oracle pins), incl. the ground-field HC = [1,0,1,0,...] closed
form. Self-cert convergence rides along (run at construction)."""
import pytest

from quiverlab import GF, Quiver, truncated_polynomial
from quiverlab.specseq.presets import hochschild_bB_ss

pytestmark = pytest.mark.oracle_crossengine


def _ground(p):
    # k as a quiver algebra: one vertex, no arrows => A = k.
    return Quiver([1], {}).algebra(relations=[], field=GF(p))


def _einf_totals(ss, top):
    Einf = ss.page(ss.convergence.e_infinity_page)
    totals = {}
    for (p, q) in Einf.spots:
        totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
    return [totals.get(n, 0) for n in range(top + 1)]


def test_ground_field_hc_closed_form():
    A = _ground(7)
    top = 6
    ss = hochschild_bB_ss(A, top)
    hc = A.cyclic_homology(top)
    assert _einf_totals(ss, top) == list(hc.dims)
    assert list(hc.dims) == [1, 0, 1, 0, 1, 0, 1]        # HC_*(k) closed form


@pytest.mark.oracle_literature
def test_dual_numbers_hc():
    # k[x]/(x^2): the (b,B) SS abuts to HC_*(k[x]/(x^2)); pin E_inf total == HC.
    A = truncated_polynomial(2, field=GF(5))
    top = 5
    ss = hochschild_bB_ss(A, top)
    hc = A.cyclic_homology(top)
    assert _einf_totals(ss, top) == list(hc.dims)


def test_max_cells_guard_is_loud():
    from quiverlab.errors import DepthLimitError
    # multi-vertex, bar basis blows up. The plan snippet's 3-cycle with relations=[]
    # is INFINITE-dimensional (refused at construction); the truncated 3-cycle is the
    # finite multi-vertex algebra whose bar basis genuinely blows up (dim 6, m*(m-1)^n).
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}).algebra(
        relations=["a*b", "b*c", "c*a"], field=GF(5))
    with pytest.raises(DepthLimitError):
        hochschild_bB_ss(A, 6, max_cells=10_000)
