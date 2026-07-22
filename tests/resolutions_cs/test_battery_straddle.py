"""Plan 12 battery: presentations OUTSIDE the old quadratic-or-monomial scope.

Three cases, each certified per instance by the two CS gates (d^2 = 0 and the
Theorem-4.1 order condition) plus degreewise agreement with the normalized bar
oracle computed LIVE (never hardcoded), the same posture as test_battery_bar:

  * straddle-monomial  k<x,y>/(xx, yy, xyx)      — straddling chains, no tails
  * qci32              k<x,y>/(x^3, y^2, yx-2xy) — CS §7.2 family with n=3: cubic
                       tip + a tailed quadratic tip (formerly refused as a set)
  * cubic-tail         k<x,y>/(xx, yy, xyx-yxy)  — cubic TIP with tail (order
                       resolves the tip as yxy -> xyx): the exact presentation the
                       Plan-04 RESTRICT boundary refused, now computing.

The two hardcoded literals at the bottom are session-verified oracle values
(bar AND minimal-A^e engines agreed, 2026-07-22) pinning the straddle repro."""
import pytest

from quiverlab import Quiver, CC, GF
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
from quiverlab.hochschild.bar import (
    hochschild_cohomology_dims as bar_coh,
    hochschild_homology_dims as bar_hom,
)
pytest.importorskip("quiverlab.groebner")


def _mk(rels, field):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=rels, field=field)


CASES = [
    ("straddle-monomial", ["x*x", "y*y", "x*y*x"]),
    ("qci32", ["x*x*x", "y*y", "y*x - 2*x*y"]),
    ("cubic-tail", ["x*x", "y*y", "x*y*x - y*x*y"]),
]
TOP = 3


@pytest.mark.parametrize("name,rels", CASES, ids=[c[0] for c in CASES])
def test_gates_close(name, rels):
    for field in (CC, GF(2), GF(3)):
        Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
        A = Q.algebra(relations=rels, field=field)
        res = ChouhySolotarResolution(A, build_reduction_system(Q, rels, field),
                                      max_degree=TOP + 1)
        res.assert_dd_zero(upto=TOP + 1, side="hom")
        res.assert_dd_zero(upto=TOP + 1, side="coh")
        res.assert_order_condition(upto=TOP + 1)


@pytest.mark.parametrize("name,rels", CASES, ids=[c[0] for c in CASES])
@pytest.mark.parametrize("field", [GF(32003), GF(2), GF(3), GF(5)])
def test_cs_matches_bar_degreewise(name, rels, field):
    A = _mk(rels, field)
    assert cs_homology_dims(A, TOP).dims == bar_hom(A, TOP).dims, name
    assert cs_cohomology_dims(A, TOP).dims == bar_coh(A, TOP).dims, name


def test_cubic_tail_matches_bar_char0():
    A = _mk(["x*x", "y*y", "x*y*x - y*x*y"], CC)       # exact char-0 domain
    assert cs_homology_dims(A, TOP).dims == bar_hom(A, TOP).dims
    assert cs_cohomology_dims(A, TOP).dims == bar_coh(A, TOP).dims


def test_straddle_monomial_literals():
    """Session-verified oracle values for k<x,y>/(xx,yy,xyx) (2026-07-22)."""
    assert list(cs_homology_dims(_mk(["x*x", "y*y", "x*y*x"], GF(32003)), 3).dims) == \
        [4, 4, 6, 9]
    assert list(cs_homology_dims(_mk(["x*x", "y*y", "x*y*x"], GF(2)), 3).dims) == \
        [4, 7, 9, 12]
