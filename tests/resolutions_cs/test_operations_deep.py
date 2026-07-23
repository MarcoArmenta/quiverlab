"""Plan 14: CS operations beyond the quadratic window.

Phase A — the comparison map Phi# is a chain map for EVERY admissible presentation
inside the bar window (homotopy-lifted expansion for non-monomial degree >= 3, which
previously raised "Skoldberg homotopy expansion (a later phase)"), asserted by the
executable gates on the Plan-12 algebras. Phase B — cap transport (unit-cap and
module identities). All gates are matrix identities over GF(p): ground truth, not
approximations."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.resolutions_cs.comparison import Comparison
pytest.importorskip("quiverlab.groebner")


def _mk(rels, p=32003):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=rels, field=GF(p))


CUBIC_TAIL = ["x*x", "y*y", "x*y*x - y*x*y"]
QCI32 = ["x*x*x", "y*y", "y*x - 2*x*y"]
STRADDLE = ["x*x", "y*y", "x*y*x"]


@pytest.mark.parametrize("rels", [CUBIC_TAIL, QCI32, STRADDLE],
                         ids=["cubic-tail", "qci32", "straddle-monomial"])
def test_chain_map_beyond_quadratic_window(rels):
    """Previously: NotImplementedError for any non-monomial presentation at n >= 3.
    Now the lifted Phi# must satisfy delta_cs Phi# = Phi# delta_bar through degree 3
    and be a two-sided iso on HH^0..3."""
    comp = Comparison(_mk(rels))
    comp.assert_chain_map(upto=3)


@pytest.mark.parametrize("rels", [CUBIC_TAIL, QCI32],
                         ids=["cubic-tail", "qci32"])
def test_transport_roundtrip_beyond_quadratic_window(rels):
    comp = Comparison(_mk(rels))
    comp.assert_transport_roundtrip_identity(upto=3)


@pytest.mark.parametrize("rels", [CUBIC_TAIL, QCI32],
                         ids=["cubic-tail", "qci32"])
def test_cup_routes_agree_and_graded_commutative(rels):
    """cup via cochain-level transport == cup via class-level transport, and cup is
    graded-commutative on HH (odd*odd anticommutes) — on a formerly refused algebra."""
    comp = Comparison(_mk(rels))
    reps1 = comp.cs_cohomology_basis(1)
    if not reps1:
        pytest.skip("HH^1 = 0 here; nothing to cup")
    u = comp.hh_class_cs(1, 0)
    lhs = comp.cup_of_cs_classes(u, u)
    rhs = comp.transport_then_bar_cup(u, u)
    assert comp.same_cohomology_class(lhs, rhs, degree=2)
    # graded commutativity: |u| odd => u∪u ~ -(u∪u) => 2(u∪u) ~ 0; over p != 2 that
    # means u∪u is a coboundary-class of order dividing 2 — assert u∪u ~ -(u∪u).
    neg = [(-x) % comp.p for x in lhs]
    assert comp.same_cohomology_class(lhs, neg, degree=2) or comp.p == 2


# ---------------------------------------------------------------------------
# Phase B: cap transport (previously: no CS wrapper at all — internals 09)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("rels", [CUBIC_TAIL, QCI32, STRADDLE],
                         ids=["cubic-tail", "qci32", "straddle-monomial"])
def test_cap_unit_identity(rels):
    """z ∩ 1 = z: capping with the CS class of the bar cup-unit is the identity on
    HH_n (n = 1, 2 where nonzero)."""
    from quiverlab.engine import tt_calculus as TT
    comp = Comparison(_mk(rels))
    unit_cs = comp.transport_class_bar_to_cs(
        [int(x) % comp.p for x in TT.unit_cochain(comp.E)], 0)
    from quiverlab.resolutions_cs.comparison import CSClass
    u0 = CSClass(0, unit_cs)
    for n in (1, 2):
        reps = comp.cs_homology_basis(n)
        if not reps:
            continue
        z = comp.hh_class_cs_hom(n, 0)
        capped = comp.cap_of_cs_classes(u0, z)
        assert comp.same_homology_class(capped, z.vec, degree=n), f"unit cap at HH_{n}"


@pytest.mark.parametrize("rels", [CUBIC_TAIL, QCI32],
                         ids=["cubic-tail", "qci32"])
def test_cap_module_identity(rels):
    """(z ∩ f) ∩ g ~ z ∩ (f ∪ g): HH_* is a module over the cup ring, checked through
    the transported operations on a formerly refused algebra."""
    from quiverlab.resolutions_cs.comparison import CSClass
    comp = Comparison(_mk(rels))
    if not comp.cs_cohomology_basis(1) or not comp.cs_homology_basis(2):
        pytest.skip("needs HH^1 != 0 and HH_2 != 0")
    f = comp.hh_class_cs(1, 0)
    g = comp.hh_class_cs(1, 0)
    z = comp.hh_class_cs_hom(2, 0)
    lhs_inner = comp.cap_of_cs_classes(f, z)                     # degree-1 chain
    lhs = comp.cap_of_cs_classes(g, CSClass(1, lhs_inner))       # degree-0 chain
    fg = comp.cup_of_cs_classes(f, g)                            # degree-2 cochain
    rhs = comp.cap_of_cs_classes(CSClass(2, fg), z)              # degree-0 chain
    assert comp.same_homology_class(lhs, rhs, degree=0)
