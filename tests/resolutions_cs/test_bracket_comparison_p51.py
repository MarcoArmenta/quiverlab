"""Plan 51 Task C: Comparison.bracket_of_cs_classes engine=native/transport/auto.

Deep bucket.  Mirrors test_native_cup's Task-4/5 selector + bridge pattern: a TINY
comparison window makes "past-window" land at LOW absolute degree, so only the small
Delta build costs.  engine='native' computes at any degree; 'transport' keeps the
window refusal; 'auto' routes native past the window (byte-identical to the old
default in-window); 'bogus' is a ValueError naming the three options.
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.resolutions_cs.comparison import Comparison
from quiverlab.resolutions_cs.bracket import native_bracket

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_selfcert]


def _kx2_gf5():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(5))


def test_bracket_engine_selector_kx2():
    """engine= forces a route: 'native' computes past the (tiny) window; 'transport'
    keeps the window refusal; an invalid engine is a ValueError naming the options."""
    comp = Comparison(_kx2_gf5(), max_cells=8)                # window 0
    assert comp.window == 0
    u1 = comp.hh_class_cs(1, 0)
    forced_native = comp.bracket_of_cs_classes(u1, u1, engine="native")   # deg 1 > 0
    assert isinstance(forced_native, list)
    with pytest.raises(NotImplementedError):
        comp.bracket_of_cs_classes(u1, u1, engine="transport")
    with pytest.raises(ValueError):
        comp.bracket_of_cs_classes(u1, u1, engine="bogus")


@pytest.mark.oracle_crossengine
def test_bracket_past_window_bridge_kx2():
    """Bridge oracle: a bracket past a TINY window computed NATIVELY equals -- mod
    coboundary -- the SAME bracket by TRANSPORT on a wider-window instance (Plan-17
    canonical CS bases make the coordinate vectors directly comparable).  k[x]/x^2
    over GF(5): [alpha, beta] at (1,2) -> HH^2 is NONZERO, so the match is not a
    vacuous coboundary coincidence."""
    small = Comparison(_kx2_gf5(), max_cells=8)              # window 0 (native past it)
    big = Comparison(_kx2_gf5())                             # wide window (transport)
    assert small.window == 0 and big.window >= 2, (small.window, big.window)

    u1s, u2s = small.hh_class_cs(1, 0), small.hh_class_cs(2, 0)
    u1b, u2b = big.hh_class_cs(1, 0), big.hh_class_cs(2, 0)

    native = small.bracket_of_cs_classes(u1s, u2s)          # max(1,2)=2 > 0 -> native
    transported = big.bracket_of_cs_classes(u1b, u2b)       # <= window -> transport

    sb = [(ch.word, j) for ch, j in small._res._basis(2, "coh")]
    bb = [(ch.word, j) for ch, j in big._res._basis(2, "coh")]
    assert sb == bb, "Plan-17 canonicalization must give element-wise identical bases"
    assert big.same_cohomology_class(native, transported, degree=2), \
        "native (tiny window) != longer transport (wide window) mod coboundary"
    zero2 = [0] * len(big._res._basis(2, "coh"))
    assert not big.same_cohomology_class(transported, zero2, degree=2), \
        "the bridged bracket class must be nonzero -- not a vacuous match"


@pytest.mark.oracle_crossengine
def test_bracket_auto_equals_native_and_transport_across_window_kx2():
    """auto == transport in-window and auto == native past-window (byte comparison
    of the routed cochains): the engine= kwarg (Task C) routes exactly as the old
    default did in-window, and native past it (mirrors cup_of_cs_classes)."""
    comp = Comparison(_kx2_gf5())                            # wide window
    comp._ensure(4)
    u1, u2 = comp.hh_class_cs(1, 0), comp.hh_class_cs(2, 0)
    # (1,1) is in-window -> auto == transport (byte-identical route).
    a11 = comp.bracket_of_cs_classes(u1, u1, engine="auto")
    t11 = comp.bracket_of_cs_classes(u1, u1, engine="transport")
    assert a11 == t11
    # force native past-window on a tiny instance; auto there == native.
    small = Comparison(_kx2_gf5(), max_cells=8)              # window 0
    a = small.bracket_of_cs_classes(u1, u2, engine="auto")   # max(1,2)>0 -> native
    n = small.bracket_of_cs_classes(u1, u2, engine="native")
    assert a == n
