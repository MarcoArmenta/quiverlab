import pytest
from quiverlab import Quiver, CC
pytest.importorskip("quiverlab.groebner")


def _square():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=CC)


def test_engine_cs_selectable():
    t = _square().hochschild_cohomology(4, engine="cs")
    assert t.dims == [1, 0, 0, 0, 0] and "Chouhy-Solotar" in (t.engine or "")


def test_default_auto_is_unchanged():
    """Compatibility: default engine='auto' gives the SAME result & label as before Plan 04.

    FLAG (transcription fix, see p04-task-13-report.md): the brief's literal assertion
    `_square().hochschild_cohomology(4) == [1, 0, 0, 0, 0]` is impossible under the FROZEN
    Plan-02 dispatch. Over CC, auto routes to the bar oracle, which walls out at d^3 for
    this 9-dim algebra (36864 x 4608 > max_cells) and raises DepthLimitError — it does NOT,
    and per Pillar-4 must NOT, silently adopt CS's [1, 0, 0, 0, 0]. Changing auto to return
    that would BE the silent redefinition the plan forbids. So this regression proves
    "auto unchanged" the only faithful way: (1) where bar can compute, auto returns bar's
    answer and label; (2) where only CS could finish (top=4), default auto still hits the
    same bar wall it did before Plan 04 rather than quietly becoming CS. Contrast
    test_opt_in_auto_cs_routes_to_cs: the SAME algebra/top reaches [1, 0, 0, 0, 0] ONLY
    when the caller explicitly opts in with auto_cs=True."""
    from quiverlab import DepthLimitError

    # (1) Same result & label as before Plan 04, where the bar oracle can reach.
    t = _square().hochschild_cohomology(2)               # engine="auto", auto_cs default False
    assert "Chouhy-Solotar" not in (t.engine or "") and t.dims == [1, 0, 0]

    # (2) Pillar-4: default auto did NOT silently route to CS — it still hits the bar wall.
    with pytest.raises(DepthLimitError):
        _square().hochschild_cohomology(4)


def test_opt_in_auto_cs_routes_to_cs():
    t = _square().hochschild_cohomology(4, auto_cs=True)
    assert "Chouhy-Solotar" in (t.engine or "")


def test_trace_claims_equal_values():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=CC)
    A.hochschild_cohomology(4, engine="cs", trace=(tr := []))
    dim_events = [e for e in tr if type(e).__name__ == "ResolutionTerm"]
    assert [e.collapsed_dim for e in dim_events][:3] == [2, 2, 2]
