"""Public product surface: routing, refusal matrix, provenance (Plan 35 §2)."""
import pytest

import quiverlab as ql
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ

pytestmark = [pytest.mark.oracle_selfcert]


def test_gfp_routes_to_bar():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    hp = A.cup_products(2)
    assert hp.basis == "bar/GF(7)"


def test_quiver_presented_qq_routes_to_cs():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=QQ)
    hp = A.cup_products(2)
    assert hp.basis.startswith("cs/")


def test_structure_constants_off_gfp_refuse():
    # a genuinely presentation-less algebra: k[x]/(x^2) from structure constants
    B = ql.Algebra.from_structure_constants(
        [[[1, 0], [0, 1]], [[0, 1], [0, 0]]], unit=[1, 0], field=QQ)
    with pytest.raises(QuiverlabError):
        B.cup_products(2)


def test_bracket_refuses_off_gfp():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=QQ)
    with pytest.raises(QuiverlabError):
        A.gerstenhaber_brackets(2)


def test_unknown_engine_refuses():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    with pytest.raises(QuiverlabError):
        A.cup_products(2, engine="fast")     # products know auto/bar/cs only


def test_explicit_cs_on_gfp_serves_cs_basis():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(7))
    hp = A.cup_products(2, engine="cs")
    assert hp.basis == "cs/GF(7)"


def test_connes_serves_both_fields():
    for F in (ql.GF(7), QQ):
        A = ql.truncated_polynomial(2, field=F)
        cb = A.connes_differentials(2)
        assert set(cb.matrices) == {0, 1}


def test_auto_gfp_falls_back_to_cs_past_the_bar_precheck():
    """engine='auto' over GF(p): when the bar cochain-basis pairing exceeds
    max_cells (the gfp route's DepthLimitError pre-check), a quiver-presented
    algebra silently falls back to the CS native route -- the Marco-2026-07-26
    dispatch amendment, extended to the product surface. The smallest cup table on
    the x^3 loop pairs the two degree-0 bases (3*3 = 9 cells), so max_cells=4 trips
    it; the result must still arrive, served by CS."""
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(7))
    hp = A.cup_products(2, max_cells=4)          # engine='auto' by default
    assert hp.basis.startswith("cs/")            # the fallback fired


def test_explicit_bar_does_not_fall_back_and_raises_at_the_precheck():
    """No silent fallback for an EXPLICIT engine='bar': it honours the wall and
    raises rather than quietly switching engines (only 'auto' falls back)."""
    from quiverlab.errors import DepthLimitError
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=ql.GF(7))
    with pytest.raises(DepthLimitError):
        A.cup_products(2, engine="bar", max_cells=4)
