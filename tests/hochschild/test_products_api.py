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
