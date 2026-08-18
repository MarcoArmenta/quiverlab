"""Plan 54 Task C -- the public Algebra.bv_operator method + engine routing + loud
typed refusals."""
import pytest

import quiverlab as ql
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.hochschild.products import BVOperator

pytestmark = pytest.mark.oracle_selfcert


def test_bv_operator_exists_and_returns_bvoperator():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    bv = A.bv_operator(3)
    assert isinstance(bv, BVOperator)
    b = bv.blocks()
    assert b["kind"] == "bv_operator"
    assert set(b) >= {"top", "hh_dims", "matrices", "ranks", "hypothesis",
                      "nakayama", "basis", "window", "references"}


def test_unknown_engine_raises():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    with pytest.raises(QuiverlabError, match="unknown engine"):
        A.bv_operator(2, engine="nope")


def test_engine_cs_not_until_p51():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    with pytest.raises(QuiverlabError, match="not available until P51"):
        A.bv_operator(2, engine="cs")


def test_engine_bar_off_gfp_raises():
    A = ql.truncated_polynomial(2, field=QQ)
    with pytest.raises(QuiverlabError, match="GF\\(p\\)"):
        A.bv_operator(2, engine="bar")


def test_non_frobenius_refused():
    A = ql.linear_path_algebra(2, field=ql.GF(5))
    with pytest.raises(QuiverlabError, match="not Frobenius"):
        A.bv_operator(2)


def test_biklz_blocked_refused():
    A = ql.NakayamaAlgebra(n=2, l=4, cyclic=True, field=ql.GF(2))
    with pytest.raises(QuiverlabError, match="is NOT implemented in v1"):
        A.bv_operator(2)
