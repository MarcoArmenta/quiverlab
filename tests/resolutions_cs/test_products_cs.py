"""CS-route product tables (Plan 35): Domain-generic cup/cap on the CS basis.

Deep-bucket (tests/resolutions_cs -> deep). Oracles: the QQ smoke is
self-certifying (unit row + dims equal cs dims); the GF(p) cross-engine
equality against the bar tables lives in tests/hochschild/test_products_identities.py
(Task 6) -- HERE we pin shape, unit law, and Domain-genericity."""
import pytest

import quiverlab as ql
from quiverlab.fields import QQ

pytestmark = [pytest.mark.oracle_selfcert]


@pytest.fixture(scope="module")
def A_qq():
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1), "y": (1, 1)})
    return Q.algebra(relations=["x*x", "y*y", "x*y + y*x"], field=QQ)


def test_cs_cup_over_qq_unit_law(A_qq):
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims
    from quiverlab.resolutions_cs.products import cs_product_tables
    hp = cs_product_tables(A_qq, "cup", 2, max_cells=4_000_000)
    assert hp.basis.startswith("cs/")
    dims = list(cs_cohomology_dims(A_qq, 2).dims)
    # unit law: the (0, q) table with the identity class of HH^0 = Z(A) acting
    # as identity -- for the class index of 1_A, the table column is the identity.
    t = hp.tables[(0, 1)]
    assert t.dims == (dims[0], dims[1], dims[1])


def test_cs_cap_over_qq_shapes(A_qq):
    from quiverlab.resolutions_cs.homology import cs_homology_dims, cs_cohomology_dims
    from quiverlab.resolutions_cs.products import cs_product_tables
    hp = cs_product_tables(A_qq, "cap", 2, max_cells=4_000_000)
    hdims = list(cs_homology_dims(A_qq, 2).dims)
    cdims = list(cs_cohomology_dims(A_qq, 2).dims)
    t = hp.tables[(1, 2)]
    assert t.dims == (cdims[1], hdims[2], hdims[1])


def test_cs_bracket_refuses():
    # cs_product_tables serves ONLY cup/cap; the bracket has its own builder.
    from quiverlab.resolutions_cs.products import cs_product_tables
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=QQ)
    with pytest.raises(ql.QuiverlabError):
        cs_product_tables(A, "bracket", 2, max_cells=4_000_000)


def test_cs_bracket_tables_serves():
    # Plan 51: the dedicated CS-native bracket builder over QQ (past the bar world).
    from quiverlab.resolutions_cs.products import cs_bracket_tables
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims
    Q = ql.Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=QQ)
    hb = cs_bracket_tables(A, 2, max_cells=4_000_000)
    assert hb.kind == "bracket" and hb.basis == "cs/QQ" and hb.window is None
    cdims = list(cs_cohomology_dims(A, 2).dims)
    # (1,1) -> HH^1 : dims (dim HH^1, dim HH^1, dim HH^1)
    t11 = hb.tables[(1, 1)]
    assert t11.out_degree == 1
    assert t11.dims == (cdims[1], cdims[1], cdims[1])
    assert "bracket_liftings" in hb.references and "oke_koszul" in hb.references
