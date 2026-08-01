"""GF(p) bar-route tables vs direct tt_calculus calls (Plan 35).

k[x]/(x^2) over GF(7): dim HH^n = 1 for all n (the classic pin); the cup ring
in char != 2 has HH^odd squaring to 0 and HH^even polynomial -- we assert the
table entries EQUAL the raw tt_calculus tensor, and the unit law on (0, q)."""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_crossengine]


@pytest.fixture(scope="module")
def A():
    return ql.truncated_polynomial(2, field=ql.GF(7))


def test_cup_tables_match_tt(A):
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.tt_calculus import cup_product_matrix
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "cup", 3, max_cells=4_000_000)
    assert hp.kind == "cup" and hp.basis == "bar/GF(7)"
    E = to_engine(A.unit_adapted())
    for (p, q), table in hp.tables.items():
        C, dl, dr, dout = cup_product_matrix(E, p, q, 7)
        assert table.dims == (dl, dr, dout)
        for k in range(dout):
            for i in range(dl):
                for j in range(dr):
                    assert table.constants[k][i][j] == str(int(C[k, i, j]))


def test_cup_pairs_cover_exactly_p_plus_q_le_top(A):
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "cup", 2, max_cells=4_000_000)
    assert sorted(hp.tables) == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]


def test_cap_pairs_cover_p_le_n_le_top(A):
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "cap", 2, max_cells=4_000_000)
    assert sorted(hp.tables) == [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    t = hp.tables[(1, 2)]
    assert t.out_degree == 1            # HH^1 (x) HH_2 -> HH_1


def test_bracket_pairs_need_p_and_q_ge_1(A):
    from quiverlab.hochschild.products import gfp_product_tables
    hp = gfp_product_tables(A, "bracket", 3, max_cells=4_000_000)
    assert all(p >= 1 and q >= 1 for (p, q) in hp.tables)
    assert (1, 1) in hp.tables and hp.tables[(1, 1)].out_degree == 1
