"""Container contract: exact-string constants, degenerate degrees present,
canonical .blocks() shape shared by both runners (Plan 35)."""
import pytest

from quiverlab.hochschild.products import ProductTable, HHProducts, ConnesB


def test_product_table_is_frozen_and_stringly_exact():
    t = ProductTable(kind="cup", degrees=(1, 1), out_degree=2,
                     dims=(2, 2, 1), constants=((("3", "0"), ("0", "1")),))
    assert t.constants[0][0][0] == "3"
    with pytest.raises(Exception):
        t.kind = "cap"


def test_hhproducts_blocks_shape():
    t = ProductTable(kind="cup", degrees=(0, 0), out_degree=0,
                     dims=(1, 1, 1), constants=((("1",),),))
    hp = HHProducts(kind="cup", top=0, tables={(0, 0): t},
                    engine="hanlab engine (F_p fast rank)",
                    basis="bar/GF(7)", window=None, references=["cup"])
    b = hp.blocks()
    assert b["kind"] == "cup" and b["top"] == 0
    assert b["basis"] == "bar/GF(7)"
    assert b["tables"][0]["degrees"] == [0, 0]
    assert b["tables"][0]["constants"] == [[["1"]]]
    assert "window" not in b            # cup carries no window key


def test_bracket_blocks_carry_window():
    hp = HHProducts(kind="bracket", top=3, tables={},
                    engine="hanlab engine (F_p fast rank)",
                    basis="bar/GF(7)", window=2, references=["bracket"])
    assert hp.blocks()["window"] == 2


def test_connesb_blocks_shape():
    cb = ConnesB(top=2, hh_dims=[1, 1, 1],
                 matrices={0: [["1"]], 1: [["0"]]}, ranks={0: 1, 1: 0},
                 engine="engine (b,B) GF(7)", references=["cyclic"])
    b = cb.blocks()
    assert b["kind"] == "connes_b" and b["ranks"] == {"0": 1, "1": 0}
    assert b["matrices"]["1"] == [["0"]]
