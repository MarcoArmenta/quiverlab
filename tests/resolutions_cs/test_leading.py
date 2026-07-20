import pytest
from quiverlab import Quiver, CC
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
pytest.importorskip("quiverlab.groebner")


def _res(field=CC, xi="2"):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", f"y*x - {xi}*x*y"]
    A = Q.algebra(relations=rels, field=field)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, field), max_degree=6)


def _norm(res, terms):
    return {(res.to_int(c), a, t, cc) for (c, a, t, cc) in terms}


def test_d2_quantum_ci_yx_matches_paper_no_swap():
    """CS §6: d_1^{CS}(1⊗yx⊗1) = y⊗x⊗1 + 1⊗y⊗x − ξ x⊗y⊗1 − ξ 1⊗x⊗y  (ξ=2). a=y is
    LEFT of the first term, c=x is RIGHT of the second: the swap-catcher."""
    res = _res()
    yx = next(c for c in res.ss.S(2) if c.word == ("y", "x"))
    assert _norm(res, res.d_terms(2, yx)) == {
        ( 1, ("y",), ("x",), ()),          #  y ⊗ x ⊗ 1
        ( 1, (),     ("y",), ("x",)),       #  1 ⊗ y ⊗ x
        (-2, ("x",), ("y",), ()),           # −ξ x ⊗ y ⊗ 1
        (-2, (),     ("x",), ("y",)),       # −ξ 1 ⊗ x ⊗ y
    }


def test_delta_leading_x3_not_mirror_reversed():
    """CS §6: δ for x³ = x⊗x²⊗1 − 1⊗x²⊗x. a=x LEFT of the +term, c=x RIGHT of the −term
    (the exact positions my earlier draft had reversed)."""
    res = _res()
    x3 = next(c for c in res.ss.S(3) if c.word == ("x", "x", "x"))
    d = _norm(res, res.delta_terms(3, x3))
    assert d == {(1, ("x",), ("x", "x"), ()), (-1, (), ("x", "x"), ("x",))}


def test_square_d2_fox_derivative():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    res = ChouhySolotarResolution(A, build_reduction_system(Q, ["a*b - c*d"], CC), max_degree=3)
    cd = res.ss.S(2)[0]
    assert _norm(res, res.d_terms(2, cd)) == {
        ( 1, (),     ("c",), ("d",)),       #  1 ⊗ c ⊗ d
        ( 1, ("c",), ("d",), ()),           #  c ⊗ d ⊗ 1
        (-1, (),     ("a",), ("b",)),       # −1 ⊗ a ⊗ b
        (-1, ("a",), ("b",), ()),           # −a ⊗ b ⊗ 1
    }
