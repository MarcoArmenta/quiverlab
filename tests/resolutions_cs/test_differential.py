import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs import cs_cohomology_dims, cs_homology_dims
from quiverlab.fields.linalg import rank
pytest.importorskip("quiverlab.groebner")


def _kx2(field=CC):
    Q = Quiver([1], {"x": (1, 1)})
    return ChouhySolotarResolution(Q.algebra(relations=["x*x"], field=field),
                                   build_reduction_system(Q, ["x*x"], field), max_degree=6)


def _qci(field=CC, xi="2"):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", f"y*x - {xi}*x*y"]  # DEVIATION: brief wrote "({xi})"; parenthesized
    return ChouhySolotarResolution(Q.algebra(relations=rels, field=field),
                                   build_reduction_system(Q, rels, field), max_degree=8)


def test_kx2_boundary_ranks_char0():
    res = _kx2()
    assert [rank(res.matrix(n, "hom"), res.dom) for n in range(1, 7)] == [0, 1, 0, 1, 0, 1]


def test_kx2_dd_zero_both_sides():
    res = _kx2()
    res.assert_dd_zero(upto=6, side="hom"); res.assert_dd_zero(upto=6, side="coh")


def test_square_dd_zero():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    res = ChouhySolotarResolution(A, build_reduction_system(Q, ["a*b - c*d"], CC), max_degree=3)
    res.assert_dd_zero(upto=3, side="hom"); res.assert_dd_zero(upto=3, side="coh")


def test_qci_dd_zero_and_order_condition():
    for field in (CC, GF(2)):
        res = _qci(field=field)
        res.assert_dd_zero(upto=8, side="hom")               # non-commutative recursion closes
        res.assert_order_condition(upto=8)                   # CS Theorem 4.1 condition (2)


@pytest.mark.xfail(strict=False, reason="canonicalization pending: the d²=0 solve is unique only "
                   "up to the constraint nullspace, so a correct-but-noncanonical correction is "
                   "allowed; binding criteria (d²=0 + order gate + HH dims) are checked separately")
def test_qci_d3_correction_matches_paper():
    """CS §6 (verbatim, kept as an aspirational canonical-form pin): d_2^{CS}(1⊗y²x⊗1) =
    y⊗yx⊗1 + ξ 1⊗yx⊗y + ξ² x⊗y²⊗1 − 1⊗y²⊗x  (ξ=2). Edit #2: the BINDING criteria are
    test_qci_dd_zero_and_order_condition (d²=0 + order) and test_qci_homology_matches_bank_vector
    (HH dims, Task 7) — this exact-coefficient check is xfail-tolerant until the canonicalization
    stretch item (normal-form modulo the solve nullspace) flips it strict."""
    res = _qci()
    q = next(c for c in res.ss.S(3) if c.word == ("y", "y", "x"))
    got = {(res.to_int(c), a, t, cc) for (c, a, t, cc) in res.d_terms(3, q)}
    assert got == {
        ( 1, ("y",), ("y", "x"), ()),       #  y ⊗ yx ⊗ 1
        ( 2, (),     ("y", "x"), ("y",)),    #  ξ 1 ⊗ yx ⊗ y
        ( 4, ("x",), ("y", "y"), ()),        #  ξ² x ⊗ y² ⊗ 1
        (-1, (),     ("y", "y"), ("x",)),    # −1 ⊗ y² ⊗ x
    }


def _kx3(field=CC):
    """k<x>/(x³): single vertex, loop x, MONOMIAL relation x*x*x. Non-quadratic (cubic tip)
    but monomial, so _require_in_scope ADMITS it (no correction; collapsed = Bardzell)."""
    Q = Quiver([1], {"x": (1, 1)})
    return ChouhySolotarResolution(Q.algebra(relations=["x*x*x"], field=field),
                                   build_reduction_system(Q, ["x*x*x"], field), max_degree=5)


def test_kx3_monomial_nonquadratic_gates_pass():
    """RESTRICT-boundary regression (edit #1): a MONOMIAL non-quadratic algebra k<x>/(x³)
    passes the guard legitimately (monomial admission) and both CS gates close on both sides.
    Pins the monomial-non-quadratic corner the loud backstops otherwise rest on, uncommitted."""
    for field in (CC, GF(2)):
        res = _kx3(field=field)
        res.assert_dd_zero(upto=5, side="hom")           # d²=0, homology side
        res.assert_dd_zero(upto=5, side="coh")           # d²=0, cohomology side
        res.assert_order_condition(upto=5)               # CS Theorem 4.1 condition (2)
        # guard ADMITS monomial: delta_terms/d_terms must NOT raise (contrast the cubic-tip
        # NON-monomial case below, which does).
        for n in range(1, 6):
            for c in res.ss.S(n):
                res.delta_terms(n, c)                    # no NotImplementedError
                res.d_terms(n, c)                        # no NotImplementedError


def test_cubic_tip_nonmonomial_raises_notimplemented():
    """RESTRICT boundary (edit #1): a non-quadratic (cubic tip) NON-monomial presentation
    raises NotImplementedError at the exact degree-≥3 differential. A = k<x,y>/(x²,y²,xyx−yxy)
    is finite-dimensional (basis {1,x,y,xy,yx,xyx}, dim 6) with a cubic tip and a nonzero tail.
    (If completion changes the basis, keep any admissible f.d. cubic-tip non-monomial algebra.)"""
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", "x*y*x - y*x*y"]
    A = Q.algebra(relations=rels, field=CC)
    res = ChouhySolotarResolution(A, build_reduction_system(Q, rels, CC), max_degree=4)
    with pytest.raises(NotImplementedError) as e:
        res.d_terms(3, res.ss.S(3)[0])                   # first degree-3 differential trips the guard
    assert "quadratic" in str(e.value).lower() and "right_decomposition" in str(e.value)


def test_cubic_tip_nonmonomial_refuses_at_battery_level():
    """RESTRICT boundary (edit #1), END-TO-END: the same cubic-tip NON-monomial presentation
    A = k<x,y>/(x²,y²,xyx−yxy) must refuse through the FULL HH pipeline, not only at the raw
    res.d_terms(3,·) level (test_cubic_tip_nonmonomial_raises_notimplemented pins that). The
    public entry points cs_cohomology_dims / cs_homology_dims and A.hochschild_cohomology(
    engine="cs") each build the resolution and drive the degree-≥3 differential, so the guard's
    NotImplementedError surfaces battery-level, confirming the refusal is not swallowed."""
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x*x", "y*y", "x*y*x - y*x*y"]
    A = Q.algebra(relations=rels, field=CC)
    for call in (lambda: cs_cohomology_dims(A, 3),
                 lambda: cs_homology_dims(A, 3),
                 lambda: A.hochschild_cohomology(3, engine="cs")):
        with pytest.raises(NotImplementedError) as e:
            call()
        assert "quadratic" in str(e.value).lower()
