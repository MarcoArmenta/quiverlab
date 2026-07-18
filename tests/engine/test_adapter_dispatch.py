import pytest
from quiverlab import CC, GF, Quiver, truncated_polynomial
from quiverlab.errors import FieldError, QuiverlabError


def _zoo(field):
    Q1 = Quiver(vertices=[1], arrows={"x": (1, 1)})
    Q2 = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    Q3 = Quiver(vertices=[1, 2], arrows={})
    return [
        Q1.algebra(relations=["x^2"], field=field),
        Q1.algebra(relations=["x^4"], field=field),
        Q2.algebra(relations=["a*b"], field=field),
        Q2.algebra(relations=[], field=field),
        Q3.algebra(relations=[], field=field),
    ]


@pytest.mark.parametrize("p", [2, 3, 5])
def test_engine_matches_bar_cohomology_and_homology(p):
    # The central correctness gate: engine (fast, numpy mod-p rank) vs the pure
    # bar complex (exact Gaussian elimination over the Domain) -- genuinely
    # independent code paths. The bar oracle is exponential, so top=4 exceeds its
    # max_cells guard (4M) for the dim-5/6 members; cap the reference depth per
    # algebra to what the bar can feasibly compute. Equality is checked over every
    # HH^/HH_ degree the independent bar path can reach.
    for A in _zoo(GF(p)):
        top = 4 if A.dim <= 4 else 3
        slow_c = A.hochschild_cohomology(top, engine="bar").dims
        fast_c = A.hochschild_cohomology(top, engine="fast").dims
        assert fast_c == slow_c, f"HH^ mismatch over GF({p}) on {A!r}"
        slow_h = A.hochschild_homology(top, engine="bar").dims
        fast_h = A.hochschild_homology(top, engine="fast").dims
        assert fast_h == slow_h, f"HH_ mismatch over GF({p}) on {A!r}"


def test_auto_dispatch_picks_fast_for_gfp_and_bar_for_cc():
    A = truncated_polynomial(2, field=GF(2))
    B = truncated_polynomial(2, field=CC)
    assert "fast" in A.hochschild_cohomology(2).engine
    assert "bar" in B.hochschild_cohomology(2).engine
    assert A.hochschild_cohomology(2).dims == [2, 2, 2]
    assert B.hochschild_cohomology(2).dims == [2, 1, 1]


def test_fast_engine_refuses_non_prime_fields_loudly():
    B = truncated_polynomial(2, field=CC)
    C = truncated_polynomial(2, field=GF(4))
    for A in (B, C):
        with pytest.raises(FieldError):
            A.hochschild_cohomology(2, engine="fast")


def test_gf4_still_works_via_bar_auto():
    A = truncated_polynomial(2, field=GF(4))
    assert A.hochschild_cohomology(2).dims == [2, 2, 2]
    assert "bar" in A.hochschild_cohomology(2).engine


def test_unknown_engine_string_raises():
    # anything outside {'auto', 'bar', 'fast'} is refused loudly on both wrappers
    A = truncated_polynomial(2, field=GF(2))
    with pytest.raises(QuiverlabError):
        A.hochschild_cohomology(2, engine="qpa")
    with pytest.raises(QuiverlabError):
        A.hochschild_homology(2, engine="qpa")


def test_fast_engine_field_error_precedes_depth_guard():
    # Ordering lock: engine='fast' on a non-prime field must raise FieldError,
    # NOT DepthLimitError, even when top is large enough to trip the bar-size
    # guard. Uses dim 3 (m-1=2) so the guard genuinely fires without the fix.
    A = truncated_polynomial(3, field=CC)
    with pytest.raises(FieldError):
        A.hochschild_cohomology(40, engine="fast")
    with pytest.raises(FieldError):
        A.hochschild_homology(40, engine="fast")
