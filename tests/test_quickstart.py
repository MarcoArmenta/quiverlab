"""The README example, verbatim: this is the Plan-01 acceptance test.
A non-coder builds a custom monomial algebra and reads off exact Hochschild
dimensions with a characteristic sweep."""
from quiverlab import CC, GF, Quiver


def test_readme_example():
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    A = Q.algebra(relations=["a*b"], field=CC)

    assert A.dim == 6            # e1, e2, e3, a, b, c
    assert A.basis_labels == ["e_1", "e_2", "e_3", "a", "b", "c"]

    hh = A.hochschild_cohomology(3)
    hh_gf2 = Q.algebra(relations=["a*b"], field=GF(2)).hochschild_cohomology(3)

    # both are exact, certified answers; printing them shows readable tables
    assert hh.dims[0] >= 1 and len(hh.dims) == 4
    assert len(hh_gf2.dims) == 4
    assert "HH^0" in repr(hh)


def test_characteristic_sweep_dual_numbers():
    from quiverlab import truncated_polynomial

    table = {}
    for name, field in [("CC", CC), ("GF(2)", GF(2)), ("GF(3)", GF(3)), ("GF(4)", GF(4))]:
        table[name] = truncated_polynomial(2, field=field).hochschild_cohomology(4).dims
    assert table["CC"] == [2, 1, 1, 1, 1]
    assert table["GF(2)"] == [2, 2, 2, 2, 2]      # the char-2 pathology, exactly
    assert table["GF(3)"] == [2, 1, 1, 1, 1]
    assert table["GF(4)"] == [2, 2, 2, 2, 2]
