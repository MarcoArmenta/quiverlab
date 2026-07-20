"""sweep(): rebuild an algebra over several fields and tabulate invariants (spec 3.9)."""
from quiverlab import CC, GF, sweep, truncated_polynomial, linear_path_algebra


def test_sweep_default_invariants_over_fields():
    tab = sweep(truncated_polynomial, 3, fields=[CC, GF(2), GF(3), GF(5)])
    # dimension is field-independent
    for f in (CC, GF(2), GF(3), GF(5)):
        assert tab.cell("dimension", f) == 3
        assert tab.cell("loewy_length", f) == 3


def test_sweep_custom_invariants_and_char_dependence():
    tab = sweep(truncated_polynomial, 2,
                fields=[GF(2), GF(3), GF(5)],
                invariants={"is_frobenius": lambda A: A.is_frobenius(),
                            "hh1": lambda A: A.hochschild_cohomology(1).dims[1]})
    assert tab.cell("is_frobenius", GF(2)) is True
    assert all(isinstance(tab.cell("hh1", f), int) for f in (GF(2), GF(3), GF(5)))


def test_sweep_records_field_errors_without_crashing():
    # an engine-backed invariant over CC must be recorded, not raised, in a sweep cell
    tab = sweep(truncated_polynomial, 2, fields=[CC, GF(5)],
                invariants={"is_frobenius": lambda A: A.is_frobenius()})
    assert tab.cell("is_frobenius", GF(5)) is True
    assert "n/a" in str(tab.cell("is_frobenius", CC)).lower()


def test_sweep_repr_and_latex_are_strings():
    tab = sweep(linear_path_algebra, 2, fields=[CC, GF(2)])
    assert isinstance(repr(tab), str) and "dimension" in repr(tab)
    assert isinstance(tab.latex(), str)
