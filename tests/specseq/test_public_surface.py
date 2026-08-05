"""Plan 42 public surface: the top-level exports and the four presets import, and
``A.hochschild_bB_ss`` round-trips to a pre-certified SpectralSequence."""
import pytest

from quiverlab import (DoubleComplex, FilteredComplex, GF, SpectralSequence,
                       truncated_polynomial)

pytestmark = pytest.mark.oracle_selfcert


def test_top_level_exports():
    assert FilteredComplex is not None
    assert DoubleComplex is not None
    assert SpectralSequence is not None


def test_preset_and_class_exports_from_specseq():
    from quiverlab.specseq import (Page, Subquotient, ConvergenceReport,
                                   cartan_eilenberg_ss, grothendieck_double_complex,
                                   hochschild_bB_ss, radical_filtration_ss)
    assert all(x is not None for x in (
        Page, Subquotient, ConvergenceReport, cartan_eilenberg_ss,
        grothendieck_double_complex, hochschild_bB_ss, radical_filtration_ss))


def test_algebra_wrapper_round_trips():
    A = truncated_polynomial(2, field=GF(5))
    ss = A.hochschild_bB_ss(4)
    assert isinstance(ss, SpectralSequence)
    # the construction is pre-certified (E_inf totals == total homology); the
    # abutment is HC_*(k[x]/(x^2)).
    Einf = ss.page(ss.convergence.e_infinity_page)
    totals = {}
    for (p, q) in Einf.spots:
        totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
    assert [totals.get(n, 0) for n in range(5)] == list(A.cyclic_homology(4).dims)
