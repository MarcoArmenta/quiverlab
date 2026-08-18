"""CS HH(A,M) == bar HH(A,M) degreewise over CC and GF(p), every named coefficient.
The resolution of A is coefficient-independent; only the collapse changes."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims

pytestmark = pytest.mark.oracle_crossengine


@pytest.mark.parametrize("field", [ql.CC, ql.GF(7)])
@pytest.mark.parametrize("mk", ["regular", "dual", "twisted_by_nakayama"])
def test_cs_matches_bar_with_coefficients(field, mk):
    A = ql.truncated_polynomial(3, field=field)     # presented, admissible (symmetric)
    M = getattr(Bimodule, mk)(A)
    top = 5
    assert cs_cohomology_dims(A, top, coefficients=M).dims == \
        A.hochschild_cohomology(top, engine="bar", coefficients=M).dims
    assert cs_homology_dims(A, top, coefficients=M).dims == \
        A.hochschild_homology(top, engine="bar", coefficients=M).dims


@pytest.mark.parametrize("mk", ["dual", "twisted_by_nakayama"])
def test_cs_matches_bar_nontrivial_nu(mk):
    # C2: bar==CS on a NON-symmetric self-injective algebra (nu != id, so the
    # twist path is genuinely exercised, not ~ regular). NakayamaAlgebra(3,3) is
    # dim 9; the reduced-bar cochain basis is dim_M*(m-1)^n = 9*8^n, so bar's d^3
    # is 169M cells (over max_cells). The bar==CS oracle therefore runs to the
    # deepest degree bar can CERTIFY here (top=2) -- deeper degrees are the
    # bar==minimal / CS==minimal cross-checks (Task 5, both non-bar engines).
    # (DEVIATION FROM PLAN top=5: bar's exponential blowup on dim 9. top=2 still
    # exercises corner_M + the twist path + the change-of-basis transform.)
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(7))
    assert A.is_symmetric() is False
    M = getattr(Bimodule, mk)(A)
    top = 2
    assert cs_cohomology_dims(A, top, coefficients=M).dims == \
        A.hochschild_cohomology(top, engine="bar", coefficients=M).dims
    assert cs_homology_dims(A, top, coefficients=M).dims == \
        A.hochschild_homology(top, engine="bar", coefficients=M).dims


def test_public_engine_cs_routes_coefficients():
    A = ql.truncated_polynomial(3, field=ql.CC)
    M = Bimodule.dual(A)
    assert A.hochschild_cohomology(4, engine="cs", coefficients=M).dims == \
        A.hochschild_cohomology(4, engine="bar", coefficients=M).dims


def test_cs_needs_presentation():
    # structure-constant-only algebra (no _quiver): CS refuses (as today), bar
    # still works. Real ctor is Algebra(domain, T, unit, ...) -- a rebuild WITHOUT
    # a presentation, so resolutions_cs cannot form its reduction system.
    from quiverlab.errors import QuiverlabError
    A = ql.truncated_polynomial(3, field=ql.CC)
    Asc = ql.Algebra(A.domain, A.T, A.unit)          # presentation-less (_quiver=None)
    assert getattr(Asc, "quiver", None) is None      # confirm no presentation
    Md = Bimodule.dual(Asc)
    with pytest.raises(QuiverlabError):
        Asc.hochschild_cohomology(3, engine="cs", coefficients=Md)
    # bar still serves the same coefficient (presentation-free reference route)
    assert Asc.hochschild_cohomology(3, engine="bar", coefficients=Md).dims == \
        A.hochschild_cohomology(3, engine="bar", coefficients=Bimodule.dual(A)).dims
