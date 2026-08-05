"""QPA as the oracle for the derived AR translate (Plan 43). qpa-marked: skips
locally without the [qpa] extra, mandatory under QUIVERLAB_REQUIRE_QPA=1.

Probes ``TauOfComplex`` / ``HomologyOfComplex`` / ``MappingCone``. The complex-object
scripting route is the documented P39 Ch.10 hazard -- ``TauOfComplex`` on a
``ProjectiveResolution`` raises inside libgap (``no method found for
DirectSumInclusions``), so ``crosscheck_tau_complex`` compares our DERIVED-category
``tau_Db`` (homology in degree 0 of the projective-resolution complex) against QPA's
MODULE AR translate ``DTr(M)`` -- a genuine cross-engine comparison, never a silent
skip (both routes are ``tau M`` for a non-projective interval module over ``kA_n``)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_derived_surface_probe():
    lg = session.libgap_handle()
    present = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
               for name in ("TauOfComplex", "HomologyOfComplex", "MappingCone")}
    # record what QPA exposes; the battery below uses the DTr module-level fallback
    # (TauOfComplex on a ProjectiveResolution does not script through libgap -- the
    # documented P39 Ch.10 hazard).
    assert present["HomologyOfComplex"]              # Ch.10 is present in QPA 1.37


@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ),
                               linear_path_algebra(4, field=QQ)])
def test_tau_complex_homology_vs_qpa(A):
    for v in A.quiver.vertices:
        M = A.simple(v)
        if M.tau().dim == 0:                         # projective end: skip
            continue
        A.crosscheck("tau_complex", M).assert_agree()
