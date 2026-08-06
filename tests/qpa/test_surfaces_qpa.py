"""Surface Jacobians are gentle -- crosscheck our is_gentle / is_special_biserial against
QPA's IsGentleAlgebra / IsSpecialBiserialAlgebra on the surface algebras (Plan 48). QPA has
NO triangulation / marked-surface surface, so the comparison is at the resulting-algebra
level, mirroring P46's recognizer parity (tests/qpa/test_recognizers_qpa.py). A standing
guard FAILS if QPA ever ships a surface constructor."""
import pytest

from quiverlab.fields import GF
from quiverlab.invariants import recognizers as rec
from quiverlab.qpa import scripts, session
from quiverlab.surfaces.qp import jacobian_of
from quiverlab.surfaces.triangulation import (annulus_triangulation, fan_triangulation,
                                             hexagon_with_internal_triangle)

pytestmark = [pytest.mark.qpa,
              pytest.mark.skipif(session.should_skip_qpa(),
                                 reason="[qpa] backend not installed")]

_CASES = [
    ("disc_fan_A3", lambda: fan_triangulation(6)),
    ("annulus_C21", lambda: annulus_triangulation(2, 1)),
    ("annulus_C22", lambda: annulus_triangulation(2, 2)),
    ("hexagon_internal", lambda: hexagon_with_internal_triangle()),
]


@pytest.mark.parametrize("name, factory", _CASES, ids=[c[0] for c in _CASES])
def test_surface_jacobian_gentle_matches_qpa(name, factory):
    J = jacobian_of(factory(), field=GF(3))
    decl = scripts.quiver_and_algebra_script(J)          # binds GAP variable `A`
    ours_g = rec.is_gentle(J)
    ours_sb = rec.is_special_biserial(J)
    qpa_g = bool(session.run(decl + "\nIsGentleAlgebra(A);"))
    qpa_sb = bool(session.run(decl + "\nIsSpecialBiserialAlgebra(A);"))
    assert ours_g is True                                # ABCP: surface Jacobians are gentle
    assert (ours_g, ours_sb) == (qpa_g, qpa_sb)


def test_qpa_has_no_surface_surface():
    # standing honest-scope guard: QPA has no triangulation / marked-surface constructor.
    for nm in ("TriangulationOfSurface", "MarkedSurface", "SurfaceAlgebra",
               "GentleAlgebraFromSurface", "AlgebraFromSurfaceTriangulation"):
        assert not bool(session.run(f'IsBoundGlobal("{nm}");')), nm
