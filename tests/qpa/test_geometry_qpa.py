"""QPA probe for a representation-geometry surface (Plan 49 / C8). qpa-marked:
skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1. QPA 1.37 has NO confirmed
orbit / canonical-decomposition / degeneration verb -- this test documents that
honestly and FAILS if one ever appears (forcing a real crosscheck). What QPA CAN
corroborate is dim End(M) = dim Hom_A(M, M), the orbit-dim factor
(dim O_M = dim GL(d) - dim End(M)), via HomOverAlgebra (the P37 hom_glue crosscheck)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.qpa import session
from quiverlab.qpa.crosscheck import crosscheck_hom_glue

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name a QPA geometry surface might plausibly use; none are expected present.
_GEOMETRY_NAMES = (
    "OrbitDimension",
    "CanonicalDecomposition",
    "DegenerationOrder",
    "DegenerateModule",
    "GenericDecomposition",
)


def test_qpa_has_no_degeneration_surface():
    lg = session.libgap_handle()
    # Scan the global name table FIRST (IsBoundGlobal REGISTERS a queried name into
    # NamesGVars, so scanning after would echo the probes back -- the Plan-35 gotcha).
    # NB: bare "orbit" is USELESS here -- GAP core has hundreds of group-theory
    # Orbit*/OrbitStabilizer* names unrelated to module degeneration; we scan only
    # for REPRESENTATION-theory compound terms a real degeneration surface would use.
    gvar_names = [str(n).lower() for n in lg.eval("NamesGVars()")]
    geometry_like = sorted(set(
        n for n in gvar_names
        if any(w in n for w in ("degeneration", "canonicaldecomposition",
                                "modulevariety", "representationvariety",
                                "genericdecomposition"))))
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _GEOMETRY_NAMES}
    present = [name for name, ok in bound.items() if ok]
    if present or geometry_like:
        pytest.fail(
            "QPA now exposes a representation-geometry surface -- wire a real "
            "crosscheck for orbit dim / canonical decomposition / degeneration "
            f"(this assert is the trip-wire that Plan 49's QPA scope note is stale). "
            f"present={present}, geometry_like={geometry_like}")


def test_end_dim_matches_qpa_hom():
    # what QPA CAN corroborate: dim End(M) = dim Hom_A(M, M), the orbit-dim factor.
    # crosscheck dim Hom against QPA's HomOverAlgebra on kA3 indecomposables.
    A = linear_path_algebra(3, field=QQ)
    for v in (1, 2, 3):
        M = A.simple(v)
        crosscheck_hom_glue(A, M, M).assert_agree()   # dim End(M) via HomOverAlgebra
    # and on a projective, whose End is likewise the orbit stabiliser dimension
    crosscheck_hom_glue(A, A.projective(1), A.projective(1)).assert_agree()
