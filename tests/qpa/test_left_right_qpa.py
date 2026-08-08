"""QPA cross-check for the left/right parts (Plan 55). qpa-marked: skips locally, mandatory
under QUIVERLAB_REQUIRE_QPA=1. QPA has NO left/right-part or support-algebra verb -- this test
documents that (fail-if-appears) and corroborates the POINTWISE ingredients (the pd<=1 / id<=1
flags that DEFINE L_A / R_A) via QPA's ProjectiveResolution / InjDimensionOfModule on the same
kA3 indecomposables. QQ-scope (M1): kA3 has distinct dim-vectors, so identification never
enters the positive-only branch."""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.qpa

_skip = pytest.mark.skipif(session.should_skip_qpa(),
                           reason="[qpa] backend not installed")


@_skip
def test_qpa_has_no_left_right_surface():
    lg = session.libgap_handle()
    for name in ("LeftPartOfModuleCategory", "RightPartOfModuleCategory",
                 "SupportAlgebra", "LeftSupport", "RightSupport"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), (
            f"QPA now exposes {name} -- wire a real crosscheck (this assert is the trip-wire "
            "that Plan 55's QPA scope note is stale)")


@_skip
def test_pd_id_flags_match_qpa_pointwise():
    # the pd<=1 / id<=1 flags that DEFINE L_A/R_A, corroborated per indecomposable by QPA:
    # our projective resolution terms agree with QPA's ProjectiveResolution through degree 2
    # (hence betti(2)==0 <=> pd<=1 matches), and our injective dimension agrees with QPA's
    # InjDimensionOfModule (hence id<=1 matches).
    from quiverlab import linear_path_algebra
    from quiverlab.fields import QQ
    A = linear_path_algebra(3, field=QQ)
    atlas = A.left_right_parts()
    assert atlas.is_complete
    for i, M in enumerate(atlas._modules):
        A.crosscheck("proj_resolution", M, 2).assert_agree()
        A.crosscheck("inj_dimension", M, 4).assert_agree()
        # the atlas flags are exactly the >=2-term / dimension predicates QPA just confirmed.
        assert atlas.pd_le_1[i] == (M.projective_resolution(2).betti(2) == 0)
        assert atlas.id_le_1[i] == (M.injective_resolution(2).betti(2) == 0)
