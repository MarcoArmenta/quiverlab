"""QPA probe for a module-category-radical / degree surface (Plan 57). qpa-marked:
skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1. QPA has NO category-radical
or degree verb (RadicalOfModule is the Jacobson radical of ONE module, rad M -- NOT
the category radical rad(X,Y)); there is no LeftDegree/RightDegree/
NilpotencyIndexOfRadical. This documents that honestly and FAILS if one ever appears;
the layer-1 dims dim rad(X,Y) = dim Hom(X,Y) for X !~ Y are the one QPA-checkable
slice, corroborated via HomOverAlgebra. The filtration / index / degrees are covered
by the theory + self-cert oracles in tests/modules/."""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a module-category radical / degree surface.
_ABSENT_NAMES = (
    "LeftDegreeOfIrreducibleMorphism",
    "RightDegreeOfIrreducibleMorphism",
    "NilpotencyIndexOfRadical",
    "RadicalOfModuleCategory",
    "DegreeOfIrreducibleMorphism",
    "RadicalPowerOfModuleCategory",
)


def test_qpa_has_no_category_radical_or_degree_surface():
    lg = session.libgap_handle()
    for name in _ABSENT_NAMES:
        assert not bool(lg.eval('IsBoundGlobal("%s")' % name)), (
            "QPA now exposes %s -- wire a real crosscheck (this assert is the "
            "trip-wire that Plan 57's QPA scope note is stale)" % name)


def test_layer1_dims_match_qpa_hom():
    # rad(X_i, X_j) = Hom(X_i, X_j) for X_i !~ X_j: dim-check the off-diagonal
    # layer-1 against QPA HomOverAlgebra on kA_3 indecomposables (GF(7) -- QPA is a
    # finite-field backend; char 7 > dim clears the trace-form radical char guard).
    from quiverlab import linear_path_algebra
    from quiverlab.fields import GF
    from quiverlab.modules.radical import radical_filtration
    A = linear_path_algebra(3, field=GF(7))
    rf = radical_filtration(A)
    assert rf.is_complete and rf.nilpotency_index == 3
    for i, Xi in enumerate(rf.indecs):
        for j, Xj in enumerate(rf.indecs):
            if i != j:
                A.crosscheck("hom_glue", Xi, Xj).assert_agree()   # P37 Hom crosscheck
                assert rf.layer_dim(i, j, 1) == A.hom(Xi, Xj)     # rad = Hom off-diagonal
