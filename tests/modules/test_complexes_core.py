"""ChainComplex core: validation, shift/truncate identities, homology.
Self-certifying (d.d=0 gates, rank identities) + cross-engine (resolution
round-trip reproduces the module)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_stalk_homology_is_the_module():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.stalk(S1, degree=0)
    assert X.homology_dims() == {0: 1}
    assert X.homology(0).dimension_vector() == S1.dimension_vector()


def test_resolution_roundtrip_homology_concentrated_in_zero():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.from_projective_resolution(S1, length=4)
    h = X.homology_dims()
    assert h.get(0) == S1.dim
    assert all(v == 0 for n, v in h.items() if n != 0)


def test_dd_nonzero_is_refused():
    A = _a3()
    P1 = A.projective(1)
    idm = P1.identity_hom().matrix
    with pytest.raises(QuiverlabError, match="d.*d|square"):
        ChainComplex({1: P1, 0: P1, -1: P1}, {1: idm, 0: idm})


def test_non_module_map_differential_refused():
    A = _a3()
    P1, S1 = A.projective(1), A.simple(1)
    bad = [[1, 1]]                        # not an intertwiner (P37 pin)
    with pytest.raises(QuiverlabError, match="module map"):
        ChainComplex({1: P1, 0: S1}, {1: bad})


def test_shift_moves_degrees_and_signs():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.from_projective_resolution(S1, length=2)
    Y = X.shift(3)
    assert Y.homology_dims().get(3) == S1.dim
    Z = Y.shift(-3)
    assert Z.homology_dims() == X.homology_dims()
    # d.d = 0 still validated after odd shift (sign flip is consistent)
    ChainComplex({n: Y.term(n) for n in Y.degrees()},
                 {n: Y.differential(n).matrix for n in Y.degrees()
                  if Y.term(n - 1).dim and Y.term(n).dim})


def test_truncate_brutal():
    A = _a3()
    S1 = A.simple(1)
    X = ChainComplex.from_projective_resolution(S1, length=4)
    T = X.truncate(0, 1)
    assert set(T.degrees()) <= {0, 1}


def test_mixed_sides_refused():
    A = _a3()
    R = A.simple(1)
    L = A.simple(1, side="left")
    with pytest.raises(QuiverlabError):
        ChainComplex({0: R, 1: L}, {1: [[0]]})


def test_toplevel_exports_and_algebra_wrapper():
    # public surface (Plan 39 Task 6): from quiverlab import ChainComplex, ChainMap
    import quiverlab
    from quiverlab import ChainComplex as CC
    from quiverlab import ChainMap as CM
    assert CC is ChainComplex
    assert CM is quiverlab.modules.complexes.ChainMap
    # A.chain_complex(...) convenience wrapper round-trips into a ChainComplex
    A = _a3()
    S1 = A.simple(1)
    X = A.chain_complex({0: S1}, {})
    assert isinstance(X, ChainComplex)
    assert X.homology_dims() == {0: 1}
    assert X.homology(0).dimension_vector() == S1.dimension_vector()
