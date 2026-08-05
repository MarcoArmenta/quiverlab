"""hyper_hom_basis reifies H^n(Hom^.(X, Y)) as chain maps X -> Y[n]. Self-cert:
every reified class is a valid chain map (check=True), the count equals
hyper_hom_dims, and degree-0 composition matches End on stalks. Cross-engine:
degree-n classes count Ext^n on a projective-resolution source."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex, identity_chain_map
from quiverlab.modules.complexes import hyper_hom_dims
from quiverlab.derived.homs import hyper_hom_basis

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


@xeng
def test_basis_count_equals_ext_on_resolution():
    A = _a3()
    for v in (1, 2, 3):
        M, N = A.simple(v), A.simple(1)
        X = ChainComplex.from_projective_resolution(M, length=5)
        Y = ChainComplex.stalk(N, 0)
        for n in range(0, 5):
            basis = hyper_hom_basis(X, Y, n)
            assert len(basis) == A.ext(M, N, n)                 # count = Ext^n
            assert len(basis) == hyper_hom_dims(X, Y, n, n)[n]  # count = dim H^n


@selfcert
def test_every_reified_class_is_a_chain_map():
    # ChainMap(check=True) inside hyper_hom_basis is the self-cert; re-assert the
    # target is exactly Y.shift(n) and each square commutes (rebuild with check).
    A = _a3()
    M, N = A.simple(2), A.simple(1)
    X = ChainComplex.from_projective_resolution(M, length=4)
    Y = ChainComplex.stalk(N, 0)
    for n in range(0, 4):
        for f in hyper_hom_basis(X, Y, n):
            assert f.tgt.degrees() == Y.shift(n).degrees()
            from quiverlab.modules.complexes import ChainMap
            ChainMap(f.src, f.tgt, {k: f.component(k) for k in f.src.degrees()},
                     check=True)                                # squares commute


@selfcert
def test_degree0_endo_composition_and_identity():
    # degree-0 hyper-Hom of a projective stalk with itself = End_A(P); ChainMap.then
    # composes them; the identity chain map is one of the classes (up to homotopy).
    A = _a3()
    P1 = A.projective(1)
    X = ChainComplex.stalk(P1, 0)
    basis = hyper_hom_basis(X, X, 0)
    assert len(basis) == A.hom(P1, P1)
    idX = identity_chain_map(X)
    # then() type-checks and stays a chain map:
    comp = basis[0].then(idX)
    assert comp.component(0) == basis[0].component(0)           # f . id = f


@selfcert
def test_nonperfect_source_refused():
    A = _a3()
    X = ChainComplex.stalk(A.simple(2), 0)      # simple: not projective
    with pytest.raises(QuiverlabError, match="perfect"):
        hyper_hom_basis(X, X, 0)
