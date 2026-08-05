"""Hyper-Hom of a perfect complex: pinned degreewise against module
Ext/Hom -- two independent engines computing the same numbers."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.complexes import ChainComplex, hyper_hom_dims

pytestmark = pytest.mark.oracle_crossengine


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_projective_stalk_hyper_hom_is_hom():
    A = _a3()
    P1 = A.projective(1)
    for v in (1, 2, 3):
        N = ChainComplex.stalk(A.simple(v), 0)
        X = ChainComplex.stalk(P1, 0)
        hh = hyper_hom_dims(X, N, 0, 3)
        assert hh[0] == A.hom(P1, A.simple(v))
        assert all(hh[n] == 0 for n in (1, 2, 3))


def test_resolution_hyper_hom_computes_ext():
    # THE pin: X = proj resolution of M as a perfect complex, Y = stalk N
    # => H^n(Hom(X, Y)) = Ext^n(M, N) for n within the resolution window.
    A = _a3()
    for v in (1, 2, 3):
        M, N = A.simple(v), A.simple(1)
        X = ChainComplex.from_projective_resolution(M, length=5)
        got = hyper_hom_dims(X, ChainComplex.stalk(N, 0), 0, 4)
        for n in range(0, 5):
            assert got[n] == A.ext(M, N, n), (v, n)


def test_shift_bookkeeping():
    A = _a3()
    M, N = A.simple(2), A.simple(1)
    X = ChainComplex.from_projective_resolution(M, length=4)
    Y = ChainComplex.stalk(N, 0)
    base = hyper_hom_dims(X, Y, 0, 3)
    shifted = hyper_hom_dims(X, Y.shift(1), -1, 2)
    assert all(shifted[n - 1] == base[n] for n in range(0, 3))


def test_nonperfect_source_refused():
    A = _a3()
    X = ChainComplex.stalk(A.simple(2), 0)     # simple: not projective
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="perfect"):
        hyper_hom_dims(X, X, 0, 1)
