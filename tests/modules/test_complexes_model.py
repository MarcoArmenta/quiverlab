"""Projective models: certified quasi-iso by construction; hyper-Ext on
stalks reproduces module Ext (cross-engine); two-term complex LES pin."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.complexes import (ChainComplex, ChainMap,
                                         hyper_ext_dims, projective_model)

pytestmark = pytest.mark.oracle_crossengine


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_model_of_stalk_is_certified_and_computes_ext():
    A = _a3()
    for v in (1, 2, 3):
        M = A.simple(v)
        X = ChainComplex.stalk(M, 0)
        P, eps = projective_model(X, length=5)
        assert P.is_perfect() and eps.is_quasi_iso()
        got = hyper_ext_dims(X, ChainComplex.stalk(A.simple(1), 0), 0, 4,
                             length=6)
        for n in range(0, 5):
            assert got[n] == A.ext(M, A.simple(1), n), (v, n)


def test_model_of_two_term_complex():
    # X = [P1 --f--> S1] (f the cover, degrees 1,0): quasi-iso to the stalk
    # of ker f shifted -- hyper-Ext must match Ext of rad P1 with a shift.
    A = _a3()
    P1, S1 = A.projective(1), A.simple(1)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)          # Q0 = P1 here
    X = ChainComplex({1: Q0, 0: S1}, {1: d0})  # NOTE: check orientation --
    # this complex has homology ker d0 = rad P1 in degree 1 and 0 in degree 0
    P, eps = projective_model(X, length=4)
    assert eps.is_quasi_iso()
    R = Q0.radical()
    got = hyper_ext_dims(X, ChainComplex.stalk(A.simple(2), 0), -1, 3)
    # SHARPENED SHIFT IDENTITY (replacing the plan's `or True` placeholder).
    # Derivation (on paper): X = [P1 --cover--> S1] in homological degrees (1, 0)
    # has H_1(X) = ker(cover) = rad P1 = R and H_0(X) = coker(cover) = 0, so X is
    # quasi-isomorphic to R placed in degree 1, i.e. X ~ R[1] (stalk(R,0) shifted
    # up by 1). The Hom total complex shifts the SOURCE contravariantly:
    #   Hom^m(X'[1], Y) = prod_p Hom(X'_{p-1}, Y_{p-m}) = Hom^{m-1}(X', Y),
    # so hyper_hom(X'[1], Y)[n] = hyper_hom(X', Y)[n-1]. With X ~ R[1]:
    #   hyper_ext(X, N)[n] = Hom_{D^b}(R[1], N[n]) = Ext^{n-1}_A(R, N).
    # Concretely Hom^n(model(X), stalk N) = Hom(Q^R_{n-1}, N), so
    # H^n = Ext^{n-1}(R, N); in particular got[0] = Ext^{-1} = 0.
    # (The plan's tentative n+1 is WRONG -- e.g. got[1] = Ext^0(R, S2) = 1, but
    #  Ext^2(R, S2) = 0; the placeholder `or True` was masking exactly this.)
    for n in range(0, 3):
        expected = A.ext(R, A.simple(2), n - 1) if n - 1 >= 0 else 0
        assert got[n] == expected, (n, got[n], expected)


def test_model_length_window_is_honest():
    # dims past the model window must not be silently wrong: the function
    # raises when [lo, hi] exceeds what `length` certifies.
    A = _a3()
    X = ChainComplex.stalk(A.simple(1), 0)
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="length|window"):
        hyper_ext_dims(X, X, 0, 30, length=3)
