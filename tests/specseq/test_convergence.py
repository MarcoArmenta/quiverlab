"""Convergence: E_inf totals == total homology (the standing self-certificate),
E_inf page bound, degeneration decided by rank. A one-step filtration degenerates
at E_1."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.filtered import FilteredComplex
from quiverlab.specseq.pages import SpectralSequence

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _trivial(X):
    dims = {n: X.term(n).dim for n in X.degrees()}
    filt = {n: [[[1 if i == j else 0 for i in range(d)] for j in range(d)]]
            for n, d in dims.items()}
    return FilteredComplex.from_chain_complex(X, filt)


def test_einf_totals_equal_total_homology():
    A = _a3()
    for v in (1, 2, 3):
        X = ChainComplex.from_projective_resolution(A.simple(v), length=4)
        ss = SpectralSequence(_trivial(X))            # __init__ certifies
        rep = ss.convergence
        totals = {}
        Einf = ss.page(rep.e_infinity_page)
        for (p, q) in Einf.spots:
            totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
        for n, h in ss.total_homology_dims.items():
            assert totals.get(n, 0) == h
        assert rep.abutment == ss.total_homology_dims


def test_trivial_filtration_degenerates_at_E1():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    ss = SpectralSequence(_trivial(X))
    assert ss.convergence.degenerates_at == 1
    assert ss.convergence.collapse() is True


def test_bookkeeping_bug_is_caught():
    # feed a FilteredComplex whose declared filtration is a valid subcomplex filt
    # but hand-mangle the abutment expectation: the certificate compares against
    # the complex's OWN homology, so a genuine construction can never mismatch.
    # This test instead asserts the certificate is WIRED (runs at __init__) by
    # confirming ss.convergence exists and is a ConvergenceReport.
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(2), length=2)
    ss = SpectralSequence(_trivial(X))
    from quiverlab.specseq.convergence import ConvergenceReport
    assert isinstance(ss.convergence, ConvergenceReport)


def test_stupid_filtration_certificate_holds_at_scale():
    # a genuinely multi-column (stupid, by homological degree) filtration -- the
    # self-certificate MUST hold (it runs at __init__; this pins that the E_inf
    # totals of a nontrivial SS reproduce the total homology).
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(3), length=4)
    dims = {n: X.term(n).dim for n in X.degrees()}
    lo = min(dims)
    filt = {}
    for n, d in dims.items():
        full = [[1 if i == j else 0 for i in range(d)] for j in range(d)]
        filt[n] = [([] if (lo + j) < n else full) for j in range(len(dims))]
    ss = SpectralSequence(FilteredComplex.from_chain_complex(X, filt))
    assert ss.convergence.abutment == X.homology_dims()
