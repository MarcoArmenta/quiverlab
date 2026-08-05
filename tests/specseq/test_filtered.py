"""FilteredComplex: validation (closed under d, exhaustive) + total homology.
Self-certifying (d.d=0 gate, filtration-subcomplex gate, rank identity)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.filtered import FilteredComplex

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _trivial_filt(terms):
    # the one-step filtration F_0 V_n = V_n (a valid, exhaustive, degenerate filt)
    return {n: [[[1 if i == j else 0 for i in range(d)] for j in range(d)]]
            for n, d in terms.items()}, 0


def test_forget_chaincomplex_total_homology_matches():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    dims = {n: X.term(n).dim for n in X.degrees()}
    filt, lo = _trivial_filt(dims)
    F = FilteredComplex.from_chain_complex(X, {n: filt[n] for n in dims})
    # the forgotten complex has the SAME homology as X (module structure ignored)
    assert F.total_homology_dims() == X.homology_dims()


def test_non_closed_filtration_is_refused():
    # a filtration piece F_0 V_1 whose image under d escapes F_0 V_0 must raise.
    A = _a3()
    P1, S1 = A.projective(1), A.simple(1)
    from quiverlab.modules.resolution import projective_cover
    Q0, d0, _ = projective_cover(S1)              # d0: Q0 -> S1, surjective
    terms = {1: Q0.dim, 0: S1.dim}
    dmats = {1: d0}
    # F_0 V_1 = a 1-dim piece mapping to a nonzero vector, but F_0 V_0 = {0}:
    bad = {1: [[[1] + [0] * (Q0.dim - 1)]], 0: [[]]}   # F_0 V_0 empty -> not closed
    with pytest.raises(QuiverlabError, match="closed|subcomplex|F_"):
        FilteredComplex(terms, dmats, bad, lo=0, dom=A.domain)


def test_non_exhaustive_filtration_is_refused():
    A = _a3()
    terms = {0: 2}
    dmats = {}
    filt = {0: [[[1, 0]]]}                          # top level spans only a line
    with pytest.raises(QuiverlabError, match="exhaustive|span"):
        FilteredComplex(terms, dmats, filt, lo=0, dom=A.domain)


def test_dd_nonzero_is_refused():
    A = _a3()
    # d_1 = d_2 = identity on a 1-dim space => d.d != 0
    terms = {2: 1, 1: 1, 0: 1}
    dmats = {2: [[1]], 1: [[1]]}
    filt, lo = _trivial_filt(terms)
    with pytest.raises(QuiverlabError, match="d.*d|differential"):
        FilteredComplex(terms, dmats, filt, lo, A.domain)
