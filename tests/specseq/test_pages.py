"""Spectral-sequence pages (Weibel 5.4). Self-certifying: d_r . d_r = 0 at every
page, E_{r+1} = H(E_r, d_r) (rank identity), and canonical reps are stable across
runs. A trivial (one-step) filtration collapses at E_1 = E_inf.

NB (Plan-42 implementation, arbitrated): a ONE-STEP filtration puts every graded
piece in the p=0 COLUMN (E^0_{0,q} = V_q, d^0 = the internal d), so the homology
lands at cells (0, n) -- the plan snippet's (n, 0) was a p/q transposition (the
STUPID filtration is what concentrates on the q=0 row). The self-cert
(E_inf totals == total homology) is layout-independent and is the load-bearing
arbiter; the corrected cell reads (0, n) are pinned here."""
import pytest

from quiverlab import GF, Quiver
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


def test_trivial_filtration_collapses_at_E1():
    # one-step filtration => E^0 = Tot, d^0 = d, E^1 = H(Tot) = E^inf immediately.
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    ss = SpectralSequence(_trivial(X))
    E1 = ss.page(1)
    # E^1 total per degree == H_n(Tot); on a single filtration level p=0, q=n
    for n, h in X.homology_dims().items():
        assert E1.dim(0, n) == h
    # already E_inf: d^1 is zero, page 2 == page 1
    assert ss.page(2).dim(0, 0) == E1.dim(0, 0)


def test_dr_squares_to_zero_and_E_next_is_homology():
    # a genuinely 2-column filtration -- pin d_r.d_r = 0 and the E_{r+1} rank identity.
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(2), length=3)
    # split the filtration by homological degree: F_p Tot_n = Tot_n if n <= p else 0
    dims = {n: X.term(n).dim for n in X.degrees()}
    lo = min(dims)
    filt = {}
    for n, d in dims.items():
        full = [[1 if i == j else 0 for i in range(d)] for j in range(d)]
        # level j (p = lo+j): full basis once p >= n, else empty
        filt[n] = [([] if (lo + j) < n else full) for j in range(len(dims))]
    F = FilteredComplex.from_chain_complex(X, filt)
    ss = SpectralSequence(F)
    for r in (0, 1, 2, 3):
        pg = ss.page(r)
        for (p, q) in pg.spots:
            dr = pg.differential(p, q)
            # d_r . d_r = 0 : compose with the differential out of the target cell
            tgt = (p - r, q + r - 1)
            dr2 = ss.page(r).differential(*tgt)
            from quiverlab.modules import linalg_mod as lm
            comp = lm.matmul(dr2, dr, A.domain) if (dr and dr2 and dr[0] and dr2[0]) else []
            assert not any(not A.domain.is_zero(x) for row in comp for x in row)


def test_E_next_is_homology_of_dr():
    # E_{r+1}^{p,q} = ker d_r / im d_r (the rank identity), across a genuine 2-column
    # filtration. Pins the induced-map coordinates, not just d_r.d_r = 0.
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(2), length=3)
    dims = {n: X.term(n).dim for n in X.degrees()}
    lo = min(dims)
    filt = {}
    for n, d in dims.items():
        full = [[1 if i == j else 0 for i in range(d)] for j in range(d)]
        filt[n] = [([] if (lo + j) < n else full) for j in range(len(dims))]
    ss = SpectralSequence(FilteredComplex.from_chain_complex(X, filt))
    from quiverlab.modules import linalg_mod as lm
    for r in (0, 1, 2):
        pr, pn = ss.page(r), ss.page(r + 1)
        cells = set(pr.spots) | set(pn.spots)
        for (p, q) in cells:
            dr_out = pr.differential(p, q)                         # E^r_{p,q} -> tgt
            dr_in = pr.differential(p + r, q - r + 1)              # src -> E^r_{p,q}
            dim_src = pr.dim(p, q)
            rk_out = lm.mat_rank(dr_out, A.domain) if (dr_out and dr_out[0]) else 0
            rk_in = lm.mat_rank(dr_in, A.domain) if (dr_in and dr_in[0]) else 0
            assert pn.dim(p, q) == dim_src - rk_out - rk_in


def test_canonical_reps_are_reproducible():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    r1 = SpectralSequence(_trivial(X)).page(1)
    r2 = SpectralSequence(_trivial(X)).page(1)
    for (p, q) in r1.spots:
        assert r1[p, q].reps == r2[p, q].reps       # byte-identical reps
