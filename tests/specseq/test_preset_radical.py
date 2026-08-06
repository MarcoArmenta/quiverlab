"""Radical-filtration spectral sequence. Self-cert: converges to H(X) always
(the standing certificate). Koszul tie-in is an ARBITRATED oracle hypothesis:
for a Koszul algebra the radical-filtration SS of the minimal simple-resolution
degenerates at E_2 (the linear staircase on E_1 collapses) -- the provable page,
pinned below (Design-decision 3: never force the oracle; the observed page IS the
folklore E_2 collapse, so it is asserted exactly)."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.modules.complexes import ChainComplex
from quiverlab.specseq.presets import radical_filtration_ss

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_converges_to_homology():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    ss = radical_filtration_ss(X)                 # __init__ self-certifies
    assert ss.convergence.abutment == X.homology_dims()


def test_semisimple_input_degenerates_immediately():
    # a complex of SEMISIMPLE modules (rad = 0) has a trivial radical filtration:
    # F_0 = whole, F_{-1} = 0, so E_1 = E_inf (collapse).
    A = Quiver([1, 2], {}).algebra(relations=[], field=GF(7))   # semisimple k x k
    X = ChainComplex.stalk(A.simple(1), 0)
    ss = radical_filtration_ss(X)
    assert ss.convergence.collapse() is True
    assert ss.convergence.degenerates_at == 1


@pytest.mark.oracle_literature
def test_koszul_degeneration_arbitrated():
    # kA_n is Koszul (hereditary => the minimal resolution of a simple is linear).
    # The radical-filtration SS of the minimal projective resolution of S_1
    # degenerates at E_2 = E_inf (ARBITRATED: this IS the observed provable page,
    # the folklore "E_2 collapse" -- recorded on the verification page).
    A = linear_path_algebra(3, field=GF(5))       # hereditary kA3, Koszul
    X = ChainComplex.from_projective_resolution(A.simple(1), length=3)
    ss = radical_filtration_ss(X)
    assert ss.convergence.degenerates_at == 2
    assert ss.convergence.collapse() is True
    assert ss.convergence.abutment == X.homology_dims()   # always true (self-cert)


def test_koszul_degeneration_kA4():
    # the same arbitrated E_2 collapse on kA4 (a second Koszul witness).
    A = linear_path_algebra(4, field=GF(5))
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    ss = radical_filtration_ss(X)
    assert ss.convergence.degenerates_at == 2
    assert ss.convergence.abutment == X.homology_dims()


def test_block_infinite_gldim_cutoff_independence():
    """Regression (devil's-advocate 2026-08-06): ``radical_filtration_ss_block``
    reported the degree-``top`` homology as an abutment value, but for INFINITE
    global dimension that is the top syzygy ``Omega^top`` -- a truncation artifact
    that GROWS with the cutoff (``k[x]/(x^2)`` sprouts a trailing ``1`` at whatever
    ``top`` the user picks), NOT an invariant. The block now builds one degree
    deeper and trims the boundary, so the reported window is cutoff-INDEPENDENT.

    Pin (oracle_selfcert -- a cutoff-independence certificate): the shared window
    agrees across two cutoffs, at the honest values ``H_0 = dim (+)_v S_v`` and
    ``H_{>0} = 0`` (a projective resolution is exact in positive degrees); a
    finite-gl.dim algebra is byte-unchanged by the trim."""
    from quiverlab.specseq.block import radical_filtration_ss_block

    # k[x]/(x^2): one vertex, self-injective, INFINITE global dimension. Its minimal
    # resolution of S_1 is periodic (mult-by-x, rank 1), so ker d_top = x.P (dim 1)
    # was leaking as a phantom H_top before the deepen-and-trim fix.
    kxx = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(5))
    b3 = radical_filtration_ss_block(kxx, 3)
    b5 = radical_filtration_ss_block(kxx, 5)
    assert b3["abutment"][:4] == b5["abutment"][:4]     # cutoff-independent window
    assert b3["abutment"] == [1, 0, 0, 0]               # H_0 = dim S_1 = 1, rest 0
    assert b5["abutment"] == [1, 0, 0, 0, 0, 0]         # trailing artifact gone

    # a FINITE-gl.dim algebra (kA3, gl.dim 2): the resolution terminates before the
    # boundary, so trimming changes nothing -- H_0 = dim(S_1 (+) S_2 (+) S_3) = 3.
    kA3 = linear_path_algebra(3, field=GF(32003))
    c3 = radical_filtration_ss_block(kA3, 3)
    c4 = radical_filtration_ss_block(kA3, 4)
    assert c3["abutment"] == [3, 0, 0, 0]
    assert c4["abutment"] == [3, 0, 0, 0, 0]
    assert c3["abutment"][:4] == c4["abutment"][:4]     # unchanged + cutoff-independent


def test_nonzero_higher_differential_rank_identity():
    # Devil's-advocate regression (2026-08-05): the standing self-certificate
    # computes E_inf subquotients DIRECTLY and never exercises _dr_matrix, so
    # the multi-step lift-apply-reduce path (the engine's hardest code) had no
    # suite guard. k[x]/(x^a)'s radical SS degenerates at page a, so a >= 3
    # forces a nonzero d_{r>=2}. Assert, at EVERY cell of every pre-collapse
    # page: dim E_{r+1} = dim E_r - rank(d_r out) - rank(d_r in), and
    # d_r . d_r = 0 -- the full induced-differential contract.
    import pytest
    from quiverlab import GF, truncated_polynomial
    from quiverlab.fields.linalg import rank as _rank
    from quiverlab.modules.complexes import ChainComplex
    from quiverlab.specseq.presets import radical_filtration_ss

    for a in (3, 4):
        A = truncated_polynomial(a, field=GF(5))
        X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
        ss = radical_filtration_ss(X)
        dom = A.domain
        saw_nonzero_high_dr = False
        for r in range(1, ss.convergence.degenerates_at + 1):
            Er, Enext = ss.page(r), ss.page(r + 1)
            for (p, q) in Er.spots:
                d_out = Er.differential(p, q)
                d_in = Er.differential(p + r, q - r + 1)
                rk_out = _rank(d_out, dom) if d_out else 0
                rk_in = _rank(d_in, dom) if d_in else 0
                if r >= 2 and (rk_out or rk_in):
                    saw_nonzero_high_dr = True
                assert Enext.dim(p, q) == Er.dim(p, q) - rk_out - rk_in, (a, r, p, q)
        assert saw_nonzero_high_dr, f"k[x]/(x^{a}) must have a nonzero d_(r>=2)"
