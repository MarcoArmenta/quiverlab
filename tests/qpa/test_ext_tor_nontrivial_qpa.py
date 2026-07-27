"""QPA (GAP) as an EXTERNAL oracle for module Ext^n(M, N) where BOTH M and N are
genuinely non-projective, non-injective, non-simple modules (Plan 34 -- Marco's
release feedback).

The standing QPA Ext bridge (``qpa/crosscheck.py::crosscheck_module_ext``) checks
only the SELF-Ext Ext^*(M, M) via ``ExtAlgebraGenerators``.  This battery instead
compares dim Ext^n(M, N) for DISTINCT explicit-matrix interior modules, computing
the right-hand side in QPA with its OWN homological machinery -- ``HomOverAlgebra``
for Ext^0, then ``ExtOverAlgebra(NthSyzygy(M, n-1), N)`` for Ext^n by dimension
shifting -- so a match is a genuine external certificate.  (QPA ships no native
Tor; the Tor <-> Ext duality route is covered in ``test_tor_qpa.py``.  Here the
NEW coverage is distinct-module Ext on interior modules, including the periodic
higher-degree Ext of self-injective algebras.)

This test does NOT edit ``qpa/crosscheck.py`` or ``qpa/scripts.py`` -- it uses the
read-only script/module helpers and the live session, mirroring ``test_tor_qpa.py``.

qpa-marked (by directory): skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.qpa_module import graded_form
from quiverlab.qpa import scripts, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


# --------------------------------------------------------------------------- #
# Interior-module constructors + a cheap, rigorous non-triviality certificate
# --------------------------------------------------------------------------- #
def _interval(A, i, j, n):
    """The uniserial interval [i,j] over the equioriented kA_n by explicit matrices."""
    dim = j - i + 1
    dimvec = {v: (1 if i <= v <= j else 0) for v in range(1, n + 1)}
    maps = {}
    for k in range(1, n):
        blk = [[0] * dim for _ in range(dim)]
        if i <= k <= j - 1:
            blk[(k + 1) - i][k - i] = 1
        maps[f"a{k}"] = blk
    return A.module(dimvec, maps, name=f"[{i},{j}]")


def _cyclic_nakayama(n, loewy, field):
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i % n + 1) for i in range(1, n + 1)}
    rels = ["*".join(f"a{((i - 1 + t) % n) + 1}" for t in range(loewy))
            for i in range(1, n + 1)]
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def _square(field):
    """The commutative square kQ/(ab-cd): non-hereditary, gl.dim 2."""
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _p_mod_soc(P):
    """P / soc P (soc = intersection of the arrow kernels)."""
    from quiverlab.modules import linalg_mod as lm
    from quiverlab.modules.radtopsoc import _intersect, quotient
    dom = P.domain
    inter = None
    for a in P.algebra.quiver.arrows:
        ker = lm.kernel_columns(P.action[a], dom)
        inter = ker if inter is None else _intersect(inter, ker, dom)
    return quotient(P, inter or [], name=f"{P.name}/soc")


def _dv(A, M):
    d = M.dimension_vector()
    return [d[v] for v in A.quiver.vertices]


def _assert_interior(A, M):
    """Certify M is non-simple, non-projective, non-injective.  Indecomposable + a
    dimension-vector distinct from every P_v and I_v is a rigorous refutation
    (isomorphism preserves the dimension vector), and cheaper than a full iso search
    in the slow QPA leg.  Over the fields used here (QQ; GF(5) with dim M < 5) the
    is_indecomposable certificate is rigorous."""
    assert M.dim > 1, f"{M.name} is simple"
    assert M.is_indecomposable(), f"{M.name} is decomposable"
    dvM = _dv(A, M)
    for v in A.quiver.vertices:
        assert dvM != _dv(A, A.projective(v)), f"{M.name} could be P_{v}"
        assert dvM != _dv(A, A.injective(v)), f"{M.name} could be I_{v}"


# --------------------------------------------------------------------------- #
# The QPA Ext oracle (dimension shifting via NthSyzygy) + our-vs-QPA crosscheck
# --------------------------------------------------------------------------- #
def _qpa_ext_dims(A, M, N, top):
    """dim Ext_A^0..top(M, N) in QPA (M, N both RIGHT A-modules): Ext^0 = dim Hom;
    Ext^n = dim Ext^1(Omega^{n-1} M, N) via NthSyzygy (mirrors test_tor_qpa.py)."""
    dvM, arrM = graded_form(M)
    dvN, arrN = graded_form(N)
    base = scripts.quiver_and_algebra_script(A)
    base += "\n" + scripts.module_decl(A, dvM, arrM, "MM")
    base += "\n" + scripts.module_decl(A, dvN, arrN, "NN")
    dims = [int(session.run(base + "\nhh := HomOverAlgebra(MM, NN);;\nLength(hh);"))]
    for n in range(1, top + 1):
        script = (base + f"\nsy := NthSyzygy(MM, {n - 1});;"
                  "\ne := ExtOverAlgebra(sy, NN);;\nLength(e[2]);")
        dims.append(int(session.run(script)))
    return dims


def _crosscheck(A, M, N, top):
    ours = [A.ext(M, N, n) for n in range(top + 1)]
    qpa = _qpa_ext_dims(A, M, N, top)
    assert ours == qpa, (M.name, N.name, ours, qpa)
    return ours


# --------------------------------------------------------------------------- #
# Hereditary kA_5: distinct interior intervals, our Ext vs QPA
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("field", [QQ, GF(5)])
@pytest.mark.parametrize("ij,kl", [((2, 3), (2, 3)), ((2, 3), (3, 4)), ((3, 4), (2, 3))])
def test_ext_interior_intervals_kA5_vs_qpa(field, ij, kl):
    """kA_5 interior intervals [2,3],[3,4],[2,4] are neither projective ([i,5]) nor
    injective ([1,j]) nor simple.  Our dim Ext^n(M,N) matches QPA's ExtOverAlgebra
    dimension-shifting computation; the [2,3]->[3,4] pair carries the AR extension
    (Ext^1 = 1)."""
    A = linear_path_algebra(5, field=field)
    M = _interval(A, *ij, 5)
    N = _interval(A, *kl, 5)
    _assert_interior(A, M)
    _assert_interior(A, N)
    _crosscheck(A, M, N, 3)


# --------------------------------------------------------------------------- #
# External anchor for the Python-battery commutative-square Tor pin
# --------------------------------------------------------------------------- #
def test_ext_square_duality_pair_vs_qpa():
    """Anchor the commutative-square Tor pin against an EXTERNAL oracle.  The Plan-34
    Python battery pins Tor_n(rad P_1, D(P_1/soc P_1)) = [2,0,0,0] on the non-hereditary
    square kQ/(ab-cd) via the duality Tor_n(M,N) = Ext^n(M, DN); that Ext side is
    Ext^n(rad P_1, P_1/soc P_1) -- BOTH RIGHT interior modules.  Since the Python
    anchor resolves only M (shared between its Tor and its Ext^*(M,DN)), the pin
    otherwise rests on a single engine path; QPA recomputes this Ext with its own
    dimension-shifting machinery and must agree.  One pair, GF(5), short session; the
    dim is field-independent here (matches the QQ Python pin, verified [2,0,0,0])."""
    S = _square(field=GF(5))
    M = S.projective(1).radical()               # rad P_1, interior (dim vector (0,1,1,1))
    N = _p_mod_soc(S.projective(1))             # P_1/soc P_1, interior ((1,1,1,0))
    _assert_interior(S, M)
    _assert_interior(S, N)
    assert _crosscheck(S, M, N, 2) == [2, 0, 0]


# --------------------------------------------------------------------------- #
# Self-injective algebras: periodic higher-degree Ext between interior modules
# --------------------------------------------------------------------------- #
def test_ext_selfinjective_kx3_higher_degree_vs_qpa():
    """k[x]/(x^3) is self-injective (projective = injective), so the length-2 nilpotent
    module M2 (x -> [[0,0],[1,0]]) is non-projective => non-injective, and non-simple.
    Its self-Ext is periodic -- dim Ext^n(M2, M2) = [2,1,1,1,1] -- exercising QPA's
    NthSyzygy on an infinite-pd module past degree 1.  Checked over GF(5) and QQ."""
    for field in (GF(5), QQ):
        B = truncated_polynomial(3, field=field)
        assert B.is_selfinjective()                      # projective = injective here
        M2 = B.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M2")
        assert M2.dim > 1 and M2.is_indecomposable()
        assert not M2.is_isomorphic(B.projective(1))     # non-projective => non-injective
        assert _crosscheck(B, M2, M2, 4) == [2, 1, 1, 1, 1]


def test_ext_selfinjective_nakayama_higher_degree_vs_qpa():
    """A self-injective cyclic Nakayama kZ_3/J^3: the interior length-2 uniserial
    rad P_1 (non-projective => non-injective, non-simple) has periodic self-Ext
    dim Ext^n = [1,0,1,0,1] to degree 4, matching QPA."""
    B = _cyclic_nakayama(3, 3, field=GF(5))
    assert B.is_selfinjective()
    M = B.projective(1).radical()
    _assert_interior(B, M)
    assert _crosscheck(B, M, M, 4) == [1, 0, 1, 0, 1]
