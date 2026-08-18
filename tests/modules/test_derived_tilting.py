"""Tilting-complex verifier + End(T). Self-cert: T=A (projective stalks) is tilting
with End(T) ~ A (corner-Cartan). Literature: the kA2 APR tilt T = P1 (+) S1 is
tilting, End(T) has the reoriented-A2 Cartan. Negative: a missing summand fails
generation; X (+) X[1] fails rigidity."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.tilting import (is_tilting_complex, end_algebra_of_complex,
                                       two_term_silting_from_presentation)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a2():
    return linear_path_algebra(2, field=QQ)      # kA2, 1->2, hereditary


@selfcert
def test_regular_module_is_tilting_and_end_is_A():
    # T = A = (+)_v stalk(P_v): the trivial tilting complex; End(T) ~ A.
    A = _a2()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_tilting_complex(T)
    assert rep.is_tilting and rep.rigid and rep.generates
    assert rep.window == (0, 0)                  # width-0 summands: Ext^{!=0}(P,P)=0
    E = end_algebra_of_complex(T)
    from quiverlab.derived.tilting import corner_cartan_of_complex
    assert corner_cartan_of_complex(T) == A.cartan_matrix()   # End(A_A) ~ A oracle


@lit
def test_ka2_apr_tilt():
    # DERIVATION (kA2, arrow a:1->2; right modules). Indecomposables: S1=(1,0),
    # S2=P2=(0,1), P1=[1,2]=(1,1). Non-projective: S1 (pd 1: 0->P2->P1->S1->0).
    # APR tilt at the sink 2: T = P1 (+) tau^{-1}(P2) = P1 (+) S1 (the unique AR
    # sequence 0->S2->P1->S1->0 gives tau^{-1}(S2)=S1). Checks: pd P1=0, pd S1=1;
    # Ext^1(S1,P1)=0 (coker(Hom(P1,P1)->Hom(P2,P1)) = coker(k->k iso)=0),
    # Ext^1(S1,S1)=0 (Hom(P2,S1)=0); Ext^1(P1,-)=0. Summands=2=#simples. So T is
    # tilting. End(T): Hom(P1,P1)=Hom(S1,S1)=k, Hom(P1,S1)=k (P1->>S1), Hom(S1,P1)=0.
    # CORNER-CARTAN ORIENTATION (arbitrated, see below): under the convention that makes
    # the theorem End(A_A)=A hold -- corner[i][j] = dim e_i End(T) e_j = dim Hom(T_j, T_i),
    # exactly P37's regular_corner_dims == cartan_matrix -- the APR corner-Cartan is
    # cartan(A^op) = [[1,0],[1,1]], the genuine REORIENTED-A2 (End(T) of the APR tilt is
    # the reflected algebra A^op, verified: End(T) has dim 3 = dim kA2 and this corner).
    # (The plan's hand-derivation wrote the transpose [[1,1],[0,1]] = A's Cartan; that is
    # A itself, not A^op -- inconsistent with the theorem-anchored End(A_A)=A pin, which
    # fixes the orientation. Systematic-debugging deviation, Plan-43 Task 3.)
    A = _a2()
    P1 = A.projective(1)
    S1 = A.simple(1)
    T = [ChainComplex.stalk(P1, 0),
         ChainComplex.from_projective_resolution(S1, length=2)]   # S1 as a perfect cx
    rep = is_tilting_complex(T)
    assert rep.is_tilting and rep.rigid and rep.generates
    from quiverlab.derived.tilting import corner_cartan_of_complex
    assert corner_cartan_of_complex(T) == [[1, 0], [1, 1]]        # = cartan(A^op), reoriented A2
    assert corner_cartan_of_complex(T) == A.opposite().cartan_matrix()


@selfcert
def test_missing_summand_fails_generation():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(1), 0)]     # one summand, two simples
    rep = is_tilting_complex(T)
    assert rep.generates is False and rep.is_tilting is False     # g-matrix not square


@selfcert
def test_shifted_copy_fails_rigidity():
    # T = X (+) X[1] with X = stalk(P1): Hom_{D^b}(X, X[1][-1]) = End(X) != 0, so
    # rigidity fails at n = -1 (and symmetrically at n = +1).
    A = _a2()
    X = ChainComplex.stalk(A.projective(1), 0)
    T = [X, X.shift(1)]
    rep = is_tilting_complex(T)
    assert rep.rigid is False and rep.is_tilting is False
    assert rep.window[0] <= -1 <= rep.window[1]


@selfcert
def test_two_term_silting_from_presentation():
    A = _a2()
    cx, rep = two_term_silting_from_presentation(A.simple(1))
    assert set(cx.degrees()) <= {0, 1} and cx.is_perfect()
    assert rep.rigid                                 # a 2-term silting object is rigid


@selfcert
def test_regular_object_is_tilting_over_nonunimodular_cartan():
    # Plan 67 Task 0 regression. The g-matrix must live in K0(K^b proj) = (+)_v Z[P_v],
    # the PROJECTIVE basis -- NOT the composition-factor basis: the regular object
    # A = (+)_v P_v is a tilting complex over EVERY algebra (AI Ex 2.2). The old _chi
    # g-matrix gave det(Cartan) and FAILED here (det C = 0 / 2) -- a systematic
    # false-negative on non-unimodular Cartan. The kA_n suite above only exercises
    # UNIMODULAR Cartan (hereditary, det C = 1), where _chi and g_proj COINCIDE
    # (det _chi = det C . det g_proj = 1 . det g_proj) -- which is exactly why the bug
    # survived P43's green suite.
    import sympy as sp
    from quiverlab import GF, NakayamaAlgebra, Quiver
    from quiverlab.fields import QQ
    # Example 2.47: 1 <=> 2 with ab = ba = 0; det Cartan = 0 (NON-unimodular).
    A = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)
    assert int(sp.Matrix(A.cartan_matrix()).det()) == 0            # non-unimodular
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_tilting_complex(T)
    assert rep.is_tilting is True and rep.generates is True and rep.det in (1, -1)
    # self-injective Nakayama kZ3/J2: det Cartan = 2.
    B = NakayamaAlgebra([2, 2, 2], cyclic=True, field=GF(32003))
    assert int(sp.Matrix(B.cartan_matrix()).det()) == 2
    TB = [ChainComplex.stalk(B.projective(v), 0) for v in B.quiver.vertices]
    repB = is_tilting_complex(TB)
    assert repB.is_tilting is True and repB.det in (1, -1)


@selfcert
def test_generation_is_certified_only_on_two_term_regular_and_apr():
    # Plan 67 fix round (MAJOR H1): the tilting verifier's generation is honest
    # THREE-VALUED. It is CERTIFIED (`is_tilting is True`) only where a completion theorem
    # reaches -- a 2-term self-orthogonal K0-basis object is 2-term silting (IJY/AIR) and,
    # two-sided rigid, tilting. That covers the regular object A (width 0) and the APR tilt
    # (width 1). Both must carry generation == "certified".
    A = _a2()
    regular = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    r = is_tilting_complex(regular)
    assert r.is_tilting is True and r.generation == "certified"
    apr = [ChainComplex.stalk(A.projective(1), 0),
           ChainComplex.from_projective_resolution(A.simple(1), length=2)]
    ra = is_tilting_complex(apr)
    assert ra.is_tilting is True and ra.generation == "certified"
    # a missing summand is not a K0-basis at all: generation == "no", is_tilting False.
    rm = is_tilting_complex([ChainComplex.stalk(A.projective(1), 0)])
    assert rm.generation == "no" and rm.generates is False and rm.is_tilting is False


@selfcert
def test_wide_rigid_k0_basis_is_unknown_not_hard_true():
    # Plan 67 fix round (MAJOR H1): a WIDE (non-2-term) two-sided-RIGID K0-basis complex is
    # Rickard-OPEN -- "rigid + (#summands = rk K0) => tilting" is exactly Rickard's rank
    # QUESTION (unresolved; thick subcategories are not K0-classified -- Krah phantom). The
    # verifier must NOT return a possibly-unsound hard True from bare K0; it returns the
    # honest "unknown" (generation == "k0_necessary_only"), even though THIS particular
    # object does in fact generate (provable by brutal-truncation triangles -- an argument
    # K0 alone cannot run). The pre-fix P43 surface returned a hard True here.
    # Λ = k[1->2->3]/(a*b), gl.dim 2; res(S1) = [P3 -> P2 -> P1] is a genuine width-2
    # perfect complex; T = res(S1) (+) P1 (+) P2 is rigid with det(g_proj) = 1.
    from quiverlab import GF
    L = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(32003))
    res_S1 = ChainComplex.from_projective_resolution(L.simple(1), length=2)
    T = [res_S1, ChainComplex.stalk(L.projective(1), 0),
         ChainComplex.stalk(L.projective(2), 0)]
    assert res_S1.degrees() == [0, 1, 2]                 # genuinely width-2 (non-2-term)
    rep = is_tilting_complex(T)
    assert rep.rigid is True and rep.generates is True and rep.det in (1, -1)
    assert rep.generation == "k0_necessary_only"
    assert rep.is_tilting == "unknown"                   # NOT a hard True
    # CONSISTENCY (the whole point of H1): the silting verifier gives the SAME "unknown"
    # on this same input -- the tilting rung and the K0-basis-only rung are reconciled.
    from quiverlab.derived.silting import is_silting_object
    srep = is_silting_object(T)
    assert srep.is_silting == "unknown"
    assert srep.is_presilting is True and srep.k0_basis is True
