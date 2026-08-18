"""Silting-object verifier (Plan 67 / AI Def 2.1). Self-cert: T = A is silting (in fact
tilting) with the positive window (1,0); T = P2 (+) P1[1] over kA2 is a GENUINE
silting-but-not-tilting object (presilting via the clean positive window, but the P43
tilting verifier fails rigidity at n = -1) -- the silting-vs-tilting separation; a
single projective P1 over kA2 is presilting but NOT silting (#summands != #simples) --
the presilting-vs-silting (generation) separation; a single projective stalk over a
local algebra is silting (Thm 2.26). SMOKE (would have caught the Task-0 K0-basis bug):
the regular object over the NON-unimodular-Cartan Example 2.47 and a self-injective
Nakayama both verify silting=True. Literature: k[x]/(x^2) silting = shifts. Honest
three-valued: the K0-basis-only case reports 'unknown', never a silent True.

WHY the kA2/local oracles ALONE mask the K0-basis bug (Task 0): kA2 is hereditary so its
Cartan is unimodular (det = 1), where g_proj and the old composition-factor _chi
coincide; a LOCAL algebra hits the #summands==1 'local' branch that never touches the
K0 g-matrix. Only a non-unimodular-Cartan, non-local input (Example 2.47, kZ3/J2)
exercises the projective-basis g_proj -- hence the two smoke tests below."""
import pytest

from quiverlab import GF, NakayamaAlgebra, Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.silting import is_silting_object, co_t_structure_of

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a2():
    return linear_path_algebra(2, field=QQ)          # kA2, 1->2, hereditary


def _comm_square():
    # 1 <=> 2 with ab = ba = 0 (AI Example 2.47) -- det Cartan = 0 (NON-unimodular).
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


@selfcert
def test_regular_module_is_silting_and_tilting():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_silting_object(T)
    assert rep.is_presilting and rep.k0_basis
    assert rep.is_silting is True
    assert rep.generation_certified_by.startswith("tilting")
    assert rep.window == (1, 0)                       # no positive shifts to scan
    assert rep.det in (1, -1)                         # g_proj identity (projective basis)


@selfcert
def test_silting_not_tilting_separation():
    # T = P2 (+) P1[1] over kA2 is a GENUINE silting object that is NOT tilting.
    # Presilting: the positive window (1,1) is clean (Hom_{D^b}(T, T[1]) = 0). It is
    # silting via the 2-term (IJY) rung. But it is NOT tilting: the P43 tilting verifier
    # fails rigidity at n = -1 (Hom_{D^b}(T, T[-1]) = End(P1)-arrow != 0). This is the
    # vanishing-condition separation: silting checks only n > 0, tilting checks n != 0.
    A = _a2()
    P1 = ChainComplex.stalk(A.projective(1), 0)
    P2 = ChainComplex.stalk(A.projective(2), 0)
    T = [P2, P1.shift(1)]
    rep = is_silting_object(T)
    assert rep.is_presilting is True                  # positive window is clean
    assert rep.window == (1, 1)
    assert rep.is_silting is True                     # IJY 2-term rung
    assert rep.generation_certified_by.startswith("2-term")
    # but it is NOT tilting -- the negative window is nonzero:
    from quiverlab.derived.tilting import is_tilting_complex
    assert is_tilting_complex(T).rigid is False       # fails rigidity at n = -1


@selfcert
def test_single_projective_is_presilting_not_silting():
    # A single indecomposable P1 over kA2: PRESILTING (positive window (1,0) is
    # vacuously clean) but NOT silting -- #summands (1) != #simples (2), so the classes
    # cannot be a K0-basis (k0_basis False). The presilting-vs-silting separation.
    A = _a2()
    rep = is_silting_object([ChainComplex.stalk(A.projective(1), 0)])
    assert rep.is_presilting is True
    assert rep.window == (1, 0)
    assert rep.k0_basis is False and rep.is_silting is False


@selfcert
def test_verifier_smoke_nonunimodular_cartan():
    # SMOKE: the ONE probe that would have caught the Task-0 K0-basis bug. Over
    # Example 2.47 (det Cartan = 0) the regular object A = P1 (+) P2 IS silting (in fact
    # tilting, AI Ex 2.2). The old _chi g-matrix gave det = det Cartan = 0 and a FALSE
    # is_silting=False.
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_silting_object(T)
    assert rep.is_silting is True and rep.det in (1, -1)
    assert rep.generation_certified_by.startswith("tilting")


@selfcert
def test_selfinjective_nakayama_regular_is_silting():
    # SMOKE: a NON-local, NON-unimodular-Cartan self-injective algebra. kZ3/J2 =
    # NakayamaAlgebra([2,2,2], cyclic): 3 simples, det Cartan = 2, self-injective. The
    # regular object is silting (tilting). Old _chi gave det = 2 -> false negative.
    B = NakayamaAlgebra([2, 2, 2], cyclic=True, field=GF(32003))
    T = [ChainComplex.stalk(B.projective(v), 0) for v in B.quiver.vertices]
    rep = is_silting_object(T)
    assert rep.is_silting is True and rep.det in (1, -1)


@lit
def test_local_algebra_silting_is_shifts():
    # k[x]/(x^2): local, indecomposable silting object A -> silt = {A[i]} (Thm 2.26).
    from quiverlab.families import truncated_polynomial
    A = truncated_polynomial(2, field=GF(32003))
    Astalk = ChainComplex.stalk(A.projective(1), 0)
    rep = is_silting_object([Astalk])
    assert rep.is_silting is True
    assert rep.generation_certified_by.startswith(("tilting", "local"))
    # a shift A[3] is also silting:
    assert is_silting_object([Astalk.shift(3)]).is_silting is True


@selfcert
def test_missing_summand_is_not_silting():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(1), 0)]      # one summand, two simples
    rep = is_silting_object(T)
    assert rep.k0_basis is False and rep.is_silting is False


@selfcert
def test_nonperfect_summand_refused():
    A = _a2()
    T = [ChainComplex.stalk(A.simple(1), 0)]          # simple: not projective
    with pytest.raises(QuiverlabError, match="perfect"):
        is_silting_object(T)


@selfcert
def test_co_t_structure_record_and_refusal():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    ct = co_t_structure_of(T)
    assert ct["bounded"] is True and len(ct["coheart"]) == 2
    assert "aihara" in " ".join(ct["references"]).lower() or ct["references"]
    # refuses on a non-silting input (does not fabricate a co-t-structure):
    with pytest.raises(QuiverlabError):
        co_t_structure_of([ChainComplex.stalk(A.projective(1), 0)])
