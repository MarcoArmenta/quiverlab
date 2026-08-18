"""Literature + cross-engine oracles for the tilted recognizer (Plan 60 / R17). Literature:
hereditary Dynkin => tilted with the projective slice + correct type (kA_n, kD_4); kA3/rad^2 =>
tilted type A3; kZ3/J^2 (= cluster-tilted A3) => not tilted; Liu's rad^2=0 sincere-not-faithful
example => not tilted (faithfulness is essential -- arXiv:1409.2054 Ex. after Thm 2.6). Cross-
engine: the gl.dim gate and the direct section search agree on rad^2=0 A5; the reconstructed
Gabriel type equals the section-graph type. QQ-scope; GF(32003) parity for kA_n only.

CORRECTION (worker, live-verified from the arXiv:1409.2054 PDF, p.9): Liu's counterexample after
Thm 2.6 is the rad^2=0 3-CYCLE (its AR quiver is a TUBE with exactly 6 indecomposables =
{P_a,P_b,P_c,S_a,S_b,S_c}, all distinct -- which FORCES every vertex to have an outgoing arrow,
i.e. an oriented cycle). That algebra IS kZ3/J^2, self-injective, so quiverlab refutes it via
Gate S (not the faithfulness prune -- the whole algebra is off the rep-finite knit). This unifies
with the record's settled correction that cluster-tilted A3 = kZ3/J^2. The plan's premise that
Liu's example is a non-self-injective acyclic quiver was a misreading of the figure's arrowheads.
Liu's PEDAGOGICAL point -- faithfulness cannot be weakened to sincereness -- is still exercised
directly on Liu's actual cut Delta = {P_b, S_b, P_a} of kZ3/J^2 below (sincere + Hom(X,tau Y)=0
but NOT faithful)."""
import pytest

from quiverlab import GF, Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.families.dynkin import dynkin_quiver
from quiverlab.modules.tilted import tilted_check, _faithful_section_search, _is_faithful

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _kZ3_J2():
    # the rad^2=0 3-cycle a->b->c->a (a=1, b=2, c=3): Liu's example AND kZ3/J^2 AND
    # cluster-tilted A3. Self-injective; 6 indecomposables (3 projective, 3 simple).
    return RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}), field=QQ)


@lit
def test_kD4_hereditary_tilted_type_D4():
    A = dynkin_quiver("D4").algebra(field=QQ)          # kD_4 (central vertex 2 + 3 leaves), acyclic
    rep = tilted_check(A)
    assert rep.verdict == "tilted" and rep.reason == "hereditary"
    assert rep.hereditary_type == "D_4"


@lit
def test_kA3_parity_gf32003_tilted():
    # kA_n has distinct dim-vectors per indec -> is_isomorphic never enters the positive-only
    # branch on a non-match -> GF(32003) is safe (M1). Verdict byte-identical to QQ.
    A = linear_path_algebra(3, field=GF(32003))
    assert tilted_check(A).verdict == "tilted"


@lit
def test_liu_counterexample_is_kZ3J2_selfinjective_not_tilted():
    # Liu's Ex. after Thm 2.6 = the rad^2=0 3-cycle = kZ3/J^2: self-injective => Gate S => not
    # tilted. (See the module docstring: the AR quiver is a 6-indec tube.)
    A = _kZ3_J2()
    assert A.is_selfinjective() is True
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted" and rep.reason == "self_injective"


@lit
def test_liu_cut_is_sincere_but_not_faithful():
    # Liu's exact cut Delta = P_b - S_b - P_a (b=2, a=1) in kZ3/J^2: SINCERE (covers every simple)
    # and Hom_A(X, tau Y) = 0, but NOT faithful -- so it is rejected. This is Liu's point that the
    # faithfulness of the cut cannot be weakened to sincereness (arXiv:1409.2054, Ex. after 2.6).
    from quiverlab.modules.hom import hom_dim
    from quiverlab.modules.morphism import direct_sum
    A = _kZ3_J2()
    Pa, Pb, Sb = A.projective(1), A.projective(2), A.simple(2)
    cut = [Pb, Sb, Pa]
    S, _, _ = direct_sum(*cut)
    # sincere: every simple is a composition factor of some cut module
    dv = S.dimension_vector()
    assert all(dv.get(v, 0) > 0 for v in A.quiver.vertices)          # SINCERE
    taus = [Y.tau() for Y in cut]
    assert all(hom_dim(X, tY) == 0 for X in cut for tY in taus)      # Hom(X, tau Y) = 0
    assert _is_faithful(A, S) is False                              # but NOT faithful => rejected


@xeng
def test_gate_and_search_agree_on_radsq_a5():
    A = RadicalSquareZero(Quiver([1, 2, 3, 4, 5],
        {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)}), field=QQ)
    # Gate G verdict:
    assert tilted_check(A).verdict == "not_tilted"
    # direct search verdict (bypassing gates):
    ar = A.ar_quiver(budget_modules=128)
    assert _faithful_section_search(ar, A, budget_sections=4096) is None
