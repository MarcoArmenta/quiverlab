"""Tilted-algebra recognizer via the Liu-Skowronski faithful-section criterion (Plan 60 / R17,
Liu Arch. Math. 61 (1993); arXiv:1409.2054 Thm 2.6; Ringel LNM 1099 (4.2) = Liu Thm 1.9(2)).
Verdict pipeline: Gate H (hereditary => tilted), Gate S (non-semisimple self-injective => not),
Gate G (gl.dim > 2 => not), then the exhaustive faithful-section search on the knitted AR quiver.
QQ-scope (the tilting count + Gabriel recovery are char 0 / char > dim)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.tilted import tilted_check

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _radsq_a5():                                   # ACLV 2.2(b) / P55 shared fixture, gl.dim 4
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


def _kA3_radsq():                                  # P55 A_lambda shape: kA3/rad^2, TILTED, type A3
    return RadicalSquareZero(Quiver([1, 2, 3], {"a1": (2, 1), "a2": (3, 2)}), field=QQ)


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_hereditary_is_tilted(n):
    A = linear_path_algebra(n, field=QQ)           # kA_n hereditary, rep-finite
    rep = tilted_check(A)
    assert rep.verdict == "tilted" and rep.reason == "hereditary"
    assert rep.hereditary_type == f"A_{n}"
    assert bool(rep) is True and A.is_tilted() is True
    # M3: rep-finite hereditary carries a CERTIFIED projective slice (A_A), n summands.
    assert len(rep.slice) == n
    from quiverlab.modules.tilted import _certify_slice
    ok, _H = _certify_slice(A, rep.slice_module)   # End_A(A_A) hereditary => a genuine slice
    assert ok is True


@lit
def test_kA3_radsq_is_tilted_type_A3():
    A = _kA3_radsq()
    rep = tilted_check(A)
    assert rep.verdict == "tilted" and rep.reason == "faithful_section_found"
    assert len(rep.slice) == 3                      # |Sigma| = n (consequence)
    assert {r["name"] for r in rep.slice} == {"S_2", "P_2", "P_3"}
    assert rep.hereditary_type == "A_3"             # End_A(S) = kA3 (verified live)


@lit
def test_selfinjective_kZ3_J2_not_tilted():
    # kZ3/J^2 = cluster-tilted A3 (3-cycle) = Jac(3-cycle, a*b*c): self-injective, knit REFUSES;
    # the non-semisimple-self-injective theorem gate gives a REAL not_tilted, not a refusal.
    A = RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}), field=QQ)
    assert A.is_selfinjective() is True
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted" and rep.reason == "self_injective"


@lit
def test_gldim_gt_2_not_tilted():
    A = _radsq_a5()                                 # gl.dim 4 (exact) -> Gate G
    assert A.global_dimension().value == 4
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted" and rep.reason == "gldim>2"


@selfcert
def test_truncated_poly_selfinjective_not_tilted():
    A = truncated_polynomial(3, field=QQ)           # k[x]/(x^3) self-injective non-semisimple
    rep = tilted_check(A)
    assert rep.verdict == "not_tilted" and rep.reason == "self_injective"


@selfcert
def test_search_refutes_radsq_a5_directly():
    # exercise the SEARCH path (bypassing Gate G): the knittable rad^2=0 A5 (9 indec, 5 orbits,
    # NOT tilted) has NO faithful section -- _faithful_section_search returns None.
    from quiverlab.modules.tilted import _faithful_section_search
    A = _radsq_a5()
    ar = A.ar_quiver(budget_modules=128)
    assert ar.is_complete and len(ar.tau_orbits) == 5   # 5 orbits, sizes [5,1,1,1,1] => 5 transversals
    assert _faithful_section_search(ar, A, budget_sections=4096) is None


@selfcert
def test_slice_module_is_tilting_and_size_n():
    from quiverlab.modules.tilting import is_tilting_module
    rep = tilted_check(_kA3_radsq())
    assert is_tilting_module(rep.slice_module).is_tilting is True   # Ringel 1.9(2), tilting half
    assert len(rep.slice) == rep.universe_size - 2 == 3             # |Sigma| = n = 3 of 5 indec


@selfcert
def test_structure_constant_algebra_refused():
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.endomorphism import end_algebra
    A = linear_path_algebra(2, field=QQ)
    E = end_algebra(A.projective(1))                # quiver is None
    with pytest.raises(QuiverlabError):
        tilted_check(E)


@selfcert
def test_delegate_matches_free_function():
    A = _kA3_radsq()
    assert A.tilted_check().verdict == tilted_check(A).verdict
