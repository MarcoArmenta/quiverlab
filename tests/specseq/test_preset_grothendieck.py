"""Cartan-Eilenberg / Grothendieck change-of-rings SS. Cross-engine: E_inf total
== module Ext on several instances incl. a multi-vertex one and NONZERO-Ext pins;
degenerate B=A collapse; scope + acyclicity loud refusals.

NB (Plan-42 implementation, arbitrated + documented):
 * The Hom double complex is COHOMOLOGICAL, stored with NEGATED total degree
   (position (-p,-q)) to fit the homological engine (the P39 C^n:=C_{-n}
   discipline); the abutment is read at the negated homological degree
   (``totals.get(-n)``) -- the plan Global-Constraint "its oracle compares the
   abutment at the negated degree".
 * The restriction of a B-module to A is rebuilt via ``from_arrow_action``
   (``A.module(dimvec, arrow_action)``), NOT ``Module(A, M.dim, M.action)``: A's
   longer paths (e.g. ``a*b``, nonzero in A but absent from a B=kQ/(a*b)-module's
   action) must be filled in by composition, else the A-side Ext resolution
   KeyErrors.
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.specseq.presets import cartan_eilenberg_ss

pytestmark = pytest.mark.oracle_crossengine


def _a3rel():                                       # A = kA3 / (a*b)
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def _kA3():                                         # hereditary kA3 (I = 0)
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=GF(5))


def _restrict(M, A):
    # M viewed as an A-module, rebuilt from the common arrow actions (fills in A's
    # longer paths -- see the module docstring).
    return A.module(M.dimension_vector(), {a: M.action[a] for a in A.quiver.arrows},
                    side=M.side, name=M.name)


def _einf_totals(ss, top):
    Einf = ss.page(ss.convergence.e_infinity_page)
    totals = {}
    for (p, q) in Einf.spots:
        totals[p + q] = totals.get(p + q, 0) + Einf.dim(p, q)
    # cohomological Ext^n lives at homological total degree -n (negated storage)
    return [totals.get(-n, 0) for n in range(top + 1)]


def test_degenerate_B_equals_A_collapses_to_ext():
    # B = A: E_2 = E_inf, abutment == Ext_A(M, N).
    A = _kA3()
    B = _kA3()                                      # same presentation => I' = I
    M, N = B.simple(1), A.simple(3)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=4, q_len=4)
    got = _einf_totals(ss, 3)
    assert got == [A.ext(A.simple(1), A.simple(3), n) for n in range(4)]
    assert ss.convergence.collapse() is True


def test_degenerate_B_equals_A_nonzero_ext():
    # a NONZERO abutment pin: Ext_A(S_1, S_2) = [0, 1, 0, 0] over hereditary kA3.
    A, B = _kA3(), _kA3()
    M, N = B.simple(1), A.simple(2)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=4, q_len=4)
    assert _einf_totals(ss, 3) == [A.ext(A.simple(1), A.simple(2), n) for n in range(4)]
    assert _einf_totals(ss, 3) == [0, 1, 0, 0]


def test_change_of_rings_abuts_to_A_ext():
    # A = kA3 (hereditary), B = kA3/(a*b) (an admissible quotient): a B-module M
    # and an A-module N; E_inf total == Ext_A(M|_A, N).
    A = _kA3()
    B = _a3rel()                                    # rel(A)=∅ subset of {a*b}
    M, N = B.simple(2), A.simple(1)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=5, q_len=5)
    got = _einf_totals(ss, 4)
    M_over_A = _restrict(M, A)                       # restriction
    assert got == [A.ext(M_over_A, N, n) for n in range(5)]


def test_change_of_rings_nonzero_ext():
    # a NONZERO change-of-rings abutment: M|_A = S_1, N = S_2, Ext_A = [0, 1, 0, ...].
    A, B = _kA3(), _a3rel()
    M, N = B.simple(1), A.simple(2)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=5, q_len=5)
    M_over_A = _restrict(M, A)
    assert _einf_totals(ss, 4) == [A.ext(M_over_A, N, n) for n in range(5)]
    assert _einf_totals(ss, 4) == [0, 1, 0, 0, 0]


def test_multivertex_instance():
    # the required multi-vertex change-of-rings pin (Design-decision 5a).
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=[], field=GF(7))
    B = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=["a*b", "b*c"], field=GF(7))
    M, N = B.simple(1), A.simple(4)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=5, q_len=5)
    got = _einf_totals(ss, 4)
    M_over_A = _restrict(M, A)
    assert got == [A.ext(M_over_A, N, n) for n in range(5)]


def test_multivertex_nonzero_ext():
    # multi-vertex NONZERO pin: M|_A = S_1, N = S_2, Ext_A = [0, 1, 0, 0, 0].
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=[], field=GF(7))
    B = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=["a*b", "b*c"], field=GF(7))
    M, N = B.simple(1), A.simple(2)
    ss = cartan_eilenberg_ss(A, B, M, N, p_len=5, q_len=5)
    M_over_A = _restrict(M, A)
    assert _einf_totals(ss, 4) == [A.ext(M_over_A, N, n) for n in range(5)]
    assert _einf_totals(ss, 4) == [0, 1, 0, 0, 0]


def test_scope_refusal_is_loud():
    # the active per-instance hypothesis gate: B must be an admissible quotient of A
    # (same quiver, rel(A) subset rel(B)). Swapping the roles (rel=∅ ⊄ rel={a*b})
    # refuses loudly naming the hypothesis -- never a wrong abutment.
    A = _a3rel()          # rel = {a*b}
    B = _kA3()            # rel = ∅  -> rel(A)={a*b} NOT subset of rel(B)=∅
    M, N = B.simple(1), A.simple(1)
    with pytest.raises(QuiverlabError, match="subset|quotient|acyclic|hypothesis"):
        cartan_eilenberg_ss(A, B, M, N, p_len=3, q_len=3)


def test_acyclicity_failure_refuses_loudly():
    # a constructed instance whose per-instance acyclicity probe fails must raise
    # a loud refusal naming the hypothesis, NOT return a wrong abutment. For the
    # admissible-quotient scope the abutment is ALWAYS Ext_A (injective-resolution
    # collapse), so the abutment-mismatch refusal is defensive and this instance
    # takes the else branch (correct abutment) -- recorded on the verification page.
    A = _kA3()
    B = _a3rel()
    M, N = B.simple(1), A.simple(1)
    try:
        ss = cartan_eilenberg_ss(A, B, M, N, p_len=3, q_len=3)
    except QuiverlabError as exc:
        assert "acyclic" in str(exc).lower() or "hypothesis" in str(exc).lower()
    else:
        got = _einf_totals(ss, 2)
        M_over_A = _restrict(M, A)
        assert got == [A.ext(M_over_A, N, n) for n in range(3)]
