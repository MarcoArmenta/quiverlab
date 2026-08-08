"""The homological string-algebra test (Plan 59 / R34, Suarez-Alvarez 2023).
Discriminating oracle: a >= 3-summand middle (verdict "not_string") is k-bar-sound and,
if is_string(A) is True, RAISES (one engine is a bug). The finite-field confirm yields
"string_over_ground_field" (inconclusive w.r.t. the k-bar theorem) and defers to the
syntactic is_string; a mismatch there is a recorded k-bar-gap, never raised. Refute
runs over QQ (the pipeline needs char 0 / char > dim); the finite-field confirm runs
over small primes on small algebras (char > dim of the middles).

Perf note (verified 2026-08-07): the AR knit that decides rep-finiteness is ~3 s/module
and grows with dim, so a rep-INFINITE input (the 2-Kronecker) is refused under a small
knit budget rather than churning the default cap for minutes; the point of that test is
only that it RAISES."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_string
from quiverlab.modules.string_homological import homological_string_test

xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _kA(n, field=QQ):
    arrows = {chr(ord("a") + i): (i + 1, i + 2) for i in range(n - 1)}
    return Quiver(list(range(1, n + 1)), arrows).algebra(relations=[], field=field)


def _kD4(field=QQ):
    # subspace orientation: three arrows into the centre 0. Rep-finite (Dynkin D4),
    # NOT special-biserial (out/in-degree 3), so NOT a string algebra.
    Q = Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)})
    return Q.algebra(relations=[], field=field)


def _gentle_a3(field=QQ):
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=field)


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_linear_A_over_QQ_never_refutes_string(n):
    # kA_n: every extension of indecomposables has <= 2 middle summands (string). Over
    # QQ (infinite field) the confirm side cannot exhaust, so the honest verdict is
    # "inconclusive" -- and NEVER "not_string" (no >= 3 witness exists).
    v = homological_string_test(_kA(n))
    assert v.is_string is True
    assert v.verdict == "inconclusive"                 # infinite field: no k-bar confirm
    assert v.witness is None


@lit
def test_kD4_is_not_string_with_the_verified_three_middle_witness():
    # THE discriminating witness (VERIFIED LIVE 2026-08-07, dev, QUIVERLAB_NO_NUMBA=1):
    # knit_ar_quiver(kD4) is complete with 12 indecomposables; the AR sequence ending at
    # the module M with dim vector {0:2,1:1,2:1,3:1} has middle E = P1 (+) P2 (+) P3, the
    # three arm-projectives {0:1,1:1} / {0:1,2:1} / {0:1,3:1} -- EXACTLY 3 summands.
    v = homological_string_test(_kD4())
    assert v.is_string is False                        # kD4 is not special-biserial
    assert v.verdict == "not_string"
    assert v.witness is not None and v.witness["summand_count"] == 3
    assert v.witness["E_dimvec"] == {0: 3, 1: 1, 2: 1, 3: 1}
    # this verdict is k-bar-sound and is_string is False, so the function RETURNS (no
    # raise); a raise would occur only if is_string had (wrongly) been True.


@xeng
@pytest.mark.parametrize("factory", [_kA, _gentle_a3])
def test_returned_verdict_is_consistent_with_is_string(factory):
    A = factory(3) if factory is _kA else factory()
    v = homological_string_test(A)
    # a returned "not_string" is only produced when is_string is False (else it raises);
    # "string_over_ground_field"/"inconclusive" defer to is_string (no assertion of equality).
    if v.verdict == "not_string":
        assert v.is_string is False


@xeng
def test_finite_field_exhaustive_confirm_on_kA3():
    # small prime, small middles (dim E <= 4 < char 5 so decompose certifies): the confirm
    # side closes -> "string_over_ground_field" (NOT a definitive k-bar "string"; that is
    # is_string's job). No k-bar gap here since is_string is True.
    v = homological_string_test(_kA(3, field=GF(5)))
    assert v.verdict == "string_over_ground_field"
    assert v.is_string is True and v.kbar_gap_note is None


@selfcert
def test_witness_middle_self_certifies():
    v = homological_string_test(_kD4())
    # the recorded witness is a genuine >= 3 decomposition of a genuine extension
    assert sum(1 for _ in v.witness["summand_dimvecs"]) == v.witness["summand_count"]
    assert v.witness["class"] in ("ar_socle", "basis")


@selfcert
def test_bounded_split_count_char_fallback_returns_three_on_kD4_middle():
    # The char-blocked refute fallback (_bounded_split_count) is a bounded Fitting-split
    # LOWER BOUND that never raises. On the kD4 AR middle P1 (+) P2 (+) P3 it returns 3 --
    # the same >= 3 witness the char-robust terminal split gives (so a char <= dim input
    # whose decompose refuses still refutes when the split reaches 3).
    from quiverlab.modules.ar import almost_split_sequence, knit_ar_quiver
    from quiverlab.modules.string_homological import _bounded_split_count
    A = _kD4()
    ar = knit_ar_quiver(A, budget_modules=64, budget_dim=32)
    M = next(v["module"] for v in ar.vertices
             if v["module"].dimension_vector() == {0: 2, 1: 1, 2: 1, 3: 1})
    E = almost_split_sequence(M).M                  # dim vector {0:3,1:1,2:1,3:1}
    cnt, dvs = _bounded_split_count(E)
    assert cnt == 3 and len(dvs) == 3              # the lower-bound fallback finds all 3


@selfcert
def test_rep_infinite_and_selfinjective_and_presentationless_refused():
    from quiverlab.core.algebra import Algebra
    from quiverlab.families import NakayamaAlgebra
    kron = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    with pytest.raises(QuiverlabError):
        # 2-Kronecker: rep-infinite (knit never completes). A small knit budget refuses
        # fast -- the rep-infinite knit is ~3 s/module and would otherwise churn the cap.
        homological_string_test(kron, knit_budget_modules=6, knit_budget_dim=12)
    sinj = NakayamaAlgebra(n=3, l=4, cyclic=True, field=QQ)   # self-injective
    with pytest.raises(QuiverlabError):
        homological_string_test(sinj)                  # knitter unsupported
    one = QQ.one()
    sc = Algebra(QQ, [[[one]]], [one])                 # presentation-less (quiver is None)
    with pytest.raises(QuiverlabError):
        homological_string_test(sc)
