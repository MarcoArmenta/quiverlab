"""ACLV Theorem B: ada over an algebraically closed field => (simply connected <=> HH^1=0),
and HH.(A) reduces to k. Plan 61's headline oracle. Over CC the verdict is emitted; over
QQ/GF(p)/QQi only dim HH^1 + the theorem note (Theorem B's alg-closed hypothesis unmet).
The gate is A.domain.is_algebraically_closed -- there is NO Algebra.field, and QQi shares
CC's SympyExactDomain class, so the flag is the ONLY sound predicate."""
import pytest

from quiverlab import GF, Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ, CC, QQi
from quiverlab.modules.recognizers_ladder import recognizer_ladder

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _radsq_nakayama_a5(field):
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=field)


@selfcert
def test_alg_closed_gate_predicate_on_constructed_algebras():
    # H1 MANDATE: assert the gate on A.domain (the constructed algebra), NEVER the CC sentinel.
    # CC declares C => True; QQ/GF(p)/QQi are not algebraically closed => False (QQi shares CC's
    # SympyExactDomain class, so the flag -- not the type -- must distinguish them).
    assert linear_path_algebra(2, field=CC).domain.is_algebraically_closed is True
    for f in (QQ, GF(7), QQi):
        assert linear_path_algebra(2, field=f).domain.is_algebraically_closed is False


@lit
def test_fixture_simply_connected_over_CC():
    # rad^2=0 A5: ada True, tree quiver => simply connected => HH^1 = 0 (Theorem B).
    L = recognizer_ladder(_radsq_nakayama_a5(CC))
    sc = L.simple_connectedness
    assert sc.applicable is True and sc.field_algebraically_closed is True
    assert sc.hh1_dim == 0 and sc.verdict is True          # simply connected


@lit
def test_theorem_B_verdict_false_on_real_ada_hh1_nonzero():
    # A REAL rep-finite ada algebra over CC with HH^1 != 0 (verified live) -- the genuine
    # negative pin, no longer only the synthetic _theorem_b_verdict branch. The "square"
    # kQ/(a*b, c*d) on 1->2->4, 1->3->4 with BOTH length-2 routes killed: quasi-tilted
    # (gl.dim 2, empty complement) hence ada, dim HH^1 = 1, pi_1 = Z (NOT simply connected),
    # so ACLV Theorem B returns verdict False. Cross-agrees with P56 (see
    # tests/modules/test_ada_hh1_p56_agreement.py, where this same algebra is the False member).
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}).algebra(
        relations=["a*b", "c*d"], field=CC)
    L = recognizer_ladder(A)
    assert L.verdict("ada") is True and L.verdict("quasi_tilted") is True
    sc = L.simple_connectedness
    assert sc.applicable is True and sc.field_algebraically_closed is True
    assert sc.hh1_dim == 1 and sc.verdict is False          # HH^1 != 0 => NOT simply connected


@selfcert
def test_theorem_B_verdict_false_when_hh1_nonzero():
    # Gate-logic coverage of the pure _theorem_b_verdict helper across all four
    # (ada x alg_closed) corners; the end-to-end REAL ada + HH^1 != 0 instance is pinned in
    # test_theorem_B_verdict_false_on_real_ada_hh1_nonzero above.
    from quiverlab.modules.recognizers_ladder import _theorem_b_verdict
    assert _theorem_b_verdict(ada=True, alg_closed=True, hh1_dim=1) is False
    assert _theorem_b_verdict(ada=True, alg_closed=True, hh1_dim=0) is True
    assert _theorem_b_verdict(ada=True, alg_closed=False, hh1_dim=0) is None   # gated off
    assert _theorem_b_verdict(ada=False, alg_closed=True, hh1_dim=0) is None   # not ada


@selfcert
def test_verdict_gated_off_non_algebraically_closed_field():
    # over QQ (not alg-closed): applicable False, hh1_dim reported, verdict None + theorem note.
    L = recognizer_ladder(_radsq_nakayama_a5(QQ))
    sc = L.simple_connectedness
    assert sc.field_algebraically_closed is False and sc.applicable is False
    assert sc.hh1_dim == 0 and sc.verdict is None and "HH" in sc.theorem


@lit
def test_not_ada_no_verdict():
    # A REAL in-scope non-ada witness (verified live): the triangle 1->2->3 + 1->3 with the
    # long route killed (I = <a*b>). Its projective P_1 lies in NEITHER part (P_1 has id >= 2
    # so P_1 not in R_A, and a bad-pd predecessor keeps it out of L_A), so ada = False; hence
    # Theorem B does not apply and the SC block is not applicable.
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3)}).algebra(
        relations=["a*b"], field=CC)
    L = recognizer_ladder(A)
    assert L.is_complete and L.verdict("ada") is False
    assert L.simple_connectedness.applicable is False and L.simple_connectedness.verdict is None
