"""The Tits-form tame/wild verdict, gated on the P56 strong-simple-connectivity
certificate over characteristic 0. rep-finite <=> weakly positive (Bongartz, simply
connected); tame <=> weakly nonnegative (BdlPS, strongly simply connected). Off scope:
verdict None, form still computed. Refs: Bongartz Math. Ann. 269 (1984); BdlPS Adv. Math.
226 (2011); Kasjan-Skowronski arXiv:1905.06028."""
import pytest

from quiverlab import CC, GF, QQi, Quiver, linear_path_algebra
from quiverlab.fields.rationals import RationalField
from quiverlab.invariants.tits import tame_wild_certificate

QQ = RationalField()

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _star(arms):
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt)
            ars[f"e{ai}_{nxt}"] = (prev, nxt)
            prev = nxt
            nxt += 1
    return Quiver(verts, ars).algebra(field=CC)


@lit
def test_dynkin_tree_rep_finite():
    c = tame_wild_certificate(linear_path_algebra(5, field=CC))   # A5 tree: ssc
    assert c.strongly_simply_connected is True
    assert c.weakly_positive is True and c.rep_type == "rep-finite"


@lit
def test_euclidean_tree_tame():
    c = tame_wild_certificate(_star([1, 2, 5]))                   # ~E8 tree: ssc, tame
    assert c.strongly_simply_connected is True
    assert c.weakly_nonnegative is True and c.weakly_positive is False
    assert c.rep_type == "tame" and c.witness is not None         # isotropic direction


@lit
def test_T237_tree_wild():
    c = tame_wild_certificate(_star([1, 2, 6]))                   # T_{2,3,7}: ssc, wild
    assert c.strongly_simply_connected is True
    assert c.weakly_nonnegative is False and c.rep_type == "wild"
    assert c.witness_value < 0                                    # the wild certificate


@selfcert
def test_kronecker_form_layer_only_not_simply_connected():
    # THE honest layer split: the m-Kronecker is NOT simply connected (parallel arrows,
    # pi1^ab = Z for m>=2), so the VERDICT is refused -- the m-Kronecker exercises the
    # FORM layer, not the algebra verdict. The form booleans are still computed.
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=CC)   # 2-Kronecker
    c = tame_wild_certificate(K)
    assert c.weakly_nonnegative is True and c.weakly_positive is False  # form: ~A1
    assert c.simply_connected is False                                 # not simply connected
    assert c.rep_type is None and "simply connected" in c.reason


@selfcert
def test_non_CC_refuses_verdict_keeps_form():
    A = linear_path_algebra(5, field=GF(7))                      # not algebraically closed
    c = tame_wild_certificate(A)
    assert c.weakly_positive is True                             # form computed (field-free)
    assert c.field_alg_closed is False and c.rep_type is None
    assert "algebraically closed" in c.reason or "characteristic" in c.reason


@selfcert
def test_P56_None_propagates():
    # a starved convex_budget forces P56's strong verdict to None; the rep-finite axis
    # (Bongartz, simple connectivity) still decides via the tree route, and the tame/wild
    # split is withheld -- never a fabricated tame/wild.
    A = linear_path_algebra(5, field=CC)
    c = tame_wild_certificate(A, convex_budget=0)
    assert c.strongly_simply_connected is None                  # P56 budget => None
    assert c.rep_type in (None, "rep-finite")                   # Bongartz on the tree
    if c.rep_type is None:
        assert "budget" in c.reason.lower() or "simply connected" in c.reason


@selfcert
def test_wild_witness_is_exact_negative():
    c = tame_wild_certificate(_star([1, 2, 6]))
    A = _star([1, 2, 6])
    from quiverlab.invariants.tits import tits_form_combinatorial as titsc
    assert titsc(A, list(c.witness)) == c.witness_value < 0


@selfcert
def test_bongartz_rep_infinite_middle_is_None_not_tame_wild():
    # DECIDE-AND-PIN (Plan 62 Task 4): the Zito example is simply connected but NOT
    # strongly (P56); its form is weakly nonnegative-not-positive. Bongartz (simple
    # connectivity) certifies representation-INFINITE, but the tame/wild split needs
    # STRONG simple connectivity -- which is not certified -- so the verdict is WITHHELD
    # (rep_type None), NEVER "tame" even though the form is weakly nonnegative. The
    # verdict never outruns its hypotheses.
    Z = Quiver([1, 2, 3, 4, 5],
               {"al": (1, 2), "be": (2, 3), "ga": (3, 5), "de": (2, 4), "ep": (4, 5)}
               ).algebra(relations=["al*be*ga - al*de*ep"], field=CC)
    c = tame_wild_certificate(Z)
    assert c.simply_connected is True and c.strongly_simply_connected is False
    assert c.weakly_positive is False and c.weakly_nonnegative is True
    assert c.rep_type is None                       # withheld -- never a guessed tame
    assert "strong" in c.reason.lower() and "infinite" in c.reason.lower()
    # ...but the rep-infinite certainty is SURFACED explicitly (ruling 4) -- Bongartz
    # (simple connectivity + not weakly positive) certifies rep-infinite; only the
    # tame/wild split is undecided, and that is what rep_type=None records.
    assert c.certified == "rep_infinite"


@selfcert
def test_reason_never_fabricates_algebraically_closed_over_non_CC_char0():
    # RULING 1(a): over a char-0 field that is NOT algebraically closed (QQ, QQ(i)),
    # the verdict is read by BASE CHANGE to the algebraic closure -- the reason string
    # must NEVER call the field itself "over an algebraically closed field" (the old
    # fabricated phrase), nor the scope note "algebraically-closed-only".
    for field in (QQ, QQi):
        c = tame_wild_certificate(linear_path_algebra(5, field=field))  # A5: ssc + wp
        assert "over an algebraically closed field" not in c.reason
        assert "algebraically-closed-only" not in c.scope_note
        # the honest base-change wording is present wherever a char-0 verdict is read
        if c.rep_type is not None:                          # base (pre-P61): char-0 fallback
            assert "base change" in c.reason.lower()
        else:                                               # post-P61: honest refusal
            assert "algebraically closed base field" in c.reason


@selfcert
def test_verdict_field_gate_is_two_level():
    # RULING 1(b): the verdict FIELD gate defers to P61's is_algebraically_closed flag
    # when the domain carries it, else falls back to characteristic 0. This pin must hold
    # on BOTH the pre-P61 base (fallback) and post-merge dev (flag) -- it asserts the
    # verdict is reached IFF the domain-level gate says so, never a hardcoded per-field
    # expectation. A5 is strongly simply connected + weakly positive, so the ONLY thing
    # that can withhold its verdict is the field gate.
    for field in (CC, QQ, QQi):                             # all characteristic 0
        A = linear_path_algebra(5, field=field)
        gate = getattr(A.domain, "is_algebraically_closed", A.domain.characteristic == 0)
        c = tame_wild_certificate(A)
        assert c.field_alg_closed is bool(gate)
        if gate:
            assert c.rep_type == "rep-finite"              # verdict reached
        else:
            assert c.rep_type is None and c.certified is None   # honest refusal
    # GF(7): the gate is False in BOTH regimes -> refusal, form still computed.
    A = linear_path_algebra(5, field=GF(7))
    gate = getattr(A.domain, "is_algebraically_closed", A.domain.characteristic == 0)
    c = tame_wild_certificate(A)
    assert gate is False and c.field_alg_closed is False and c.rep_type is None
    assert c.weakly_positive is True                        # field-free form still computed
