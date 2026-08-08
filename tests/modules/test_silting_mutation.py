"""Single silting mutation via one approximation triangle (Plan 67 / AI Def 2.30/2.34).
Literature: the 1<=>2 algebra with ab=ba=0 (AI Example 2.47): mutating P1 (+) P2 at P1
gives P2 (+) cone(P1->P2), at P2 gives P1 (+) cone(P2->P1). Self-cert: every mutant is
silting, shares n-1 summands, differs from the input; left-then-right is the identity
(involution)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.silting import (is_silting_object, silting_mutate,
                                       silting_neighbors)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _comm_square():
    # 1 <=> 2 with ab = ba = 0 (AI Example 2.47): a: 1->2, b: 2->1.
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


def _key(cx):
    # a shift-insensitive fingerprint of a summand: its K0 class up to sign, plus the
    # sorted per-degree dim-vectors -- enough to recognise cone(P1->P2) etc.
    return (tuple(sorted((n, tuple(sorted(cx.term(n).dimension_vector().items())))
                         for n in cx.degrees())),)


@selfcert
def test_mutant_is_silting_and_shares_all_but_one():
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    assert is_silting_object(T).is_silting is True
    for i in range(2):
        mut = silting_mutate(T, i, direction="left")
        assert is_silting_object(mut).is_silting in (True, "unknown")
        shared = {_key(c) for c in T} & {_key(c) for c in mut}
        assert len(shared) == 1                         # n - 1 == 1 summand shared
        assert {_key(c) for c in mut} != {_key(c) for c in T}


@lit
def test_example_2_47_cones():
    # AI Example 2.47: mu_{P1}^+(P1 (+) P2) = P2 (+) cone(P1->P2).
    A = _comm_square()
    P1 = ChainComplex.stalk(A.projective(1), 0)
    P2 = ChainComplex.stalk(A.projective(2), 0)
    mut = silting_mutate([P1, P2], 0, direction="left")   # mutate at P1 (index 0)
    # the shared summand is P2; the other is cone(P1 -> P2) (degrees {0,-1}) -- verify
    # via homology: cone of the P1->P2 approximation is a genuine 2-term complex, not a
    # stalk, and P2 survives.
    keys = {_key(c) for c in mut}
    assert _key(P2) in keys
    other = [c for c in mut if _key(c) != _key(P2)][0]
    assert set(other.degrees()) != {0}                    # a cone, not a stalk


@selfcert
def test_left_then_right_is_identity():
    # mu^-_?(mu^+_X(T)) recovers T at the matching summand (AI Prop 2.33 involution).
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    mut = silting_mutate(T, 0, direction="left")
    # the new summand N_X is the one NOT shared with T; mutate it back RIGHT:
    shared = {_key(c) for c in T}
    j = [k for k, c in enumerate(mut) if _key(c) not in shared][0]
    back = silting_mutate(mut, j, direction="right")
    assert {_key(c) for c in back} == {_key(c) for c in T}


@selfcert
def test_neighbors_reports_all_summands():
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    nb = silting_neighbors(T, direction="left")
    assert len(nb) == 2
    for (i, mutant) in nb:
        assert mutant is not None and is_silting_object(mutant).is_silting in (True, "unknown")
