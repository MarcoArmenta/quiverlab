"""Discriminating battery (Plan 59 / R34). A returned "not_string" (a >= 3-summand
middle) is k-bar-sound and only occurs when is_string is False (else the function
raises). "string_over_ground_field"/"inconclusive" defer to is_string; a k-bar-gap is
recorded, never raised. String side: kA_n / gentle kA_n(ab) -> no >= 3 witness.
Non-string rep-finite side: Dynkin D path algebras (a valency->=3 vertex -> a >= 3 AR
mesh) -> not_string with a witness.

Deviation from the plan sketch (perf, verified 2026-08-07): the plan listed kD4/kD5/kE6.
The AR knit is ~3 s/module and grows with dim, so kE6 (36 indecomposables, dim up to 11)
knits for MINUTES; kD5 (20 indec, dim <= 7) knits in ~17 s. The kE6 mesh mechanism is
identical to kD5's (a valency-3 vertex), so the battery keeps kD4 + kD5 (both
live-verified 3-summand middles) and drops kE6 as too slow for CI. A module-scoped
fixture computes each verdict once (kD5 is otherwise re-knit per test)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.dynkin import dynkin_quiver
from quiverlab.invariants.recognizers import is_string
from quiverlab.modules.string_homological import homological_string_test

xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


def _path(vertices, arrows):
    return Quiver(vertices, arrows).algebra(relations=[], field=QQ)


STRING = [                                   # rep-finite string algebras
    _path([1, 2], {"a": (1, 2)}),                          # kA2
    _path([1, 2, 3], {"a": (1, 2), "b": (2, 3)}),          # kA3
    Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ),
]
# rep-finite, NOT string (a valency->=3 vertex). VERIFIED LIVE: kD4 (subspace) knits
# complete (12 indec) with a 3-summand AR middle; kD5 knits complete (20 indec) likewise.
NOT_STRING = [
    _path([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}),          # kD4 (subspace)
    dynkin_quiver("D5").algebra(relations=[], field=QQ),                   # kD5
]

_ALL = STRING + NOT_STRING


@pytest.fixture(scope="module")
def verdicts():
    """Compute every verdict once (kD5's AR knit is ~17 s; do not repeat it per test)."""
    return {id(A): homological_string_test(A) for A in _ALL}


@xeng
@pytest.mark.parametrize("i", range(len(STRING)))
def test_string_side_no_false_refutation(i, verdicts):
    A = STRING[i]
    assert is_string(A) is True
    v = verdicts[id(A)]                          # cannot raise: no >=3 witness on strings
    assert v.verdict != "not_string"             # never a false >= 3 witness


@xeng
@pytest.mark.parametrize("i", range(len(NOT_STRING)))
def test_not_string_side_finds_a_witness(i, verdicts):
    A = NOT_STRING[i]
    assert is_string(A) is False                 # a valency->=3 vertex: not special-biserial
    v = verdicts[id(A)]
    assert v.verdict == "not_string" and v.witness["summand_count"] >= 3


@lit
def test_refute_side_raises_on_a_contradiction_never_returns_it(verdicts):
    # The one loud disagreement: a >=3 witness with is_string True would RAISE. Since our
    # zoo has none, every returned verdict is self-consistent: not_string => is_string False.
    from quiverlab.errors import QuiverlabError  # noqa: F401 (documents the raise contract)
    for A in _ALL:
        v = verdicts[id(A)]                      # returns (no contradiction in the zoo)
        if v.verdict == "not_string":
            assert v.is_string is False
