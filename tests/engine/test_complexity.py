"""Complexity-estimator behaviour, including the counterexample-detector regression.

cx = 1 + (polynomial growth degree of dim HH_n);  cx = 0  <=>  HH_n eventually vanishes.
A cx-0 algebra of infinite global dimension would be a counterexample to Han's
conjecture, so the homology detector MUST recognise an eventually-vanishing sequence.
"""
# hanlab sys.path shim dropped: quiverlab uses absolute package imports.
import pytest
from quiverlab.engine.scan2 import complexity_diagnostic
from quiverlab.engine.scan3 import complexity_of


# ---- known-good classifications (oracle; expected to pass) ----
@pytest.mark.parametrize("seq,expected_cx", [
    ([2, 1, 1, 1, 1], 1),        # eventually constant nonzero  -> cx 1
    ([4, 4, 5, 6, 7, 8, 9], 2),  # HH_0 then linear growth      -> cx 2
    ([3, 3, 3, 3], 1),           # constant                     -> cx 1
])
def test_complexity_diagnostic_known(seq, expected_cx):
    assert complexity_diagnostic(seq)["complexity"] == expected_cx


# ---- BUG (expected RED until fixed): eventual vanishing must be cx 0 ----
# A homology sequence decaying to zero through a non-constant tail is the Han
# counterexample signature; the detector must not mislabel it as "growing".
@pytest.mark.parametrize("seq", [
    [4, 2, 1, 0, 0],
    [3, 2, 1, 0, 0, 0],
    [5, 3, 2, 1, 0, 0],
])
def test_complexity_diagnostic_detects_eventual_vanishing(seq):
    assert complexity_diagnostic(seq)["complexity"] == 0


# ---- the two estimators must agree on the complexity value ----
@pytest.mark.parametrize("seq", [
    [2, 1, 1, 1, 1],
    [4, 4, 5, 6, 7, 8, 9],
    [4, 2, 1, 0, 0],
    [3, 3, 3, 3],
    [2, 2, 1, 0, 0],
])
def test_estimators_agree(seq):
    assert complexity_diagnostic(seq)["complexity"] == complexity_of(seq)
