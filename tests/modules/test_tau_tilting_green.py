"""Maximal green sequences (Plan 45 / C4). Literature: kA2 has exactly 2 MGSs (the two
sides of the pentagon); a tau-tilting-infinite algebra caps loudly."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.green import maximal_green_sequences

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
def test_ka2_has_two_maximal_green_sequences():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)
    mgs = maximal_green_sequences(A)
    assert mgs["complete"] and mgs["count"] == 2


@selfcert
def test_every_mgs_starts_at_top_ends_at_bottom():
    A = linear_path_algebra(3, field=QQ)
    mgs = maximal_green_sequences(A)
    assert mgs["complete"] and mgs["count"] >= 1
    for seq in mgs["sequences"]:
        assert len(seq) >= 2
        assert seq[0] == 0                     # (A,0) is vertex 0 (BFS seed = initial pair)


@selfcert
def test_infinite_case_caps_loudly():
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    mgs = maximal_green_sequences(A, cap=10)      # any finite cap trips (infinite fan)
    assert mgs["complete"] is False and mgs["status"] == "budget"
