"""Plan 74 Task 5: free-action detection + the honest covering slice.

The flagship (Z/2 on dual numbers) is NOT free -- sigma fixes the single vertex;
its orbit-Nakayama A|xG arises from the kG-idempotent splitting k[Z/2] = k x k,
NOT from a free vertex orbit. A genuine FREE action swaps two vertices."""
import pathlib

import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.families.basic import truncated_polynomial
from quiverlab.families.skew_group import (
    GroupAction, QuiverAutomorphism, is_free_action)
from quiverlab.fields import QQ

pytestmark = pytest.mark.oracle_selfcert

_BACKLOG = (pathlib.Path(__file__).resolve().parents[2]
            / "docs" / "plans" / "DEEPER-ENGINES-BACKLOG.md")


def test_flagship_not_free_kg_splitting():
    A = truncated_polynomial(2, field=QQ)
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": -1})
    ga = GroupAction.cyclic(sig, 2)
    # sigma fixes the single vertex -> NOT a free action (the splitting is from kG)
    assert is_free_action(A, ga) is False


def test_vertex_swap_is_free():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2), "b": (2, 1)})
    A = Q.algebra(relations=["a*b", "b*a"], field=QQ)
    swap = QuiverAutomorphism({1: 2, 2: 1}, {"a": "b", "b": "a"})
    ga = GroupAction.cyclic(swap, 2)
    assert is_free_action(A, ga) is True


def test_covering_transport_is_v1_deferred():
    # The general Galois-covering HH transport is NOT claimed in v1 -- only the
    # detection + the free-orbit oracle. The deferral is honestly ledgered.
    text = _BACKLOG.read_text(encoding="utf-8")
    assert "P74 skew-group covering-reduction HH transport" in text
