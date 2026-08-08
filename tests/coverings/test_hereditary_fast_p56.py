"""Hereditary short-circuit (adjudicated Fix 1): a relation-free (hereditary) algebra
has NO minimal relations, so fundamental_group / is_simply_connected must skip
loewy_length + the block linear algebra entirely and return pi1 = free on the Betti
generators immediately. Without the short-circuit, A_40 profiled at ~7.5s (and A_16
at ~7.25s) purely in loewy_length + path enumeration. Unmarked: a performance /
contract test, not an oracle (Plan-32 ruling)."""
import time

from quiverlab import GF, linear_path_algebra
from quiverlab.invariants.coverings import fundamental_group, is_simply_connected


def test_hereditary_path_is_fast_and_free():
    A = linear_path_algebra(40, field=GF(7))       # kA_40: a tree, relations == []
    assert A.relations == []
    t0 = time.perf_counter()
    g = fundamental_group(A)
    r = is_simply_connected(A)
    elapsed = time.perf_counter() - t0
    assert g.free_rank == 0 and g.invariant_factors == () and g.components == 1
    assert r.verdict is True                       # a tree is simply connected (R1)
    assert elapsed < 2.0, f"hereditary short-circuit regressed: {elapsed:.2f}s"
