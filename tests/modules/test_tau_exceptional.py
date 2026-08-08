"""tau-exceptional sequences via the Jasso tau-perpendicular reduction (Plan 65 / R27;
Buan-Marsh J. Algebra 585 (2021); Jasso IMRN 2015; DIJ 2019). Enumeration via the ordered
support-tau-tilt bijection (count = n! * #sTt); NO mutation-BFS (transitivity only rank 2,
Buan-Hanson-Marsh 2402.10301). tau-tilting-finite + char-0/char>dim scope, QQ.

Also pins Plan 65 Task 0 (M1): the P45 `mutate` D_4-star defect fix -- `exchange_graph`
must return status="complete" with 50 vertices on ALL D_4-star orientations (the mixed
star tripped a spurious status="error" while still discovering all 50 vertices)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


# --------------------------------------------------------------------------- #
# Task 0 (M1): the P45 mutate D_4-star fix -- regression pins
# --------------------------------------------------------------------------- #
_D4_STARS = {
    "subspace": {"a": (2, 1), "b": (3, 1), "c": (4, 1)},
    "source": {"a": (1, 2), "b": (1, 3), "c": (1, 4)},
    "mixed": {"a": (1, 2), "b": (3, 1), "c": (4, 1)},
}


@selfcert
@pytest.mark.parametrize("orient", ["subspace", "source", "mixed"])
def test_d4_star_exchange_graph_complete_all_orientations(orient):
    """M1 / Task 0: #sTt(D_4) = 50 with status='complete' on every star orientation.
    Before the fix the mixed star returned status='error' (a spurious P45 mutate miss on
    summands 0/3 of one full-support tau-tilting module) even though all 50 vertices were
    found; the fix arbitrates the decomposable-cone case and restores status='complete'."""
    from quiverlab.tautilting.mutation import exchange_graph
    A = Quiver([1, 2, 3, 4], _D4_STARS[orient]).algebra(field=QQ)
    eg = exchange_graph(A, budget_pairs=2000)
    assert eg.status == "complete"
    assert eg.is_complete is True
    assert len(eg.vertices) == 50
    assert eg.n_regular is True


def _matching_index(p1, p0):
    """The exchangeable summand of ``p1`` whose removal recovers the almost-complete pair
    shared with ``p0`` -- the g-column of ``p1`` NOT present in ``p0``."""
    from quiverlab.tautilting.rigid import g_columns
    g0 = {tuple(c) for c in g_columns(p0)}
    for k, c in enumerate(g_columns(p1)):
        if tuple(c) not in g0:
            return k
    raise AssertionError("no matching exchange index")


@selfcert
def test_d4_mixed_star_previously_failing_pair_certifies():
    """The mutate fix must certify (not merely relabel): the previously-failing full-support
    tau-tilting module {(0,1,0,0),(1,1,1,0),(1,1,0,1),(1,1,1,1)} raised at summands 0 and 3
    ('no valid exchange found'); now every summand exchanges and the involution holds. (The
    n_regular=True check above already certifies all 50 vertices; this pins the exact pair.)"""
    from quiverlab.tautilting.mutation import exchange_graph, mutate
    A = Quiver([1, 2, 3, 4], _D4_STARS["mixed"]).algebra(field=QQ)
    eg = exchange_graph(A, budget_pairs=2000)
    target = {(0, 1, 0, 0), (1, 1, 1, 0), (1, 1, 0, 1), (1, 1, 1, 1)}
    pair = None
    for rec in eg.vertices:
        dvs = {tuple(dv[w] for w in sorted(dv)) for dv in rec["summand_dimvecs"]}
        if not rec["support"] and dvs == target:
            pair = rec["pair"]
            break
    assert pair is not None, "the previously-failing full-support pair was not found"
    for k in range(4):                      # k=0 and k=3 were the failing exchanges
        nb = mutate(pair, k)                 # certifies internally; raises if it cannot
        back = mutate(nb, _matching_index(nb, pair))   # involution
        assert back.g_key() == pair.g_key()
