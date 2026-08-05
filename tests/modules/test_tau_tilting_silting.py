"""2-term silting bridge (Plan 45 / C4). Self-cert: each pair maps to a 2-term silting
object (P_1->P_0 summands + P_v[1] shifts); the count == #pairs completes the AIR four-way
identity. Soft P43: the ChainComplex wrapper is present iff modules.complexes imports."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.pairs import initial_pair, terminal_pair
from quiverlab.tautilting.silting import silting_count, two_term_silting

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_initial_pair_silting_is_the_projectives():
    A = linear_path_algebra(2, field=QQ)
    s = two_term_silting(initial_pair(A))
    # (A,0): each P_v presents as [0 -> P_v], i.e. P1 empty, P0 = {v}.
    p0s = sorted(tuple(sd["P0"]) for sd in s["summands"])
    assert p0s == [(1,), (2,)] and all(sd["P1"] == [] for sd in s["summands"])


@selfcert
def test_terminal_pair_silting_is_shifted_projectives():
    A = linear_path_algebra(2, field=QQ)
    s = two_term_silting(terminal_pair(A))
    # (0,A): each killed P_v presents as P_v[1], i.e. P1 = {v}, P0 empty.
    p1s = sorted(tuple(sd["P1"]) for sd in s["summands"])
    assert p1s == [(1,), (2,)] and all(sd["P0"] == [] for sd in s["summands"])


@lit
@pytest.mark.parametrize("n, count", [(2, 5), (3, 14)])
def test_silting_leg_of_four_way_identity(n, count):
    A = linear_path_algebra(n, field=QQ)
    sc = silting_count(A)
    assert sc["complete"] and sc["count"] == count
