"""tau-tilting mutation + the exchange-graph BFS (Plan 45 / C4). Literature: the number
of support tau-tilting modules of the LINEAR kA_n is Catalan(n+1) (A1->2, A2->5, A3->14);
the exchange graph is n-regular; a tau-tilting-INFINITE algebra (2-Kronecker) trips the
budget loudly. Self-cert: mutation is an involution and swaps exactly one g-column."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.mutation import exchange_graph, mutate
from quiverlab.tautilting.pairs import initial_pair

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _catalan(m):
    from math import comb
    return comb(2 * m, m) // (m + 1)


@lit
@pytest.mark.parametrize("n, count", [(1, 2), (2, 5), (3, 14)])
def test_support_tau_tilting_count_is_catalan(n, count):
    A = linear_path_algebra(n, field=QQ)
    eg = exchange_graph(A, budget_pairs=512)
    assert eg.is_complete and eg.status == "complete"
    assert len(eg.vertices) == count == _catalan(n + 1)


@lit
def test_exchange_graph_is_n_regular():
    A = linear_path_algebra(3, field=QQ)          # n = 3
    eg = exchange_graph(A)
    assert eg.is_complete and eg.n_regular
    deg = {i: 0 for i in range(len(eg.vertices))}
    for (i, j) in eg.arrows:
        deg[i] += 1
        deg[j] += 1                                # undirected regularity
    assert all(d == 3 for d in deg.values())


@selfcert
def test_mutation_is_an_involution_and_swaps_one_column():
    A = linear_path_algebra(2, field=QQ)          # n = 2, pentagon exchange graph
    p0 = initial_pair(A)
    for k in range(2):
        p1 = mutate(p0, k)
        assert p1.g_key() != p0.g_key()
        # exactly one g-column changed
        shared = p0.g_key() & p1.g_key()
        assert len(shared) == 1                    # n - 1 = 1
        # involution: mutating back at the matching summand returns p0
        back = mutate(p1, _matching_index(p1, p0))
        assert back.g_key() == p0.g_key()


@selfcert
def test_wild_algebra_trips_budget_loudly():
    # the 2-Kronecker is tau-tilting-INFINITE: the BFS cannot close.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    eg = exchange_graph(A, budget_pairs=10)      # any finite budget trips (infinite fan)
    assert eg.is_complete is False and eg.status == "budget"


def _matching_index(p1, p0):
    # the exchangeable summand of p1 whose removal recovers the shared almost-complete
    # pair with p0: the g-column of p1 NOT present in p0.
    from quiverlab.tautilting.rigid import g_columns
    g0 = {tuple(c) for c in g_columns(p0)}
    for k, c in enumerate(g_columns(p1)):
        if tuple(c) not in g0:
            return k
    raise AssertionError("pairs are not adjacent")
