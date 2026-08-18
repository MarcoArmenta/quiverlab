"""Pure-combinatorics graph primitives underpinning pi1 (spanning forest,
predecessors, induced subquiver, undirected components). Self-certifying."""
import pytest
from quiverlab import Quiver

pytestmark = pytest.mark.oracle_selfcert


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_predecessors_successors():
    Q = _square()
    assert Q.predecessors(4) == {2, 3}
    assert Q.successors(1) == {2, 3}
    assert Q.predecessors(1) == set()


def test_spanning_forest_betti_number():
    Q = _square()
    tree, _ = Q.spanning_forest()
    # connected: |tree| = |V| - 1; non-tree count = first Betti number = |E|-|V|+c
    assert len(tree) == 3
    non_tree = set(Q.arrows) - tree
    assert len(non_tree) == 4 - 4 + 1                  # Betti = 1


def test_undirected_components_disconnected():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)})
    comps = Q.undirected_components()
    assert sorted(map(sorted, comps)) == [[1, 2], [3, 4]]


def test_induced_subquiver_drops_boundary_arrows():
    Q = Quiver([1, 2, 3, 4, 5],
               {"al": (1, 2), "be": (2, 3), "ga": (3, 5), "de": (2, 4), "ep": (4, 5)})
    sub = Q.induced_subquiver({2, 3, 4, 5})            # drop vertex 1 and arrow "al"
    assert set(sub.vertices) == {2, 3, 4, 5}
    assert set(sub.arrows) == {"be", "ga", "de", "ep"}


def test_spanning_forest_one_tree_per_component():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)})
    tree, _ = Q.spanning_forest()
    assert len(tree) == 2                              # one edge per component, c=2
