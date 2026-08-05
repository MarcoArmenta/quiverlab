"""Jacobian algebra Jac(Q, W) = kQ/(d_a W) (Plan 44 / C7). The 3-cycle triangle is the
hand-derived pin; Jac(triangle, abc) = cyclic Nakayama kZ_3/J^2 is the cross-engine oracle;
an under-constrained potential is the honest NotFiniteDimensionalError refusal.
Derksen-Weyman-Zelevinsky 2008; Labardini-Fragoso 2009.

NOTE (Plan-44 deviation, documented): the plan's proposed 'preprojective algebra =
Jac(double Q, sum_a (a a* - a* a))' cross-engine oracle is mathematically FALSE. The
cyclic derivative of a commutator 2-cycle is identically 0 (d_a(a a* - a* a) = a* - a* =
0), and the double of a Dynkin (tree) quiver is bipartite, hence has NO oriented 3-cycles,
so no cubic potential exists whose Jacobian could reproduce the quadratic mesh relations.
A Dynkin preprojective algebra is therefore NOT a Jacobian algebra of its double. The
genuinely-true cross-engine oracle used instead: Jac(3-cycle triangle, abc) is the cyclic
Nakayama algebra kZ_3/J^2 (both dim 6, identical Cartan)."""
import pytest

from quiverlab import NakayamaAlgebra, Quiver
from quiverlab.errors import NotFiniteDimensionalError
from quiverlab.families.jacobian import JacobianAlgebra, Potential, cyclic_derivative
from quiverlab.fields import QQ


def _triangle():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)})


@pytest.mark.oracle_selfcert
def test_potential_rejects_non_cyclic_word():
    Q = _triangle()
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError):
        Potential(Q, [(1, ("a", "b"))])          # a*b : 1 -> 3, NOT a cycle


@pytest.mark.oracle_selfcert
def test_cyclic_derivatives_of_abc():
    Q = _triangle()
    W = Potential(Q, [(1, ("a", "b", "c"))])     # the 3-cycle a*b*c
    assert cyclic_derivative(W, "a") == "b*c"    # d_a(abc) = bc
    assert cyclic_derivative(W, "b") == "c*a"    # d_b(abc) = ca
    assert cyclic_derivative(W, "c") == "a*b"    # d_c(abc) = ab


@pytest.mark.oracle_literature
def test_triangle_jacobian_dimension():
    # Jac = kQ/(bc, ca, ab): all length-2 paths die => rad^2 = 0 => dim = 3 vertices +
    # 3 arrows = 6. Hand-derived.
    Q = _triangle()
    W = Potential(Q, [(1, ("a", "b", "c"))])
    J = JacobianAlgebra(Q, W, field=QQ)
    assert J.dim == 6


@pytest.mark.oracle_crossengine
def test_jacobian_of_triangle_is_cyclic_nakayama():
    # Jac(1->2->3->1, abc) = kZ_3 / J^2 (all length-2 paths die) = the cyclic Nakayama
    # algebra NakayamaAlgebra(n=3, l=2). Cross-check the DWZ cyclic-derivative engine
    # against the independent Kupisch-series builder: same dim + same Cartan matrix.
    Q = _triangle()
    W = Potential(Q, [(1, ("a", "b", "c"))])
    J = JacobianAlgebra(Q, W, field=QQ)
    N = NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ)   # kZ_3/J^2, dim 6
    assert J.dim == N.dim == 6
    assert J.cartan_matrix() == N.cartan_matrix()


@pytest.mark.oracle_selfcert
def test_infinite_jacobian_refused_loudly():
    # one vertex, two loops, W = 0: relations empty => free algebra k<x,y> => infinite.
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    W = Potential(Q, [])                          # empty potential
    with pytest.raises(NotFiniteDimensionalError):
        JacobianAlgebra(Q, W, field=QQ)
