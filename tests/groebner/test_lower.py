"""General kQ/I lowering + Quiver.algebra dispatch (spec §3.3, §5 components 3, 4)."""
from quiverlab.combinat import Quiver
from quiverlab.core.monomial import build_monomial_algebra
from quiverlab.combinat.relations import parse_relations
from quiverlab.fields import CC
from quiverlab.groebner import Dispatch
from quiverlab.groebner.lower import groebner_algebra


def _loop():
    return Quiver([1], {"x": (1, 1)})


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def _same_algebra(A, B):
    return (A.dim == B.dim and A.basis_labels == B.basis_labels
            and A.T == B.T and A.unit == B.unit)


def test_fixture1_monomial_route_equivalence():
    """x^3: the Groebner route and the Plan-01 monomial route agree elementwise."""
    Q = _loop()
    rels = parse_relations(["x^3"], Q)
    A_mono = build_monomial_algebra(Q, rels, CC)
    A_grob = groebner_algebra(Q, rels, CC)
    assert _same_algebra(A_mono, A_grob)
    assert A_grob.dim == 3
    assert A_grob.basis_labels == ["e_1", "x", "x*x"]


def test_fixture2_square_dim_nine():
    Q = _square()
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    assert A.dim == 9
    assert A.basis_labels == ["e_1", "e_2", "e_3", "e_4", "a", "b", "c", "d", "a*b"]


def test_square_realizes_commutativity_in_structure_constants():
    """In the algebra, c*d = a*b (the reduction cd -> ab)."""
    Q = _square()
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    labels = A.basis_labels

    def e(name):
        v = [A.domain.zero()] * A.dim
        v[labels.index(name)] = A.domain.one()
        return v

    cd = A.multiply(e("c"), e("d"))
    ab = A.multiply(e("a"), e("b"))
    assert cd == ab
    assert cd[labels.index("a*b")] == A.domain.one()


def test_dispatch_records_route():
    Q = _square()
    trace = []
    Q.algebra(relations=["a*b - c*d"], field=CC, trace=trace)
    routes = [ev.route for ev in trace if isinstance(ev, Dispatch)]
    assert routes == ["groebner"]

    trace2 = []
    Q2 = _loop()
    Q2.algebra(relations=["x^3"], field=CC, trace=trace2)
    assert [ev.route for ev in trace2 if isinstance(ev, Dispatch)] == ["monomial"]


def test_monomial_dispatch_still_uses_plan01_path():
    """A monomial input through Quiver.algebra equals the direct monomial build."""
    Q = _loop()
    A_dispatch = Q.algebra(relations=["x^3"], field=CC)
    A_mono = build_monomial_algebra(Q, parse_relations(["x^3"], Q), CC)
    assert _same_algebra(A_dispatch, A_mono)


def test_fixture6_quantum_ci_dim_four():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    A = Q.algebra(relations=["x^2", "y^2", "y*x - x*y"], field=CC)
    assert A.dim == 4
    assert A.basis_labels == ["e_1", "x", "y", "x*y"]
