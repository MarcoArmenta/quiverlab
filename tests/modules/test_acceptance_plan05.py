"""Plan-05 acceptance: the modules + invariants surface, end to end (spec 3.5/3.6/3.9)."""
import sympy as sp
from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial, sweep

t = sp.Symbol("t")


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def test_module_surface():
    A = _square()
    S1, S4 = A.simple(1), A.simple(4)
    assert A.projective(1).dimension_vector() == {1: 1, 2: 1, 3: 1, 4: 1}
    assert A.injective(4).dimension_vector() == {1: 1, 2: 1, 3: 1, 4: 1}
    assert A.projective(1).radical().radical().dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}
    assert A.hom(S1, S1) == 1 and A.hom(S1, S4) == 0
    assert A.ext(S1, S4, 2) == 1
    assert int(A.global_dimension()) == 2
    res = S1.projective_resolution(4)
    assert res.pd() == 2 and res.betti(0) == 1


def test_invariant_surface():
    A = _square()
    assert A.loewy_length() == 3
    assert A.center()[0] == 1
    assert A.coxeter_polynomial().domain == sp.ZZ            # unimodular
    from quiverlab.invariants.spectral import spectral_radius
    assert spectral_radius(A.coxeter_polynomial().as_expr()) is not None
    assert truncated_polynomial(2, field=GF(32003)).complexity(6) == 1


def test_sweep_surface():
    tab = sweep(truncated_polynomial, 3, fields=[CC, GF(2), GF(3)])
    assert tab.cell("dimension", GF(2)) == 3
