"""Plan 03 acceptance: a general (non-monomial) kQ/I built by Quiver.algebra, with
Hochschild cohomology computed through the existing Plan-01 bar path, exact and
characteristic-independent. Also asserts the ReductionSystem shape Plan 04 consumes."""
from dataclasses import fields as dc_fields

from quiverlab.combinat import Quiver
from quiverlab.fields import CC, GF
from quiverlab.groebner import ReductionSystem, build_reduction_system


def _commutative_square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_readme_general_kqi_hochschild():
    """The README example: commutative square kQ/(a*b - c*d), a genuine non-monomial
    presentation, lowered through Groebner completion and fed to the bar complex."""
    Q = _commutative_square()
    A = Q.algebra(relations=["a*b - c*d"], field=CC)      # non-monomial -> Groebner route
    assert A.dim == 9
    # HH^0 = center = 1 ; HH^1 = 0 ; HH^2 = 0  (kA_2 x kA_2 by Kunneth; both factors
    # are trees, HH^{>=1}(kA_2)=0). HH^2 via the bar oracle on a 9-dim algebra is a
    # ~4608x576 exact rank: ~1.5s over GF(32003), ~17s over exact CC -- the sparse
    # differential lets rref terminate early; well within the max_cells guard.
    hh = A.hochschild_cohomology(2)
    assert hh.dims == [1, 0, 0]

    # exact over any field, same dimensions:
    Ap = Q.algebra(relations=["a*b - c*d"], field=GF(32003))
    assert Ap.dim == 9
    assert Ap.hochschild_cohomology(2).dims == hh.dims


def test_reduction_system_is_the_plan04_contract():
    rs = build_reduction_system(_commutative_square(), ["a*b - c*d"], CC)
    # exact frozen field set + methods Chouhy-Solotar (Plan 04) depends on
    assert {f.name for f in dc_fields(ReductionSystem)} == {
        "quiver", "domain", "order", "rules", "irreducibles", "degree_bound", "is_confluent"}
    assert rs.leading_words() == (("c", "d"),)
    assert rs.irreducibles == (("a",), ("b",), ("c",), ("d",), ("a", "b"))
    # every ambiguity of the completed system resolves to zero (CS reduction-finiteness)
    from quiverlab.groebner.complete import s_polynomial
    for amb in rs.ambiguities():
        assert rs.reduce(s_polynomial(amb, rs.domain)) == {}
