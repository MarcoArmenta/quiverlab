"""Property/invariance tests for the Groebner engine (spec §5 c.3; §8 ring 4)."""
import itertools
import random

from quiverlab.combinat import Quiver
from quiverlab.combinat.relations import parse_relations
from quiverlab.core.monomial import build_monomial_algebra
from quiverlab.fields import CC, GF, QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import first_factor, lc_add
from quiverlab.groebner.system import build_reduction_system
from quiverlab.groebner.lower import groebner_algebra


def _random_reduce(comb, rules, order, dom, rng):
    """Reduce by applying rules at RANDOMLY chosen reducible occurrences. For a
    confluent system this must reach the same normal form as reduce_comb."""
    work = {w: c for w, c in comb.items() if not dom.is_zero(c)}
    while True:
        reducible = [w for w in work if first_factor(w, rules) is not None]
        if not reducible:
            return work
        w = rng.choice(reducible)
        rule, i = first_factor(w, rules)
        coeff = work.pop(w)
        u, v = w[:i], w[i + len(rule.lead):]
        for tc, tw in rule.tail:
            lc_add(work, u + tw + v, dom.mul(coeff, tc), dom)


def test_unique_normal_form_random_orders_agree():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rs = build_reduction_system(Q, ["x^2", "y^2", "y*x - x*y"], QQ)   # Fixture 6, confluent
    rng = random.Random(20260718)
    letters = ["x", "y"]
    words = [tuple(p) for n in range(1, 6) for p in itertools.product(letters, repeat=n)]
    for w in words:
        deterministic = rs.reduce({w: QQ.coerce(1)})
        for _ in range(5):
            assert _random_reduce({w: QQ.coerce(1)}, rs.rules, rs.order, rs.domain, rng) \
                == deterministic


def test_monomial_route_equivalence_family():
    cases = [
        (Quiver([1], {"x": (1, 1)}), ["x^2"]),
        (Quiver([1], {"x": (1, 1)}), ["x^4"]),
        (Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}), ["a*b"]),
        # Two loops: needs y^2 too, else y^n is irreducible for all n and kQ/I is
        # infinite-dimensional (both routes correctly reject that). With y^2 the
        # algebra is finite (basis e_1, x, y) and the two routes must agree.
        (Quiver([1], {"x": (1, 1), "y": (1, 1)}), ["x^2", "y^2", "x*y", "y*x"]),
    ]
    for Q, rels in cases:
        parsed = parse_relations(rels, Q)
        A_mono = build_monomial_algebra(Q, parsed, CC)
        A_grob = groebner_algebra(Q, parsed, CC)
        assert A_mono.dim == A_grob.dim
        assert A_mono.basis_labels == A_grob.basis_labels
        assert A_mono.T == A_grob.T
        assert A_mono.unit == A_grob.unit


def test_characteristic_independence_square():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    dims = {}
    for field, key in [(CC, "CC"), (GF(32003), "p")]:
        A = Q.algebra(relations=["a*b - c*d"], field=field)
        assert A.dim == 9
        dims[key] = A.hochschild_cohomology(1).dims
    assert dims["CC"] == dims["p"]
    assert dims["CC"] == [1, 0]         # HH^0 = 1 (center), HH^1 = 0 (kA_2 x kA_2, Kunneth)
