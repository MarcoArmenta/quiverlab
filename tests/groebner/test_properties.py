"""Property/invariance tests for the Groebner engine (spec §5 c.3; §8 ring 4)."""
import itertools
import random

from quiverlab.combinat import Quiver
from quiverlab.combinat.relations import parse_relations
from quiverlab.core.monomial import build_monomial_algebra
from quiverlab.fields import CC, GF, QQ
from quiverlab.groebner.reduction import lc_add
from quiverlab.groebner.system import build_reduction_system
from quiverlab.groebner.lower import groebner_algebra


def _reduction_sites(work, rules):
    """Every reduction SITE of the combination: each (word, rule, position)
    triple where a rule's lead occurs at that position in that word. Scans ALL
    positions and ALL rules (not just the leftmost, as first_factor would) so the
    caller can pick uniformly among the genuine branch points of the rewriting
    system rather than following one forced path."""
    sites = []
    for w in work:
        for rule in rules:
            lead = rule.lead
            n = len(lead)
            for i in range(len(w) - n + 1):
                if w[i:i + n] == lead:
                    sites.append((w, rule, i))
    return sites


def _random_reduce(comb, rules, dom, rng):
    """Reduce by rewriting at a RANDOMLY chosen site among ALL available
    (word, rule, position) triples. For a confluent system every such sequence of
    choices reaches the same normal form as the deterministic reducer. Returns
    (normal_form, max_sites), where max_sites is the largest number of
    simultaneously-available sites seen during the run — a witness that the walk
    actually faced multi-way order choices rather than a single forced move."""
    work = {w: c for w, c in comb.items() if not dom.is_zero(c)}
    max_sites = 0
    while True:
        sites = _reduction_sites(work, rules)
        max_sites = max(max_sites, len(sites))
        if not sites:
            return work, max_sites
        w, rule, i = rng.choice(sites)
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
    max_sites_seen = 0
    for w in words:
        deterministic = rs.reduce({w: QQ.coerce(1)})
        for _ in range(5):
            nf, max_sites = _random_reduce({w: QQ.coerce(1)}, rs.rules, rs.domain, rng)
            assert nf == deterministic
            max_sites_seen = max(max_sites_seen, max_sites)
    # Permanent anti-tautology guard: the random walk must have faced a genuine
    # branch point (>= 2 simultaneously-available reduction sites) at least once.
    # If it never did, every move was forced and the test proves nothing about
    # confluence -- fail loudly rather than pass vacuously.
    print(f"\nmax simultaneous reduction sites seen: {max_sites_seen}")
    assert max_sites_seen >= 2, (
        f"confluence test degenerated into a tautology: never saw >= 2 reduction "
        f"sites (max={max_sites_seen}); it is not exercising reduction-order choices")


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
