"""Plan 75 Task I4 -- the RING side of the Gerstenhaber-Schack theorem, honestly scoped.

Gerstenhaber-Schack prove an isomorphism of graded RINGS `HH^*(kP) = H^*(Delta(P); k)`.
What v1 ships and pins here is:

  * the cup product on an incidence algebra satisfies the Gerstenhaber-algebra identities
    (graded commutativity + associativity) -- a re-use of the shipped Plan-35 surface on a
    new input class, so it is a `selfcert`, not a new claim; and
  * the CS cup engine's `HH^p` DIMENSIONS agree with the order-complex route degreewise --
    an independent second computation of the same graded vector space (the ring engine
    resolves `kP` over `(kP)^e`; the order-complex route never builds a resolution at all).

What v1 does NOT ship is a pin of the cup PAIRING RANK against the simplicial cup pairing.
That needs a poset whose order complex has a nonzero product `H^1 (x) H^1 -> H^2` -- i.e. a
closed surface -- and the smallest ones are out of the test budget (measurements in
`test_cup_rank_against_a_surface` below). So the ring isomorphism is recorded on the
verification page as **dimension-level in v1**, with the cup-rank check deferred and
runnable on demand. No cup pin is fabricated.

Note the two quantities kept apart throughout: the order-complex route computes `HH^*`
(cohomology), which is NOT `HH_*` (Hochschild homology) -- e.g. the diamond incidence
algebra has `HH^* = [1,0,0]` but `HH_* = [4,0,0,0,0]`.
"""
import itertools
import os
from fractions import Fraction

import pytest

import quiverlab as ql
from quiverlab.fields import QQ

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert

# C(3,3): 3 minimal + 3 maximal, Delta ~ S^1 (a hexagon). Not a lattice.
CROWN33 = [("a", "X"), ("b", "X"), ("b", "Y"), ("c", "Y"), ("c", "Z"), ("a", "Z")]
# C(2,2): the 4-cycle, Delta ~ S^1.
CROWN22 = [[1, 3], [1, 4], [2, 3], [2, 4]]
# The Boolean lattice B_3: a global bound, so Delta is a cone.
B3 = [[x, x | (1 << b)] for x in range(8) for b in range(3) if not (x >> b) & 1]

_POSETS = [("crown33", CROWN33), ("crown22", CROWN22), ("b3", B3)]


def _tensor(table):
    """``constants`` as an exact rational tensor ``[k][i][j]`` (the tables are EXACT
    strings; over QQ they may be fractions, so never ``int()`` them)."""
    return [[[Fraction(c) for c in row] for row in mat] for mat in table.constants]


@selfcert
@pytest.mark.parametrize("name,covers", _POSETS)
def test_cup_is_graded_commutative_on_an_incidence_algebra(name, covers):
    """`f (cup) g = (-1)^{pq} g (cup) f` on `HH^*(kP)` -- the shipped Plan-35 identity, on
    a new input class. Table-level, so it is basis-independent as an assertion about the
    ring even though the constants themselves are basis-dependent."""
    A = ql.IncidenceAlgebra(covers, field=QQ)
    hp = A.cup_products(2)
    checked = 0
    for (p, q), t in hp.tables.items():
        if (q, p) not in hp.tables:
            continue
        M, N = _tensor(t), _tensor(hp.tables[(q, p)])
        sign = 1 if (p * q) % 2 == 0 else -1
        for k in range(t.dims[2]):
            for i in range(t.dims[0]):
                for j in range(t.dims[1]):
                    assert M[k][i][j] == sign * N[k][j][i], \
                        f"graded commutativity fails at {(p, q)} on {name}"
                    checked += 1
    assert checked, "the identity must be exercised, not vacuously satisfied"


@selfcert
@pytest.mark.parametrize("name,covers", _POSETS)
def test_cup_is_associative_on_an_incidence_algebra(name, covers):
    """`(f (cup) g) (cup) h = f (cup) (g (cup) h)`, contracted at table level."""
    A = ql.IncidenceAlgebra(covers, field=QQ)
    hp = A.cup_products(2)
    for p, q, r in [(0, 0, 1), (0, 1, 1), (1, 1, 0)]:
        need = [(p, q), (q, r), (p + q, r), (p, q + r)]
        if any(k not in hp.tables for k in need):
            continue
        Tpq, Tqr = hp.tables[(p, q)], hp.tables[(q, r)]
        Tpq_r, Tp_qr = hp.tables[(p + q, r)], hp.tables[(p, q + r)]
        A_, B_, C_, D_ = (_tensor(Tpq), _tensor(Tpq_r), _tensor(Tqr), _tensor(Tp_qr))
        for i, j, l in itertools.product(range(Tpq.dims[0]), range(Tpq.dims[1]),
                                         range(Tqr.dims[1])):
            for k in range(Tpq_r.dims[2]):
                left = sum(A_[m][i][j] * B_[k][m][l] for m in range(Tpq.dims[2]))
                right = sum(C_[m][j][l] * D_[k][i][m] for m in range(Tqr.dims[2]))
                assert left == right, f"associativity fails at {(p, q, r)} on {name}"


@xeng
@pytest.mark.parametrize("name,covers,dims", [("crown33", CROWN33, [1, 1, 0]),
                                              ("crown22", CROWN22, [1, 1, 0]),
                                              ("b3", B3, [1, 0, 0])])
def test_cup_engine_dims_equal_the_order_complex_dims(name, covers, dims):
    """The dimension-level half of the ring isomorphism, computed TWICE by genuinely
    independent routes: the Plan-35 cup engine resolves `kP` over `(kP)^e` (Chouhy-Solotar)
    and reports `dim HH^p` as its table dims; the Plan-75 route builds the order complex and
    never resolves anything. They must agree degreewise."""
    A = ql.IncidenceAlgebra(covers, field=QQ)
    assert A.incidence_cohomology(2).dims == dims
    hp = A.cup_products(2)
    # Each table (p, q) reports (dim HH^p, dim HH^q, dim HH^{p+q}); read the degrees off
    # whichever tables exist and demand every reading agrees with the order complex.
    seen = {}
    for (p, q), t in hp.tables.items():
        for deg, d in ((p, t.dims[0]), (q, t.dims[1]), (p + q, t.dims[2])):
            if deg < len(dims):
                assert seen.setdefault(deg, d) == d, f"cup tables disagree on HH^{deg}"
                assert d == dims[deg], (
                    f"{name}: cup engine says dim HH^{deg} = {d}, order complex says "
                    f"{dims[deg]}")
    assert set(seen) >= {0, 1, 2}


# --------------------------------------------------------------------------- #
# the deferred half: the cup PAIRING RANK against a surface
# --------------------------------------------------------------------------- #

def _surface_face_poset(facets):
    simp = set()
    for f in facets:
        for k in range(1, len(f) + 1):
            for c in itertools.combinations(sorted(f), k):
                simp.add(c)
    covers = [(t[:i] + t[i + 1:], t) for t in simp if len(t) >= 2
              for i in range(len(t)) if (t[:i] + t[i + 1:]) in simp]
    from quiverlab.families.poset import Poset
    return Poset(covers, elements=sorted(simp, key=lambda s: (len(s), s)))


# The 7-vertex (Csaszar) torus: triangles {i, i+1, i+3} and {i, i+2, i+3} mod 7. Verified
# live: every edge lies in exactly 2 facets, f = (7, 21, 14), chi = 0, and the order-complex
# route gives H^* = [1, 2, 1] -- the torus. Its incidence algebra has dim 168.
TORUS7 = sorted({tuple(sorted((i % 7, (i + 1) % 7, (i + 3) % 7))) for i in range(7)}
                | {tuple(sorted((i % 7, (i + 2) % 7, (i + 3) % 7))) for i in range(7)})


@xeng
@pytest.mark.skipif(
    os.environ.get("QUIVERLAB_INCIDENCE_CUP_RANK") != "1",
    reason="DEFERRED, budget-MEASURED (Plan 75 Task I4 step 2, the plan's permitted "
           "fallback): the cup PAIRING RANK needs a poset whose order complex is a closed "
           "surface, and both smallest candidates were run live and produced nothing. "
           "(i) The 7-vertex Csaszar torus face poset over QQ (dim kP = 168, the CS route) "
           "ran past a 25-minute box on A.cup_products(2). (ii) RP^2_6 over GF(2) "
           "(dim kP = 121) routes through the BAR cup instead, which already costs ~10s at "
           "dim 12 (crown C(3,3), measured over GF(2)/GF(3)/GF(32003)) and returned no "
           "result at dim 121. Meanwhile the order-complex route answers BOTH in 0.02s -- "
           "which is the whole point of the fast path, and also why the ring check is the "
           "expensive half. The ring isomorphism is therefore recorded on the verification "
           "page as DIMENSION-LEVEL in v1 (pinned by the two tests above); no cup value is "
           "fabricated. Set QUIVERLAB_INCIDENCE_CUP_RANK=1 to run it anyway.")
def test_cup_rank_against_a_surface():
    """On the torus the cup pairing `H^1 (x) H^1 -> H^2` is the intersection form: a
    perfect skew pairing `k^2 x k^2 -> k`, so its rank is 2 -- a BASIS-INDEPENDENT
    invariant (the Plan-35 cross-engine discipline: never compare raw constants)."""
    A = ql.IncidenceAlgebra(_surface_face_poset(TORUS7), field=QQ)
    assert A.incidence_cohomology(2).dims == [1, 2, 1]
    t = A.cup_products(2).tables[(1, 1)]
    assert t.dims == (2, 2, 1)
    M = [[_tensor(t)[0][i][j] for j in range(2)] for i in range(2)]
    # rank of the 2x2 pairing matrix over QQ
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    assert det != 0, f"the torus cup pairing must be nondegenerate; got {M}"
