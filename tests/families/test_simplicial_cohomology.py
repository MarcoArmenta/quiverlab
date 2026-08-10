"""The order-complex simplicial engine: exact field cohomology (rank) + integral homology
(Smith normal form). The SECOND, independent HH oracle for incidence algebras (Plan 75 / R9).

Char-sensitivity via universal coefficients: RP^2 has H^*(GF2) = [1,1,1] but H^*(QQ) = [1,0,0],
and the INTEGER computation explains why (H_1(RP^2; ZZ) = Z/2) in one pass, for every field.
"""
import itertools

import pytest

from quiverlab.families.poset import Poset
from quiverlab.fields import GF, QQ
from quiverlab.hochschild.simplicial import (
    OrderComplex, cohomology_dims_from_integral, integral_homology,
    simplicial_cohomology_dims)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert

CROWN = [("a", "X"), ("b", "X"), ("b", "Y"), ("c", "Y"), ("c", "Z"), ("a", "Z")]

# The minimal 6-vertex RP^2 triangulation. Independently validated in the plan: every edge
# lies in exactly 2 facets, f = (6, 15, 10), Euler characteristic 1 = chi(RP^2).
RP2_FACETS = [(1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 6), (1, 2, 6),
              (2, 3, 5), (2, 4, 5), (2, 4, 6), (3, 4, 6), (3, 5, 6)]


def _face_poset(facets):
    """The face poset of a simplicial complex: elements are the faces, covers are
    inclusions of codimension 1. This is Gerstenhaber-Schack's ORIGINAL setting."""
    simp = set()
    for f in facets:
        for k in range(1, len(f) + 1):
            for c in itertools.combinations(sorted(f), k):
                simp.add(c)
    covers = [(t[:i] + t[i + 1:], t) for t in simp if len(t) >= 2
              for i in range(len(t)) if (t[:i] + t[i + 1:]) in simp]
    return Poset(covers, elements=sorted(simp, key=lambda s: (len(s), s)))


# ------------------------------------------------------------------ field cohomology
@lit
def test_crown_is_circle():
    oc = OrderComplex.of(Poset(CROWN))
    assert simplicial_cohomology_dims(oc, 3, QQ) == [1, 1, 0, 0]        # S^1


@lit
def test_rp2_char_sensitive():
    oc = OrderComplex.of(_face_poset(RP2_FACETS))
    assert simplicial_cohomology_dims(oc, 2, QQ) == [1, 0, 0]           # LIVE-VERIFIED
    assert simplicial_cohomology_dims(oc, 2, GF(2)) == [1, 1, 1]        # torsion in H_1
    assert simplicial_cohomology_dims(oc, 2, GF(3)) == [1, 0, 0]


@lit
def test_b3_is_contractible():
    """B_3 is bounded, so Delta(B_3) is a cone: H^{>=1} = 0 -- the Plan-33 pin, now a
    theorem instance rather than a coincidence."""
    P = Poset([(x, x | (1 << b)) for x in range(8) for b in range(3) if not (x >> b) & 1])
    assert simplicial_cohomology_dims(P.order_complex(), 2, QQ) == [1, 0, 0]


# ------------------------------------------------------------------ integral homology
@selfcert
def test_integral_homology_torsion_certificate():
    oc = OrderComplex.of(_face_poset(RP2_FACETS))
    ih = integral_homology(oc, 2)               # [(free_rank, invariant_factors), ...]
    assert ih[0] == (1, ())                     # H_0 = Z
    assert ih[1] == (0, (2,))                   # H_1 = Z/2  -> the char-2 sensitivity
    assert ih[2] == (0, ())                     # H_2 = 0


@selfcert
@pytest.mark.parametrize("dom,char", [(QQ, 0), (GF(2), 2), (GF(3), 3)])
def test_universal_coefficients_reproduces_the_field_answer(dom, char):
    """The payoff: ONE integer Smith normal form must reproduce the field answer in EVERY
    characteristic. If this ever disagrees with the direct field rank, the fast path's
    whole claim (every characteristic from one computation) is void."""
    oc = OrderComplex.of(_face_poset(RP2_FACETS))
    ih = integral_homology(oc, 2)
    assert cohomology_dims_from_integral(ih, 2, char) == simplicial_cohomology_dims(oc, 2, dom)


@selfcert
def test_boundary_squares_to_zero():
    for P in (Poset(CROWN), _face_poset(RP2_FACETS)):
        oc = OrderComplex.of(P)
        for k in range(1, oc.top_dimension):
            A, B = oc.boundary(k), oc.boundary(k + 1)
            if not A or not B:
                continue
            rows, cols, inner = len(A), len(B[0]), len(B)
            assert inner == len(A[0])
            for i in range(rows):
                for j in range(cols):
                    assert sum(A[i][t] * B[t][j] for t in range(inner)) == 0


@selfcert
def test_euler_characteristic_is_the_alternating_betti_sum():
    """chi = sum (-1)^p f_p = sum (-1)^p dim H_p, over any field -- an independent
    consistency tie between the face vector and the computed ranks."""
    for P in (Poset(CROWN), _face_poset(RP2_FACETS)):
        oc = OrderComplex.of(P)
        dims = simplicial_cohomology_dims(oc, oc.top_dimension, QQ)
        assert oc.euler_characteristic() == sum((-1) ** p * d for p, d in enumerate(dims))


@selfcert
def test_rp2_triangulation_is_valid():
    """Guards the fixture itself: every edge in exactly 2 facets, f = (6, 15, 10),
    chi = 1. A mistyped facet list would otherwise silently pin the wrong space."""
    edges = {}
    for f in RP2_FACETS:
        for e in itertools.combinations(sorted(f), 2):
            edges[e] = edges.get(e, 0) + 1
    assert set(edges.values()) == {2}
    verts = {v for f in RP2_FACETS for v in f}
    assert (len(verts), len(edges), len(RP2_FACETS)) == (6, 15, 10)
    assert len(verts) - len(edges) + len(RP2_FACETS) == 1        # chi(RP^2) = 1
