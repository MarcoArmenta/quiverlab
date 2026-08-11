"""The order complex Delta(P) of a Poset: chains-by-dimension, face vector, lattice /
boundedness predicates; and the incidence-algebra provenance stash A._poset (Plan 75 / R9)."""
import pytest

from quiverlab.families import IncidenceAlgebra
from quiverlab.families.poset import Poset
from quiverlab.fields import CC

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _b3():
    return [(x, x | (1 << b)) for x in range(8) for b in range(3) if not (x >> b) & 1]


CROWN_C33 = [("a", "X"), ("b", "X"), ("b", "Y"),
             ("c", "Y"), ("c", "Z"), ("a", "Z")]


@lit
def test_b3_face_vector():
    P = Poset(_b3())
    assert P.order_complex().face_vector() == (8, 19, 18, 6)      # LIVE-VERIFIED
    assert P.is_bounded() and P.is_lattice()


@lit
def test_crown_c33_is_s1_face_vector():
    P = Poset(CROWN_C33)
    # S^1: six vertices, six edges, and NO 2-chains -- a crown has height 1.
    assert P.order_complex().face_vector() == (6, 6)
    assert not P.is_lattice() and not P.is_bounded()


@lit
def test_crown_c22_face_vector():
    P = Poset([("a", "X"), ("a", "Y"), ("b", "X"), ("b", "Y")])
    assert P.order_complex().face_vector() == (4, 4)
    assert not P.is_lattice()


@selfcert
def test_chains_are_strictly_increasing_and_unique():
    oc = Poset(_b3()).order_complex()
    P = oc.poset
    for level in oc.chains:
        assert len(set(level)) == len(level)              # no chain enumerated twice
        for ch in level:
            for i in range(len(ch) - 1):
                assert P.leq(ch[i], ch[i + 1]) and ch[i] != ch[i + 1]


@selfcert
def test_euler_characteristic_matches_face_vector():
    oc = Poset(CROWN_C33).order_complex()
    assert oc.euler_characteristic() == 6 - 6 == 0        # chi(S^1) = 0


@selfcert
def test_global_bound_is_weaker_than_bounded():
    """A poset with only a global MINIMUM still has a contractible (cone) order complex,
    so `has_global_bound` -- not `is_bounded` -- is what the vanishing theorem needs."""
    P = Poset([("w", "x"), ("w", "y")])                   # a global min, no global max
    assert P.has_global_bound() and not P.is_bounded()


@selfcert
def test_incidence_provenance_stashed():
    A = IncidenceAlgebra([(1, 2), (2, 3)], field=CC)
    assert A._poset is not None and A._poset.leq(1, 3)     # provenance for the fast path
