"""HH^*(kP) via the Gerstenhaber-Schack / Cibils theorem = simplicial H^* of the order
complex, and the cross-check that it EQUALS the general engine (Plan 75 / R9).

Extends the Plan-33 B_3 pin from a coincidence to a theorem instance: B_3 is bounded, so
its order complex is a cone, so HH^{>=1} = 0. Crowns are the non-lattice posets whose order
complexes are circles (nonvanishing H^1). RP^2 is the characteristic split.
"""
import itertools

import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.families import IncidenceAlgebra
from quiverlab.fields import CC, GF, QQ
from quiverlab.resolutions_cs.homology import cs_cohomology_dims

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


def _b3():
    return [(x, x | (1 << b)) for x in range(8) for b in range(3) if not (x >> b) & 1]


CROWN = [("a", "X"), ("b", "X"), ("b", "Y"), ("c", "Y"), ("c", "Z"), ("a", "Z")]
# the thickened crown: one crown edge replaced by a diamond -- a NON-LATTICE,
# NON-MONOMIAL (genuine commutativity relations) poset that still has H^1 != 0.
THICK = [("a", "m"), ("a", "n"), ("m", "A"), ("n", "A"), ("b", "A"),
         ("b", "B"), ("c", "B"), ("c", "C"), ("a", "C")]
RP2_FACETS = [(1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 6), (1, 2, 6),
              (2, 3, 5), (2, 4, 5), (2, 4, 6), (3, 4, 6), (3, 5, 6)]


def _rp2_face_poset():
    simp = set()
    for f in RP2_FACETS:
        for k in range(1, len(f) + 1):
            for c in itertools.combinations(sorted(f), k):
                simp.add(c)
    covers = [(t[:i] + t[i + 1:], t) for t in simp if len(t) >= 2
              for i in range(len(t)) if (t[:i] + t[i + 1:]) in simp]
    return covers, sorted(simp, key=lambda s: (len(s), s))


@lit
def test_b3_contractible_theorem():
    A = IncidenceAlgebra(_b3(), field=CC)
    r = A.incidence_cohomology(4)
    assert r.dims == [1, 0, 0, 0, 0]
    assert r.contractible is True                  # PROVED (global bound => cone)
    assert "cone" in r.contractible_reason


@lit
def test_crown_nonvanishing_h1():
    r = IncidenceAlgebra(CROWN, field=CC).incidence_cohomology(4)
    assert r.dims == [1, 1, 0, 0, 0]
    assert r.contractible is False


@xeng
def test_incidence_fast_equals_general_engine():
    """The fast simplicial path == the general CS engine, degreewise, on a NON-MONOMIAL
    non-lattice poset with H^1 != 0."""
    A = IncidenceAlgebra(THICK, field=CC)
    fast = A.incidence_cohomology(3).dims
    assert fast == list(cs_cohomology_dims(A, 3).dims) == [1, 1, 0, 0]


@xeng
@pytest.mark.parametrize("covers,expected", [(CROWN, [1, 1, 0]), (_b3(), [1, 0, 0])])
def test_fast_equals_general_across_posets(covers, expected):
    A = IncidenceAlgebra(covers, field=QQ)
    assert A.incidence_cohomology(2).dims == expected
    assert list(A.hochschild_cohomology(2, verbose=False).dims) == expected


@lit
def test_rp2_char_dependent_hh():
    covers, elems = _rp2_face_poset()
    A2 = IncidenceAlgebra(covers, elements=elems, field=GF(2))
    Aq = IncidenceAlgebra(covers, elements=elems, field=QQ)
    assert A2.dim == 121                                   # LIVE-VERIFIED
    assert A2.incidence_cohomology(2).dims == [1, 1, 1]
    assert Aq.incidence_cohomology(2).dims == [1, 0, 0]
    assert A2.incidence_cohomology(2).char_dependent is True
    assert Aq.incidence_cohomology(2).char_dependent is False   # char 0 -> never


@selfcert
def test_torsion_is_reported_and_explains_the_split():
    """The reported integral torsion is not decoration: it is the REASON the GF(2) answer
    differs, and `char_dependent` must be True exactly when a factor divides the char."""
    covers, elems = _rp2_face_poset()
    r2 = IncidenceAlgebra(covers, elements=elems, field=GF(2)).incidence_cohomology(2)
    r3 = IncidenceAlgebra(covers, elements=elems, field=GF(3)).incidence_cohomology(2)
    assert 2 in r2.torsion
    assert r2.char_dependent is True and r3.char_dependent is False
    assert r3.dims == [1, 0, 0]                            # 3 does not divide 2


@selfcert
def test_no_provenance_refuses():
    from quiverlab import Quiver
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=CC)   # kA2, no poset
    with pytest.raises(QuiverlabError, match="poset"):
        A.incidence_cohomology(2)


@selfcert
def test_contractible_is_three_valued_not_guessed():
    """A poset with no global bound whose H^n happens to vanish in range must NOT be
    claimed contractible -- the report says `None` and says why."""
    # two disjoint 2-chains joined into a zigzag with no global bound but contractible-ish
    A = IncidenceAlgebra([("p", "q"), ("r", "q"), ("r", "s")], field=CC)
    r = A.incidence_cohomology(3)
    assert r.dims == [1, 0, 0, 0]
    assert r.contractible is None                  # NOT True -- never proved
    assert "not a proof" in r.contractible_reason
