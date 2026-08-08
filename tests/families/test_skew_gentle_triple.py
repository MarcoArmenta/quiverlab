"""The skew-gentle triple, the recognizer, and the associated gentle algebra (Plan 68).
Self-cert: a valid triple validates; the associated gentle algebra IS gentle; dim law
is set up (checked against the split in Task 2). Literature: HZZ Def 1.1/1.3 -- the
recognizer accepts the gentle-pair triples and refuses non-gentle ones."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_gentle
from quiverlab.skewgentle.triple import (SkewGentleTriple, associated_gentle,
                                         is_skew_gentle_triple)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _Q_arrow():                      # Q: 1 --a--> 2
    return Quiver([1, 2], {"a": (1, 2)})


def test_valid_triple_special_vertex_2():
    assert is_skew_gentle_triple(_Q_arrow(), relations=[], special={2})


@selfcert
def test_associated_gentle_is_gentle_and_admissible():
    t = SkewGentleTriple.make(_Q_arrow(), relations=[], special={2})
    Ag = associated_gentle(t, field=QQ)
    assert is_gentle(Ag)                              # A^g genuinely gentle (nilpotent loop)
    assert Ag.dim == 5                                # e1,e2,a,eps,a*eps  (HZZ Lemma 1.5)


@lit
def test_sp_empty_is_the_plain_gentle():
    # Sp = empty: the associated gentle algebra IS the plain gentle A(Q,I).
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    t = SkewGentleTriple.make(Q, relations=["a*b"], special=set())
    Ag = associated_gentle(t, field=QQ)
    plain = Q.algebra(relations=["a*b"], field=QQ)
    assert Ag.dim == plain.dim == 5


@selfcert
def test_special_vertex_outside_Q0_refused():
    with pytest.raises(QuiverlabError):
        SkewGentleTriple.make(_Q_arrow(), relations=[], special={99}).validate()


@selfcert
def test_non_gentle_associated_pair_refused():
    # a vertex with 3 arrows out cannot be a gentle pair -> not a skew-gentle triple.
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (1, 4)})
    assert is_skew_gentle_triple(Q, relations=[], special=set()) is False


@selfcert
def test_length3_relation_refused():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)})
    with pytest.raises(QuiverlabError):
        SkewGentleTriple.make(Q, relations=["a*b*c"], special=set()).validate()


@selfcert
def test_special_vertex_with_existing_loop_refused():
    # HZZ Lemma 1.6: at most one loop per vertex; a special vertex may not already
    # carry a loop in Q.
    Q = Quiver([1, 2], {"a": (1, 2), "L": (2, 2)})
    with pytest.raises(QuiverlabError):
        SkewGentleTriple.make(Q, relations=[], special={2}).validate()
