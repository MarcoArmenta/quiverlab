"""Ringel slice-theorem reconstruction (Plan 60 / R17, Liu Thm 1.9(2)): H = End_A(S) presented
via P44 is hereditary of the reported type; A = End_H(D(S)) is theorem-guaranteed and the
checkable invariants (dim, Cartan, section-graph type) are reported. The tilting-but-not-slice
discriminator: A_A of kA3/rad^2 is a tilting module whose End has a relation (NOT hereditary), so
it is NOT a slice -- Thm 1.9(2) rejects it. QQ-scope (presented_form is char 0 / char > dim)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.tilted import tilted_check, _certify_slice

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA3_radsq():
    return RadicalSquareZero(Quiver([1, 2, 3], {"a1": (2, 1), "a2": (3, 2)}), field=QQ)


@lit
def test_reconstructed_H_is_hereditary_kA3():
    rep = tilted_check(_kA3_radsq())
    H = rep.hereditary_algebra
    from quiverlab.invariants.recognizers import is_hereditary
    assert is_hereditary(H)                              # End_A(S) presented = hereditary kA3
    assert len(list(H.quiver.vertices)) == 3 and not H.relations
    assert rep.reconstruction["type"] == "A_3"
    assert rep.reconstruction["dim_H"] == H.dim


@selfcert
def test_tilting_but_not_slice_rejected():
    # A_A = S_1 (+) P_2 (+) P_3 of kA3/rad^2 IS a tilting module (pd 0) but End_A(A_A) = A has a
    # relation => NOT hereditary => NOT a slice (Ringel 1.9(2)). _certify_slice must return False.
    A = _kA3_radsq()
    from quiverlab.modules.morphism import direct_sum
    S = direct_sum(A.simple(1), A.projective(2), A.projective(3))[0]   # = A_A (S_1 = P_1)
    from quiverlab.modules.tilting import is_tilting_module
    assert is_tilting_module(S).is_tilting is True       # tilting ...
    ok, _H = _certify_slice(A, S)
    assert ok is False                                   # ... but not a slice


@xeng
def test_reconstruction_roundtrip_invariants():
    # A = End_H(D(S)) is theorem-guaranteed; check the checkable invariants agree.
    A = _kA3_radsq()
    rep = tilted_check(A)
    assert rep.reconstruction["dim_A"] == A.dim
    # section type (via _section_quiver + dynkin_type) == recovered-Gabriel type (routes agree)
    assert rep.reconstruction["type"] == rep.hereditary_type == "A_3"


@lit
def test_hereditary_reconstruction_is_self():
    A = linear_path_algebra(3, field=QQ)                 # H = A
    rep = tilted_check(A)
    assert rep.hereditary_type == "A_3"
