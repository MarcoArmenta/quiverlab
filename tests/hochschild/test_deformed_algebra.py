"""Plan 78 Task 3: the presented deformed algebra A_alpha + the Ext-algebra handoff.

The radical direction QuantumCI(0) -> QuantumCI(-1) is the flat two-ways build (Pin 4);
the unit direction (k[x]/x^2, f(x,x)=1) is refused at PARSE level by RelationError (Pin 1),
and an infinite direction relays NotFiniteDimensionalError.
"""
import pytest

import quiverlab as ql
from quiverlab.fields import QQ
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild.deformations import (
    deformed_algebra, _canonical_radical_direction, _ext_summary)

pytestmark = pytest.mark.deep


@pytest.mark.oracle_crossengine
@pytest.mark.oracle_literature
def test_deformed_quantumci0_to_minus1_is_flat_and_matches():
    """deformed_algebra(QuantumCI(0), {'x*y':'y*x'}) builds k<x,y>/(x^2,y^2,xy-yx) =
    QuantumCI(-1): flat (dim 4) and its HH equals QuantumCI(-1) built directly (Pin 4)."""
    A = ql.QuantumCI(0, field=QQ)
    A_alpha = deformed_algebra(A, {"x*y": "y*x"}, t="1")
    assert A_alpha.dim == 4 == A.dim
    direct = ql.QuantumCI(-1, field=QQ)
    got = [A_alpha.hochschild_cohomology(3, engine="cs")[i] for i in range(4)]
    want = [direct.hochschild_cohomology(3, engine="cs")[i] for i in range(4)]
    assert got == want == [4, 4, 5, 6]


@pytest.mark.oracle_selfcert
def test_unit_direction_refused_at_parse_level():
    """A unit direction on k[x]/x^2 (f(x,x)=1) -> loud QuiverlabError relaying the
    parse-level RelationError (the length-0 term 'has no arrows', Pin 1)."""
    kx2 = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=QQ)
    with pytest.raises(QuiverlabError) as ei:
        deformed_algebra(kx2, {"x*x": "1"}, t="1")
    msg = str(ei.value)
    assert "RelationError" in msg
    assert "no arrows" in msg
    assert "radical" in msg.lower()


@pytest.mark.oracle_selfcert
def test_nonadmissible_direction_relays_loudly():
    """A parseable-but-non-admissible direction (a length-1 perturbation -> the ideal
    leaves rad^2) relays AdmissibilityError as a loud QuiverlabError (the distinct case
    from the parse-level unit RelationError; the NotFiniteDimensionalError infinite case
    goes through the SAME relay except-clause)."""
    A = ql.QuantumCI(0, field=QQ)
    with pytest.raises(QuiverlabError) as ei:
        deformed_algebra(A, {"x*y": "x"}, t="1")     # length-1 term -> not admissible
    assert "AdmissibilityError" in str(ei.value)
    assert "radical" in str(ei.value).lower()


@pytest.mark.oracle_selfcert
def test_nonflat_direction_refused():
    """A radical direction that changes the dimension (non-flat) is refused loudly --
    it is not a flat formal deformation."""
    A = ql.QuantumCI(0, field=QQ)
    with pytest.raises(QuiverlabError) as ei:
        deformed_algebra(A, {"x*x": "x*y"}, t="1")   # dim jumps 4 -> 6
    assert "not flat" in str(ei.value).lower()


@pytest.mark.oracle_literature
def test_ext_algebra_handoff_on_deformed():
    """Pin 5: deformed_algebra(...).ext_algebra(3) returns a YonedaPresentation (the RRR
    handoff), and _ext_summary describes it."""
    A = ql.QuantumCI(0, field=QQ)
    A_alpha = deformed_algebra(A, {"x*y": "y*x"}, t="1")
    E = A_alpha.ext_algebra(3)
    assert type(E).__name__ == "YonedaPresentation"
    summary = _ext_summary(E)
    assert summary.startswith("E(A_alpha):")
    assert "koszul=" in summary


@pytest.mark.oracle_selfcert
def test_canonical_radical_direction_found_for_quantumci0():
    """The report's canonical radical direction search finds a flat perturbation of a
    monomial relation for QuantumCI(0)."""
    A = ql.QuantumCI(0, field=QQ)
    direction, A_alpha = _canonical_radical_direction(A)
    assert direction is not None
    assert A_alpha is not None and A_alpha.dim == A.dim
    # it perturbs a base monomial relation by a parallel path
    (r, p), = direction.items()
    assert r in [str(x) for x in A.relations]


@pytest.mark.oracle_selfcert
def test_direction_validation_rejects_unknown_relation():
    A = ql.QuantumCI(0, field=QQ)
    with pytest.raises(QuiverlabError):
        deformed_algebra(A, {"z*z": "y*x"}, t="1")   # not a relation of A
    with pytest.raises(QuiverlabError):
        deformed_algebra(A, {}, t="1")               # empty direction
