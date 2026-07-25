"""QPA (GAP) as the oracle for LEFT modules (Plan 24, Tier 1b).

QPA is right-module native. A left A-module IS a right A^op-module, so left-side
quantities are crosschecked by FEEDING QPA THE OPPOSITE ALGEBRA A.opposite():
the left module's underlying right-A^op representation is exactly QPA's input, and
the left surface (A.simple(v, side="left").tau() etc.) is tied to that
QPA-validated A^op right-module computation.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import Quiver, GF, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _cn(n, field=QQ):
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i % n + 1) for i in range(1, n + 1)}
    rels = [f"a{i}*a{i % n + 1}" for i in range(1, n + 1)]
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def _right_twins(Aop):
    return ([Aop.simple(v) for v in Aop.quiver.vertices]
            + [Aop.projective(v) for v in Aop.quiver.vertices]
            + [Aop.injective(v) for v in Aop.quiver.vertices])


@pytest.mark.parametrize("A", [linear_path_algebra(2, field=QQ),
                               linear_path_algebra(3, field=QQ), _cn(3)])
def test_left_tau_crosscheck_via_opposite(A):
    # The left surface routes to the right-module tau over A^op; QPA validates that
    # A^op computation (DTr/TrD run on the opposite algebra).
    Aop = A.opposite()
    for v in A.quiver.vertices:
        left_S = A.simple(v, side="left")            # left A-module
        right_twin = Aop.simple(v)                   # same representation, right A^op-module
        assert left_S.tau().dimension_vector() == right_twin.tau().dimension_vector()
        assert (left_S.tau_minus().dimension_vector()
                == right_twin.tau_minus().dimension_vector())
    for M in _right_twins(Aop):                       # QPA on the opposite algebra
        Aop.crosscheck("tau", M).assert_agree()
        Aop.crosscheck("tau_minus", M).assert_agree()


@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ), _cn(3)])
def test_left_injdim_crosscheck_via_opposite(A):
    # inj.dim of a left A-module = inj.dim of its underlying right A^op-module (a
    # categorical invariant); QPA's InjDimensionOfModule over A^op is the oracle.
    Aop = A.opposite()
    for v in A.quiver.vertices:
        left_S = A.simple(v, side="left")
        right_twin = Aop.simple(v)
        assert left_S.injective_dimension(bound=8) == right_twin.injective_dimension(bound=8)
        Aop.crosscheck("inj_dimension", right_twin, 8).assert_agree()


def test_left_proj_and_inj_resolution_crosscheck_via_opposite():
    A = linear_path_algebra(3, field=QQ)
    Aop = A.opposite()
    for v in A.quiver.vertices:
        left_S = A.simple(v, side="left")
        right_twin = Aop.simple(v)
        # left projective/injective resolution term dim-vectors = the A^op right ones.
        assert (left_S.projective_resolution(4).dimension_vectors()
                == right_twin.projective_resolution(4).dimension_vectors())
        Aop.crosscheck("proj_resolution", right_twin, 4).assert_agree()
        Aop.crosscheck("inj_resolution", right_twin, 4).assert_agree()


@pytest.mark.parametrize("field", [GF(2), GF(3)])
def test_left_tau_dimvec_over_gfp(field):
    A = linear_path_algebra(3, field=field)
    Aop = A.opposite()
    for v in A.quiver.vertices:
        left_S = A.simple(v, side="left")
        right_twin = Aop.simple(v)
        assert left_S.tau().dimension_vector() == right_twin.tau().dimension_vector()
        rep = Aop.crosscheck("tau", right_twin)
        assert rep.ours == rep.qpa, rep              # A^op right tau dim vectors vs QPA
