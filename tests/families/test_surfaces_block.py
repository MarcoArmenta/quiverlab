"""surface_block descriptor + the AG (gentle-dictionary) tie (Plan 48, deep).
Self-cert: the block's invariants match the surface; gentle Jacobians carry an AG
invariant equal to strings.ag.ag_invariant of the SAME algebra."""
import pytest

from quiverlab.fields import QQ
from quiverlab.strings.ag import ag_invariant
from quiverlab.surfaces.block import surface_block
from quiverlab.surfaces.qp import jacobian_of
from quiverlab.surfaces.triangulation import annulus_triangulation, fan_triangulation

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
def test_block_invariants_match_surface():
    T = fan_triangulation(6)                              # A_3
    b = surface_block(T, field=QQ)
    assert b["kind"] == "surface" and b["arc_count"] == 3
    assert b["euler_characteristic"] == 1 and b["is_gentle"] is True
    assert b["dim"] == jacobian_of(T, field=QQ).dim


@xeng
@pytest.mark.parametrize("factory", [lambda: fan_triangulation(7),
                                     lambda: annulus_triangulation(2, 2)])
def test_block_ag_invariant_matches_engine(factory):
    T = factory()
    b = surface_block(T, field=QQ)
    J = jacobian_of(T, field=QQ)
    engine_ag = [list(pair) for pair in ag_invariant(J).pairs]
    assert b["ag_invariant"] == engine_ag
