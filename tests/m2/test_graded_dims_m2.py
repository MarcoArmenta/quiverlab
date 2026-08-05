"""Live M2 battery: graded dims of single-vertex kQ/I over GF(p) --
an independent nc-Groebner engine recomputes our Hilbert data."""
import pytest

from quiverlab import GF, Quiver, truncated_polynomial
from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.skipif(session.should_skip_m2(),
                                reason="Macaulay2 not installed")


CASES = [
    ("truncated_x5", lambda: truncated_polynomial(5, field=GF(32003)), 8),
    ("dual_numbers", lambda: Quiver([1], {"x": (1, 1)}).algebra(
        relations=["x*x"], field=GF(7)), 6),
    ("two_loops_radsq", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "x*y", "y*x", "y*y"], field=GF(5)), 6),
    ("exterior_2", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y+y*x"], field=GF(32003)), 6),
    ("quantum_ci_q2", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x-2*x*y"], field=GF(32003)), 6),
    ("straddle_monomial", lambda: Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y*x"], field=GF(3)), 8),
]


@pytest.mark.parametrize("name,build,top", CASES, ids=[c[0] for c in CASES])
def test_graded_dims_match_m2(name, build, top):
    cc.crosscheck_graded_dims(build(), top=top).assert_agree()
