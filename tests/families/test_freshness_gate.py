"""Plan 06 freshness gate. Every assumption Plan 06 builds on is probed here.
If this file is not green, the ground truth drifted -- STOP and re-baseline the
plan before writing any family code."""
import importlib

import pytest

from quiverlab.combinat import Quiver
from quiverlab.core import Algebra
from quiverlab.fields import CC, GF
from quiverlab.fields.primefield import PrimeField
from quiverlab.hochschild.table import HHTable


def test_general_quiver_algebra_is_delivered():
    """Plan 03 must be on disk: a non-monomial relation lowers to an Algebra."""
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)   # general route
    assert isinstance(A, Algebra)
    assert A.dim == 9
    assert A.hochschild_cohomology(1).dims == [1, 0]


def test_monomial_route_and_starter_builders_unchanged():
    from quiverlab.families.basic import linear_path_algebra, truncated_polynomial
    assert linear_path_algebra(3).dim == 6          # 1->2->3, dim n(n+1)/2
    assert truncated_polynomial(3).dim == 3         # k[x]/(x^3): e, x, x^2


def test_cartan_convention_is_paths_source_to_target():
    from quiverlab.invariants.cartan import cartan_matrix
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=CC)  # kA_2
    assert cartan_matrix(A) == [[1, 1], [0, 1]]     # C[i][j] = #paths i->j


def test_domains_and_hhtable_shapes():
    assert isinstance(GF(5), PrimeField)
    assert CC.parse_entry("i") is not None          # non-rational exact token
    t = HHTable([1, 0, 0], "HH^", "x")
    assert t.dims == [1, 0, 0] and t[0] == 1


def test_groebner_reduction_system_surface_present():
    g = importlib.import_module("quiverlab.groebner")
    assert hasattr(g, "build_reduction_system") and hasattr(g, "ReductionSystem")


@pytest.mark.skipif(
    importlib.util.find_spec("quiverlab.engine.resolutions_cs") is None,
    reason="Plan 04 CS backend not yet delivered",
)
def test_cs_periodic_backend_present_for_zoo_and_quantum():
    from quiverlab.engine import resolutions_periodic as rp
    assert hasattr(rp, "QuantumCIResolution")
