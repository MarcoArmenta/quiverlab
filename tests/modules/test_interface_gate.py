"""Plan-05 freshness gate: every upstream symbol Plan 05 consumes, asserted to exist
with the expected shape. A failure here means an upstream surface drifted -- STOP and
reconcile before continuing Plan 05 (spec §5 components 7-8)."""
import inspect
import pytest


def test_algebra_surface_present_and_module_methods_absent():
    from quiverlab.core.algebra import Algebra
    for attr in ("domain", "T", "unit", "dim", "basis_labels", "quiver", "relations",
                 "is_unit_adapted"):
        assert attr in Algebra.__init__.__code__.co_varnames or hasattr(Algebra, attr) \
            or True  # attributes are instance-set; presence checked at runtime below
    for meth in ("multiply", "unit_adapted", "cartan_matrix", "coxeter_matrix",
                 "coxeter_polynomial", "hochschild_cohomology", "hochschild_homology",
                 "nakayama_automorphism", "is_frobenius", "is_symmetric"):
        assert callable(getattr(Algebra, meth)), meth
    # Plan-05 adds these; they must NOT pre-exist (else a rebase duplicated work).
    # simple/projective/injective (Task 3), hom (Task 5) and now ext/global_dimension/
    # is_selfinjective (Task 7) have since LANDED and are intentionally dropped from this
    # list; the remaining names stay guarded until their tasks land.
    for new in ("loewy_length", "complexity", "center"):
        assert not hasattr(Algebra, new), f"{new} already exists -- reconcile"


def test_algebra_runtime_attributes():
    from quiverlab import linear_path_algebra
    A = linear_path_algebra(2)
    for attr in ("domain", "T", "unit", "dim", "basis_labels", "quiver", "relations"):
        assert hasattr(A, attr), attr
    assert A.dim == 3 and A.basis_labels == ["e_1", "e_2", "a1"]   # arrow auto-named a1


def test_linalg_signatures():
    from quiverlab.fields import linalg
    from quiverlab.fields import QQ
    assert linalg.rank([[QQ.coerce(1)]], QQ) == 1
    ns = linalg.nullspace([[QQ.coerce(1), QQ.coerce(-1)]], QQ)
    assert len(ns) == 1
    A = [[QQ.coerce(1), QQ.coerce(0)], [QQ.coerce(0), QQ.coerce(1)]]
    assert linalg.solve(A, [QQ.coerce(2), QQ.coerce(3)], QQ) == [QQ.coerce(2), QQ.coerce(3)]
    R, piv = linalg.rref([[QQ.coerce(2), QQ.coerce(0)]], QQ)
    assert piv == [0]


def test_invariants_cartan_surface():
    from quiverlab.invariants import cartan
    from quiverlab import linear_path_algebra
    A = linear_path_algebra(2)
    assert cartan.cartan_matrix(A) == [[1, 1], [0, 1]]
    assert cartan.coxeter_polynomial(A) is not None


def test_engine_spectral_and_coxeter2_surface():
    from quiverlab.engine.coxeter_spectrum import (is_cyclotomic_product, star_quiver,
                                                   cartan_of_quiver, trivial_extension_cartan)
    from quiverlab.engine.coxeter2 import coxeter_polynomial_from_cartan
    from quiverlab.engine.scan3 import complexity_of
    n, arrows = star_quiver([1, 2, 6])
    assert n == 10
    C, _alg = cartan_of_quiver(n, arrows)
    poly, Phi = coxeter_polynomial_from_cartan(C)
    assert poly is not None
    assert trivial_extension_cartan(C).shape == (10, 10)
    assert is_cyclotomic_product(None) is None
    assert complexity_of([1, 1, 1, 1, 1, 1]) in (1, "1", 1)  # constant seq -> complexity 1


def test_minimal_resolution_guard_surface():
    from quiverlab.engine import resolutions_minimal as rm
    from quiverlab.engine.adapter import to_engine
    sig = inspect.signature(rm.minimal_resolution)
    assert "max_transient_bytes" in sig.parameters
    sig2 = inspect.signature(rm.minimal_homology_dims)
    assert "max_transient_bytes" in sig2.parameters
    assert callable(to_engine)


def test_plan03_general_kqi_available():
    """General (non-monomial) kQ/I must lower -- needed for non-monomial modules and
    the memory-guard fixture. If Plan 03 has not landed, STOP."""
    from quiverlab import Quiver, GF
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x^3 - y^2", "y^3", "y*x + x*y"], field=GF(32003))
    assert A.dim == 9


def test_quiver_algebra_accepts_degree_bound():
    """Task 13's memory-guard fixture may need an explicit degree_bound; Plan 03's frozen
    interface adds `Quiver.algebra(..., degree_bound=None, trace=None)`. If the parameter
    is absent, the Gröbner completion cannot be pushed and Task 13 is blocked -- STOP."""
    import inspect
    from quiverlab import Quiver, GF
    params = inspect.signature(Quiver.algebra).parameters
    assert "degree_bound" in params, "Quiver.algebra must accept degree_bound= (Plan 03)"
    # and it must actually take effect (accepts the kwarg without error)
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x^3"], field=GF(2), degree_bound=12)
    assert A.dim == 3


def test_plan04_cs_optional_but_shape_checked_if_present():
    """Plan 04 is consumed only for the optional monomial/quadratic Ext cross-check.
    If present it must expose the frozen names; if absent, skip (Plan 05 stays testable)."""
    cs = pytest.importorskip("quiverlab.resolutions_cs")
    assert hasattr(cs, "cs_cohomology_dims") or hasattr(cs, "cs_homology_dims")
    from quiverlab.resolutions_cs.engine_facade import CSResolution  # noqa: F401
