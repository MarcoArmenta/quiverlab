"""Item 6 (TODO P3): correctness of the Search II structural builders.

Previously these were only PRINTED in scan2.__main__, never asserted. The strongest
free check is that a trivial extension is a SYMMETRIC algebra, so HH^n == HH_n in every
degree -- a one-line oracle that catches sign/index bugs in the dual-bimodule action.
The one-point (triangular) extension is NOT self-injective and is pinned against a frozen
HH sequence; the module builders are checked against the module axioms.
"""
import numpy as np
import pytest

from quiverlab.engine.hh_engine import (
    truncated_polynomial,
    check_associative,
    hochschild_homology_dims,
)
from quiverlab.engine.scan2 import (
    trivial_extension, triangular_extension, module_simple, module_regular,
    module_semisimple, kA,
)

# hanlab __init__ alias, reproduced locally. hochschild_homology_dims routes through
# quiverlab.engine.resolutions (Task 9); cohomology_dims lives in scan3 (Task 4). The
# HH-computing tests below self-heal per-test once those land; the module-builder
# axiom tests run now.
homology_dims = hochschild_homology_dims
PRIME = 32003
P = PRIME


# ---- trivial extension is symmetric: HH^n == HH_n in every degree ----
@pytest.mark.parametrize("Bfn,N", [
    (lambda: truncated_polynomial(2), 6),
    (lambda: truncated_polynomial(3), 3),
    (lambda: kA(2), 3),
])
def test_trivial_extension_is_symmetric(Bfn, N):
    scan3 = pytest.importorskip("quiverlab.engine.scan3")   # cohomology_dims (Task 4)
    pytest.importorskip("quiverlab.engine.resolutions")     # homology backend (Task 9)
    cohomology_dims = scan3.hochschild_cohomology_dims
    TB = trivial_extension(Bfn(), "T(B)")
    ok, _ = check_associative(TB)
    assert ok
    assert homology_dims(TB, N)[P] == cohomology_dims(TB, N)[P]


def test_trivial_extension_of_kx2_matches_commutative_ci():
    pytest.importorskip("quiverlab.engine.resolutions")     # homology backend (Task 9)
    # T(k[x]/(x^2)) is the dim-4 symmetric algebra with the same HH as k[x,y]/(x^2,y^2).
    TB = trivial_extension(truncated_polynomial(2), "T(k[x]/(x^2))")
    assert homology_dims(TB, 6)[P] == [4, 4, 5, 6, 7, 8, 9]


# ---- one-point (triangular) extension: NOT self-injective, frozen oracle ----
def test_triangular_extension_frontier_homology():
    pytest.importorskip("quiverlab.engine.resolutions")     # homology backend (Task 9)
    # k[x]/(x^3)[k]: the Search II frontier (no 2-truncated cycle, non-self-injective).
    a3 = truncated_polynomial(3)
    acts, d = module_simple(3)
    T = triangular_extension(a3, acts, d, "k[x]/(x^3)[k]")
    assert T.m == 5
    ok, _ = check_associative(T)
    assert ok
    assert homology_dims(T, 5)[P] == [4, 2, 2, 2, 2, 2]


def test_triangular_extension_kx2_simple():
    pytest.importorskip("quiverlab.engine.resolutions")     # homology backend (Task 9)
    a2 = truncated_polynomial(2)
    acts, d = module_simple(2)
    T = triangular_extension(a2, acts, d, "k[x]/(x^2)[k]")
    assert T.m == 4
    assert homology_dims(T, 6)[P] == [3, 1, 1, 1, 1, 1, 1]


# ---- module builders satisfy the module axioms ----
def test_module_simple_unit_acts_as_identity():
    acts, dim = module_simple(3)
    assert dim == 1
    assert np.array_equal(acts[0], np.eye(1, dtype=np.int64))   # 1_A acts as id
    assert np.array_equal(acts[1], np.zeros((1, 1), dtype=np.int64))  # x acts as 0


def test_module_regular_is_left_multiplication():
    # regular module of k[x]/(x^3): action of x^p is the shift x^p . x^j = x^{p+j}.
    acts, dim = module_regular(3)
    assert dim == 3
    assert np.array_equal(acts[0], np.eye(3, dtype=np.int64))       # 1 acts as id
    # x acts as the nilpotent shift e_j -> e_{j+1}
    x = acts[1]
    assert x[1, 0] == 1 and x[2, 1] == 1 and x.sum() == 2
    # module axiom: action of x composed twice == action of x^2
    assert np.array_equal((x @ x) % P, acts[2] % P)


def test_module_semisimple_identity_action():
    acts, dim = module_semisimple(2, 3)
    assert dim == 3
    assert np.array_equal(acts[0], np.eye(3, dtype=np.int64))
    assert np.array_equal(acts[1], np.zeros((3, 3), dtype=np.int64))
