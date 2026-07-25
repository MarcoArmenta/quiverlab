import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.groebner import build_reduction_system

# Reduction systems are the CS entry currency (there is NO A.reduction_system()).
@pytest.fixture
def kx2_rs():
    Q = Quiver([1], {"x": (1, 1)})
    return build_reduction_system(Q, ["x*x"], CC)


@pytest.fixture
def square_rs():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return build_reduction_system(Q, ["a*b - c*d"], CC)


@pytest.fixture
def qci_rs():
    def build(xi="2", field=CC):
        Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
        # NOTE: coefficient is written WITHOUT parentheses. The Plan-03 relation
        # grammar (quiverlab.combinat.relations) rejects parenthesized coefficients
        # -- "y*x - (2)*x*y" raises RelationError "unknown arrow '(2)'". The bare
        # form "y*x - 2*x*y" parses. Grammar extension is Plan-06 territory; do not
        # reintroduce the parens here.
        return build_reduction_system(Q, ["x*x", "y*y", f"y*x - {xi}*x*y"], field)
    return build


@pytest.fixture(scope="session")
def qci_gf5_diag4():
    """A shared quantum-CI / GF(5) `Comparison` with the comparison window forced to 0
    and the lifted diagonal Δ built to degree 4 (the ~28s cost).  SESSION-scoped so the
    cup and cap deep pins reuse the ONE heavy Δ_4 build instead of rebuilding it per
    test (Plan 21 review item -- recover the +2:18 battery cost of per-test Δ rebuilds).

    window=0 forces every cup/cap onto the native route, so "past-window" lands at low
    absolute degree.  The resolution is grown to max_degree 5 (Δ_4 needs S up to 4);
    callers must stay at total degree <= 4 so no test triggers a resolution REBUILD
    (which would invalidate the shared Δ cache).  Every consumer is read-mostly (append
    caches) -- safe to share."""
    from quiverlab.resolutions_cs.comparison import Comparison
    from quiverlab.resolutions_cs.diagonal import diagonal
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=GF(5))
    comp = Comparison(A, window=0)
    comp._ensure(4)                     # CS resolution to max_degree 5
    diagonal(comp._res, 4)              # the ~28s build; cached on comp._res._tensor_complex
    return comp
