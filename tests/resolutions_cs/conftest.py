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
