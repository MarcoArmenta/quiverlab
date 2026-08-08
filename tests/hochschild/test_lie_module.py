"""HH^* as a graded Lie module over HH^1 (Plan 71 / R12; Gerstenhaber deg-1 = Lie derivative).
Field-general on the normalized bar cochain complex; char-0 weights/decomposition behind P70's gate."""
import pytest
from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.families import truncated_polynomial
from quiverlab.hochschild.lie_module import lie_module_action, lie_derivative_on_hh

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _rad2loops(field):
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "x*y", "y*x", "y*y"], field=field)


def _kron(field):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=field)


# ---------------------------------------------------------------------------
# Task 1: field-general action rho_n(D), dims, module axiom, inner-acts-zero
# ---------------------------------------------------------------------------
@xeng
@pytest.mark.parametrize("n,dims", [(2, [2, 1, 1]), (3, [3, 2, 2])])   # k[x]/x^n
def test_hh_dims_match_bar_engine(n, dims):
    A = truncated_polynomial(n, field=QQ)
    L = lie_module_action(A, top=2)
    assert L.hh_dims == dims == A.hochschild_cohomology(top=2).dims


@xeng
def test_kronecker_hh_dims():
    A = _kron(QQ)
    L = lie_module_action(A, top=2)
    assert L.hh_dims == [1, 3, 0] and L.hh1_dim == 3


@selfcert
@pytest.mark.parametrize("build", ["trunc3", "rad2loops", "kron"])
def test_module_axiom_and_inner_zero(build):
    A = ({"trunc3": lambda: truncated_polynomial(3, field=QQ),
          "rad2loops": lambda: _rad2loops(QQ),
          "kron": lambda: _kron(QQ)}[build])()
    L = lie_module_action(A, top=2)
    assert L.module_axiom_ok is True and L.inner_acts_zero is True


@selfcert
def test_lie_derivative_primitive_shape():
    """lie_derivative_on_hh returns the d_n x d_n action matrix on HH^n."""
    A = truncated_polynomial(3, field=QQ)
    from quiverlab.invariants.hh1_lie import derivations
    B = A.unit_adapted()
    Der = derivations(B)
    R = lie_derivative_on_hh(A, Der[0], 1)
    assert len(R) == 2 and all(len(row) == 2 for row in R)      # dim HH^1(k[x]/x^3) = 2


@selfcert
def test_basis_provenance_tag():
    A = truncated_polynomial(2, field=QQ)
    L = lie_module_action(A, top=2)
    assert L.basis == "der_inn/bar"
