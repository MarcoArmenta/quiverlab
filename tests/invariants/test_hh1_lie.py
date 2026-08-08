"""HH^1 as a Lie algebra = Der(A)/Inn(A) (Plan 70 / R11; RSS 1903.12145; Strametz 2006).
Field-general via the algebra's structure constants; char-0 classification behind a hard gate."""
import pytest

from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.families import truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.hh1_lie import (
    derivations, inner_derivations, hh1_lie_structure,
)

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _kron(field=QQ):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=field)   # kK2, dim 4


# -------------------------------------------------------------------- Task 1
@xeng
@pytest.mark.parametrize("n,expect", [(2, 1), (3, 2), (4, 3), (5, 4)])   # dim Der(k[x]/x^n)=n-1 (QQ)
def test_dim_hh1_matches_bar_engine_truncpoly(n, expect):
    A = truncated_polynomial(n, field=QQ)
    L = hh1_lie_structure(A)
    assert L.dim == expect
    assert L.dim == A.hochschild_cohomology(top=2).dims[1]               # cross-engine anchor
    assert L.dim_inn == 0                                                # commutative


@xeng
def test_dim_hh1_kronecker():
    A = _kron()
    L = hh1_lie_structure(A)
    assert L.dim == 3 == A.hochschild_cohomology(top=2).dims[1]


@selfcert
def test_inn_is_ideal_of_der():
    """[Der, Inn] subset Inn (so HH^1 = Der/Inn is a Lie algebra)."""
    A = _kron()
    L = hh1_lie_structure(A)
    assert L.dim_der - L.dim_inn == L.dim
    assert len(derivations(A)) == L.dim_der and len(inner_derivations(A)) == L.dim_inn


@selfcert
def test_der_inn_basis_provenance():
    A = _kron()
    L = hh1_lie_structure(A)
    assert L.basis == "der_inn"
    assert isinstance(L.constants, list) and len(L.constants) == L.dim


# -------------------------------------------------------------------- Task 2
@lit
@pytest.mark.parametrize("n,solv,nilp,ds", [
    (2, True, True,  [1, 0]),
    (3, True, False, [2, 1, 0]),
    (4, True, False, [3, 2, 0]),
    (5, True, False, [4, 3, 1, 0]),
])
def test_truncpoly_char0_solvable(n, solv, nilp, ds):
    A = truncated_polynomial(n, field=QQ)
    L = hh1_lie_structure(A)
    assert L.solvable is solv and L.nilpotent is nilp
    assert L.derived_series_dims == ds


@lit
@pytest.mark.parametrize("p,dim,solvable,perfect", [(2, 2, True, False), (3, 3, False, True),
                                                    (5, 5, False, True), (7, 7, False, True)])
def test_witt_char_p(p, dim, solvable, perfect):
    """n = char = p: HH^1 = W_1 (Jacobson-Witt), simple for p>=3 (perfect, not solvable),
    2-dim solvable for p=2. THE char-p side of the dichotomy."""
    A = truncated_polynomial(p, field=GF(p))
    L = hh1_lie_structure(A)                     # over char p: field-general block only
    assert L.dim == dim and L.solvable is solvable and L.perfect is perfect


@lit
@pytest.mark.parametrize("n,p,solvable", [(4, 2, True), (6, 2, True), (6, 3, False)])
def test_char_divides_n_but_not_prime(n, p, solvable):
    """Card refinement: char|n with n!=p does NOT force the simple Witt case."""
    A = truncated_polynomial(n, field=GF(p))
    assert hh1_lie_structure(A).solvable is solvable


@selfcert
def test_abelian_and_perfect_flags():
    assert hh1_lie_structure(truncated_polynomial(2, field=QQ)).abelian is True   # k[x]/x^2
    assert hh1_lie_structure(_kron(GF(3))).perfect is True                        # sl2 perfect


@selfcert
def test_is_solvable_nilpotent_thin_verdicts():
    from quiverlab.invariants.hh1_lie import is_solvable_hh1, is_nilpotent_hh1
    A = truncated_polynomial(3, field=QQ)
    assert is_solvable_hh1(A) is True and is_nilpotent_hh1(A) is False
    assert _kron().is_solvable_hh1() is False


@selfcert
def test_oversize_budget_refused_loudly():
    A = truncated_polynomial(3, field=QQ)
    with pytest.raises(QuiverlabError):
        hh1_lie_structure(A, budget=2)               # dim 3 > budget 2
