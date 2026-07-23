"""Plan 17: the CS correction solve is canonicalized modulo its nullspace.

reduce_mod_nullspace(x, A, dom) is the unique element of x + Null(A) with zero
coordinates at every free (non-pivot) column of A's RREF.  solve() already
returns that representative (free variables set to 0), so canonicalization is a
NO-OP today -- the point is that it is now an explicit, tested guarantee instead
of a solver-convention accident (see the WARNING block that used to live in
test_battery_bank_oracle.py).  The adversarial test below proves the CS
differential bytes no longer depend on WHICH solution the solver returns."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace, reduce_mod_nullspace, solve


def _dom_cases():
    return (GF(5), QQ)


def test_coset_invariance_and_idempotence():
    """Every solution of A y = b canonicalizes to the SAME vector; applying the
    reduction twice equals applying it once."""
    for dom in _dom_cases():
        i = dom.coerce
        # rank-2 system with a 2-dim nullspace (4 unknowns)
        A = [[i(1), i(2), i(0), i(1)],
             [i(0), i(0), i(1), i(3)],
             [i(1), i(2), i(1), i(4)]]          # row3 = row1 + row2 (dependent)
        b = [i(1), i(2), i(3)]
        x0 = solve(A, b, dom)
        assert x0 is not None
        canon = reduce_mod_nullspace(x0, A, dom)
        assert canon == reduce_mod_nullspace(canon, A, dom)      # idempotent
        for v in nullspace(A, dom):
            shifted = [dom.add(a_, b_) for a_, b_ in zip(x0, v)]
            assert reduce_mod_nullspace(shifted, A, dom) == canon  # coset-invariant


def test_solver_output_is_already_canonical():
    """solve()'s free-variables-zero particular solution IS the canonical
    representative -- the no-op property that keeps every byte pin passing."""
    for dom in _dom_cases():
        i = dom.coerce
        A = [[i(1), i(2), i(0), i(1)],
             [i(0), i(0), i(1), i(3)]]
        b = [i(4), i(1)]
        x0 = solve(A, b, dom)
        assert reduce_mod_nullspace(x0, A, dom) == x0


def test_full_rank_is_untouched():
    """No nullspace -> the reduction returns the input unchanged."""
    dom = GF(7)
    i = dom.coerce
    A = [[i(1), i(1)], [i(0), i(1)]]
    b = [i(3), i(2)]
    x0 = solve(A, b, dom)
    assert reduce_mod_nullspace(x0, A, dom) == x0
