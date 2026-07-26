"""Deep companion to ``tests/invariants/test_identity_oracles.py`` (Plan 29,
Part 2 item 1): the one piece measurement forced out of the fast bucket.

Happel's trace identity (Happel, Linear Algebra Appl. 258 (1997) 169-177,
registry key ``happel_trace``) on A_4. A_4 has dimension 10, so its degree-2 bar
differential is 7290 x 810 = 5.9M cells -- past the default 4M guard. The fast
file checks A_4 at top=1 (its Euler characteristic is complete by gl.dim 1);
here we raise the guard and reach top=2 so the vanishing tail HH^2 = 0 is
directly OBSERVED rather than inferred, closing the identity end-to-end past the
exponential bar blow-up.

tests/engine/ auto-assigns to the deep bucket (tests/conftest.py).
"""
from quiverlab import GF, Quiver

F = GF(32003)


def _trace(A):
    M = A.coxeter_matrix()
    return sum(M[i][i] for i in range(len(M)))


def _euler(dims):
    return sum((-1) ** i * d for i, d in enumerate(dims))


def test_happel_trace_A4_top2_observes_vanishing_tail():
    A4 = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=[], field=F)
    # d^2 is 7290 x 810 = 5.9M cells; raise the guard just enough to reach it.
    hh = A4.hochschild_cohomology(2, max_cells=8_000_000).dims
    assert hh == [1, 0, 0]
    assert hh[2] == 0                        # vanishing tail OBSERVED at top=2
    assert _trace(A4) == -_euler(hh)         # Happel's identity, tail included
