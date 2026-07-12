"""Starter algebra families (building blocks; the full catalog is Plan 06)."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError


def truncated_polynomial(n: int, field=None):
    """k[x]/(x^n) as the one-loop quiver with the monomial relation x^n, n >= 2."""
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 2):
        raise QuiverlabError(f"truncated_polynomial({n!r}): need an integer n >= 2",
                             hint="n = 1 is the ground field; use n >= 2")
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=[f"x^{n}"], field=field)


def linear_path_algebra(n: int, field=None):
    """The hereditary path algebra of 1 -> 2 -> ... -> n (dimension n(n+1)/2)."""
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
        raise QuiverlabError(f"linear_path_algebra({n!r}): need an integer n >= 1",
                             hint="e.g. linear_path_algebra(3)")
    arrows = {f"a{i}": (i, i + 1) for i in range(1, n)}
    Q = Quiver(vertices=list(range(1, n + 1)), arrows=arrows)
    return Q.algebra(relations=[], field=field)
