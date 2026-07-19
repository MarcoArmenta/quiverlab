"""Exterior algebra Lambda(k^n) = k<x_1..x_n>/(x_i^2, x_i x_j + x_j x_i). General route."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError


def ExteriorAlgebra(n, field=None):
    if not (isinstance(n, int) and n >= 1):
        raise QuiverlabError(f"ExteriorAlgebra({n!r}): need integer n >= 1", hint="e.g. 3")
    names = [f"x{i}" for i in range(1, n + 1)]
    Q = Quiver([1], {a: (1, 1) for a in names})
    rels = [f"{a}^2" for a in names]
    rels += [f"{names[i]}*{names[j]} + {names[j]}*{names[i]}"
             for i in range(n) for j in range(i + 1, n)]
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("quantum_ci", "chouhy_solotar")
    return A
