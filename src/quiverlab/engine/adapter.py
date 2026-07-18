"""Bridge between quiverlab's domain-generic Algebra and the hanlab engine.

The engine computes over prime fields F_p with numpy int64 structure
constants. `to_engine` converts a quiverlab Algebra over GF(p) (PrimeField);
every other domain is refused loudly (spec: the fast engine is a GF(p)
accelerator, the pure bar path serves all fields)."""
import numpy as np

from quiverlab.errors import DepthLimitError, FieldError
from quiverlab.fields.primefield import PrimeField

_ENGINE_GUARD_HINT = (
    "the fast engine still builds the exponential normalized-bar basis "
    "(m*(m-1)^n cells per cochain degree); it accelerates the rank, not the "
    "basis size -- raise max_cells only if you know what you are doing"
)


def _guard(m, sizes, kind, max_cells):
    """Refuse (loudly, like the bar oracle) when a differential matrix would be
    larger than max_cells. `sizes` is a list of (label, rows, cols)."""
    for label, rows, cols in sizes:
        if rows * cols > max_cells:
            raise DepthLimitError(
                f"fast engine {kind} {label}: differential matrix would have "
                f"{rows} x {cols} entries (> max_cells = {max_cells})",
                hint=_ENGINE_GUARD_HINT,
            )


def to_engine(A):
    dom = A.domain
    if not isinstance(dom, PrimeField):
        raise FieldError(
            f"the fast engine computes over prime fields only; this algebra is over {dom.name}",
            hint="use engine='bar' (any field), or construct the algebra over GF(p)",
        )
    from quiverlab.engine.hh_engine import Algebra as EngineAlgebra

    m = A.dim
    T = np.zeros((m, m, m), dtype=np.int64)
    for i in range(m):
        for j in range(m):
            vec = A.T[i][j]
            for t in range(m):
                T[i, j, t] = int(vec[t])
    unit = np.array([int(c) for c in A.unit], dtype=np.int64)
    return EngineAlgebra(m, T, unit)


def engine_cohomology_dims(A, top, max_cells=4_000_000):
    """HH^0..HH^top dimensions over GF(p) via the engine, as a plain list[int].

    Guards the exponential bar-basis size against max_cells with the same
    semantics as the pure bar oracle (coboundary d^n: C^n -> C^{n+1})."""
    from quiverlab.engine.scan3 import hochschild_cohomology_dims
    m = A.dim
    _guard(m, [(f"d^{n}", m * (m - 1) ** (n + 1), m * (m - 1) ** n)
               for n in range(top + 1)], "coboundary", max_cells)
    E = to_engine(A.unit_adapted())
    p = A.domain.p
    out = hochschild_cohomology_dims(E, top, primes=(p,))
    return [int(d) for d in out[p]]


def engine_homology_dims(A, top, max_cells=4_000_000):
    """HH_0..HH_top dimensions over GF(p) via the engine, as a plain list[int].

    Guards the exponential bar-basis size against max_cells with the same
    semantics as the pure bar oracle (boundary b_n: C_n -> C_{n-1})."""
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    m = A.dim
    _guard(m, [(f"b_{n}", m * (m - 1) ** (n - 1), m * (m - 1) ** n)
               for n in range(1, top + 2)], "boundary", max_cells)
    E = to_engine(A.unit_adapted())
    p = A.domain.p
    out = hochschild_homology_dims(E, top, primes=(p,))
    return [int(d) for d in out[p]]
