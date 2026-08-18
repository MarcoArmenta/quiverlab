"""Plan 54 Task E -- the P52 twisted-homology adapter (the cross-plan contract, spec s3).

P54's semisimple-nu route consumes P52's chain-level twisted Hochschild HOMOLOGY
``HH_*(A, {}_1A_sigma)`` -- the accessor
:func:`quiverlab.hochschild.coefficients.twisted_homology_classes` (per-degree
``{basis, boundary, classes}``). This module is the thin adapter so the rest of
P54 is P52-API-agnostic:

- :func:`twisted_homology_quotient` wraps the accessor into numpy int64 arrays over
  ``GF(p)`` (class reps as columns in the bar chain coordinates, per-degree boundary
  matrices, dims);
- :func:`bar_to_engine_perm` is the ONE reindexing between P52's bar chain ordering
  ``(s, J)`` (M-slot ``s`` outermost) and the engine's cochain ordering ``(J, s)``
  (M-slot innermost) -- pinned by a round-trip test (spec s3 point 5);
- :func:`nu_inverse` inverts the Nakayama matrix exactly over the Domain (the second
  twist direction the arbiter discriminates).

CONVENTION (mirrors P52's block VERBATIM, spec s3): the twist lives in the RIGHT
action and enters on the **FIRST face** (``m . a_1 = m . sigma(a_1)``); the last face
(left action ``a_n . m``) is untwisted; reduced-bar normalized; basis ordering
``(s, J)`` with ``s`` outermost. P54's twisted Connes ``B_sigma`` (``twisted_connes``)
is built to descend against exactly this boundary.
"""
import numpy as np

from quiverlab.errors import QuiverlabError
from quiverlab.fields.primefield import PrimeField


def nu_inverse(nu_matrix, domain):
    """The exact inverse of the Nakayama automorphism matrix over ``domain``
    (columns = images). ``nu`` is an automorphism, hence invertible; a singular
    matrix raises loudly (would be an internal invariant failure)."""
    from quiverlab.fields.linalg import solve
    dom = domain
    m = len(nu_matrix)
    N = [[dom.coerce(nu_matrix[i][j]) for j in range(m)] for i in range(m)]
    cols = []
    for j in range(m):
        e = [dom.one() if i == j else dom.zero() for i in range(m)]
        x = solve(N, e, dom)
        if x is None:
            raise QuiverlabError(
                "BV twisted route: the Nakayama automorphism matrix is singular "
                "(internal invariant) -- please report this algebra")
        cols.append(x)
    # cols[j] is column j of nu^{-1}; assemble row-major
    return [[cols[j][i] for j in range(m)] for i in range(m)]


def bar_to_engine_perm(E, n, m):
    """The permutation ``perm`` with ``perm[e] = b`` mapping the engine cochain
    index ``e`` (basis ``(J, s)``) to the P52 bar chain index ``b`` (basis
    ``(s, J)``) of the SAME degenerate-free tensor. Round-trip pinned in the tests."""
    from quiverlab.engine.scan3 import cochain_basis
    from quiverlab.hochschild import bar
    eb = list(cochain_basis(E, n))                     # [(J, s), ...]
    bb = bar._cochain_basis(m, n, m)                   # [(s, J), ...]
    idx = {b: i for i, b in enumerate(bb)}
    return [idx[(s, J)] for (J, s) in eb]


def twisted_homology_quotient(A, sigma_matrix, top, max_cells=4_000_000):
    """P52 twisted Hochschild homology ``HH_*(A, {}_1A_sigma)`` as GF(p) numpy data.

    ``sigma_matrix`` is an m x m Domain-int matrix (columns = images) IN THE
    UNIT-ADAPTED basis (``A`` must already be unit-adapted). Returns a list of
    ``top + 2`` per-degree dicts: indices ``0..top`` carry ``reps`` (class-rep
    columns as int64 vectors over the bar chain basis), ``boundary`` (the int64
    ``b_n`` matrix, ``None`` at ``n = 0``) and ``dim`` (# classes); the trailing
    index ``top + 1`` carries ONLY its ``boundary`` (``b_{top+1}``, needed by the
    class-level Connes projection at degree ``top - 1``) with ``reps = []`` --
    computing the degree-``top+1`` homology CLASSES (an extra exact nullspace) is
    deliberately avoided. GF(p) only."""
    from quiverlab.hochschild import bar
    from quiverlab.hochschild.coefficients import Bimodule, twisted_homology_classes
    dom = A.domain
    if not isinstance(dom, PrimeField):
        raise QuiverlabError("BV twisted route is GF(p) only")
    p = dom.p
    if not A.is_unit_adapted:
        raise QuiverlabError(
            "twisted_homology_quotient expects a unit-adapted algebra "
            "(internal invariant)")
    psi = [[int(sigma_matrix[i][j]) for j in range(A.dim)] for i in range(A.dim)]
    M = Bimodule.twisted(A, psi=psi, name="twisted _1A_sigma")
    data = twisted_homology_classes(A, M, top, max_cells=max_cells)

    def to_np(bn):
        return (np.array([[int(x) % p for x in row] for row in bn], dtype=np.int64)
                if bn else None)

    out = []
    for n in range(top + 1):
        d = data[n]
        reps = [np.array([int(x) % p for x in cls], dtype=np.int64)
                for cls in d["classes"]]
        out.append({"reps": reps, "boundary": to_np(d["boundary"]), "dim": len(reps)})
    # b_{top+1} only (no degree-(top+1) nullspace): the matrix from the bar complex.
    bt1, _, _ = bar.boundary_matrix(A, top + 1, max_cells, coefficients=M)
    out.append({"reps": [], "boundary": to_np(bt1), "dim": 0})
    return out
