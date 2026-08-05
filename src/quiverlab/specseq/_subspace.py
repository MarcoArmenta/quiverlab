"""Exact subspace linear algebra on column-span bases over any Domain.

A subspace of ``V = dom^D`` is held as a ``list`` of coordinate columns (each a
length-``D`` list over ``dom``) spanning it -- the internal currency of the whole
spectral-sequence engine. Every operation is exact Gaussian elimination via
``fields.linalg`` / ``modules.linalg_mod`` (no floats; the ``src/`` AST gate scans
this file). Shared by ``filtered.py`` (filtration validation) and ``pages.py``
(the Weibel subquotient pages)."""
from quiverlab.fields import linalg as flinalg
from quiverlab.modules import linalg_mod as lm


def colmat(cols):
    """The matrix whose columns are ``cols`` (``[]`` for the empty span)."""
    return lm.cols_to_matrix(cols) if cols else []


def span_dim(cols, dom):
    """``dim span(cols)`` = rank of the column matrix (0 for the empty span)."""
    M = colmat(cols)
    return flinalg.rank(M, dom) if M else 0


def dim_sum(U, W, dom):
    """``dim(span(U) + span(W))`` -- the rank of the concatenated columns."""
    both = list(U) + list(W)
    return span_dim(both, dom)


def in_span(vec, basis_cols, dom):
    """Is ``vec`` in the column-span of ``basis_cols``? (Exact solve; the empty
    span contains only the zero vector.)"""
    if all(dom.is_zero(x) for x in vec):
        return True
    B = colmat(basis_cols)
    if not B:
        return False
    return flinalg.solve(B, list(vec), dom) is not None


def span_contains(sub_cols, super_cols, dom):
    """Is ``span(sub_cols) <= span(super_cols)``? (rank-based, exact)."""
    if not sub_cols:
        return True
    r_super = span_dim(super_cols, dom)
    r_both = dim_sum(super_cols, sub_cols, dom)
    return r_both == r_super


def reduce_to_independent(cols, dom):
    """A maximal linearly-independent subset of ``cols`` (a basis of their span),
    picked in column order (deterministic)."""
    if not cols:
        return []
    G = colmat(cols)
    piv = lm.column_space_pivots(G, dom)
    return [lm.col(G, j) for j in piv]


def image(dmat, cols, dom):
    """Columns of ``dmat @ [cols]`` (the image of ``span(cols)`` under ``dmat``)."""
    if not cols or not dmat or not dmat[0]:
        return []
    B = colmat(cols)
    P = lm.matmul(dmat, B, dom)
    return [lm.col(P, j) for j in range(len(P[0]))] if (P and P[0]) else []


def intersect(U, W, dom):
    """A basis of ``span(U) ∩ span(W)`` (both in the same ambient ``dom^D``), via
    the kernel of ``[U | -W]`` projected to the ``U`` side (the radtopsoc pattern)."""
    if not U or not W:
        return []
    BU, BW = colmat(U), colmat(W)
    D = len(BU)
    stacked = [BU[i] + [dom.neg(x) for x in BW[i]] for i in range(D)]
    ker = lm.kernel_columns(stacked, dom)
    ku = len(U)
    out = []
    for z in ker:
        a = z[:ku]
        out.append(lm.matvec(BU, a, dom))
    return reduce_to_independent(out, dom)


def preimage_selecting(A, target_cols, dom):
    """Basis of ``{y : A y in span(target_cols)}`` as ``y``-columns.

    Solves ``A y - T z = 0`` (``T = colmat(target_cols)``): the nullspace of the
    stacked ``[A | -T]`` projected to the ``y`` (first) coordinates, reduced to an
    independent set. ``A`` maps ``dom^c -> dom^D``; the empty target span forces
    ``A y = 0`` (i.e. ``ker A``)."""
    c = len(A[0]) if (A and A[0]) else 0
    if c == 0:
        return []
    T = colmat(target_cols)
    if not T:                                   # target span = {0}: y in ker A
        ker = lm.kernel_columns(A, dom)
        return ker
    D = len(A)
    stacked = [list(A[i]) + [dom.neg(x) for x in T[i]] for i in range(D)]
    ker = lm.kernel_columns(stacked, dom)
    out = [z[:c] for z in ker]
    return reduce_to_independent(out, dom)
