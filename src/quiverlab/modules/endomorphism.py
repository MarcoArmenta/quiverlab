"""End_A(M) as a structure-constant Algebra (Plan 37 / C1).

Basis = ``hom_basis(M, M)``; the product of basis elements ``b_i * b_j`` is
``b_i.then(b_j)`` (left-to-right, the house path convention -- this makes
``End(A_A) ~ A`` as a k-algebra, and the regular-module oracle
``E.loewy_length() == A.loewy_length()`` accepts this order: both ``dim`` and
Loewy length are opposite-invariant, so the oracle passes with either order and
we keep the house left-to-right ``.then``). The unit is ``id_M``.

The returned :class:`~quiverlab.core.algebra.Algebra` is presentation-less
(no quiver): arithmetic, ``center``, ``loewy_length`` and decompose-style analysis
work; ``.simple`` / ``.projective`` need a quiver presentation and refuse loudly.

Loewy length / center read the Jacobson radical of ``End(M)``. That radical is
computed EXACTLY as the radical of the trace form ``tr_M(H_i H_j)``
(``rad(End) = End^{perp}`` -- the Dickson / Cohen--Ivanyos--Wales identity, valid
when ``char k = 0`` or ``char k > dim M``, the same bound
:mod:`quiverlab.modules.decompose` relies on) and the algebra is returned in a
radical-adapted basis (the radical basis vectors labeled non-``e_`` so the
label-heuristic Loewy walk reads the true radical). When ``char <= dim M`` the
trace-form radical is unreliable, so we do NOT fabricate labels: the algebra is
returned presentation-less and ``loewy_length`` is honestly unavailable there
rather than risk a wrong nilpotency index. Float-free.
"""
from __future__ import annotations

from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import solve
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.morphism import hom_basis


def _vec(mat, d):
    """Column-major vectorization of a d x d matrix."""
    return [mat[i][j] for j in range(d) for i in range(d)]


def end_algebra(M):
    """``End_A(M)`` as a structure-constant :class:`Algebra` (Plan 37 / C1)."""
    if M.dim == 0:
        raise QuiverlabError(
            "End of the zero module is the zero ring -- not an Algebra with 1 here; "
            "refused.")
    basis = hom_basis(M, M)
    dom = M.domain
    d = M.dim
    B = lm.cols_to_matrix([_vec(f.matrix, d) for f in basis])   # (d*d) x r
    T = []
    for bi in basis:
        row = []
        for bj in basis:
            prod = bi.then(bj)                                   # b_i then b_j
            x = solve(B, _vec(prod.matrix, d), dom)
            if x is None:
                raise QuiverlabError(
                    "End(M) product left the Hom basis -- hom_space is inconsistent "
                    "(bug)")
            row.append(x)
        T.append(row)
    unit = solve(B, _vec(M.identity_hom().matrix, d), dom)
    E = Algebra.from_structure_constants(T, list(unit), field=dom, check=True)
    relabeled = _with_radical_labels(E, [f.matrix for f in basis], dom, d)
    return relabeled if relabeled is not None else E


def _with_radical_labels(E, hmats, dom, d):
    """An isomorphic copy of ``E`` in a radical-adapted basis (the last ``dim rad E``
    basis vectors span ``rad(E)``, labeled non-``e_``) so the Loewy/center label
    heuristic reads the true radical -- or ``None`` when ``rad(E)`` is not rigorously
    computable over this field (``char p <= dim M``: the trace-form radical is
    unreliable, e.g. ``k[x]/(x^p)``). rad(E) = radical of the trace form
    ``tr_M(H_i H_j)`` (Dickson / Cohen--Ivanyos--Wales, ``char 0`` or ``char > dim M``)."""
    char = dom.characteristic
    n = d
    if not (char == 0 or char > n):
        return None
    r = len(hmats)
    if r <= 1:
        E.basis_labels = ["e_0"] if r == 1 else []      # E = k*id: rad = 0
        return E
    # trace-form Gram matrix T[i][j] = tr_M(H_i H_j)
    T = lm.zeros(r, r, dom)
    for i in range(r):
        for j in range(r):
            prod = lm.matmul(hmats[i], hmats[j], dom)
            s = dom.zero()
            for t in range(n):
                s = dom.add(s, prod[t][t])
            T[i][j] = s
    rad_coords = lm.kernel_columns(T, dom)               # nullspace = rad(E) coordinates
    dim_rad = len(rad_coords)
    if dim_rad == 0:
        E.basis_labels = [f"e_{i}" for i in range(r)]    # semisimple: rad = 0
        return E
    ident = lm.identity(r, dom)
    std = [lm.col(ident, j) for j in range(r)]
    comp_idx = lm.independent_modulo(std, rad_coords, dom)   # complement of rad
    comp = [std[i] for i in comp_idx]
    P = lm.cols_to_matrix(comp + [list(c) for c in rad_coords])
    E2 = E.change_of_basis(P)                            # complement first, radical last
    E2.basis_labels = ([f"e_{i}" for i in range(len(comp))]
                       + [f"r_{i}" for i in range(dim_rad)])
    return E2
