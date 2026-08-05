"""End_A(M) as a structure-constant Algebra (Plan 37 / C1).

Basis = ``hom_basis(M, M)``; the product of basis elements ``b_i * b_j`` is
the COMPOSITE ``b_i o b_j`` (apply ``b_j`` first -- function-composition
order, i.e. ``b_j.then(b_i)``). ARBITRATED (devil's-advocate round,
2026-08-05): on the non-self-opposite source quiver ``1->2, 1->3`` the
corner dimensions of ``End((+) P_v)`` equal the CARTAN MATRIX of ``A``
under this order (``End(A_A) ~ A``) and its transpose under the path
order ``b_i.then(b_j)`` (``~ A^op``); the original kA2 regular-module
oracle was blind to the difference (kA2 is self-opposite). Pinned by
``regular_corner_dims`` and ``test_end_of_regular_is_A_not_Aop``. The
unit is ``id_M``.

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
    T, unit, basis, _B = _structure_constants(M)
    dom = M.domain
    d = M.dim
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


def _structure_constants(M):
    """(T, unit, basis, B) of End_A(M) in the hom_basis, product
    ``b_i * b_j = b_i o b_j`` (function composition = ``b_j.then(b_i)`` --
    see the module docstring for the arbitration)."""
    basis = hom_basis(M, M)
    dom = M.domain
    d = M.dim
    B = lm.cols_to_matrix([_vec(f.matrix, d) for f in basis])   # (d*d) x r
    T = []
    for bi in basis:
        row = []
        for bj in basis:
            prod = bj.then(bi)                    # b_i o b_j: apply b_j first
            x = solve(B, _vec(prod.matrix, d), dom)
            if x is None:
                raise QuiverlabError(
                    "End(M) product left the Hom basis -- hom_space is "
                    "inconsistent (bug)")
            row.append(x)
        T.append(row)
    unit = solve(B, _vec(M.identity_hom().matrix, d), dom)
    return T, unit, basis, B


def regular_corner_dims(A):
    """dim(e_v . End(A_A) . e_w) for the summand projectors e_v of the
    regular module (+)_v P_v, computed THROUGH the structure constants --
    the sided oracle: equals ``cartan_matrix(A)`` iff End(A_A) ~ A."""
    from quiverlab.modules.morphism import direct_sum
    verts = list(A.quiver.vertices)
    D, incls, projs = direct_sum(*[A.projective(v) for v in verts])
    T, unit, basis, B = _structure_constants(D)
    dom = D.domain
    d = D.dim
    E0 = Algebra.from_structure_constants(T, list(unit), field=dom, check=False)
    eps = [solve(B, _vec(projs[i].then(incls[i]).matrix, d), dom)
           for i in range(len(verts))]
    n = len(basis)
    out = []
    for i in range(len(verts)):
        row = []
        for j in range(len(verts)):
            vecs = []
            for k in range(n):
                bk = [dom.one() if m == k else dom.zero() for m in range(n)]
                vecs.append(E0.multiply(eps[i], E0.multiply(bk, eps[j])))
            row.append(lm.mat_rank(vecs, dom))
        out.append(row)
    return out
