"""Radical, top and socle of a right A-module (spec §3.6).

rad M = M*(rad A) = sum over arrows of image(action[arrow])  (a submodule);
top M = M / rad M  (semisimple);
soc M = {m : m*(rad A) = 0} = intersection over arrows of ker(action[arrow]).
Submodule/quotient modules inherit the action by restriction/co-restriction, solved
exactly over the Domain via fields.linalg (solve_columns)."""
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module


def _rad_image_cols(M):
    """Column-vectors spanning M*(rad A) = sum_arrows image(action[arrow])."""
    dom = M.domain
    gens = []
    for aname in M.algebra.quiver.arrows:
        Aa = M.action[aname]
        for j in range(M.dim):
            gens.append(lm.col(Aa, j))
    if not gens:
        return []
    G = lm.cols_to_matrix(gens)
    piv = lm.column_space_pivots(G, dom)
    return [lm.col(G, j) for j in piv]


def submodule(M, basis_cols, name="sub"):
    """The submodule of M spanned by basis_cols (assumed A-stable). Its action[b] is
    the coordinates of action[b] applied to each basis column, expressed back in
    basis_cols (solved over the Domain)."""
    dom = M.domain
    B = lm.cols_to_matrix(basis_cols) if basis_cols else lm.zeros(M.dim, 0, dom)
    n = len(basis_cols)
    action = {}
    for label, Ab in M.action.items():
        if n == 0:
            action[label] = lm.zeros(0, 0, dom)
            continue
        images = [lm.matvec(Ab, c, dom) for c in basis_cols]   # b acts on each generator
        V = lm.cols_to_matrix(images)
        coeffs = lm.solve_columns(B, V, dom)                   # express in basis_cols
        assert coeffs is not None, f"submodule not A-stable under {label}"
        action[label] = lm.cols_to_matrix(coeffs)
    return Module(M.algebra, n, action, name=name)


def quotient(M, sub_cols, name="quot"):
    """The quotient module M / <sub_cols>. Pick coset representatives = a basis of M
    completing sub_cols; action[b] is read on the representatives modulo the submodule."""
    dom = M.domain
    sub_piveb = sub_cols
    ident = lm.identity(M.dim, dom)
    std = [lm.col(ident, j) for j in range(M.dim)]
    # representatives: standard vectors independent modulo the submodule
    rep_idx = lm.independent_modulo(std, sub_piveb, dom)
    reps = [std[i] for i in rep_idx]
    n = len(reps)
    # basis of the WHOLE space: submodule columns then representatives
    whole = [list(c) for c in sub_cols] + reps
    W = lm.cols_to_matrix(whole)
    action = {}
    s = len(sub_cols)
    for label, Ab in M.action.items():
        cols = []
        for r in reps:
            img = lm.matvec(Ab, r, dom)
            x = lm.solve_columns(W, lm.cols_to_matrix([img]), dom)[0]
            cols.append(x[s:])                 # drop the submodule part -> class in quotient
        action[label] = lm.cols_to_matrix(cols) if cols else lm.zeros(n, n, dom)
    return Module(M.algebra, n, action, name=name)


def radical(M):
    return submodule(M, _rad_image_cols(M), name=f"rad {M.name}")


def top(M):
    return quotient(M, _rad_image_cols(M), name=f"top {M.name}")


def socle(M):
    """soc M = intersection over arrows of ker(action[arrow])."""
    dom = M.domain
    inter = None
    arrows = list(M.algebra.quiver.arrows)
    if not arrows:                              # semisimple: soc = M
        return submodule(M, [lm.col(lm.identity(M.dim, dom), j) for j in range(M.dim)],
                         name=f"soc {M.name}")
    for aname in arrows:
        ker = lm.kernel_columns(M.action[aname], dom)
        if inter is None:
            inter = ker
        else:
            inter = _intersect(inter, ker, dom)
    return submodule(M, inter or [], name=f"soc {M.name}")


def _intersect(cols1, cols2, dom):
    """Basis of the intersection of two subspaces given by column bases (via the
    kernel of [B1 | -B2] projected to the B1 side)."""
    if not cols1 or not cols2:
        return []
    B1, B2 = lm.cols_to_matrix(cols1), lm.cols_to_matrix(cols2)
    stacked = [row1 + [dom.neg(x) for x in row2] for row1, row2 in zip(B1, B2)]
    ker = lm.kernel_columns(stacked, dom)
    k1 = len(cols1)
    out = []
    for z in ker:
        coeff = z[:k1]
        vec = lm.matvec(B1, coeff, dom)
        out.append(vec)
    # reduce to an independent set
    if not out:
        return []
    G = lm.cols_to_matrix(out)
    piv = lm.column_space_pivots(G, dom)
    return [lm.col(G, j) for j in piv]
