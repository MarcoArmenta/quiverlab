"""Serialize a quiverlab right-module into QPA's representation form (Plan 23,
item 4 oracle). QPA's ``RightModuleOverPathAlgebra(A, dimvec, [[arrow, mat], ...])``
wants, per arrow ``a: s -> t``, a ``dim_s x dim_t`` matrix acting on ROW vectors
(``m -> m * mat``). Our modules use the column/anti-homomorphism convention
(``m*b = action[b] @ m``, ``m`` a column), so the QPA matrix is the transpose of
the target-in-source-basis block.

The same graded form is the natural serialization for the future no-code module
schema (Tier 1b item 3: ``{dims: {v: n_v}, maps: {arrow: [[...]]}}``)."""
from quiverlab.modules import linalg_mod as lm


def graded_form(M):
    """(dimvec_list, arrow_matrices) in vertex order / QPA's row convention.

    dimvec_list[k] = dim (M e_{v_k}); arrow_matrices[a] = the dim_s x dim_t matrix
    (row convention) for each arrow a with a NONZERO action between positive-dim
    vertices (arrows omitted default to the zero map in QPA)."""
    A = M.algebra
    dom = A.domain
    verts = list(A.quiver.vertices)
    Ev, dv = {}, {}
    for v in verts:
        Pv = M.vertex_projection(v)
        piv = lm.column_space_pivots(Pv, dom)
        Ev[v] = [lm.col(Pv, j) for j in piv]        # a k-basis of M e_v
        dv[v] = len(piv)
    dimvec = [dv[v] for v in verts]
    arrows = {}
    for a, (s, t) in A.quiver.arrows.items():
        if dv[s] == 0 or dv[t] == 0:
            continue
        Es = lm.cols_to_matrix(Ev[s])               # M.dim x dim_s
        Et = lm.cols_to_matrix(Ev[t])               # M.dim x dim_t
        AEs = lm.matmul(M.action[a], Es, dom)       # a applied to each source basis vector
        Ca = lm.solve_columns(Et, AEs, dom)         # coords in the target basis (dim_t x dim_s)
        assert Ca is not None, (
            f"graded_form: arrow {a!r} image escapes the target vertex space "
            "(module is not vertex-graded)")
        Ca_mat = lm.cols_to_matrix(Ca)              # dim_t x dim_s
        if all(dom.is_zero(x) for row in Ca_mat for x in row):
            continue
        arrows[a] = lm.transpose(Ca_mat)            # dim_s x dim_t (QPA row convention)
    return dimvec, arrows
