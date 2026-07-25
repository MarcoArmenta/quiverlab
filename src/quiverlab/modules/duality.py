"""Duality D, transpose Tr, and the Auslander-Reiten translates tau / tau^-
(Plan 23, Tier 1b item 1), on the A^op engine (``modules/opposite.py``).

  D  = Hom_k(-, k):   right A-mod  ->  right A^op-mod   (contravariant, D.D = id)
  Tr = coker(Hom_A(P_0,A) -> Hom_A(P_1,A)) from the minimal presentation
       P_1 -> P_0 -> M -> 0:  right A-mod -> right A^op-mod
  tau  M = D (Tr M)          tau^- M = Tr (D M)          (both right A-modules)

Right-module (anti-homomorphism) convention throughout: m*b = action[b] @ m,
action[x*y] = action[y] @ action[x].
"""
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import _other_side
from quiverlab.modules.opposite import opposite_algebra, reverse_label


def _zero_module(A, side="right"):
    from quiverlab.modules.module import Module
    action = {lab: lm.zeros(0, 0, A.domain) for lab in A.basis_labels}
    return Module(A, 0, action, name="0", side=side)


def _offsets(dims):
    offs, o = [], 0
    for d in dims:
        offs.append(o)
        o += d
    return offs


def dualize(M):
    """D M = Hom_k(M, k): transpose every action matrix and reverse every basis label.
    Preserves dimension vectors. Contravariant and D.D = id. At the representation
    level the algebra flips M.algebra -> M.algebra^op; the CATEGORICAL side flips too
    (Plan 24), so D exchanges the two sides over the SAME base algebra -- D of a right
    A-module is a LEFT A-module and vice versa."""
    from quiverlab.modules.module import Module
    Rop = opposite_algebra(M.algebra)
    action = {reverse_label(lab): lm.transpose(mat) for lab, mat in M.action.items()}
    return Module(Rop, M.dim, action, name=f"D({M.name})", side=_other_side(M.side))


def transpose_module(M):
    """Tr M as a right (M.algebra)^op-module, from the minimal projective
    presentation P_1 -> P_0 -> M -> 0 (minimality => no spurious projective
    summands). Tr M = coker(d_1^*), d_1^* the corner-transpose of d_1."""
    from quiverlab.modules.builders import projective
    from quiverlab.modules.radtopsoc import quotient
    from quiverlab.modules.resolution import _direct_sum, minimal_resolution

    R = M.algebra
    Rop = opposite_algebra(R)
    dom = R.domain

    out_side = _other_side(M.side)      # Tr is contravariant: it flips the side
    terms, dmats = minimal_resolution(M, 1)
    v_list = terms[0].vertices          # P_0 summand vertices
    w_list = terms[1].vertices          # P_1 summand vertices
    if not w_list:                      # M projective (or zero) => Tr M = 0
        return _zero_module(Rop, side=out_side)
    d1 = dmats[1]                       # P_1 -> P_0 in k-bases (P0.dim x P1.dim)

    # Target module N = P_1^{op} = (+)_j projective(Rop, w_j). It carries Tr's flipped
    # side; the cokernel below inherits it.
    Sop = [projective(Rop, w) for w in w_list]
    N, off1op = _direct_sum(Sop, name=f"Tr({M.name})_cover", side=out_side)
    posmap = [{lab: k for k, lab in enumerate(s._pv_basis_labels)} for s in Sop]
    off1op_start = [s for (s, _d) in off1op]

    # A-side P_0 / P_1 summands to slice corner elements out of d1.
    S0 = [projective(R, v) for v in v_list]
    S1 = [projective(R, w) for w in w_list]
    off0 = _offsets([s.dim for s in S0])
    off1 = _offsets([s.dim for s in S1])

    # Source summands over Rop (only their basis labels are needed).
    S0op = [projective(Rop, v) for v in v_list]

    cols_dstar = []
    for i, S0i in enumerate(S0):
        # h_i in N = image of the source generator g_i = e_{v_i}: sum_j ybar_{ij}.
        h_i = [dom.zero()] * N.dim
        for j in range(len(w_list)):
            gen_col = off1[j]                       # generator (e_{w_j}) column of P_1 summ. j
            for k in range(S0i.dim):
                val = d1[off0[i] + k][gen_col]
                if dom.is_zero(val):
                    continue
                rlab = reverse_label(S0i._pv_basis_labels[k])   # A^op-label at w_j
                loc = posmap[j].get(rlab)
                assert loc is not None, (
                    "transpose_module: corner element off the vertex grading "
                    "(module-map violation)")
                h_i[off1op_start[j] + loc] = val
        # d_1^* columns for source summand i: (g_i .^op p) |-> N.action[p] @ h_i
        for p in S0op[i]._pv_basis_labels:
            cols_dstar.append(lm.matvec(N.action[p], h_i, dom))

    d1star = lm.cols_to_matrix(cols_dstar) if cols_dstar else lm.zeros(N.dim, 0, dom)
    piv = lm.column_space_pivots(d1star, dom) if (d1star and d1star[0]) else []
    image_cols = [lm.col(d1star, j) for j in piv]
    return quotient(N, image_cols, name=f"Tr({M.name})")


def tau(M):
    """AR translate tau M = D(Tr M), a right (M.algebra)-module. tau(projective)=0."""
    out = dualize(transpose_module(M))
    out.name = f"tau({M.name})"
    return out


def tau_minus(M):
    """inverse AR translate tau^- M = Tr(D M), a right (M.algebra)-module.
    tau^-(injective)=0."""
    out = transpose_module(dualize(M))
    out.name = f"tau^-({M.name})"
    return out
