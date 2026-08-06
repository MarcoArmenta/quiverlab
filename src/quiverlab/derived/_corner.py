"""The corner-transpose ``d^*`` of a map between projective modules (Plan 43 shared
home). P41 did NOT factor this out of ``duality._presentation_transpose``; this is the
one shared implementation the derived AR translate (``tau.py``) consumes and, after a
byte-stability check, ``duality._presentation_transpose`` delegates to.

Given a module map ``d: (+)_j P_{to_verts[j]} -> (+)_i P_{from_verts[i]}`` between
projective right ``A``-modules (matrix ``d`` in the resolution layout: rows =
``(+)P_{from}`` = target, cols = ``(+)P_{to}`` = source, and within each source block
the canonical generator ``e_v`` is the FIRST column), the corner-transpose is
``d^* = Hom_A(d, A): Hom_A((+)P_from, A) -> Hom_A((+)P_to, A)``, i.e.
``(+)_i A e_{from_verts[i]} -> (+)_j A e_{to_verts[j]}`` -- a right ``A^op``-module map,
matrix rows = ``(+)Ae_to`` (target of ``d^*``), cols = ``(+)Ae_from`` (source of ``d^*``).

Reconstructed from the generator images alone (``d`` is a module map, so ``d`` is
determined by where it sends each source generator ``e_v``); the slicing is byte-for-byte
the ``d1star`` construction of ``duality._presentation_transpose`` (``duality.py``).
Float-free."""
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.builders import projective
from quiverlab.modules.duality import _offsets, _zero_module
from quiverlab.modules.module import _other_side
from quiverlab.modules.opposite import opposite_algebra, reverse_label
from quiverlab.modules.resolution import _direct_sum


def corner_transpose(d, from_verts, to_verts, A, side="right"):
    """Full corner-transpose data ``(Rop, N, M0op, d1star, out_side)`` (mirrors
    ``duality._presentation_transpose``'s return). ``N = (+)Ae_to`` (target of ``d^*``),
    ``M0op = (+)Ae_from`` (source of ``d^*``); ``coker(d1star) = Tr`` of the presented
    module. ``out_side`` is the contravariantly-flipped side of ``side``."""
    R = A
    Rop = opposite_algebra(R)
    dom = R.domain
    out_side = _other_side(side)

    v_list = list(from_verts)                 # P_0 / target-of-d summand vertices
    w_list = list(to_verts)                   # P_1 / source-of-d summand vertices

    # M0op = Hom(P_0, A) = (+)_i projective(Rop, v_i) -- the ordered column basis of d1star.
    S0op = [projective(Rop, v) for v in v_list]
    M0op, _ = (_direct_sum(S0op, name="corner_0", side=out_side) if S0op
               else (_zero_module(Rop, side=out_side), []))

    if not w_list:                            # projective source => Tr = 0, d1star empty
        N = _zero_module(Rop, side=out_side)
        return Rop, N, M0op, lm.zeros(0, M0op.dim, dom), out_side

    # Target module N = (+)_j projective(Rop, w_j).
    Sop = [projective(Rop, w) for w in w_list]
    N, off1op = _direct_sum(Sop, name="corner_cover", side=out_side)
    posmap = [{lab: k for k, lab in enumerate(s._pv_basis_labels)} for s in Sop]
    off1op_start = [s for (s, _d) in off1op]

    # A-side P_0 / P_1 summands to slice corner elements out of d.
    S0 = [projective(R, v) for v in v_list]
    off0 = _offsets([s.dim for s in S0])
    off1 = _offsets([projective(R, w).dim for w in w_list])

    cols_dstar = []
    for i, S0i in enumerate(S0):
        # h_i in N = image of the source generator g_i = e_{v_i}: sum_j ybar_{ij}.
        h_i = [dom.zero()] * N.dim
        for j in range(len(w_list)):
            gen_col = off1[j]                       # generator (e_{w_j}) column of P_1 summ. j
            for k in range(S0i.dim):
                val = d[off0[i] + k][gen_col]
                if dom.is_zero(val):
                    continue
                rlab = reverse_label(S0i._pv_basis_labels[k])   # A^op-label at w_j
                loc = posmap[j].get(rlab)
                assert loc is not None, (
                    "corner_transpose: corner element off the vertex grading "
                    "(module-map violation)")
                h_i[off1op_start[j] + loc] = val
        # d^* columns for source summand i: (g_i .^op p) |-> N.action[p] @ h_i
        for p in S0op[i]._pv_basis_labels:
            cols_dstar.append(lm.matvec(N.action[p], h_i, dom))

    d1star = lm.cols_to_matrix(cols_dstar) if cols_dstar else lm.zeros(N.dim, 0, dom)
    return Rop, N, M0op, d1star, out_side


def _corner_transpose_matrix(d, from_verts, to_verts, A, side="right"):
    """Just the corner-transpose matrix ``d^*`` (the ``tau.py`` accessor)."""
    return corner_transpose(d, from_verts, to_verts, A, side)[3]
