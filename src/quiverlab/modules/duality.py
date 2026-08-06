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


def _presentation_transpose(M):
    """The corner-transpose ``d_1^*`` of the minimal projective presentation
    ``P_1 --d_1--> P_0 -> M -> 0`` as a genuine ``R^op``-module map
    ``M0op = Hom_A(P_0, A) --d_1^*--> N = Hom_A(P_1, A)`` between the ``R^op``-projective
    covers. Returns ``(Rop, N, M0op, d1star, out_side)`` where ``coker(d1star) = Tr M``
    and, dualising, ``ker(D d1star) = tau M`` / ``coker(D d1star) = nu M`` (Plan 41).

    Single source of the corner-slicing: :func:`transpose_module` and
    :func:`quiverlab.modules.ar.nakayama_functor` both consume this. Since Plan 43 the
    slicing lives once in :func:`quiverlab.derived._corner.corner_transpose` (the shared
    home the derived AR translate also imports); this presents the minimal projective
    presentation of ``M`` to it. ``N`` (and thus ``d1star``) is empty exactly when ``M``
    is projective (``P_1 = 0``)."""
    from quiverlab.derived._corner import corner_transpose
    from quiverlab.modules.resolution import minimal_resolution

    R = M.algebra
    terms, dmats = minimal_resolution(M, 1)
    v_list = terms[0].vertices          # P_0 summand vertices (target of d_1)
    w_list = terms[1].vertices          # P_1 summand vertices (source of d_1)
    d1 = dmats[1]                       # P_1 -> P_0 in k-bases (P0.dim x P1.dim)
    return corner_transpose(d1, from_verts=v_list, to_verts=w_list, A=R, side=M.side)


def transpose_module(M):
    """Tr M as a right (M.algebra)^op-module, from the minimal projective
    presentation P_1 -> P_0 -> M -> 0 (minimality => no spurious projective
    summands). Tr M = coker(d_1^*), d_1^* the corner-transpose of d_1."""
    from quiverlab.modules.radtopsoc import quotient
    _Rop, N, _M0op, d1star, out_side = _presentation_transpose(M)
    if N.dim == 0:                      # M projective (or zero) => Tr M = 0
        return _zero_module(_Rop, side=out_side)
    piv = lm.column_space_pivots(d1star, dom=M.domain) if (d1star and d1star[0]) else []
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
