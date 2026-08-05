"""Minimal left/right add(M)-approximations (Plan 44 / C7).

A right add(M)-approximation ``f: X -> C`` (``X in add M``) means every ``M -> C``
factors through ``f`` (``Hom_A(M, -)`` surjectivity). The MINIMAL one is unique up to
isomorphism (ASS I.2, I.5); its source has no superfluous summand. Dually a left
add(M)-approximation ``g: C -> X`` (``X in add M``) has every ``C -> M`` factor through
``g`` (``Hom_A(-, M)`` surjectivity).

We build the minimal approximation by RETRACT PRUNING. This is the fix for the
DECOMPOSABLE-``M`` bug of the old "generators mod rad End(M)" count: that count
multiplied the WHOLE module ``M`` by the generator number, so it could never drop a
superfluous indecomposable summand of ``M`` -- e.g. ``right_add_approximation(A_A, S_v)``
returned the regular module ``A_A`` instead of the projective cover ``P_v``, and
``right_add_approximation(A_A, A_A)`` returned ``A_A^2`` instead of ``A_A``. Pruning:

  1. decompose ``M`` into its distinct indecomposable summand TYPES ``M_1, ..., M_r``
     (``add M = add((+) M_i)``; multiplicities inside ``M`` are irrelevant to ``add M``);
  2. start from the FULL approximation -- one source block per basis element of every
     ``Hom(M_i, C)`` (right) resp. ``Hom(C, M_i)`` (left). This trivially surjects: the
     identity of each ``M_i`` recovers all of ``Hom(M_i, C)`` resp. ``Hom(C, M_i)``;
  3. PRUNE: while some single indecomposable block can be dropped with the map staying
     an approximation (the surjectivity test), drop it.

Pruning provably terminates (blocks strictly decrease) and lands on the minimal
approximation: an irredundant right approximation is right-minimal, hence isomorphic to
the (unique) minimal one. Sketch: if ``f`` is not right-minimal then ``X = X_m (+) X_s``
with ``f|_{X_s} = 0`` and ``X_s != 0`` (ARS); a Krull--Schmidt exchange replaces an
indecomposable summand ``N`` of ``X_s`` by a physical block ``B`` of the same type with
``X = N (+) ker(pi_B)``, and since ``f|_N = 0`` the map factors through ``ker(pi_B) =``
(X with block B removed), so ``B`` is removable -- contradicting irredundancy. The left
case is exactly dual. The right-minimal / left-minimal approximation is unique, so the
result is independent of the pruning order.

The distinct summand types come from :func:`quiverlab.modules.decompose.decompose`, so
both constructors inherit its char refusal loudly (rigorous over ``char 0`` or
``char > dim M``; never a silently wrong add(M) generator set). Float-free, exact linear
algebra over the Domain (ASS I.2, I.5).
"""
from __future__ import annotations

from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import _assert_comparable
from quiverlab.modules.morphism import ModuleHom, direct_sum, hom_basis


def _add_summand_types(M, budget):
    """The distinct indecomposable summand TYPES of ``M`` -- one representative per
    isomorphism class (``add M = add((+) M_i)``; multiplicities inside ``M`` do not
    matter for ``add M``). Raises loudly off the ``decompose`` char scope (never a
    silently wrong add(M))."""
    from quiverlab.modules.decompose import decompose
    return [Mi for (Mi, _mult) in decompose(M, budget=budget)]


def _vec_hom(g):
    """Column-major flatten of a :class:`ModuleHom`'s ``tgt.dim x src.dim`` matrix. Homs
    with the SAME (src, tgt) flatten into the SAME coordinate space, so ranks compare
    like-with-like."""
    m = g.matrix
    tgt = len(m)
    src = len(m[0]) if m else 0
    return [m[i][j] for j in range(src) for i in range(tgt)]


def _surjects(image_cols, target_count, dom):
    """Does the span of ``image_cols`` (columns in a fixed hom space ``Hom(U, V)``) fill
    that whole hom space? The columns already lie in ``Hom(U, V)`` (dim ``target_count``),
    so this is exactly ``rank(image_cols) == target_count``. An empty target is filled
    vacuously; a nonempty target with no image columns is not."""
    if target_count == 0:
        return True
    if not image_cols:
        return False
    return lm.mat_rank(lm.cols_to_matrix(image_cols), dom) == target_count


def _right_approx_blocks(M, C, budget):
    """Ingredients of the right ``add(M)``-approximation of ``C``: the distinct summand
    ``types``, the FULL (trivially surjective) starting ``blocks`` -- one ``(type index t,
    hom g: M_t -> C)`` per basis element of every ``Hom(M_t, C)`` -- and the ``is_approx``
    predicate on a block list. Exposed so the retract prune is testable in isolation (the
    idempotence self-cert)."""
    dom = C.domain
    types = _add_summand_types(M, budget)
    to_C = [hom_basis(Mi, C) for Mi in types]                 # to_C[t]: [M_t -> C]
    target_count = [len(h) for h in to_C]                     # dim Hom(M_t, C)
    between = [[hom_basis(types[s], types[t]) for t in range(len(types))]
               for s in range(len(types))]                    # between[s][t]: [M_s -> M_t]
    blocks = [(t, g) for t in range(len(types)) for g in to_C[t]]

    def is_approx(bl):
        # for each SOURCE type M_s the induced map Hom(M_s, X) -> Hom(M_s, C), whose image
        # is the span of { chi.then(g) : block (t, g), chi: M_s -> M_t }, must be onto.
        for s in range(len(types)):
            if target_count[s] == 0:
                continue
            cols = [_vec_hom(chi.then(g)) for (t, g) in bl for chi in between[s][t]]
            if not _surjects(cols, target_count[s], dom):
                return False
        return True

    return types, blocks, is_approx


def _left_approx_blocks(M, C, budget):
    """Dual of :func:`_right_approx_blocks`: the full starting ``blocks`` are ``(type index
    t, hom g: C -> M_t)`` per basis element of every ``Hom(C, M_t)``, and ``is_approx``
    tests ``Hom(-, M_s)``-surjectivity (blocks live in the TARGET ``M^k``)."""
    dom = C.domain
    types = _add_summand_types(M, budget)
    from_C = [hom_basis(C, Mi) for Mi in types]               # from_C[t]: [C -> M_t]
    target_count = [len(h) for h in from_C]                   # dim Hom(C, M_t)
    between = [[hom_basis(types[s], types[t]) for t in range(len(types))]
               for s in range(len(types))]                    # between[s][t]: [M_s -> M_t]
    blocks = [(t, g) for t in range(len(types)) for g in from_C[t]]

    def is_approx(bl):
        # for each TARGET type M_s the induced map Hom(X, M_s) -> Hom(C, M_s), whose image
        # is the span of { g.then(chi) : block (t, g), chi: M_t -> M_s }, must be onto.
        for s in range(len(types)):
            if target_count[s] == 0:
                continue
            cols = [_vec_hom(g.then(chi)) for (t, g) in bl for chi in between[t][s]]
            if not _surjects(cols, target_count[s], dom):
                return False
        return True

    return types, blocks, is_approx


def right_add_approximation(M, C, budget=512):
    """The minimal right ``add(M)``-approximation ``f: X ->> C``: every ``M -> C`` factors
    through ``f`` and ``X in add M`` has no superfluous summand (unique up to iso). Built
    by retract pruning (see the module docstring); ``check=True`` re-certifies ``f`` IS an
    A-module map. Inherits the ``decompose`` char refusal loudly."""
    _assert_comparable(M, C, "add-approximation")
    types, blocks, is_approx = _right_approx_blocks(M, C, budget)
    _prune(blocks, is_approx)
    return _assemble_right(types, blocks, M, C, C.domain)


def left_add_approximation(M, C, budget=512):
    """The minimal left ``add(M)``-approximation ``g: C >-> X``: every ``C -> M`` factors
    through ``g`` and ``X in add M`` has no superfluous summand (unique up to iso). Built
    DIRECTLY (not by dualizing) so ``g.src`` IS the passed ``C`` (dual of the right case:
    ``Hom(-, M_s)``-surjectivity, blocks pruned from the TARGET). Inherits the
    ``decompose`` char refusal loudly."""
    _assert_comparable(M, C, "add-approximation")
    types, blocks, is_approx = _left_approx_blocks(M, C, budget)
    _prune(blocks, is_approx)
    return _assemble_left(types, blocks, M, C, C.domain)


def _prune(blocks, is_approx):
    """Retract-prune ``blocks`` IN PLACE: drop indecomposable blocks one at a time while
    the reduced map is still an approximation, until none is droppable. Terminates (blocks
    strictly decrease) at the irredundant = minimal approximation (unique up to iso, so
    order-independent). The starting ``blocks`` MUST be an approximation."""
    changed = True
    while changed:
        changed = False
        for idx in range(len(blocks)):
            if is_approx(blocks[:idx] + blocks[idx + 1:]):
                del blocks[idx]
                changed = True
                break


def _assemble_right(types, blocks, M, C, dom):
    """Assemble ``f = [g_0 | ... ]: X = (+) M_{t_i} ->> C`` from the surviving right
    blocks; the empty case (``Hom(M, C) = 0``) is the zero map out of the zero module."""
    if not blocks:
        from quiverlab.modules.duality import _zero_module
        Z = _zero_module(M.algebra, side=M.side)
        return ModuleHom(Z, C, lm.zeros(C.dim, 0, dom), check=False)
    D, _incls, _projs = direct_sum(*[types[t] for (t, _g) in blocks])
    fmat = lm.zeros(C.dim, D.dim, dom)
    off = 0
    for (t, g) in blocks:
        Mt = types[t]
        for r in range(C.dim):
            gr = g.matrix[r]
            for j in range(Mt.dim):
                fmat[r][off + j] = gr[j]
        off += Mt.dim
    return ModuleHom(D, C, fmat, check=True)


def _assemble_left(types, blocks, M, C, dom):
    """Assemble ``g: C >-> X = (+) M_{t_i}`` stacking the surviving left blocks; the empty
    case (``Hom(C, M) = 0``) is the zero map into the zero module."""
    if not blocks:
        from quiverlab.modules.duality import _zero_module
        Z = _zero_module(M.algebra, side=M.side)
        return ModuleHom(C, Z, lm.zeros(0, C.dim, dom), check=False)
    D, _incls, _projs = direct_sum(*[types[t] for (t, _g) in blocks])
    gmat = lm.zeros(D.dim, C.dim, dom)
    off = 0
    for (t, g) in blocks:
        Mt = types[t]
        for r in range(Mt.dim):
            gr = g.matrix[r]
            for j in range(C.dim):
                gmat[off + r][j] = gr[j]
        off += Mt.dim
    return ModuleHom(C, D, gmat, check=True)
