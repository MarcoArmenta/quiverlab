"""Auslander-Reiten completion (Plan 41 / C3): the general chain-map lift, the
Nakayama functor, stable Hom, the End(M)-action on Ext^1, almost-split
sequences, irreducible maps, and AR-quiver knitting. Right modules; exact over
any Domain. Every constructed object self-certifies -- correctness never rests
on trusting a construction."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.resolution import minimal_resolution


# --------------------------------------------------------------------------- #
# Vec / unvec: the single column-stacked flatten shared by Tasks 3-6. A matrix
# (rows x cols) flattens column-major to a length rows*cols vector, and unvec
# rebuilds a (rows x cols) matrix.
# --------------------------------------------------------------------------- #
def _vec(mat):
    """Column-major flatten of a matrix (list[list])."""
    if not mat:
        return []
    rows, cols = len(mat), len(mat[0])
    return [mat[i][j] for j in range(cols) for i in range(rows)]


def _unvec(vec, rows, cols):
    """Inverse of _vec: a length rows*cols column-major vector -> rows x cols matrix."""
    return [[vec[j * rows + i] for j in range(cols)] for i in range(rows)]


def _trace(mat, dom):
    s = dom.zero()
    for i in range(len(mat)):
        s = dom.add(s, mat[i][i])
    return s


# --------------------------------------------------------------------------- #
# The general chain-map lift (Task 1).
# --------------------------------------------------------------------------- #
def _summand_layout(term, algebra):
    """[(vertex, offset, dim, basis_labels)] for the direct-sum blocks of a
    resolution term P_n = (+)_s P_{v_s} -- rebuilt from projective(A, v), whose
    ordered path basis matches the block layout the resolution used (the exact
    idiom of ext_algebra._SimpleResolution.info / complex_reps._terms_info)."""
    from quiverlab.modules.builders import projective
    out, off = [], 0
    for v in term.vertices:
        Pv = projective(algebra, v)
        out.append((v, off, Pv.dim, Pv._pv_basis_labels))
        off += Pv.dim
    return out


def _solve_generator_image(dtarget, source_mod, vertex, rhs, dom):
    """Solve ``dtarget . y = rhs`` for a generator image ``y`` in the vertex-``vertex``
    component of the DOMAIN ``source_mod`` of ``dtarget`` (the exact
    ``ext_algebra._solve_in_component`` routine: restrict to the idempotent-``e_vertex``
    columns of the domain, solve, canonicalise with ``reduce_mod_nullspace``, scatter
    back). Empty component => rhs must be 0 (the only image is 0)."""
    ev = source_mod.action[f"e_{vertex}"]
    dim_src = source_mod.dim
    cols_v = [c for c in range(dim_src) if not dom.is_zero(ev[c][c])]
    if not cols_v:
        return [dom.zero()] * dim_src
    Dred = [[dtarget[r][c] for c in cols_v] for r in range(len(dtarget))]
    z = linalg.solve(Dred, rhs, dom)
    if z is None:
        raise QuiverlabError(
            "ar: chain-map lift solve is inconsistent",
            hint="on a genuine minimal resolution the obstruction is a boundary; "
                 "report the presentation")
    z = linalg.reduce_mod_nullspace(z, Dred, dom)
    y = [dom.zero()] * dim_src
    for c, zc in zip(cols_v, z):
        y[c] = zc
    return y


def _map_from_generator_images(layout, target_mod, images, dom):
    """The full matrix P_src -> target sending each summand generator (the block's
    first basis vector) to images[s], extended A-linearly on the path basis:
    (generator_s . path) |-> images[s] . path = action[path] @ images[s]."""
    G = lm.zeros(target_mod.dim, sum(dv for _, _, dv, _ in layout), dom)
    for (v, off, dv, labels), y in zip(layout, images):
        for k, plabel in enumerate(labels):
            col = lm.matvec(target_mod.action[plabel], y, dom)
            for r in range(target_mod.dim):
                G[r][off + k] = col[r]
    return G


def _term_exhausted(terms, n):
    """True once the resolution has run out at degree ``n`` (the module has finite
    projective dimension ``< n``): the higher chain-map components are maps between
    zero modules, so the lift is COMPLETE and truncates here."""
    return n >= len(terms) or terms[n].module is None or terms[n].dim == 0


def lift_endomorphism_along_resolution(phi, res=None, degrees=1):
    """Lift a module endomorphism ``phi: M -> M`` to a chain map over ``M``'s minimal
    projective resolution.

    Returns ``[phi_0, ..., phi_degrees]`` with ``phi_n`` the ``dim P_n x dim P_n``
    matrix of the degree-``n`` component: ``eps . phi_0 == phi . eps`` (``phi_0``
    covers ``phi``; ``eps = dmats[0]`` the augmentation) and ``d_n . phi_n ==
    phi_{n-1} . d_n`` (a chain map). Canonical (``reduce_mod_nullspace``), so the lift
    is byte-reproducible. Pass a precomputed ``res = minimal_resolution(M, degrees+1)``
    to avoid recomputation."""
    from quiverlab.modules.morphism import ModuleHom
    if isinstance(phi, ModuleHom):
        M, phimat = phi.src, phi.matrix
        if phi.tgt is not phi.src:
            raise QuiverlabError("ar: lift needs an endomorphism M -> M")
    else:                                         # (matrix, M) -- read from res owner
        raise QuiverlabError("ar: pass a ModuleHom endomorphism")
    dom = M.domain
    if res is None:
        res = minimal_resolution(M, degrees + 1)
    terms, dmats = res
    if _term_exhausted(terms, 0):
        # M = 0: the empty chain map. (Nonzero M always has a nonzero cover P_0.)
        return [lm.zeros(0, 0, dom)]
    layout0 = _summand_layout(terms[0], M.algebra)
    P0 = terms[0].module
    eps = dmats[0]                                # P_0 ->> M augmentation
    # phi_0 covers phi: eps . phi_0 = phi . eps, solved per P_0 generator against eps.
    imgs0 = []
    for (v, off, dv, labels) in layout0:
        gen_col = [eps[r][off] for r in range(len(eps))]        # eps(generator_s) in M
        rhs = lm.matvec(phimat, gen_col, dom)                   # phi(eps(gen_s))
        imgs0.append(_solve_generator_image(eps, P0, v, rhs, dom))
    phis = [_map_from_generator_images(layout0, P0, imgs0, dom)]
    # phi_{n}: d_n . phi_n = phi_{n-1} . d_n, per generator via the component solve.
    for n in range(1, degrees + 1):
        if _term_exhausted(terms, n):
            break                             # pd M < n: the chain map is complete
        layout = _summand_layout(terms[n], M.algebra)
        Pn = terms[n].module
        dn = dmats[n]
        imgs = []
        for (v, off, dv, labels) in layout:
            gencol = [dn[r][off] for r in range(len(dn))]       # d_n(generator_s)
            rhs = lm.matvec(phis[n - 1], gencol, dom)           # phi_{n-1}(d_n gen)
            imgs.append(_solve_generator_image(dn, Pn, v, rhs, dom))
        phis.append(_map_from_generator_images(layout, Pn, imgs, dom))
    return phis
