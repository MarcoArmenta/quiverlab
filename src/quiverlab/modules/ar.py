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


# --------------------------------------------------------------------------- #
# The Nakayama functor nu / nu^- (Task 2).
# --------------------------------------------------------------------------- #
def _dense_transpose(mat, rows, cols, dom):
    """The exact ``cols x rows`` transpose of a ``rows x cols`` matrix (empty-shape
    safe: ``lm.transpose`` collapses a 0-row matrix to ``[]``, which would break a
    ModuleHom of nonzero target dim)."""
    out = lm.zeros(cols, rows, dom)
    for i in range(rows):
        for j in range(cols):
            out[j][i] = mat[i][j]
    return out


def nakayama_functor(M, _return_ker=False):
    r"""The Nakayama functor value ``nu M = D Hom_A(M, A)``.

    For ``M`` with minimal presentation ``P_1 --d_1--> P_0 -> M -> 0``, applying
    ``nu = D Hom_A(-, A)`` yields the induced injective map ``g: nu P_1 -> nu P_0`` with
    ``ker g = tau M`` and ``coker g = nu M`` (the two ends of
    ``0 -> tau M -> nu P_1 -> nu P_0 -> nu M -> 0``). ``g = D(d_1^*)`` is the k-dual of the
    corner-transpose ``d_1^*`` that :func:`duality._presentation_transpose` builds; since
    ``D`` is exact and contravariant, ``ker(D d_1^*) = D(coker d_1^*) = D(Tr M) = tau M``
    and ``coker(D d_1^*) = D(ker d_1^*) = D Hom_A(M, A) = nu M`` -- exactly, up to the
    per-instance isomorphism the certificate below asserts.

    Self-certified: ``g`` is a validated module map (``check=True``) and
    ``is_isomorphic(ker g, tau M)`` ties ``nu`` to the trusted translate. ``nu(P_v) = I_v``
    (a projective input has ``tau M = 0`` and ``nu M = coker(0 -> nu P_0) = I_v``)."""
    from quiverlab.modules.duality import _presentation_transpose, dualize, tau as tau_of
    from quiverlab.modules.hom import is_isomorphic
    from quiverlab.modules.morphism import ModuleHom
    dom = M.domain
    _Rop, N, M0op, d1star, _out_side = _presentation_transpose(M)
    nuP1 = dualize(N)                    # D Hom(P_1, A) = (+) I_{w_j}
    nuP0 = dualize(M0op)                 # D Hom(P_0, A) = (+) I_{v_i}
    g = _dense_transpose(d1star, N.dim, M0op.dim, dom)   # D(d_1^*): nu P_1 -> nu P_0
    # Validate the intertwining when it is meaningful; a map out of / into the zero
    # module (projective M => nu P_1 = 0) is trivially a module map, but the empty-shape
    # matmul in _is_module_map returns a false negative, so skip check there. The
    # ker g ~ tau M certificate below covers every case regardless.
    check = nuP1.dim > 0 and nuP0.dim > 0
    gmap = ModuleHom(nuP1, nuP0, g, check=check)         # self-cert: intertwines D-actions
    K, _iota = gmap.kernel()
    tauM = tau_of(M)
    if not is_isomorphic(K, tauM):
        raise QuiverlabError(
            "nakayama: ker g not isomorphic to tau M -- the induced injective map "
            "is inconsistent with the trusted AR translate (bug)")
    nuM, _proj = gmap.cokernel()
    nuM.name = f"nu({M.name})"
    if _return_ker:
        return nuM, K
    return nuM


def nakayama_functor_minus(M):
    r"""The inverse Nakayama functor ``nu^- M = Hom_A(DA, M)``, implemented as the
    opposite-algebra dual ``nu^-_A M = D( nu_{A^op}( D M ) )`` (``D`` contravariant on
    each side of the ``A^op`` Nakayama functor -- no separate injective-copresentation
    code). ``nu^-(I_v) = P_v`` self-certifies it (asserted in the battery)."""
    from quiverlab.modules.duality import dualize
    out = dualize(nakayama_functor(dualize(M)))
    out.name = f"nu^-({M.name})"
    return out


# --------------------------------------------------------------------------- #
# Stable Hom (mod projectives) (Task 3).
# --------------------------------------------------------------------------- #
def _projective_maps_columns(M, N):
    """Columns (over the vec of the tgt.dim x src.dim matrix space) spanning
    P(M,N) = { (compose with cover) . g : g in Hom(M, P(N)) } -- the maps M -> N that
    factor through the projective cover P(N) ->> N. (Any map factoring through a
    projective factors through the cover of N.)"""
    from quiverlab.modules.morphism import hom_basis
    pi = N.projective_cover_hom()                 # P(N) ->> N  (P37 Task 6)
    cols = []
    for g in hom_basis(M, pi.src):                # g: M -> P(N)
        comp = g.then(pi).matrix                  # M -> N through the projective
        cols.append(_vec(comp))
    return cols


def stable_hom_dim(M, N):
    """dim underline{Hom}(M, N) = dim Hom(M,N) - dim P(M,N), where P(M,N) is the
    subspace of maps that factor through a projective (= through the projective cover
    of N). Right modules over one fixed algebra/side."""
    from quiverlab.modules.hom import _assert_comparable
    from quiverlab.modules.morphism import hom_basis
    _assert_comparable(M, N, "stable Hom")
    dom = M.domain
    homs = hom_basis(M, N)
    hom_cols = [_vec(f.matrix) for f in homs]
    proj_cols = _projective_maps_columns(M, N)
    total = lm.mat_rank(lm.cols_to_matrix(hom_cols), dom) if hom_cols else 0
    if not proj_cols:
        return total
    both = lm.mat_rank(lm.cols_to_matrix(hom_cols + proj_cols), dom)
    # P(M,N) subset Hom(M,N), so rank(hom_cols + proj_cols) == rank(hom_cols):
    #   dim P(M,N) = rank(proj_cols); underline dim = total - dim P(M,N).
    dim_proj = lm.mat_rank(lm.cols_to_matrix(proj_cols), dom)
    assert both == total, "stable_hom_dim: P(M,N) is not inside Hom(M,N) (bug)"
    return total - dim_proj


def hom_factors_through_projective(f):
    """True iff the module map ``f`` factors through a projective module -- i.e. ``f``
    lies in ``P(src, tgt)`` (it factors through the projective cover of ``tgt``)."""
    dom = f.domain
    proj_cols = _projective_maps_columns(f.src, f.tgt)
    fvec = [_vec(f.matrix)]
    if not proj_cols:
        return f.is_zero()
    sol = lm.solve_columns(lm.cols_to_matrix(proj_cols),
                           lm.cols_to_matrix(fvec), dom)
    return sol is not None


# --------------------------------------------------------------------------- #
# The End(M)-action on Ext^1(M, N) (Task 4).
# --------------------------------------------------------------------------- #
def _ext1_data(A, M, N):
    """``(cocycle_mats, coboundary_mats, terms, dmats)`` for ``Ext^1(M, N)``.
    ``cocycle_mats[j]`` = ``f_j: P_1 -> N`` (``N.dim x P_1.dim``), the reconstructed class
    reps over the SAME ambient Hom basis the Ext dims use (``complex_reps`` accessor,
    no second Ext complex); ``coboundary_mats = { psi . d_1 : psi in Hom(P_0, N) }`` (the
    ``delta^0`` image)."""
    from quiverlab.modules.complex_reps import _reconstruct_cocycle, ext_cocycle_data
    from quiverlab.modules.morphism import hom_basis
    dom = A.domain
    terms, dmats, homs, cols_by_deg = ext_cocycle_data(A, M, N, 1)
    homs1 = homs[1] if len(homs) > 1 else []
    d1 = dmats[1] if len(dmats) > 1 else []
    width = len(d1[0]) if (d1 and d1[0]) else 0
    cocycle_mats = [_reconstruct_cocycle(col, homs1, N.dim, width, dom)
                    for col in cols_by_deg.get(1, [])]
    coboundary_mats = [lm.matmul(psi.matrix, d1, dom)
                       for psi in hom_basis(terms[0].module, N)] if (d1 and d1[0]) else []
    return cocycle_mats, coboundary_mats, terms, dmats


def _ext1_action_matrix(phi, cocycle_mats, coboundary_mats, res, dom):
    """The ``e x e`` matrix of the pullback action ``[f] |-> [f . phi_1]`` on the
    ``Ext^1`` basis, ``phi`` an endomorphism of ``M`` and ``phi_1`` its degree-1 lift
    (Task 1). Column ``j`` = coordinates of ``[f_j . phi_1]`` over the cocycle basis (mod
    coboundaries). This is the RIGHT End(M)-module structure by precomposition; Task 5's
    non-split guard arbitrates whether this order carries the almost-split socle."""
    e = len(cocycle_mats)
    if e == 0:
        return lm.zeros(0, 0, dom)
    phis = lift_endomorphism_along_resolution(phi, res, degrees=1)
    if len(phis) < 2:                          # M projective => Ext^1 = 0, unreachable here
        raise QuiverlabError("ar: no degree-1 lift (M has projective dimension 0)")
    phi1 = phis[1]
    basis_cols = [_vec(f) for f in cocycle_mats]
    bnd_cols = [_vec(b) for b in coboundary_mats]
    B = lm.cols_to_matrix(basis_cols + bnd_cols)
    out = lm.zeros(e, e, dom)
    for j, f in enumerate(cocycle_mats):
        composed = lm.matmul(f, phi1, dom)     # f_j . phi_1 : P_1 -> N (a cocycle)
        sol = lm.solve_columns(B, lm.cols_to_matrix([_vec(composed)]), dom)
        if sol is None:
            raise QuiverlabError(
                "ar: f.phi_1 is not a cocycle mod coboundaries "
                "(the lift or the class basis is inconsistent)")
        for i in range(e):
            out[i][j] = sol[0][i]              # cocycle-basis part only (drop coboundaries)
    return out


def end_action_on_ext1(M, N):
    """``(basis, action)`` for the End(M)-action on ``Ext^1(M, N)``.

    ``basis`` is the ``End(M)`` basis (``hom_basis(M, M)``); ``action[i]`` is the
    ``e x e`` matrix (``e = dim Ext^1(M, N)``) of ``[f] |-> [f . (phi_i)_1]`` -- the
    left-to-right ("precompose the degree-1 lift") action in the ``ext_reps`` cohomology
    basis of ``Ext^1(M, N)``. Representation axioms are certified in the battery; whether
    this order carries the almost-split socle is decided by Task 5's non-split guard."""
    from quiverlab.modules.morphism import hom_basis
    A = M.algebra
    dom = M.domain
    basis = hom_basis(M, M)
    cocycle_mats, coboundary_mats, terms, dmats = _ext1_data(A, M, N)
    res = (terms, dmats)
    action = [_ext1_action_matrix(phi, cocycle_mats, coboundary_mats, res, dom)
              for phi in basis]
    return basis, action


# --------------------------------------------------------------------------- #
# Almost-split sequences (Task 5).
# --------------------------------------------------------------------------- #
def _rad_end_basis(M):
    """``(H, rad_coords)``: ``H = hom_space(M, M)`` (the End(M) basis as matrices) and
    ``rad_coords`` a basis of ``rad End(M)`` as coordinate columns over ``H``. For ``M``
    indecomposable, End is local and ``rad End = End^{perp}`` w.r.t. the trace form
    ``tr_M(H_i H_j)`` (Dickson / Cohen--Ivanyos--Wales: ``char 0`` or ``char > dim M`` --
    the same scope :mod:`quiverlab.modules.decompose` relies on). Refuses loudly outside
    that char scope."""
    from quiverlab.modules.hom import hom_space
    dom = M.domain
    char = dom.characteristic
    if not (char == 0 or char > M.dim):
        raise QuiverlabError(
            f"almost_split: rad End(M) unreliable in char {char} <= dim M = {M.dim} "
            "(the trace-form radical over-counts, e.g. on k[x]/(x^p))",
            hint="run over QQ or a characteristic > dim M (e.g. GF(32003))")
    H = hom_space(M, M)
    r = len(H)
    T = lm.zeros(r, r, dom)
    for i in range(r):
        for j in range(r):
            T[i][j] = _trace(lm.matmul(H[i], H[j], dom), dom)
    return H, lm.kernel_columns(T, dom)


def _joint_kernel(mats, e, dom):
    """The joint kernel ``cap_r ker(mats[r])`` as coordinate columns (each ``mats[r]`` an
    ``e x e`` matrix). Empty ``mats`` (rad End = 0) => the whole space is annihilated."""
    if not mats:
        return [lm.col(lm.identity(e, dom), j) for j in range(e)] if e else []
    stacked = lm.vstack(mats)
    return lm.kernel_columns(stacked, dom)


def _combine_matrices(mats, coeffs, rows, cols, dom):
    """The linear combination ``sum coeffs[j] * mats[j]`` as a ``rows x cols`` matrix."""
    out = lm.zeros(rows, cols, dom)
    for c, mat in zip(coeffs, mats):
        if dom.is_zero(c):
            continue
        for i in range(rows):
            oi, mi = out[i], mat[i]
            for j in range(cols):
                oi[j] = dom.add(oi[j], dom.mul(c, mi[j]))
    return out


def _rad_action_matrices(action, rad_coords, e, dom):
    """The action of each ``rad End(M)`` basis element on ``Ext^1``: the same End-coord
    combination ``sum_i rc[i] * action[i]`` of the ``action`` matrices."""
    out = []
    for rc in rad_coords:
        Arad = lm.zeros(e, e, dom)
        for i, c in enumerate(rc):
            if dom.is_zero(c):
                continue
            ai = action[i]
            for a in range(e):
                Aa, ra = Arad[a], ai[a]
                for b in range(e):
                    Aa[b] = dom.add(Aa[b], dom.mul(c, ra[b]))
        out.append(Arad)
    return out


def almost_split_sequence(M):
    r"""The almost-split (Auslander-Reiten) sequence ``0 -> tau M -> E -> M -> 0`` for
    ``M`` certified-indecomposable and NON-projective.

    Algorithm (ARS IV.1-IV.3): ``End(M)`` is local (Fitting), so ``rad End(M)`` is the
    trace-form radical (char 0 or char > dim M); it acts on ``Ext^1(M, tau M)`` (Task 4)
    and the classes annihilated by ``rad End(M)`` form the SOCLE, which is nonzero and
    (over ``End/rad = k``) 1-dimensional. ANY nonzero socle class ``xi`` gives THE
    almost-split sequence: reconstruct its cocycle ``f: P_1 -> tau M`` and push out.

    Certification basis: by ARS IV.1-IV.3 a sequence ``0 -> tau M -> E -> M -> 0`` whose
    class lies in ``soc Ext^1(M, tau M)`` with both ends indecomposable IS almost split.
    The returned :class:`ShortExactSequence` additionally carries the computational
    self-certs -- exact (Yoneda ``assert_exact`` + P37 SES), non-split (P37 ``is_split``
    False), ends indecomposable (Plan 30) -- so nothing rests on trusting the socle pick:
    a wrong pick splits and is refused loudly.

    Refuses loudly for a projective, decomposable, or undecidable (budget/char) input."""
    from quiverlab.modules.decompose import is_indecomposable
    from quiverlab.modules.duality import tau as tau_of
    from quiverlab.modules.morphism import ModuleHom
    from quiverlab.modules.ses import ShortExactSequence
    from quiverlab.modules.yoneda import baer_extension
    A = M.algebra
    dom = M.domain
    if not is_indecomposable(M):                 # raises loudly if undecidable
        raise QuiverlabError("almost_split: M must be indecomposable")
    tauM = tau_of(M)
    if tauM.dim == 0:                            # M indecomposable & tau M = 0 <=> M projective
        raise QuiverlabError(
            "almost_split: no AR sequence ends at a projective M (tau M = 0)")
    basis, action = end_action_on_ext1(M, tauM)
    e = len(action[0]) if action else 0
    if e == 0:
        raise QuiverlabError(
            "almost_split: Ext^1(M, tau M) = 0 -- no extension to realize (bug: a "
            "non-projective indecomposable has a nonzero almost-split class)")
    _H, rad_coords = _rad_end_basis(M)
    rad_mats = _rad_action_matrices(action, rad_coords, e, dom)
    socle_cols = _joint_kernel(rad_mats, e, dom)  # cap_r ker(A_rad) = soc Ext^1(M, tau M)
    if not socle_cols:
        raise QuiverlabError(
            "almost_split: empty socle of Ext^1(M, tau M) -- the End(M)-action or the "
            "class basis is inconsistent")
    xi = socle_cols[0]                            # any nonzero socle class
    cocycle_mats, _cob, terms, dmats = _ext1_data(A, M, tauM)
    width = len(dmats[1][0]) if (len(dmats) > 1 and dmats[1] and dmats[1][0]) else 0
    f = _combine_matrices(cocycle_mats, xi, tauM.dim, width, dom)
    seq = baer_extension(M, tauM, f, terms, dmats)   # 0 -> tau M -> E -> M -> 0
    seq.assert_exact()
    E = seq.modules[1]
    ses = ShortExactSequence(ModuleHom(seq.modules[0], E, seq.maps[0], check=False),
                             ModuleHom(E, seq.modules[2], seq.maps[1], check=False))
    if ses.is_split():
        raise QuiverlabError(
            "almost_split: constructed sequence splits -- the socle class was not the "
            "almost-split class (bug/convention)")
    if not (is_indecomposable(seq.modules[0]) and is_indecomposable(M)):
        raise QuiverlabError("almost_split: an end is decomposable (bug)")
    return ses


# --------------------------------------------------------------------------- #
# Irreducible maps: dim rad(M,N)/rad^2(M,N) (Task 6).
# --------------------------------------------------------------------------- #
def _rad_basis(M, N):
    """A spanning column set (vec) of rad(M, N) in the ``tgt.dim x src.dim`` vec space:
    ALL of Hom(M, N) when ``M`` is not iso to ``N`` (a map between non-isomorphic
    indecomposables is never invertible, so it lies in the radical); ``rad End(M)`` when
    ``M ~ N``. ``M``, ``N`` indecomposable."""
    from quiverlab.modules.hom import is_isomorphic
    from quiverlab.modules.morphism import hom_basis
    if not is_isomorphic(M, N):
        return [_vec(f.matrix) for f in hom_basis(M, N)]
    H, rad_coords = _rad_end_basis(M)
    dom = M.domain
    return [_vec(_combine_matrices(H, rc, M.dim, M.dim, dom)) for rc in rad_coords]


def irreducible_maps(M, N, within):
    """The AR-quiver arrow multiplicity ``dim rad(M, N) / rad^2(M, N)`` for
    indecomposable ``M``, ``N``, with ``rad^2`` computed by composing radical maps
    through the finite indecomposable set ``within`` (the knitted component's modules):
    ``rad^2(M, N) = sum_{L in within} span{ h.then(g) : h in rad(M, L), g in rad(L, N) }``.
    Exact once ``within`` contains all indecomposables (rep-finite)."""
    dom = M.domain
    rad_MN = _rad_basis(M, N)
    if not rad_MN:
        return 0
    dim_rad = lm.mat_rank(lm.cols_to_matrix(rad_MN), dom)
    rad2 = []
    for L in within:
        left = _rad_basis(M, L)                        # vecs of rad(M, L)  (L.dim x M.dim)
        right = _rad_basis(L, N)                       # vecs of rad(L, N)  (N.dim x L.dim)
        for hvec in left:
            h = _unvec(hvec, L.dim, M.dim)             # M -> L
            for gvec in right:
                g = _unvec(gvec, N.dim, L.dim)         # L -> N
                rad2.append(_vec(lm.matmul(g, h, dom)))  # g . h : M -> N
    if not rad2:
        return dim_rad
    both = lm.mat_rank(lm.cols_to_matrix(rad_MN + rad2), dom)
    assert both == dim_rad, "irreducible_maps: rad^2 not inside rad (bug)"
    dim_rad2 = lm.mat_rank(lm.cols_to_matrix(rad2), dom)
    return dim_rad - dim_rad2
