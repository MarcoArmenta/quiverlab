"""Tilting-complex verifier + End(T) as an algebra (Plan 43 / Rickard).

``is_tilting_complex`` reports two conditions on a list of perfect summands.
(1) **rigidity** -- ``Hom_{D^b}(T, T[n]) = 0`` for all ``n != 0`` -- is DECIDED
(never semi-decided) on the EXACT window outside which hyper-Hom is provably the zero
cochain group (``n in [min lo_i - max hi_j, max hi_i - min lo_j]``), reported honestly.
(2) **generation** of ``K^b(proj)`` -- ``thick(T) = K^b(proj A)`` -- is honest
**three-valued** (``TiltingReport.generation``), because the necessary K0 datum does NOT
decide it in general. The K0 g-matrix (rows = summand classes in the **projective** basis
``K0(K^b proj A) = (+)_v Z[P_v]`` via :func:`g_proj`, NOT the composition-factor basis)
being square (``#summands = #simples``) and unimodular (``det = +-1``) means the classes
are a ``Z``-basis of ``K0`` (AI Thm 2.27) -- this is **NECESSARY** for generation but
**NOT sufficient**: whether ``rigid + (#summands = rk K0) => tilting`` holds is exactly
**Rickard's rank QUESTION, still OPEN** (thick subcategories are not classified by K0 --
Krah phantom, arXiv:2302.12502; partial answer Zhang for self-orthogonal tau-tilting
modules of finite pd). So generation is CERTIFIED (``"certified"``) only where a
completion theorem reaches: a **2-term** self-orthogonal K0-basis object is 2-term silting
(Iyama-Joergensen-Yang / AIR) and, being two-sided rigid, tilting -- this covers the
regular object ``A = (+)_v P_v`` (width 0) and every 2-term / APR tilt. A **wider**
(non-2-term) rigid K0-basis object is Rickard-open (``"k0_necessary_only"``): the verifier
does NOT independently build ``thick(T)``, so it honestly returns ``is_tilting = "unknown"``
rather than a possibly-unsound hard ``True`` (Plan 67 fix round; the P43 shipped surface
previously returned ``True`` here -- an affirmative answer to Rickard's open question).
This mirrors ``derived/silting.py``'s K0-basis-only rung EXACTLY (both refuse to certify
generation from bare K0). The projective basis is load-bearing: ``_chi`` (the
composition-factor Euler characteristic) gives ``det(Cartan . g_proj)`` and is a systematic
false-negative on non-unimodular Cartan (self-injective/symmetric; Plan 67 Task 0).

``end_algebra_of_complex`` builds ``End_{D^b}(T) = (+)_{i,j} Hom_{D^b}(T_i, T_j)`` as a
structure-constant :class:`Algebra` -- the derived-equivalent algebra (Rickard) --
composing degree-0 hyper-Hom classes with :meth:`ChainMap.then` reduced to canonical
homotopy representatives, exactly the ``endomorphism._structure_constants`` template
(product ``b_a * b_b = b_a o b_b = b_b.then(b_a)``). ``corner_cartan_of_complex`` is its
corner-Cartan; for ``T = A`` it equals ``cartan_matrix(A)`` (the ``End(A_A) ~ A`` oracle).

Float-free; every constructed algebra is validated (``from_structure_constants(check=True)``
-- associativity + unit), every window reported."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.complexes import (ChainComplex, identity_chain_map,
                                         _block_diag_module, _hom_total_blocks,
                                         _delta_total, _flatten)


@dataclass
class TiltingReport:
    is_tilting: object          # True | False | "unknown" (three-valued -- see the module
                                # docstring / the `generation` field; "unknown" == rigid but
                                # generation only K0-necessary, Rickard's rank question OPEN)
    rigid: bool                 # hyper_hom_dims(T, T)[n] == 0 for all n != 0 in window
    generates: bool             # K0 g-matrix square & unimodular (det +-1) -- the Z-basis
                                # of K0. NECESSARY for generation, NOT sufficient in general
                                # (Rickard's rank question, OPEN); read `generation` for the
                                # honest three-valued verdict.
    generation: str             # "certified" (provable -- 2-term/regular, IJY completion) |
                                # "k0_necessary_only" (rigid K0-basis but non-2-term:
                                # Rickard-open, thick(T) not independently built) | "no"
                                # (not a K0-basis at all)
    window: tuple               # (n_min, n_max): the EXACT rigidity-check range
    g_matrix: list              # rows = summand K0 classes in the PROJECTIVE basis
                                # K0(K^b proj A) = (+)_v Z[P_v] via g_proj (Plan 67 Task 0)
                                # -- NOT the composition-factor basis (see g_proj below)
    det: int                    # det(g_proj) -- unimodularity in the projective K0 basis


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _span(cx):
    ds = cx.degrees()
    return (ds[0], ds[-1]) if ds else (0, 0)


def _is_two_term(summands):
    """The whole object is concentrated in two consecutive degrees up to ONE common shift
    (2-term is defined up to shift): every summand's span has width <= 1 AND all their tops
    fit a common 2-window. True for the empty/stalk cases. This is the hypothesis of the
    Iyama-Joergensen-Yang / AIR 2-term completion theorem (a 2-term presilting K0-basis
    object is silting), so it is exactly the boundary between CERTIFIED and Rickard-open
    generation. Shared by ``derived/silting.py`` (imported there -- one source of truth)."""
    spans = [_span(T) for T in summands if T.degrees()]
    if not spans:
        return True
    widths = {hi - lo for lo, hi in spans}
    if not (widths <= {0, 1}):
        return False
    tops = {hi for _, hi in spans}
    return max(tops) - min(tops) <= 1          # a common 2-window covers them


def _chi(cx, verts):
    """The K0 class chi(cx) = sum_n (-1)^n dim-vec(cx_n) in the COMPOSITION-FACTOR basis
    K0(mod A) = (+)_v Z[S_v] (dim-vec = the class in that basis). This is NOT the basis in
    which generation of K^b(proj A) is decided -- use :func:`g_proj` (the PROJECTIVE basis)
    for that. Kept only for reference / debugging; it has no consumer in ``src`` (Plan 67
    Task 0 replaced its use in ``is_tilting_complex``). ``chi = C . g_proj`` where ``C`` is
    the Cartan matrix, so ``det(chi) = det(C) . det(g_proj)`` -- equal to ``det(g_proj)``
    ONLY when ``det C = +-1`` (hereditary etc.), which is why the P43 bug survived a green
    kA_n suite but was a systematic false-negative on non-unimodular Cartan."""
    out = [0] * len(verts)
    for n in cx.degrees():
        dv = cx.term(n).dimension_vector()
        for i, v in enumerate(verts):
            out[i] += (-1) ** n * dv.get(v, 0)
    return out


def g_proj(cx, verts):
    """K0 class of a perfect complex in the PROJECTIVE basis of
    ``K0(K^b proj A) = (+)_v Z[P_v]``:
    ``g_proj(cx)[v] = sum_n (-1)^n . (multiplicity of P_v in cx.term(n))``.

    The multiplicity of ``P_v`` in a projective ``Q`` is ``dim_k (top Q)_v``
    (``top(P_v^{m}) = S_v^{m}``). Uses the ``_proj_vertices`` provenance when present
    (``from_projective_resolution`` / ``_direct_sum_complex`` / cone-mutants carry it);
    else it reads the multiplicity from the TOP of each term -- so it is correct for a bare
    ``ChainComplex.stalk(A.projective(v))`` (which carries NO provenance).

    This is the correct g-matrix for deciding generation of ``K^b(proj A)`` (AI Thm 2.27:
    the summand K0 classes are a Z-basis iff ``det g_proj = +-1``). ``_chi`` (the
    composition-factor Euler characteristic) equals ``Cartan . g_proj`` and is WRONG on
    non-unimodular Cartan (self-injective/symmetric; AI Example 2.47, ``det C = 0``) --
    Plan 67 Task 0."""
    idx = {v: i for i, v in enumerate(verts)}
    out = [0] * len(verts)
    prov = getattr(cx, "_proj_vertices", None)
    for n in cx.degrees():
        if prov is not None and n in prov:
            mult = {}
            for v in prov[n]:
                mult[v] = mult.get(v, 0) + 1
        else:
            mult = cx.term(n).top().dimension_vector()
        for v, m in mult.items():
            if v in idx:
                out[idx[v]] += (-1) ** n * m
    return out


def _direct_sum_complex(summands):
    """The perfect complex ``(+)_i T_i`` (block-diagonal terms and differentials),
    reusing ``complexes._block_diag_module`` degreewise. Perfect (each summand is);
    carries the concatenated projective-summand vertex provenance."""
    dom = summands[0].domain
    degs = set()
    for T in summands:
        degs |= set(T.degrees())
    terms, dmats, proj = {}, {}, {}
    for n in sorted(degs):
        M = None
        for T in summands:
            Tn = T.term(n)
            M = Tn if M is None else _block_diag_module(M, Tn, dom, name=f"Tsum_{n}")
        if M is not None and M.dim:
            terms[n] = M
        vs = []
        for T in summands:
            vs += list(getattr(T, "_proj_vertices", {}).get(n, []))
        if vs:
            proj[n] = vs
    for n in sorted(degs):
        tgt = sum(T.term(n - 1).dim for T in summands)
        src = sum(T.term(n).dim for T in summands)
        if tgt == 0 or src == 0:
            continue
        M = lm.zeros(tgt, src, dom)
        ro = co = 0
        for T in summands:
            d = T._dmats.get(n)
            tr, tc = T.term(n - 1).dim, T.term(n).dim
            if d and d[0]:
                for i in range(tr):
                    di, Mi = d[i], M[ro + i]
                    for j in range(tc):
                        Mi[co + j] = di[j]
            ro += tr
            co += tc
        dmats[n] = M
    C = ChainComplex(terms, dmats, check=False)      # block-diag of complexes is a complex
    C._perfect = True                                # each summand perfect
    if proj:
        C._proj_vertices = proj
    return C


def _cochain_vec(f, X, Y, dom):
    """Coordinate vector of a degree-0 chain map ``f: X -> Y`` in
    ``Hom^0(X, Y) = (+)_p Hom(X_p, Y_p)``, in the block/basis order of
    ``_hom_total_blocks(X, Y, 0)`` (the ``_place_hom`` idiom)."""
    blocks, cdim = _hom_total_blocks(X, Y, 0, dom)
    vec = [dom.zero()] * cdim
    for b in blocks:
        p = b["p"]                                   # q = p - 0 = p
        fp = f.component(p)                           # Y_p.dim x X_p.dim
        flat = lm.cols_to_matrix([_flatten(h) for h in b["homs"]])
        coords = lm.solve_columns(flat, lm.cols_to_matrix([_flatten(fp)]), dom)
        if coords is None:
            raise QuiverlabError(
                "end_algebra_of_complex: a chain-map component left the Hom-block "
                "basis -- not a module map (bug)")
        for k, c in enumerate(coords[0]):
            vec[b["offset"] + k] = c
    return vec


# --------------------------------------------------------------------------- #
# the verifier
# --------------------------------------------------------------------------- #
def is_tilting_complex(summands):
    """Decide whether ``summands`` assemble a tilting complex (Rickard). Every summand
    must be a certified perfect complex (loud otherwise). Returns a
    :class:`TiltingReport`. Rigidity is DECIDED on the EXACT reported window. Generation
    is honest **three-valued** (``.generation`` / ``.is_tilting``): the ``Z``-unimodular
    K0 g-matrix (``.generates``) is only the NECESSARY K0-basis condition, so
    ``is_tilting`` is ``True`` only when generation is CERTIFIED (2-term / regular, via the
    IJY completion theorem), ``"unknown"`` when the object is rigid with a K0-basis but is
    non-2-term (Rickard's rank question -- OPEN -- so generation is not independently
    certified; never a possibly-unsound hard ``True``), and ``False`` otherwise. A single
    object has ``generates=False`` -- a 2-term *silting* object need not be *tilting*; read
    ``.rigid`` + the g-vector."""
    import sympy as sp
    from quiverlab.modules.complexes import hyper_hom_dims
    if not summands:
        raise QuiverlabError("is_tilting_complex: need at least one summand")
    for T in summands:
        if not T.is_perfect():
            raise QuiverlabError("is_tilting_complex: every summand must be a "
                                 "certified perfect complex")
    A = summands[0].algebra
    verts = list(A.quiver.vertices)
    spans = [_span(T) for T in summands]
    n_min = min(lo for lo, _ in spans) - max(hi for _, hi in spans)
    n_max = max(hi for _, hi in spans) - min(lo for lo, _ in spans)
    Tsum = _direct_sum_complex(summands)
    hh = hyper_hom_dims(Tsum, Tsum, n_min, n_max)
    rigid = all(hh.get(n, 0) == 0 for n in range(n_min, n_max + 1) if n != 0)
    # generation is decided in the PROJECTIVE K0 basis (+)_v Z[P_v] via g_proj -- NOT the
    # composition-factor basis _chi, which is det(Cartan . g_proj) and a systematic
    # false-negative on non-unimodular Cartan (Plan 67 Task 0; AI Ex 2.2 / Thm 2.27).
    g = [g_proj(T, verts) for T in summands]
    det = int(sp.Matrix(g).det()) if len(g) == len(verts) else 0
    generates = (len(g) == len(verts)) and det in (1, -1)   # K0-necessary (a Z-basis of K0)
    # Generation is three-valued (honest). det = +-1 is NECESSARY but NOT sufficient in
    # general -- "rigid + (#summands = rk K0) => tilting" is exactly Rickard's rank
    # QUESTION, still OPEN (thick subcategories are not K0-classified; Krah phantom). It is
    # CERTIFIED only where a completion theorem reaches THIS engine: a 2-term self-orthogonal
    # K0-basis object is 2-term silting (IJY / AIR) and, being two-sided rigid, tilting --
    # covering the regular object A (width 0) and every 2-term / APR tilt. A wider (non-2-
    # term) rigid K0-basis object is Rickard-open ("k0_necessary_only") -- the verifier does
    # NOT build thick(T), so is_tilting returns "unknown", NOT a possibly-unsound hard True.
    # This mirrors derived/silting.py's K0-basis-only rung exactly (both refuse to certify
    # generation from bare K0). Not a K0-basis at all => "no".
    if not generates:
        generation = "no"
    elif _is_two_term(summands):
        generation = "certified"
    else:
        generation = "k0_necessary_only"
    is_tilting = (True if (rigid and generation == "certified")
                  else ("unknown" if (rigid and generation == "k0_necessary_only")
                        else False))
    return TiltingReport(is_tilting=is_tilting, rigid=rigid, generates=generates,
                         generation=generation, window=(n_min, n_max),
                         g_matrix=g, det=det)


def corner_cartan_of_complex(summands):
    """``[i][j] = dim_k Hom_{D^b}(T_j, T_i)`` (degree 0) as an integer matrix -- the
    corner-Cartan of ``End(T)`` (corner ``e_i End(T) e_j = Hom(T_j, T_i)``). Equals
    ``cartan_matrix(A)`` when ``T = A`` (the ``regular_corner_dims`` analogue). The
    ``Hom(T_j, T_i)`` orientation (transpose of the naive ``[i][j]=Hom(T_i,T_j)``) is
    pinned by that ``T = A`` oracle -- ``dim Hom_A(P_j, P_i) = C[i][j]`` (P37)."""
    from quiverlab.derived.homs import hyper_hom_basis
    return [[len(hyper_hom_basis(Tj, Ti, 0)) for Tj in summands] for Ti in summands]


def end_algebra_of_complex(summands):
    """``End_{D^b}(T)`` as a structure-constant :class:`Algebra` (Rickard -- the
    derived-equivalent algebra). Degree-0 hyper-Hom classes between all summand pairs,
    composed with :meth:`ChainMap.then` and reduced to canonical homotopy reps, then
    handed to ``Algebra.from_structure_constants`` (validated: associativity + unit)."""
    from quiverlab.core.algebra import Algebra
    from quiverlab.derived.homs import hyper_hom_basis
    from quiverlab.fields.linalg import reduce_mod_nullspace, solve
    dom = summands[0].domain
    N = len(summands)
    basis, block_of = [], []
    block_maps, block_gi, coboundary = {}, {}, {}
    for i in range(N):
        for j in range(N):
            cls = hyper_hom_basis(summands[i], summands[j], 0)
            gis = []
            for f in cls:
                gis.append(len(basis))
                basis.append(f)
                block_of.append((i, j))
            block_maps[(i, j)] = cls
            block_gi[(i, j)] = gis
            dn1, _s, _t = _delta_total(summands[i], summands[j], -1, dom)
            coboundary[(i, j)] = ([lm.col(dn1, c) for c in range(len(dn1[0]))]
                                  if (dn1 and dn1[0]) else [])
    n = len(basis)
    if n == 0:
        raise QuiverlabError("end_algebra_of_complex: T has no degree-0 endomorphisms")

    def _coords_in_block(comp, i, j):
        """Coordinates of a degree-0 chain map ``comp: T_i -> T_j`` over the block's
        hyper-Hom classes (the homotopy/coboundary freedom is solved and discarded --
        the class is well-defined; reduce_mod_nullspace canonicalises)."""
        bmaps = block_maps[(i, j)]
        cols = [_cochain_vec(f, summands[i], summands[j], dom) for f in bmaps]
        cols = cols + coboundary[(i, j)]
        if not cols:
            return []
        target = _cochain_vec(comp, summands[i], summands[j], dom)
        Amat = lm.cols_to_matrix(cols)
        x = solve(Amat, target, dom)
        if x is None:
            raise QuiverlabError(
                "end_algebra_of_complex: a composite is not in the class-span + "
                "coboundary (bug)")
        x = reduce_mod_nullspace(x, Amat, dom)
        return x[:len(bmaps)]                          # drop the coboundary coordinates

    Tcon = []
    for a in range(n):
        ia, _ja = block_of[a]
        row = []
        for b in range(n):
            ib, jb = block_of[b]
            vec = [dom.zero()] * n
            if jb == ia:                               # b_a o b_b composable
                comp = basis[b].then(basis[a])         # T_ib -> T_ja
                coords = _coords_in_block(comp, ib, _ja)
                for k, gi in enumerate(block_gi[(ib, _ja)]):
                    vec[gi] = coords[k]
            row.append(vec)
        Tcon.append(row)
    unit = [dom.zero()] * n
    for i in range(N):
        coords = _coords_in_block(identity_chain_map(summands[i]), i, i)
        for k, gi in enumerate(block_gi[(i, i)]):
            unit[gi] = dom.add(unit[gi], coords[k])
    return Algebra.from_structure_constants(Tcon, unit, field=dom, check=True)


def two_term_silting_from_presentation(M):
    """The 2-term complex ``[P_1 --d_1--> P_0]`` (degrees 1, 0) of ``M``'s minimal
    projective presentation, with the rigidity report (the AIR bridge P45 consumes).
    ``TiltingReport.generates`` for a single object is ``False`` -- correct: a 2-term
    *silting* object need not be *tilting*; the consumer reads ``.rigid`` + the g-vector,
    not ``.is_tilting``."""
    from quiverlab.modules.resolution import minimal_resolution
    terms, dmats = minimal_resolution(M, 1)            # P_1 --d_1--> P_0 --> M
    t = {0: terms[0].module}
    d = {}
    prov = {0: list(terms[0].vertices)}
    if terms[1].module is not None and terms[1].dim:
        t[1] = terms[1].module
        prov[1] = list(terms[1].vertices)
        if dmats[1] and dmats[1][0]:
            d[1] = dmats[1]
    cx = ChainComplex(t, d, check=True)
    cx._perfect = True
    cx._proj_vertices = prov
    rep = is_tilting_complex([cx])                      # rigidity report (generation n/a)
    return cx, rep
