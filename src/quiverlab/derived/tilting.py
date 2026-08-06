"""Tilting-complex verifier + End(T) as an algebra (Plan 43 / Rickard).

``is_tilting_complex`` DECIDES (never semi-decides) two conditions on a list of
perfect summands: (1) **rigidity** -- ``Hom_{D^b}(T, T[n]) = 0`` for all ``n != 0`` --
scanned on the EXACT window outside which hyper-Hom is provably the zero cochain group
(``n in [min lo_i - max hi_j, max hi_i - min lo_j]``), reported honestly; (2)
**generation** of ``K^b(proj)`` -- the K0 g-matrix (rows = summand Euler characteristics
in ``K0 = Z^{#vertices}``) is square (``#summands = #simples``) and unimodular
(``det = +-1``), so the classes are a ``Z``-basis of ``K0``. ``T`` is a tilting complex
iff both hold (Rickard).

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
    is_tilting: bool
    rigid: bool                 # hyper_hom_dims(T, T)[n] == 0 for all n != 0 in window
    generates: bool             # K0 g-matrix square & unimodular (det +-1)
    window: tuple               # (n_min, n_max): the EXACT rigidity-check range
    g_matrix: list              # rows = summand K0 classes chi (Euler char per vertex)
    det: int


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _span(cx):
    ds = cx.degrees()
    return (ds[0], ds[-1]) if ds else (0, 0)


def _chi(cx, verts):
    """The K0 class chi(cx) = sum_n (-1)^n dim-vec(cx_n) in Z^{#vertices}."""
    out = [0] * len(verts)
    for n in cx.degrees():
        dv = cx.term(n).dimension_vector()
        for i, v in enumerate(verts):
            out[i] += (-1) ** n * dv.get(v, 0)
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
    :class:`TiltingReport`; rigidity is decided on the EXACT reported window, generation
    by ``Z``-unimodularity of the K0 g-matrix (a single object has ``generates=False`` --
    a 2-term *silting* object need not be *tilting*; read ``.rigid`` + the g-vector)."""
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
    g = [_chi(T, verts) for T in summands]
    det = int(sp.Matrix(g).det()) if len(g) == len(verts) else 0
    generates = (len(g) == len(verts)) and det in (1, -1)
    return TiltingReport(is_tilting=rigid and generates, rigid=rigid,
                         generates=generates, window=(n_min, n_max),
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
