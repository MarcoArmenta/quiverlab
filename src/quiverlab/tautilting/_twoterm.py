"""2-term silting complexes over an exact Domain -- the engine behind AIR mutation
(Plan 45 / C4, internal). A support tau-tilting pair (M, P) IS the 2-term silting object
``T = (+)_i [P_1(M_i) -> P_0(M_i)] (+) (+)_{v in P} P_v[1]`` in K^b(proj A) (AIR 2014
Thm 3.2). Mutation is the uniform silting exchange (Aihara-Iyama): at summand k, the cone
of the minimal left add(U)-approximation (or the cocone of the right one), U = T minus
summand k -- whichever lands 2-term is THE tau-tilting neighbour. Every crossing
(module<->support) is handled by this one construction; there is NO module-level shortcut
(a killed projective re-enters as a genuine module -- e.g. over kA2 un-killing vertex 1 of
(S2,{1}) brings back P1, which no approximation of P1 against S2 produces; only the cocone
+ minimal-complex reduction does it).

A complex is stored cohomologically: ``terms`` = ``{degree: [vertex, ...]}`` (each vertex a
projective summand P_v) and ``diffs`` = ``{i: matrix}`` with ``d^i: C^i -> C^{i+1}`` a dense
``dim C^{i+1} x dim C^i`` Domain matrix in the concatenated projective k-bases. A 2-term
presentation P_1 -> P_0 is degrees -1, 0. Float-free; all arithmetic is module-map linear
algebra, so ``coker`` of a reduced 2-term complex is a genuine module."""
from __future__ import annotations

from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import hom_space
from quiverlab.modules.morphism import ModuleHom, direct_sum


# --------------------------------------------------------------------------- #
# projective cache + block geometry
# --------------------------------------------------------------------------- #
_PCACHE = {}


def _proj(A, v):
    key = (id(A), v)
    P = _PCACHE.get(key)
    if P is None:
        from quiverlab.modules.builders import projective
        P = projective(A, v)
        _PCACHE[key] = P
    return P


def _pdim(A, v):
    return _proj(A, v).dim


def _proj_sum(A, verts):
    """The block direct sum ``(+) P_v`` (a Module), the zero module for an empty list."""
    if not verts:
        from quiverlab.modules.duality import _zero_module
        return _zero_module(A, side="right")
    D, _, _ = direct_sum(*[_proj(A, v) for v in verts])
    return D


def _offsets(A, verts):
    """``[(offset, dim), ...]`` of each summand block in the concatenated k-basis."""
    out = []
    o = 0
    for v in verts:
        d = _pdim(A, v)
        out.append((o, d))
        o += d
    return out


class PComplex:
    """A bounded complex of projective right A-modules (cohomological)."""

    def __init__(self, algebra, terms, diffs):
        self.algebra = algebra
        self.domain = algebra.domain
        self.terms = {i: list(vs) for i, vs in terms.items() if vs}
        self.diffs = {i: d for i, d in diffs.items() if d and d[0]}

    def degrees(self):
        return sorted(self.terms)

    def dim_at(self, i):
        return sum(_pdim(self.algebra, v) for v in self.terms.get(i, []))

    def is_zero(self):
        return not self.terms


# --------------------------------------------------------------------------- #
# pair <-> complexes
# --------------------------------------------------------------------------- #
def complex_of_module(M):
    """The minimal projective presentation ``P_1 -> P_0`` of ``M`` (degrees -1, 0)."""
    from quiverlab.modules.resolution import minimal_resolution
    A = M.algebra
    terms_r, dmats = minimal_resolution(M, 1)
    p0 = list(terms_r[0].vertices)
    p1 = list(terms_r[1].vertices) if len(terms_r) > 1 else []
    terms, diffs = {}, {}
    if p0:
        terms[0] = p0
    if p1:
        terms[-1] = p1
        diffs[-1] = dmats[1]
    return PComplex(A, terms, diffs)


def complex_of_support(A, v):
    """The shifted projective ``P_v[1]`` (degree -1 only) -- a killed projective."""
    return PComplex(A, {-1: [v]}, {})


def to_complexes(pair):
    """The indecomposable 2-term silting summands of a pair: module presentations (in
    order) then support shifts ``P_v[1]`` (vertex order)."""
    A = pair.algebra
    verts = list(A.quiver.vertices)
    out = [complex_of_module(Mi) for Mi in pair.summands]
    for v in sorted(pair.support, key=verts.index):
        out.append(complex_of_support(A, v))
    return out


def summand_class(C):
    """Classify a 2-term summand ``C`` (after reduction): ``("module", Module)``,
    ``("support", v)``, or ``("bad", reason)`` when ``C`` is not concentrated in degrees
    {-1, 0} -- the signal that a mutation branch is the WRONG facet."""
    R = reduce_complex(C)
    degs = R.degrees()
    if any(i not in (-1, 0) for i in degs):
        return ("bad", f"non-2-term: degrees {degs}")
    p0 = R.terms.get(0, [])
    pm = R.terms.get(-1, [])
    if not p0 and not pm:
        return ("bad", "reduced to zero")
    if not pm:                                         # [0 -> P0]: projective module
        return ("module", _proj_sum(R.algebra, p0))
    if not p0:                                         # P_v[1]: shifted projective
        if len(pm) != 1:
            return ("bad", f"shifted projective decomposable: {pm}")
        return ("support", pm[0])
    P1 = _proj_sum(R.algebra, pm)
    P0 = _proj_sum(R.algebra, p0)
    d = R.diffs.get(-1) or lm.zeros(P0.dim, P1.dim, R.domain)
    Y, _pi = ModuleHom(P1, P0, d, check=False).cokernel()
    return ("module", Y)


# --------------------------------------------------------------------------- #
# direct sum, shift
# --------------------------------------------------------------------------- #
def direct_sum_complexes(comps):
    """Block direct sum; returns ``(sum, placements)`` where ``placements[i][deg] =
    (offset, dim)`` locates input complex ``i`` at each degree of the sum."""
    A = comps[0].algebra
    dom = A.domain
    all_degs = sorted({i for C in comps for i in C.terms})
    terms = {i: [] for i in all_degs}
    placements = [dict() for _ in comps]
    off = {i: 0 for i in all_degs}
    for ci, C in enumerate(comps):
        for i in all_degs:
            d_here = C.dim_at(i)
            placements[ci][i] = (off[i], d_here)
            terms[i].extend(C.terms.get(i, []))
            off[i] += d_here
    diffs = {}
    for i in all_degs:
        rows = sum(C.dim_at(i + 1) for C in comps)
        cols = sum(C.dim_at(i) for C in comps)
        if rows == 0 or cols == 0:
            continue
        D = lm.zeros(rows, cols, dom)
        roff = coff = 0
        for C in comps:
            r, c = C.dim_at(i + 1), C.dim_at(i)
            di = C.diffs.get(i)
            if di and r and c:
                for a in range(r):
                    Da, dia = D[roff + a], di[a]
                    for b in range(c):
                        Da[coff + b] = dia[b]
            roff += r
            coff += c
        diffs[i] = D
    return PComplex(A, terms, diffs), placements


def shift(C, s):
    """``C[s]``: ``(C[s])^i = C^{i+s}`` with differential ``(-1)^s d`` (keeps d.d = 0)."""
    dom = C.domain
    terms = {i - s: vs for i, vs in C.terms.items()}
    if s % 2 == 0:
        diffs = {i - s: d for i, d in C.diffs.items()}
    else:
        diffs = {i - s: [[dom.neg(x) for x in row] for row in d]
                 for i, d in C.diffs.items()}
    return PComplex(C.algebra, terms, diffs)


# --------------------------------------------------------------------------- #
# Hom in K^b(proj) and composition
# --------------------------------------------------------------------------- #
def _hom_proj(A, srcverts, tgtverts):
    if not srcverts or not tgtverts:
        return []
    return hom_space(_proj_sum(A, srcverts), _proj_sum(A, tgtverts))


def _flatten_map(fmap, degs, shapes, dom):
    """Column-major flatten of a chain map over a FIXED (degree, shape) list, so maps with
    the same (X, Y) flatten into one coordinate space."""
    out = []
    for i in degs:
        rows, cols = shapes[i]
        m = fmap.get(i)
        for b in range(cols):
            for a in range(rows):
                out.append(m[a][b] if (m and rows and cols) else dom.zero())
    return out


def hom_kb(X, Y):
    """A basis of ``Hom_{K^b(proj A)}(X, Y)`` as chain maps (each ``{degree: matrix}``),
    i.e. chain maps modulo null-homotopies. Pure Domain linear algebra."""
    A = X.algebra
    dom = A.domain
    degs = sorted(set(X.terms) | set(Y.terms))
    homs = {i: _hom_proj(A, X.terms.get(i, []), Y.terms.get(i, [])) for i in degs}
    sizes = {i: len(homs[i]) for i in degs}
    order = [i for i in degs if sizes[i] > 0]
    total = sum(sizes[i] for i in order)
    if total == 0:
        return []
    starts, off = {}, 0
    for i in order:
        starts[i] = off
        off += sizes[i]

    def combine(i, coeffs):
        rows, cols = Y.dim_at(i), X.dim_at(i)
        out = lm.zeros(rows, cols, dom)
        for c, mat in zip(coeffs, homs[i]):
            if dom.is_zero(c):
                continue
            for a in range(rows):
                oa, ma = out[a], mat[a]
                for b in range(cols):
                    oa[b] = dom.add(oa[b], dom.mul(c, ma[b]))
        return out

    # chain-map constraint d_Y^i . f^i - f^{i+1} . d_X^i = 0 (in Hom(X^i, Y^{i+1}))
    eqs = []
    for i in degs:
        dX, dY = X.diffs.get(i), Y.diffs.get(i)
        eq_rows, eq_cols = Y.dim_at(i + 1), X.dim_at(i)
        if eq_rows == 0 or eq_cols == 0:
            continue
        block = [[dom.zero()] * total for _ in range(eq_rows * eq_cols)]
        if dY and sizes.get(i, 0):
            for bi, mat in enumerate(homs[i]):
                prod = lm.matmul(dY, mat, dom)          # Y^{i+1} x X^i
                col = starts[i] + bi
                for a in range(eq_rows):
                    pa = prod[a]
                    for b in range(eq_cols):
                        block[a * eq_cols + b][col] = dom.add(
                            block[a * eq_cols + b][col], pa[b])
        if dX and sizes.get(i + 1, 0):
            for bi, mat in enumerate(homs[i + 1]):
                prod = lm.matmul(mat, dX, dom)          # Y^{i+1} x X^i
                col = starts[i + 1] + bi
                for a in range(eq_rows):
                    pa = prod[a]
                    for b in range(eq_cols):
                        block[a * eq_cols + b][col] = dom.sub(
                            block[a * eq_cols + b][col], pa[b])
        eqs.extend(block)
    chain = lm.kernel_columns(eqs, dom) if eqs else \
        [lm.col(lm.identity(total, dom), j) for j in range(total)]
    # null-homotopies: h^i: X^i -> Y^{i-1} induces f^i = d_Y^{i-1} h^i + h^{i+1} d_X^i,
    # so a single h^i contributes to f^i (via d_Y^{i-1} h^i) AND f^{i-1} (via h^i d_X^{i-1}).
    hbases = {i: _hom_proj(A, X.terms.get(i, []), Y.terms.get(i - 1, [])) for i in degs}
    nulls = []
    for i in degs:
        for hmat in hbases.get(i, []):
            vec = [dom.zero()] * total
            dYprev = Y.diffs.get(i - 1)
            if dYprev is not None and sizes.get(i, 0):
                _accum(vec, lm.matmul(dYprev, hmat, dom), homs[i], starts.get(i), dom)
            dXprev = X.diffs.get(i - 1)
            if dXprev is not None and sizes.get(i - 1, 0):
                _accum(vec, lm.matmul(hmat, dXprev, dom), homs[i - 1],
                       starts.get(i - 1), dom)
            nulls.append(vec)
    chosen = lm.independent_modulo(chain, nulls, dom) if chain else []
    out = []
    for j in chosen:
        coords = chain[j]
        out.append({i: combine(i, coords[starts[i]:starts[i] + sizes[i]]) for i in order})
    return out


def _accum(vec, target, basis, start, dom):
    """Express ``target`` in ``basis`` coords and add into ``vec`` at ``start``."""
    if not basis or start is None:
        return
    rows = len(target)
    cols = len(target[0]) if rows else 0
    tvec = [target[a][b] for b in range(cols) for a in range(rows)]
    B = lm.cols_to_matrix([[m[a][b] for b in range(cols) for a in range(rows)]
                           for m in basis])
    x = linalg.solve(B, tvec, dom)
    if x is None:
        return
    for k, c in enumerate(x):
        vec[start + k] = dom.add(vec[start + k], c)


def compose_kb(f, g, dom):
    """Chain-map composition ``f`` then ``g`` (``f: X->Y``, ``g: Y->Z``): degreewise
    ``g^i . f^i``."""
    out = {}
    for i in set(f) & set(g):
        fi, gi = f[i], g[i]
        if fi and fi[0] and gi and gi[0]:
            out[i] = lm.matmul(gi, fi, dom)
    return out


# --------------------------------------------------------------------------- #
# mapping cone
# --------------------------------------------------------------------------- #
def cone(fmap, X, Y):
    """The mapping cone of a chain map ``fmap: X -> Y``: ``cone^n = X^{n+1} (+) Y^n``,
    ``d^n = [[-d_X^{n+1}, 0], [f^{n+1}, d_Y^n]]``."""
    A = X.algebra
    dom = A.domain
    all_n = sorted(set(list(Y.terms) + [i - 1 for i in X.terms]))
    terms = {}
    for n in all_n:
        vs = list(X.terms.get(n + 1, [])) + list(Y.terms.get(n, []))
        if vs:
            terms[n] = vs
    diffs = {}
    ns = sorted(terms)
    for n in ns:
        rows = _dim(A, terms.get(n + 1, []))
        cols = _dim(A, terms.get(n, []))
        if rows == 0 or cols == 0:
            continue
        D = lm.zeros(rows, cols, dom)
        # source split: X^{n+1} (width xw) | Y^n (width yw)
        xw = _dim(A, X.terms.get(n + 1, []))
        yw = _dim(A, Y.terms.get(n, []))
        # target split: X^{n+2} (height xh) | Y^{n+1} (height yh)
        xh = _dim(A, X.terms.get(n + 2, []))
        yh = _dim(A, Y.terms.get(n + 1, []))
        # top-left: -d_X^{n+1} : X^{n+1} -> X^{n+2}
        dX = X.diffs.get(n + 1)
        if dX and xh and xw:
            for a in range(xh):
                for b in range(xw):
                    D[a][b] = dom.neg(dX[a][b])
        # bottom-left: f^{n+1} : X^{n+1} -> Y^{n+1}
        f1 = fmap.get(n + 1)
        if f1 and yh and xw:
            for a in range(yh):
                for b in range(xw):
                    D[xh + a][b] = f1[a][b]
        # bottom-right: d_Y^n : Y^n -> Y^{n+1}
        dY = Y.diffs.get(n)
        if dY and yh and yw:
            for a in range(yh):
                for b in range(yw):
                    D[xh + a][xw + b] = dY[a][b]
        diffs[n] = D
    return PComplex(A, terms, diffs)


def _dim(A, verts):
    return sum(_pdim(A, v) for v in verts)


# --------------------------------------------------------------------------- #
# minimal-complex reduction (Gaussian elimination of iso blocks)
# --------------------------------------------------------------------------- #
def _inv(mat, dom):
    """Inverse of a square Domain matrix (solve against identity columns)."""
    n = len(mat)
    ident = lm.identity(n, dom)
    cols = []
    for j in range(n):
        x = linalg.solve(mat, lm.col(ident, j), dom)
        if x is None:
            return None
        cols.append(x)
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def _sub_block(mat, roff, rdim, coff, cdim):
    return [[mat[roff + a][coff + b] for b in range(cdim)] for a in range(rdim)]


def reduce_complex(C):
    """The minimal model of ``C``: repeatedly cancel an iso block ``phi: P_v -> P_v`` in
    some differential (Gaussian elimination / delooping lemma), until every differential
    has its entries in the radical. Homotopy-equivalent to ``C``, so the summand's iso type
    (hence :func:`summand_class`) is preserved."""
    A = C.algebra
    dom = C.domain
    terms = {i: list(vs) for i, vs in C.terms.items()}
    diffs = {i: [row[:] for row in d] for i, d in C.diffs.items()}
    changed = True
    while changed:
        changed = False
        for i in sorted(diffs):
            d = diffs.get(i)
            if not d or not d[0]:
                continue
            src = terms.get(i, [])
            tgt = terms.get(i + 1, [])
            soff = _offsets(A, src)
            toff = _offsets(A, tgt)
            found = None
            for bt, (ro, rd) in enumerate(toff):        # target summand
                for bs, (co, cd) in enumerate(soff):    # source summand
                    if tgt[bt] != src[bs]:
                        continue
                    blk = _sub_block(d, ro, rd, co, cd)
                    if lm.mat_rank(blk, dom) == rd:      # iso block
                        found = (i, bt, bs, ro, rd, co, cd, blk)
                        break
                if found:
                    break
            if found:
                _cancel(A, dom, terms, diffs, *found)
                diffs = {k: v for k, v in diffs.items() if v and v[0]}
                changed = True
                break
    return PComplex(A, terms, diffs)


def _cancel(A, dom, terms, diffs, i, bt, bs, ro, rd, co, cd, blk):
    """Cancel the iso block at (target summand ``bt`` of C^{i+1}, source summand ``bs`` of
    C^i) in ``diffs[i]``: Schur-complement ``d^i``, then drop the summand from C^i / C^{i+1}
    and the matching row/column from the neighbouring differentials."""
    phi_inv = _inv(blk, dom)
    d = diffs[i]
    src = terms[i]
    tgt = terms[i + 1]
    soff = _offsets(A, src)
    toff = _offsets(A, tgt)
    # Schur complement: d[a'][b'] -= d[a'][bs] phi^{-1} d[bt][b']  (block algebra)
    row_bt = _sub_block(d, ro, rd, 0, _dim(A, src))         # d[bt][*]  (rd x cols)
    col_bs = _sub_block(d, 0, _dim(A, tgt), co, cd)         # d[*][bs]  (rows x cd)
    corr_mid = lm.matmul(phi_inv, row_bt, dom)             # cd x cols
    corr = lm.matmul(col_bs, corr_mid, dom)                # rows x cols
    newd = [[dom.sub(d[a][b], corr[a][b]) for b in range(_dim(A, src))]
            for a in range(_dim(A, tgt))]
    # drop source block bs (cols co..co+cd) and target block bt (rows ro..ro+rd)
    keep_rows = [r for r in range(_dim(A, tgt)) if not (ro <= r < ro + rd)]
    keep_cols = [c for c in range(_dim(A, src)) if not (co <= c < co + cd)]
    diffs[i] = [[newd[r][c] for c in keep_cols] for r in keep_rows]
    # neighbour differentials: d^{i-1}: C^{i-1}->C^i loses target-rows of bs;
    # d^{i+1}: C^{i+1}->C^{i+2} loses source-cols of bt.
    d_prev = diffs.get(i - 1)
    if d_prev and d_prev[0]:
        diffs[i - 1] = [[d_prev[r][c] for c in range(len(d_prev[0]))]
                        for r in keep_cols]               # rows of d_prev index C^i = src
    d_next = diffs.get(i + 1)
    if d_next and d_next[0]:
        diffs[i + 1] = [[d_next[r][c] for c in keep_rows]
                        for r in range(len(d_next))]      # cols of d_next index C^{i+1}=tgt
    terms[i] = [src[b] for b in range(len(src)) if b != bs]
    terms[i + 1] = [tgt[b] for b in range(len(tgt)) if b != bt]
    if not terms[i]:
        terms.pop(i, None)
    if not terms[i + 1]:
        terms.pop(i + 1, None)


# --------------------------------------------------------------------------- #
# minimal left / right add(U)-approximations in K^b (retract pruning)
# --------------------------------------------------------------------------- #
def _prune(blocks, is_approx):
    changed = True
    while changed:
        changed = False
        for idx in range(len(blocks)):
            if is_approx(blocks[:idx] + blocks[idx + 1:]):
                del blocks[idx]
                changed = True
                break


def _surjects(cols, target, dom):
    if target == 0:
        return True
    if not cols:
        return False
    return lm.mat_rank(lm.cols_to_matrix(cols), dom) == target


def min_left_approx(X, types):
    """Minimal left ``add(types)``-approximation ``f: X -> E`` (E in add types): every
    ``X -> T`` (T in add types) factors through ``f``. Returns ``(fmap, E)``."""
    A = X.algebra
    dom = A.domain
    from_X = [hom_kb(X, T) for T in types]
    target = [len(h) for h in from_X]
    between = [[hom_kb(types[t], types[s]) for s in range(len(types))]
               for t in range(len(types))]
    # flatten shapes for Hom(X, types[s])
    sshapes = _shapes(X, types)
    blocks = [(t, g) for t in range(len(types)) for g in from_X[t]]

    def is_approx(bl):
        for s in range(len(types)):
            if target[s] == 0:
                continue
            cols = [_flat(compose_kb(g, chi, dom), types[s], X, sshapes[s], dom)
                    for (t, g) in bl for chi in between[t][s]]
            if not _surjects(cols, target[s], dom):
                return False
        return True

    _prune(blocks, is_approx)
    return _assemble_left(X, blocks, types)


def min_right_approx(X, types):
    """Minimal right ``add(types)``-approximation ``g: E -> X`` (E in add types): every
    ``T -> X`` factors through ``g``. Returns ``(gmap, E)``."""
    A = X.algebra
    dom = A.domain
    to_X = [hom_kb(T, X) for T in types]
    target = [len(h) for h in to_X]
    between = [[hom_kb(types[s], types[t]) for t in range(len(types))]
               for s in range(len(types))]
    sshapes = _shapes_to(types, X)
    blocks = [(t, g) for t in range(len(types)) for g in to_X[t]]

    def is_approx(bl):
        for s in range(len(types)):
            if target[s] == 0:
                continue
            cols = [_flat(compose_kb(chi, g, dom), types[s], X, sshapes[s], dom)
                    for (t, g) in bl for chi in between[s][t]]
            if not _surjects(cols, target[s], dom):
                return False
        return True

    _prune(blocks, is_approx)
    return _assemble_right(X, blocks, types)


def _shapes(X, types):
    """Per target type s, the ``{deg: (rows, cols)}`` of a chain map X -> types[s]."""
    out = []
    for T in types:
        degs = sorted(set(X.terms) | set(T.terms))
        out.append({i: (T.dim_at(i), X.dim_at(i)) for i in degs})
    return out


def _shapes_to(types, X):
    """Per source type s, the ``{deg: (rows, cols)}`` of a chain map types[s] -> X."""
    out = []
    for T in types:
        degs = sorted(set(T.terms) | set(X.terms))
        out.append({i: (X.dim_at(i), T.dim_at(i)) for i in degs})
    return out


def _flat(fmap, src, tgt, shape, dom):
    degs = sorted(shape)
    out = []
    for i in degs:
        rows, cols = shape[i]
        m = fmap.get(i)
        for b in range(cols):
            for a in range(rows):
                out.append(m[a][b] if (m and rows and cols) else dom.zero())
    return out


def _assemble_left(X, blocks, types):
    A = X.algebra
    dom = A.domain
    if not blocks:
        E = PComplex(A, {}, {})
        return {}, E
    E, placements = direct_sum_complexes([types[t] for (t, _g) in blocks])
    fmap = {}
    degs = sorted(E.terms)
    for i in degs:
        rows, cols = E.dim_at(i), X.dim_at(i)
        M = lm.zeros(rows, cols, dom)
        for j, (t, g) in enumerate(blocks):
            off, dh = placements[j][i]
            gi = g.get(i)
            if gi and dh and cols:
                for a in range(dh):
                    for b in range(cols):
                        M[off + a][b] = gi[a][b]
        fmap[i] = M
    return fmap, E


def _assemble_right(X, blocks, types):
    A = X.algebra
    dom = A.domain
    if not blocks:
        E = PComplex(A, {}, {})
        return {}, E
    E, placements = direct_sum_complexes([types[t] for (t, _g) in blocks])
    gmap = {}
    degs = sorted(E.terms)
    for i in degs:
        rows, cols = X.dim_at(i), E.dim_at(i)
        M = lm.zeros(rows, cols, dom)
        for j, (t, g) in enumerate(blocks):
            off, dh = placements[j][i]
            gi = g.get(i)
            if gi and rows and dh:
                for a in range(rows):
                    for b in range(dh):
                        M[a][off + b] = gi[a][b]
        gmap[i] = M
    return gmap, E


# --------------------------------------------------------------------------- #
# the two mutation branches at a summand
# --------------------------------------------------------------------------- #
def left_mutation_summand(Xk, U):
    """The left-mutation new summand ``cone(f)`` for the minimal left add(U)-approximation
    ``f: Xk -> E``."""
    fmap, E = min_left_approx(Xk, U)
    return cone(fmap, Xk, E)


def right_mutation_summand(Xk, U):
    """The right-mutation new summand ``cocone(g) = cone(g)[-1]`` for the minimal right
    add(U)-approximation ``g: E -> Xk``."""
    gmap, E = min_right_approx(Xk, U)
    return shift(cone(gmap, E, Xk), -1)
