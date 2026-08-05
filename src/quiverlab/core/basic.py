"""Basic-ization + Gabriel-quiver recovery of a structure-constant Algebra (Plan 44 / C7).

Exact Wedderburn over char 0 or char > dim A ONLY (the trace-form radical -- the sole
exact radical route in the library -- is rigorous exactly there; Dickson / Cohen-Ivanyos-
Wales, the same bound ``decompose.py`` uses). Pipeline:

  1. ``rad A`` = nullspace of the REGULAR-representation trace form ``G[i][j] = tr(L_i L_j)``.
  2. ``S = A / rad``, a semisimple structure-constant algebra (rad is a two-sided
     nilpotent ideal; ``check=True`` on ``S`` self-certifies the quotient).
  3. a complete set of primitive orthogonal idempotents of ``S``, by splitting-element
     refinement (min-poly + coprime CRT idempotent on each corner). A non-split
     division-algebra block (corner dim > 1 with no splitting element) refuses loudly.
  4. lift each idempotent to ``A``, orthogonally, by the Newton iteration
     ``e <- 3e^2 - 2e^3`` (terminates: rad nilpotent, the defect ``e^2 - e`` is squared
     each step) within the shrinking corner ``g A g``, ``g <- g - e``.
  5. group by the Schur corner-dim test ``dim_k(ebar_i S ebar_j) >= 1`` (== ``Ae_i ~ Ae_j``).

Then ``basic_algebra`` = ``eAe`` (one idempotent per class); ``gabriel_quiver`` reads
arrows off ``rad/rad^2`` of the basic algebra (``i -> j`` from ``e_i (rad/rad^2) e_j`` --
path composition is LEFT-TO-RIGHT, so an element of ``e_i A e_j`` is a path ``i -> j``,
ASS); ``presented_form`` extracts ``kQ/I`` by the length-lex kernel enumeration (the
``families/_present.py`` idiom) and certifies ``dim(kQ/I) == dim(basic)`` PLUS a
multiplicativity check of the reconstructed iso ``kQ/I -> B``. Vertices are ordered by a
topological sort of the Ext-quiver (sources first) so an acyclic recovery matches the
conventionally-ordered original Cartan matrix exactly (the ``kA2`` / ``kA3`` arbiters).

``recognizers.is_basic`` stays unchanged (it returns True for every ``kQ/I`` -- a presented
algebra IS basic); this surface handles the presentation-LESS / non-basic
structure-constant case (``End(M)``, ``M_n(k)``, ...). Float-free; loud off char-scope.
"""
from __future__ import annotations

from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import solve
from quiverlab.modules import linalg_mod as lm


# --------------------------------------------------------------------------- #
# Small exact polynomial toolkit over the Domain (ascending coefficient lists).
# Reuses decompose's _poly_mul / _poly_pow; adds Euclidean division + ext-gcd for
# the CRT splitting idempotent.
# --------------------------------------------------------------------------- #
def _pdeg(p, dom):
    d = len(p) - 1
    while d >= 0 and dom.is_zero(p[d]):
        d -= 1
    return d


def _ptrim(p, dom):
    d = _pdeg(p, dom)
    return p[:d + 1] if d >= 0 else []


def _pscale(p, c, dom):
    return [dom.mul(c, x) for x in p]


def _psub(a, b, dom):
    n = max(len(a), len(b))
    z = dom.zero()
    out = [dom.sub(a[i] if i < len(a) else z, b[i] if i < len(b) else z)
           for i in range(n)]
    return _ptrim(out, dom)


def _pdivmod(a, b, dom):
    from quiverlab.modules.decompose import _poly_mul
    a = _ptrim(a, dom)
    b = _ptrim(b, dom)
    db = _pdeg(b, dom)
    if db < 0:
        raise ZeroDivisionError("polynomial division by zero")
    binv = dom.inv(b[db])
    q = []
    r = list(a)
    while True:
        dr = _pdeg(r, dom)
        if dr < db:
            break
        coef = dom.mul(r[dr], binv)
        shift = dr - db
        if shift >= len(q):
            q += [dom.zero()] * (shift - len(q) + 1)
        q[shift] = coef
        term = [dom.zero()] * shift + _pscale(b, coef, dom)
        r = _psub(r, term, dom)
    return _ptrim(q, dom), _ptrim(r, dom)


def _pgcdext(a, b, dom):
    """(g, u, v) with u*a + v*b = g, g monic (extended Euclid over k[x])."""
    from quiverlab.modules.decompose import _poly_mul
    old_r, r = _ptrim(a, dom), _ptrim(b, dom)
    old_s, s = [dom.one()], []
    old_t, t = [], [dom.one()]
    while _pdeg(r, dom) >= 0:
        q, _rem = _pdivmod(old_r, r, dom)
        old_r, r = r, _psub(old_r, _poly_mul(q, r, dom), dom)
        old_s, s = s, _psub(old_s, _poly_mul(q, s, dom), dom)
        old_t, t = t, _psub(old_t, _poly_mul(q, t, dom), dom)
    dg = _pdeg(old_r, dom)
    lead_inv = dom.inv(old_r[dg])
    return (_pscale(old_r, lead_inv, dom),
            _pscale(old_s, lead_inv, dom),
            _pscale(old_t, lead_inv, dom))


# --------------------------------------------------------------------------- #
# 1. radical of the regular trace form
# --------------------------------------------------------------------------- #
def _left_mult(A, i):
    """The matrix ``L_i`` of left multiplication by basis element ``b_i``: column ``k`` is
    ``b_i * b_k = A.T[i][k]``."""
    n = A.dim
    return [[A.T[i][k][t] for k in range(n)] for t in range(n)]


def _radical_coords(A):
    """``rad A`` as coordinate column-vectors: nullspace of the regular-rep trace form
    ``G[i][j] = tr(L_i L_j)``. Rigorous iff ``char 0`` or ``char > dim A`` -- else loud."""
    dom = A.domain
    n = A.dim
    char = dom.characteristic
    if not (char == 0 or char > n):
        raise QuiverlabError(
            f"basic-ization: the trace-form radical is unreliable in characteristic "
            f"{char} <= dim A = {n}", hint="work over QQ or a prime > dim A")
    L = [_left_mult(A, i) for i in range(n)]
    G = lm.zeros(n, n, dom)
    for i in range(n):
        for j in range(n):
            prod = lm.matmul(L[i], L[j], dom)
            s = dom.zero()
            for t in range(n):
                s = dom.add(s, prod[t][t])
            G[i][j] = s
    return lm.kernel_columns(G, dom)                     # [] when A is semisimple


# --------------------------------------------------------------------------- #
# 2. semisimple quotient S = A/rad
# --------------------------------------------------------------------------- #
def _combine(vectors, coeffs, dom):
    n = len(vectors[0]) if vectors else 0
    out = [dom.zero()] * n
    for c, v in zip(coeffs, vectors):
        if dom.is_zero(c):
            continue
        for t in range(n):
            out[t] = dom.add(out[t], dom.mul(c, v[t]))
    return out


def _semisimple_quotient(A, rad_cols):
    """``S = A/rad`` as a structure-constant Algebra, plus the complement basis ``comp``
    (an ``A``-representative of each ``S``-basis vector) and the projection ``proj`` (an
    ``A``-vector -> ``S``-coordinates). ``S`` is built with ``check=True`` -- its
    associativity/unit validation IS the certificate that ``rad`` is a two-sided ideal."""
    dom = A.domain
    n = A.dim
    ident = lm.identity(n, dom)
    std = [lm.col(ident, j) for j in range(n)]
    comp_idx = lm.independent_modulo(std, rad_cols, dom)         # complement of rad
    comp = [std[j] for j in comp_idx]
    m = len(comp)
    P = lm.cols_to_matrix(comp + [list(c) for c in rad_cols])   # [comp | rad] basis of A

    def proj(vec):
        x = solve(P, list(vec), dom)
        if x is None:
            raise QuiverlabError("basic-ization: A-basis [comp | rad] is singular (bug)")
        return x[:m]

    T = [[[proj(A.multiply(comp[a], comp[b]))[c] for c in range(m)]
          for b in range(m)] for a in range(m)]
    unit_S = proj(list(A.unit))
    S = Algebra.from_structure_constants(T, list(unit_S), field=dom, check=True)
    return S, comp, proj


# --------------------------------------------------------------------------- #
# 3. primitive idempotents of the semisimple quotient (splitting-element refinement)
# --------------------------------------------------------------------------- #
def _corner_space(S, e, dom):
    """A basis (list of independent ``S``-vectors) of the corner ``e S e``."""
    gens = [S.multiply(e, S.multiply(S._basis_vec(k), e)) for k in range(S.dim)]
    G = lm.cols_to_matrix(gens)
    piv = lm.column_space_pivots(G, dom)
    return [gens[j] for j in piv]


def _element_min_poly(S, s, e, dom):
    """Minimal polynomial of ``s`` in the corner algebra with unit ``e`` (ascending monic
    coeffs; ``s^0 = e``)."""
    powers = [list(e)]                                   # s^0 = e
    cur = list(e)
    for _ in range(S.dim + 1):
        cur = S.multiply(s, cur)                         # s^{k+1}
        B = lm.cols_to_matrix(powers)
        x = solve(B, list(cur), dom)
        if x is not None:                                # s^d = sum x_i s^i
            return [dom.neg(c) for c in x] + [dom.one()]
        powers.append(list(cur))
    raise QuiverlabError("basic-ization: corner min poly exceeded the corner dimension "
                         "(bug in the exact linear-algebra kernel)")


def _eval_poly_at(S, coeffs, s, e, dom):
    """``p(s)`` for the ascending polynomial ``coeffs``, evaluated in the corner with unit
    ``e`` (Horner: ``acc <- s*acc + c*e``)."""
    n = S.dim
    acc = [dom.zero()] * n
    for c in reversed(coeffs):
        acc = S.multiply(s, acc)
        acc = [dom.add(acc[t], dom.mul(c, e[t])) for t in range(n)]
    return acc


def _split_idempotent(S, e, dom):
    """Split ``e = e1 + e2`` into orthogonal idempotents via a splitting element of the
    corner ``eSe`` whose corner min poly has >= 2 coprime factors (the CRT idempotent), or
    ``None`` when the corner is a division algebra (dim 1). A corner of dim > 1 with no
    splitting element among the corner basis + pairwise sums is a NON-split division-algebra
    block -- refused loudly (honest scope)."""
    from quiverlab.modules.decompose import (_factor_min_poly, _factoring_supported,
                                             _poly_mul, _poly_pow)
    if not _factoring_supported(dom):
        raise QuiverlabError(
            "basic-ization: cannot factor over this domain within the certified budget",
            hint="supported: GF(p), QQ, and CC number fields (as in decompose)")
    corner = _corner_space(S, e, dom)
    cdim = len(corner)
    if cdim <= 1:
        return None                                      # e primitive (division algebra dim 1)
    cands = list(corner)
    for i in range(cdim):                                # pairwise sums reach non-basis splitters
        for j in range(i + 1, cdim):
            cands.append([dom.add(corner[i][t], corner[j][t]) for t in range(S.dim)])
    for s in cands:
        m = _element_min_poly(S, s, e, dom)
        if len(m) <= 2:                                  # degree <= 1: scalar * e, no split
            continue
        factors = _factor_min_poly(m, dom)
        if len(factors) < 2:                             # prime power: no coprime split
            continue
        f = _poly_pow(factors[0][0], factors[0][1], dom)
        g = [dom.one()]
        for fc, mult in factors[1:]:
            g = _poly_mul(g, _poly_pow(fc, mult, dom), dom)
        _gg, _u, v = _pgcdext(f, g, dom)                 # u f + v g = 1 (f, g coprime)
        e1 = _eval_poly_at(S, _poly_mul(v, g, dom), s, e, dom)   # (v g)(s): f-primary idempotent
        e2 = [dom.sub(e[t], e1[t]) for t in range(S.dim)]
        return e1, e2
    raise QuiverlabError(
        f"basic-ization: a Wedderburn block is a non-split division algebra (corner "
        f"dim {cdim} > 1, no splitting element) -- out of scope",
        hint="basic-ization certifies split semisimple quotients; a genuine "
             "division-algebra block (e.g. a quaternion corner over QQ) is refused")


def _semisimple_primitive_idempotents(S):
    """A complete set of primitive orthogonal idempotents of the SEMISIMPLE ``S`` (as
    ``S``-coordinate vectors), by refining ``1_S`` via splitting elements. Each split
    strictly lowers the corner dimension, so this terminates."""
    from quiverlab.modules.decompose import _factoring_supported
    if not _factoring_supported(S.domain):
        raise QuiverlabError(
            "basic-ization: cannot factor over this domain within the certified budget",
            hint="supported: GF(p), QQ, and CC number fields (as in decompose)")
    out = []
    stack = [list(S.unit)]
    while stack:
        e = stack.pop()
        split = _split_idempotent(S, e, S.domain)
        if split is None:
            out.append(e)
        else:
            stack.extend(split)
    return out


# --------------------------------------------------------------------------- #
# 4. orthogonal lift to A + completeness certificate
# --------------------------------------------------------------------------- #
def _newton_lift(A, x):
    """Lift a near-idempotent ``x`` (``x^2 == x`` mod rad) to an exact idempotent of ``A``
    via ``e <- 3e^2 - 2e^3``, iterated until ``A.multiply(e, e) == e``. Terminates: rad is
    nilpotent, and the defect ``e^2 - e`` is squared each step."""
    dom = A.domain
    three, two = dom.coerce(3), dom.coerce(2)
    e = list(x)
    for _ in range(A.dim + 1):                           # rad nilpotency index <= dim
        if A.multiply(e, e) == e:
            return e
        e2 = A.multiply(e, e)
        e3 = A.multiply(e2, e)
        e = [dom.sub(dom.mul(three, e2[t]), dom.mul(two, e3[t])) for t in range(A.dim)]
    raise QuiverlabError("basic-ization: Newton idempotent lift did not converge",
                         hint="the radical was not nilpotent -- please report this algebra")


def _assert_complete_orthogonal(A, idems):
    dom = A.domain
    total = [dom.zero()] * A.dim
    for i, ei in enumerate(idems):
        if A.multiply(ei, ei) != ei:
            raise QuiverlabError("basic-ization: a lifted element is not idempotent (bug)")
        for j, ej in enumerate(idems):
            if i != j and any(not dom.is_zero(x) for x in A.multiply(ei, ej)):
                raise QuiverlabError(
                    "basic-ization: lifted idempotents are not orthogonal (bug)")
        total = [dom.add(total[t], ei[t]) for t in range(A.dim)]
    if total != list(A.unit):
        raise QuiverlabError(
            "basic-ization: lifted idempotents do not sum to the unit (bug)")


def primitive_idempotents(A):
    """A complete set of orthogonal PRIMITIVE idempotents of ``A`` (coordinate vectors),
    summing to ``A.unit``. Char 0 / char > dim only, else loud (the trace-form radical)."""
    rad = _radical_coords(A)                             # loud off char-scope
    S, comp, _proj = _semisimple_quotient(A, rad)
    prim_S = _semisimple_primitive_idempotents(S)
    dom = A.domain
    collected = []
    g = list(A.unit)
    for sc in prim_S:
        x = _combine(comp, sc, dom)                      # S-coords -> A-representative
        y = A.multiply(g, A.multiply(x, g))              # corner preimage g x g
        e = _newton_lift(A, y)
        collected.append(e)
        g = [dom.sub(g[t], e[t]) for t in range(A.dim)]  # remaining corner
    _assert_complete_orthogonal(A, collected)
    return collected


# --------------------------------------------------------------------------- #
# 5. iso classes (Schur corner-dim test)
# --------------------------------------------------------------------------- #
def _corner_nonzero(S, bi, bj, dom):
    for k in range(S.dim):
        v = S.multiply(bi, S.multiply(S._basis_vec(k), bj))
        if any(not dom.is_zero(x) for x in v):
            return True
    return False


def idempotent_classes(A):
    """Partition ``primitive_idempotents(A)`` by ``Ae_i ~ Ae_j`` (returned as index groups),
    certified by the Schur corner-dim test ``dim_k(ebar_i S ebar_j) >= 1`` in the semisimple
    quotient ``S = A/rad`` (same Wedderburn block <=> isomorphic projectives)."""
    dom = A.domain
    rad = _radical_coords(A)
    S, _comp, proj = _semisimple_quotient(A, rad)
    idems = primitive_idempotents(A)
    bars = [proj(e) for e in idems]
    n = len(idems)
    groups = []
    assigned = [False] * n
    for i in range(n):
        if assigned[i]:
            continue
        cls = [i]
        assigned[i] = True
        for j in range(i + 1, n):
            if not assigned[j] and _corner_nonzero(S, bars[i], bars[j], dom):
                cls.append(j)
                assigned[j] = True
        groups.append(cls)
    return groups


# --------------------------------------------------------------------------- #
# 6. corner algebra eAe
# --------------------------------------------------------------------------- #
def _corner_algebra(A, e):
    """``e A e`` as a structure-constant Algebra (unit ``e``; ``check=True`` self-certifies)."""
    dom = A.domain
    gens = [A.multiply(e, A.multiply(A._basis_vec(k), e)) for k in range(A.dim)]
    piv = lm.column_space_pivots(lm.cols_to_matrix(gens), dom)
    basis = [gens[j] for j in piv]
    m = len(basis)
    B = lm.cols_to_matrix(basis) if basis else lm.zeros(A.dim, 0, dom)
    T = []
    for a in range(m):
        row = []
        for b in range(m):
            x = solve(B, A.multiply(basis[a], basis[b]), dom)
            if x is None:
                raise QuiverlabError("basic-ization: corner product left eAe (bug)")
            row.append(x)
        T.append(row)
    unit_e = solve(B, list(e), dom)
    return Algebra.from_structure_constants(T, list(unit_e), field=dom, check=True)


def basic_algebra(A):
    """The basic algebra ``eAe`` of ``A`` (``e`` = one primitive idempotent per iso class),
    as a structure-constant Algebra. Morita-equivalent to ``A``, presentation-less."""
    idems = primitive_idempotents(A)
    classes = idempotent_classes(A)
    dom = A.domain
    e = [dom.zero()] * A.dim
    for cls in classes:
        rep = idems[cls[0]]
        e = [dom.add(e[t], rep[t]) for t in range(A.dim)]
    return _corner_algebra(A, e)


# --------------------------------------------------------------------------- #
# 7. Gabriel quiver + presented form
# --------------------------------------------------------------------------- #
def _rad_square(B, radB):
    """A basis of ``rad^2 B`` = span of ``{r_i r_j : r_i, r_j in rad}``."""
    dom = B.domain
    prods = [B.multiply(ri, rj) for ri in radB for rj in radB]
    if not prods:
        return []
    piv = lm.column_space_pivots(lm.cols_to_matrix(prods), dom)
    return [prods[j] for j in piv]


def _loewy_length(B, radB):
    """Loewy length of ``B`` = nilpotency index of ``rad B`` (smallest ``L`` with
    ``rad^L = 0``), from the radical filtration -- ``B`` may be presentation-less (no
    ``e_``-labels), so ``invariants.scalar.loewy_length`` cannot serve here."""
    dom = B.domain

    def reduce(vecs):
        vs = [v for v in vecs if any(not dom.is_zero(x) for x in v)]
        if not vs:
            return []
        piv = lm.column_space_pivots(lm.cols_to_matrix(vs), dom)
        return [vs[j] for j in piv]

    cur = reduce(radB)
    L = 1
    while cur:
        nxt = reduce([B.multiply(c, r) for c in cur for r in radB])
        if not nxt:
            return L + 1
        cur = nxt
        L += 1
        if L > B.dim + 1:                                # safety: no chain longer than dim
            return L + 1
    return 1                                             # semisimple: rad = 0


def _corner_arrow_reps(B, ei, ej, radB, rad2, dom):
    """Representatives (``B``-vectors) of a basis of ``e_i (rad/rad^2) e_j`` -- the arrows
    ``i -> j`` of the Gabriel quiver, each an ``e_i (rad) e_j`` element independent modulo
    ``e_i (rad^2) e_j``."""
    rc = [B.multiply(ei, B.multiply(r, ej)) for r in radB]
    r2c = [B.multiply(ei, B.multiply(r, ej)) for r in rad2]
    keep = lm.independent_modulo(rc, r2c, dom)
    return [rc[k] for k in keep]


def _topological_order(n, edges):
    """A deterministic topological order of ``0..n-1`` under ``edges`` (list of ``(i, j)``
    meaning ``i -> j``), sources first; identity fallback if there is a cycle (loops
    ignored). Makes an acyclic recovery's Cartan upper-triangular = the conventionally
    ordered original."""
    import heapq
    adj = {i: [] for i in range(n)}
    indeg = [0] * n
    seen = set()
    for (i, j) in edges:
        if i != j and (i, j) not in seen:
            seen.add((i, j))
            adj[i].append(j)
            indeg[j] += 1
    avail = [v for v in range(n) if indeg[v] == 0]
    heapq.heapify(avail)
    order = []
    while avail:
        v = heapq.heappop(avail)
        order.append(v)
        for w in adj[v]:
            indeg[w] -= 1
            if indeg[w] == 0:
                heapq.heappush(avail, w)
    return order if len(order) == n else list(range(n))


def _gabriel_data(A):
    """Shared basic-algebra + arrow data: ``(B, fi, arrow_reps, order)`` where ``fi`` are
    the primitive idempotents of the basic algebra ``B``, ``arrow_reps`` a list of
    ``(i, j, rep)`` (fi-index endpoints + a ``rad/rad^2`` representative), and ``order`` a
    topological vertex permutation of ``0..len(fi)-1``."""
    B = basic_algebra(A)
    dom = B.domain
    fi = primitive_idempotents(B)                        # basic => each in its own class
    radB = _radical_coords(B)
    rad2 = _rad_square(B, radB)
    arrow_reps = []
    for i in range(len(fi)):
        for j in range(len(fi)):
            for rep in _corner_arrow_reps(B, fi[i], fi[j], radB, rad2, dom):
                arrow_reps.append((i, j, rep))
    order = _topological_order(len(fi), [(i, j) for (i, j, _r) in arrow_reps])
    return B, fi, arrow_reps, order, radB


def gabriel_quiver(A):
    """The Gabriel (Ext) quiver of ``A``: vertices ``1..r`` (topologically ordered), arrows
    ``i -> j`` one per basis vector of ``e_i (rad B / rad^2 B) e_j`` in the basic algebra
    ``B`` (path composition left-to-right)."""
    from quiverlab.combinat.quiver import Quiver
    _B, fi, arrow_reps, order, _radB = _gabriel_data(A)
    pos = {old: new for new, old in enumerate(order)}
    verts = list(range(1, len(fi) + 1))
    arrows = {}
    a = 0
    for (i, j, _rep) in arrow_reps:
        a += 1
        arrows[f"a{a}"] = (pos[i] + 1, pos[j] + 1)
    return Quiver(verts, arrows)


def _verify_iso(Bnew, B, imgs, dom):
    """Certify the reconstructed ``pi: Bnew=kQ/I -> B`` (basis label -> ``imgs`` B-vector) is
    an algebra ISO: surjective (image spans ``B``) and multiplicative on every basis pair.
    Combined with ``dim Bnew == dim B`` this is the recovery certificate beyond dimension."""
    if lm.mat_rank(lm.cols_to_matrix(imgs), dom) != B.dim:
        raise QuiverlabError(
            "Gabriel recovery: the reconstructed pi is not surjective onto B (not an iso)")
    for a in range(Bnew.dim):
        for b in range(Bnew.dim):
            prod = Bnew.multiply(Bnew._basis_vec(a), Bnew._basis_vec(b))
            pi_prod = [dom.zero()] * B.dim
            for k, c in enumerate(prod):
                if dom.is_zero(c):
                    continue
                ik = imgs[k]
                for t in range(B.dim):
                    pi_prod[t] = dom.add(pi_prod[t], dom.mul(c, ik[t]))
            direct = B.multiply(imgs[a], imgs[b])
            if pi_prod != direct:
                raise QuiverlabError(
                    "Gabriel recovery: the reconstructed pi is not multiplicative "
                    "(kQ/I -> B is not an algebra map)")


def presented_form(A):
    """A genuine ``kQ/I`` presentation of the basic algebra of ``A`` (``Q = gabriel_quiver``),
    certified per instance by ``dim(kQ/I) == dim(basic_algebra(A))`` PLUS a multiplicativity
    check of the reconstructed iso ``kQ/I -> B``. Loud refusal off char-scope or on a
    non-split division-algebra block."""
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.families._present import present_from_pi
    B, fi, arrow_reps, order, radB = _gabriel_data(A)
    dom = B.domain
    pos = {old: new for new, old in enumerate(order)}
    verts = list(range(1, len(fi) + 1))
    arrows = {}
    img = {}
    vertex_idem = {}
    for new, old in enumerate(order):
        vertex_idem[new + 1] = fi[old]                   # synthetic vertex -> basic idempotent
    a = 0
    for (i, j, rep) in arrow_reps:
        a += 1
        name = f"a{a}"
        arrows[name] = (pos[i] + 1, pos[j] + 1)
        img[name] = rep
    Q = Quiver(verts, arrows)
    base_bound = _loewy_length(B, radB) + 2
    Bnew = present_from_pi(Q, img, B, dom, B.dim, base_bound, citations=("assem_book",))

    # multiplicativity certificate: reconstruct pi on Bnew's normal-form basis and verify.
    def pimap(label):
        if label.startswith("e_"):
            return vertex_idem[int(label[2:])]
        acc = None
        for a_lbl in label.split("*"):
            acc = img[a_lbl] if acc is None else B.multiply(acc, img[a_lbl])
        return acc

    imgs = [pimap(lab) for lab in Bnew.basis_labels]
    _verify_iso(Bnew, B, imgs, dom)
    return Bnew
