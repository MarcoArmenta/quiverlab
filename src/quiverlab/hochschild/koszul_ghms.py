"""Koszul kernels ``K_n`` and the GHMS comultiplicative minimal bimodule resolution
(Plan 75 / R10; Green-Hartman-Marcos-Solberg 2005, arXiv:math/0508177).

For a Koszul ``A = kQ/I`` write ``S = k^{Q_0}``, ``V = kQ_1`` (arrows) and ``R subset
V (x)_S V`` for the quadratic relation space. The **Koszul kernels** are

    K_0 = S,  K_1 = V,  K_n = intersection over i of V^{(x)i} (x) R (x) V^{(x)(n-2-i)}.

GHMS: the MINIMAL graded projective ``A^e``-resolution of ``A`` is
``P_n = A (x)_S K_n (x)_S A``, so ``rank_{A^e} P_n = dim_k K_n``, and that dimension is
the ``n``-th coefficient of the Koszul-dual Hilbert series -- equivalently the total
graded Betti number ``sum_{i,j} dim Ext^n(S_i, S_j)``. Closed-form terms, no syzygy search,
and it runs over ANY exact ``Domain`` (the syzygy engine is GF(p)-only).

**The recursion actually used.** Rather than intersecting all ``n-1`` summands, this
computes

    K_n = (V (x) K_{n-1})  intersect  (K_{n-1} (x) V),

which is EXACT and much cheaper: ``V (x) K_{n-1}`` is the intersection over ``i >= 1`` and
``K_{n-1} (x) V`` the intersection over ``i <= n-3``, so their meet is the full
intersection over ``0 <= i <= n-2``. One subspace meet per degree instead of ``n-1``.

Coordinates are on the composable paths of length ``n`` (a path basis already encodes the
corner grading ``e_v ... e_w``, so nothing extra is needed to respect it).
"""
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import nullspace, rref, solve

_GHMS_CITATIONS = ("green_hartman_marcos_solberg", "priddy", "froberg_koszul")


def _paths_of_length(quiver, n):
    """Composable arrow words of length ``n``, in a deterministic order."""
    arrows = quiver.arrows
    names = sorted(arrows)
    if n == 0:
        return [()]
    words = [(a,) for a in names]
    for _ in range(n - 1):
        nxt = []
        for w in words:
            tgt = arrows[w[-1]][1]
            for a in names:
                if arrows[a][0] == tgt:
                    nxt.append(w + (a,))
        words = nxt
        if not words:
            break
    return words


def _relation_space(A):
    """A basis of ``R subset V (x)_S V`` in coordinates on the length-2 path basis."""
    dom = A.domain
    basis = _paths_of_length(A.quiver, 2)
    index = {w: i for i, w in enumerate(basis)}
    rows = []
    for rel in (A.relations or []):
        vec = [dom.zero()] * len(basis)
        hit = False
        for c, w in rel.terms:
            w = tuple(w)
            if w not in index:
                raise QuiverlabError(
                    f"relation term {w!r} is not a length-2 path: the GHMS route needs a "
                    "QUADRATIC presentation",
                    hint="check is_quadratic(A) before calling the Koszul kernels")
            vec[index[w]] = dom.coerce(c) if not hasattr(c, "parent") else c
            hit = True
        if hit:
            rows.append(vec)
    return rref(rows, dom)[0] if rows else []


def _intersect(spanA, spanB, dom, ambient):
    """A basis of ``span(A) intersect span(B)`` inside ``k^ambient``.

    Solves ``sum_j c_j A_j = sum_k d_k B_k`` by taking the nullspace of the ``ambient x
    (p+q)`` matrix whose columns are the ``A_j`` and ``-B_k``; each nullspace vector's
    ``A``-part reconstructs a common element.
    """
    if not spanA or not spanB:
        return []
    p, q = len(spanA), len(spanB)
    rows = []
    for r in range(ambient):
        rows.append([spanA[j][r] for j in range(p)]
                    + [dom.neg(spanB[k][r]) for k in range(q)])
    out = []
    for c in nullspace(rows, dom):
        vec = [dom.zero()] * ambient
        for j in range(p):
            if not dom.is_zero(c[j]):
                for r in range(ambient):
                    vec[r] = dom.add(vec[r], dom.mul(c[j], spanA[j][r]))
        if any(not dom.is_zero(x) for x in vec):
            out.append(vec)
    return rref(out, dom)[0] if out else []


def _tensor_left(quiver, dom, K, paths_prev, paths_now):
    """``V (x) K``: prepend each arrow to each basis vector of ``K``."""
    arrows = quiver.arrows
    idx_prev = {w: i for i, w in enumerate(paths_prev)}
    idx_now = {w: i for i, w in enumerate(paths_now)}
    out = []
    for a in sorted(arrows):
        tgt = arrows[a][1]
        for k in K:
            vec = [dom.zero()] * len(paths_now)
            hit = False
            for w, i in idx_prev.items():
                if dom.is_zero(k[i]):
                    continue
                if (arrows[w[0]][0] if w else None) != tgt:
                    continue
                vec[idx_now[(a,) + w]] = k[i]
                hit = True
            if hit:
                out.append(vec)
    return out


def _tensor_right(quiver, dom, K, paths_prev, paths_now):
    """``K (x) V``: append each arrow to each basis vector of ``K``."""
    arrows = quiver.arrows
    idx_prev = {w: i for i, w in enumerate(paths_prev)}
    idx_now = {w: i for i, w in enumerate(paths_now)}
    out = []
    for a in sorted(arrows):
        src = arrows[a][0]
        for k in K:
            vec = [dom.zero()] * len(paths_now)
            hit = False
            for w, i in idx_prev.items():
                if dom.is_zero(k[i]):
                    continue
                if (arrows[w[-1]][1] if w else None) != src:
                    continue
                vec[idx_now[w + (a,)]] = k[i]
                hit = True
            if hit:
                out.append(vec)
    return out


def koszul_kernels(A, top):
    """``[K_0, K_1, ..., K_top]`` as lists of coordinate vectors.

    ``K_0`` is indexed by the VERTICES (``dim K_0 = |Q_0|``, not 1 -- the multi-vertex
    point), ``K_1`` by the arrows, and ``K_n`` (``n >= 2``) by the composable length-``n``
    paths. ``dim K_n`` is the ``n``-th Koszul-dual Hilbert coefficient.

    Requires a quadratic presentation; the Koszulity VERDICT itself is not asserted here
    (the kernels are defined for any quadratic algebra) -- ``GHMSResolution`` is where the
    Koszul gate belongs, because it is the RESOLUTION claim that needs it.
    """
    from quiverlab.modules.koszul import _require_presentation, is_quadratic
    _require_presentation(A, "koszul_kernels")
    if not is_quadratic(A):
        raise QuiverlabError(
            "the Koszul kernels K_n are defined for a QUADRATIC presentation (ideal "
            "generated in degree 2)",
            hint="is_quadratic(A) is False -- a minimal relation has length != 2")
    dom = A.domain
    Q = A.quiver
    verts = sorted(Q.vertices, key=str)
    out = [[[dom.one() if i == j else dom.zero() for j in range(len(verts))]
            for i in range(len(verts))]]
    if top == 0:
        return out
    arrows_basis = _paths_of_length(Q, 1)
    out.append([[dom.one() if i == j else dom.zero() for j in range(len(arrows_basis))]
                for i in range(len(arrows_basis))])
    if top == 1:
        return out
    K = _relation_space(A)
    out.append(K)
    paths_prev = _paths_of_length(Q, 2)
    for n in range(3, top + 1):
        paths_now = _paths_of_length(Q, n)
        if not paths_now or not K:
            out.append([])
            K, paths_prev = [], paths_now
            continue
        left = _tensor_left(Q, dom, K, paths_prev, paths_now)
        right = _tensor_right(Q, dom, K, paths_prev, paths_now)
        K = _intersect(left, right, dom, len(paths_now))
        out.append(K)
        paths_prev = paths_now
    return out


def koszul_betti(A, top):
    """``[dim K_0, ..., dim K_top]`` -- the minimal bimodule ranks, i.e. the graded Betti
    numbers, i.e. the Koszul-dual Hilbert coefficients."""
    return [len(k) for k in koszul_kernels(A, top)]


# --------------------------------------------------------------------------- #
# the GHMS resolution
# --------------------------------------------------------------------------- #
def _koszul_gate(A, window=4):
    """The THREE-VALUED Koszulity gate (Plan-75 Scope gate 2).

    Order matters. The cheap ``g_quadratic_certificate`` (Priddy PBW) is tried FIRST and,
    when True, PROVES Koszulity -- build. Otherwise we must consult the full
    ``ext_algebra`` verdict, because a False from the PBW certificate is INCONCLUSIVE, not
    a disproof. Gating on ``g_quadratic`` alone would mislabel the genuinely-not-Koszul
    preprojective ``A_3`` as merely "uncertified", which is precisely the error the plan's
    critic caught: ``A_3``'s Ext algebra grows a NEW GENERATOR in degree 3, so it is not
    Koszul at all, and the refusal must say so.
    """
    from quiverlab.modules.koszul import g_quadratic_certificate
    try:
        if g_quadratic_certificate(A):
            return True, "G-quadratic (Priddy PBW): a quadratic Grobner basis proves Koszul"
    except QuiverlabError:
        pass
    E = A.ext_algebra(window)
    verdict = getattr(E, "koszul", None)
    if verdict is True:
        return True, "ext_algebra certifies Koszul"
    if verdict is False:
        obs = getattr(E, "koszul_obstruction", None)
        raise QuiverlabError(
            f"this algebra is NOT Koszul, so the GHMS resolution does not apply: the "
            f"Ext-algebra obstruction is {obs}",
            hint="the minimal bimodule resolution of a non-Koszul algebra is NOT "
                 "A (x)_S K_n (x)_S A -- the Koszul kernels can die while the true "
                 "resolution continues (preprojective A_3: kernels [3,4,3,0,0,0] vs true "
                 "ranks [3,4,3,3,4,3,3]). Use engine='auto'/'cs'/'bar' instead")
    raise QuiverlabError(
        "Koszulity is UNDECIDED for this algebra within the Ext window, so the GHMS "
        "resolution cannot be certified (the PBW certificate was inconclusive and the "
        "Ext-algebra verdict is None -- not a disproof, but not a licence either)",
        hint="raise the ext_algebra window, or use engine='auto'/'cs'/'bar'")


def _corner_of(A, word):
    """``(source, target)`` of an arrow word (a vertex pair)."""
    arrows = A.quiver.arrows
    return (arrows[word[0]][0], arrows[word[-1]][1])


def _split_by_corner(A, basis, paths):
    """Re-express a kernel basis as CORNER-HOMOGENEOUS vectors.

    ``K_n`` is corner-graded (it lives inside ``V^{(x)n}``, whose path basis is), but an
    ``rref`` basis can straddle corners. Restricting each vector to one corner's
    coordinates and re-reducing recovers a corner-homogeneous basis of the SAME space --
    which the bimodule terms ``e_v K_n e_w`` need.
    """
    dom = A.domain
    if not basis or not paths:
        return []
    by_corner = {}
    for idx, w in enumerate(paths):
        by_corner.setdefault(_corner_of(A, w), []).append(idx)
    out = []
    for corner, idxs in sorted(by_corner.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        rows = []
        for vec in basis:
            sub = [vec[i] for i in idxs]
            if any(not dom.is_zero(x) for x in sub):
                rows.append(sub)
        if not rows:
            continue
        for red in rref(rows, dom)[0]:
            full = [dom.zero()] * len(paths)
            for pos, i in enumerate(idxs):
                full[i] = red[pos]
            out.append((full, corner))
    return out


class GHMSResolution:
    """The GHMS comultiplicative minimal ``A^e``-resolution of a KOSZUL algebra.

    ``P_n = A (x)_S K_n (x)_S A`` with ``rank_{A^e} P_n = dim K_n``, and

        d_n(1 (x) w (x) 1) = sum_a x_a (x) w'_a (x) 1  -  (-1)^n sum_b 1 (x) w''_b (x) y_b

    where ``w = sum_a x_a (x) w'_a`` under ``K_n subset V (x) K_{n-1}`` and
    ``w = sum_b w''_b (x) y_b`` under ``K_n subset K_{n-1} (x) V``. In the path
    coordinates used here both splittings are READ OFF directly -- ``w'_a`` is the slice of
    ``w`` on words beginning with ``a``, ``w''_b`` the slice on words ending with ``b`` --
    so no solve is needed to find them; a solve is used only to express those slices in the
    ``K_{n-1}`` basis, and that solve is unique because a basis is independent (hence the
    differential is canonical/byte-reproducible by construction, with nothing to
    canonicalize away).

    THE SIGN IS ARBITRATED, NOT ASSUMED: ``assert_dd_zero`` pins it here, and Task II3's
    cross-engine anchors (GHMS == minimal syzygy engine over GF(p); GHMS == bar in the bar
    window) pin it again independently. A wrong sign fails at least one.
    """

    def __init__(self, A, top=6, window=4, _skip_gate=False):
        self.algebra = A
        self.top = top
        self.koszul_reason = None
        if not _skip_gate:
            _ok, self.koszul_reason = _koszul_gate(A, window=window)
        self.kernels = koszul_kernels(A, top)
        Q = A.quiver
        self._paths = [_paths_of_length(Q, n) for n in range(top + 1)]
        # corner-homogeneous bases; degree 0 and 1 are already corner-homogeneous
        self.basis = []
        for n, K in enumerate(self.kernels):
            if n <= 1:
                verts = sorted(Q.vertices, key=str)
                if n == 0:
                    self.basis.append([(vec, (v, v)) for vec, v in zip(K, verts)])
                else:
                    names = sorted(Q.arrows)
                    self.basis.append([(vec, (Q.arrows[a][0], Q.arrows[a][1]))
                                       for vec, a in zip(K, names)])
            else:
                self.basis.append(_split_by_corner(A, K, self._paths[n]))
        self._diffs = {}

    def rank(self, n):
        """``dim K_n`` -- the ``A^e``-rank of ``P_n``."""
        return len(self.basis[n]) if n < len(self.basis) else 0

    def term(self, n):
        """The corner types ``[(v, w), ...]`` of ``P_n``'s generators."""
        return [corner for _v, corner in self.basis[n]] if n < len(self.basis) else []

    # ------------------------------------------------------------ differential
    def differential(self, n):
        """``d_n`` as a matrix of ``A^e`` entries.

        ``D[i][j]`` is a list of ``(coeff, u, v)`` terms meaning the bimodule map
        ``a (x) w_j (x) b  |->  coeff * (u a) (x) w_i (x) (b v)``; ``u`` and ``v`` are
        coordinate vectors in ``A`` (an arrow, or the unit).
        """
        if n in self._diffs:
            return self._diffs[n]
        A = self.algebra
        dom = A.domain
        if n <= 0 or n > self.top:
            return []
        src, tgt = self.basis[n], self.basis[n - 1]
        if not src or not tgt:
            self._diffs[n] = []
            return []
        names = sorted(A.quiver.arrows)
        arrow_vec = {a: A._basis_vec(A.basis_labels.index(a)) for a in names}
        one = list(A.unit)
        if n == 1:
            # DEGREE 1 IS SPECIAL: the two splittings land in K_0 = S, which is indexed by
            # VERTICES, not by the (single, empty) length-0 path. Solving against the
            # path-coordinate basis here would be degenerate -- and it silently was, until
            # the multi-vertex cohomology exposed it (HH^0 of the diamond came out 3
            # instead of Z(A) = 1; homology had masked it because that quiver is acyclic,
            # so C_1 = 0 and the bad d_1 never contributed a rank).
            #   d_1(1 (x) a (x) 1) = a (x) [t(a)] (x) 1  -  1 (x) [s(a)] (x) a
            # (the minus is (-1)^1, the same augmentation-pinned sign as everywhere else,
            # and it is exactly what makes mu . d_1 = 0).
            verts = [c[0] for _v, c in self.basis[0]]
            D = [[[] for _ in range(len(src))] for _ in range(len(tgt))]
            for j, (_w, (sv, tv)) in enumerate(src):
                a_name = sorted(A.quiver.arrows)[j]
                av = arrow_vec[a_name]
                D[verts.index(tv)][j].append((dom.one(), av, one))
                D[verts.index(sv)][j].append((dom.neg(dom.one()), one, av))
            self._diffs[n] = D
            return D
        prev_paths, now_paths = self._paths[n - 1], self._paths[n]
        idx_prev = {w: i for i, w in enumerate(prev_paths)}
        tgt_vecs = [vec for vec, _c in tgt]
        # columns of the K_{n-1} basis, for the coordinate solve
        Mprev = [[tgt_vecs[j][r] for j in range(len(tgt_vecs))] for r in range(len(prev_paths))]
        names = sorted(A.quiver.arrows)
        arrow_vec = {a: A._basis_vec(A.basis_labels.index(a)) for a in names}
        one = list(A.unit)
        # THE SIGN IS (-1)^n on the RIGHT-splitting term, and it is pinned by the
        # AUGMENTATION, not by d.d = 0. At n = 1 the resolution must satisfy
        # mu . d_1 = 0, i.e. d_1(1(x)a(x)1) = a(x)1 - 1(x)a; (-1)^1 = -1 delivers that,
        # while -(-1)^n would give +1 and a d_1 with mu . d_1 = 2a != 0 (visible as
        # HH_0 = 2 instead of 3 on the exterior algebra). Both conventions ALTERNATE, so
        # d.d = 0 cannot tell them apart -- see the sign tests.
        sign = dom.one() if (n % 2 == 0) else dom.neg(dom.one())    # (-1)^n
        D = [[[] for _ in range(len(src))] for _ in range(len(tgt))]
        for j, (w, _corner) in enumerate(src):
            # LEFT splitting: slice by FIRST arrow
            for a in names:
                slice_vec = [dom.zero()] * len(prev_paths)
                hit = False
                for wi, word in enumerate(now_paths):
                    if word[0] == a and not dom.is_zero(w[wi]):
                        slice_vec[idx_prev[word[1:]]] = w[wi]
                        hit = True
                if not hit:
                    continue
                c = solve(Mprev, slice_vec, dom)
                if c is None:
                    raise QuiverlabError(
                        "K_n is not contained in V (x) K_{n-1}: the comultiplicative "
                        "structure fails, so this presentation is not Koszul in the way "
                        "the GHMS theorem requires",
                        hint="internal invariant -- please report this algebra")
                for i in range(len(tgt)):
                    if not dom.is_zero(c[i]):
                        D[i][j].append((c[i], arrow_vec[a], one))
            # RIGHT splitting: slice by LAST arrow
            for b in names:
                slice_vec = [dom.zero()] * len(prev_paths)
                hit = False
                for wi, word in enumerate(now_paths):
                    if word[-1] == b and not dom.is_zero(w[wi]):
                        slice_vec[idx_prev[word[:-1]]] = w[wi]
                        hit = True
                if not hit:
                    continue
                c = solve(Mprev, slice_vec, dom)
                if c is None:
                    raise QuiverlabError(
                        "K_n is not contained in K_{n-1} (x) V: the comultiplicative "
                        "structure fails (see the left-splitting note)",
                        hint="internal invariant -- please report this algebra")
                for i in range(len(tgt)):
                    if not dom.is_zero(c[i]):
                        D[i][j].append((dom.mul(sign, c[i]), one, arrow_vec[b]))
        self._diffs[n] = D
        return D

    # ------------------------------------------------------------- self-certs
    def _expand(self, terms):
        """Expand a list of ``(coeff, u, v)`` into ``A (x) A^op`` coordinates."""
        A, dom = self.algebra, self.algebra.domain
        m = A.dim
        acc = {}
        for c, u, v in terms:
            for i in range(m):
                if dom.is_zero(u[i]):
                    continue
                for j in range(m):
                    if dom.is_zero(v[j]):
                        continue
                    key = (i, j)
                    acc[key] = dom.add(acc.get(key, dom.zero()),
                                       dom.mul(c, dom.mul(u[i], v[j])))
        return {k: x for k, x in acc.items() if not dom.is_zero(x)}

    def assert_dd_zero(self, top=None):
        """``d_{n-1} . d_n = 0`` for every ``n`` in range -- the sign arbiter.

        Composition of ``(u, v)`` then ``(u', v')`` is ``(u' u, v v')``: the first map
        left-multiplies by ``u`` and right-multiplies by ``v``, the second then applies
        ``u'`` on the left and ``v'`` on the right, and left/right actions commute.
        """
        A, dom = self.algebra, self.algebra.domain
        top = self.top if top is None else top
        return self._dd_zero_range(A, dom, top)

    def _dd_zero_range(self, A, dom, top):
        for n in range(2, top + 1):
            dn, dn1 = self.differential(n), self.differential(n - 1)
            if not dn or not dn1:
                continue
            for k in range(len(dn[0])):
                for i in range(len(dn1)):
                    terms = []
                    for j in range(len(dn)):
                        for (c1, u1, v1) in dn[j][k]:
                            for (c2, u2, v2) in dn1[i][j]:
                                terms.append((dom.mul(c2, c1),
                                              A.multiply(u2, u1),
                                              A.multiply(v1, v2)))
                    bad = self._expand(terms)
                    if bad:
                        raise QuiverlabError(
                            f"GHMS d_{n-1} . d_{n} != 0 at entry ({i}, {k}): {len(bad)} "
                            "nonzero A (x) A^op coordinates -- the comultiplicative sign "
                            "convention is wrong",
                            hint="the sign is -(-1)^n on the RIGHT splitting term")
        return True


# --------------------------------------------------------------------------- #
# the HH collapse: engine="ghms"
# --------------------------------------------------------------------------- #
def _basis_corners(A):
    """``(source, target)`` of every basis element of ``A``, read off its path label."""
    arrows = A.quiver.arrows
    out = []
    for lab in A.basis_labels:
        s = str(lab)
        if s.startswith("e_"):
            v = s[2:]
            match = [x for x in A.quiver.vertices if str(x) == v]
            out.append((match[0], match[0]) if match else (v, v))
        else:
            parts = s.split("*")
            out.append((arrows[parts[0]][0], arrows[parts[-1]][1]))
    return out


def _corner_indices(corners, src, tgt):
    """Indices of the ``A``-basis spanning ``e_src A e_tgt``."""
    return [i for i, (s, t) in enumerate(corners) if s == src and t == tgt]


def _collapse_complex(res, top, *, side):
    """The corner-typed collapse of ``P_.`` (Plan-16 convention).

    ``side="hom"``: ``A (x)_{A^e} P_n``. A generator ``w`` with corner ``(v, w)`` gives the
    free bimodule ``A e_v (x) e_w A``, and ``A (x)_{A^e} (A e_v (x) e_w A) = e_w A e_v``, so
    a differential term ``(c, u, v)`` -- meaning ``a (x) w (x) b |-> c (u a) (x) w' (x)
    (b v)`` -- sends ``m |-> c * v m u``. That is the ``b . w . a`` order.

    ``side="coh"``: ``Hom_{A^e}(P_n, A) = e_v A e_w`` (the SWAPPED tag), and the same term
    contributes ``f |-> c * u f v`` -- the ``a . w . b`` order. The two orders are the
    covariance flip, not a free choice.

    Returns ``(dims, matrices)`` where ``matrices[n]`` is ``d_n : C_n -> C_{n-1}`` for
    ``hom`` and ``delta^n : C^n -> C^{n+1}`` for ``coh``.
    """
    A, dom = res.algebra, res.algebra.domain
    corners = _basis_corners(A)
    layout, dims = [], []
    for n in range(top + 2):
        slots, off = [], 0
        for j, (_vec, (v, w)) in enumerate(res.basis[n] if n < len(res.basis) else []):
            idxs = _corner_indices(corners, w, v) if side == "hom" else \
                _corner_indices(corners, v, w)
            slots.append((j, idxs, off))
            off += len(idxs)
        layout.append(slots)
        dims.append(off)

    mats = {}
    for n in range(1, top + 2):
        D = res.differential(n)
        if not D:
            mats[n] = []
            continue
        if side == "hom":
            # d_n : C_n -> C_{n-1}.  m lives in generator j's corner (degree n) and is
            # carried to generator i's corner (degree n-1) by  m |-> c * v m u.
            rows, cols = dims[n - 1], dims[n]
            row_slots = {i: (idxs, off) for i, idxs, off in layout[n - 1]}
            col_slots = {j: (idxs, off) for j, idxs, off in layout[n]}
        else:
            # delta^{n-1} : C^{n-1} -> C^n.  m lives in generator i's corner (degree n-1)
            # and is carried to generator j's corner (degree n) by  m |-> c * u m v.
            # NOTE the direction: cohomology runs UP, so the roles of i and j swap
            # relative to the homology assembly. Using the homology layout with the
            # a.w.b action would build a different linear map whose rank is NOT rank
            # delta (it silently gave HH^0 = 3 instead of Z(A) = 2 on the exterior algebra).
            rows, cols = dims[n], dims[n - 1]
            row_slots = {j: (idxs, off) for j, idxs, off in layout[n]}
            col_slots = {i: (idxs, off) for i, idxs, off in layout[n - 1]}
        if rows == 0 or cols == 0:
            mats[n] = []
            continue
        M = [[dom.zero()] * cols for _ in range(rows)]
        for i in range(len(D)):
            for j in range(len(D[0])):
                if not D[i][j]:
                    continue
                if side == "hom":
                    cidx, coff = col_slots[j]
                    ridx, roff = row_slots[i]
                else:
                    cidx, coff = col_slots[i]
                    ridx, roff = row_slots[j]
                rpos = {b: r for r, b in enumerate(ridx)}
                for col, b in enumerate(cidx):
                    m = A._basis_vec(b)
                    for (c, u, v) in D[i][j]:
                        img = A.multiply(v, A.multiply(m, u)) if side == "hom" else \
                            A.multiply(u, A.multiply(m, v))
                        for r, coeff in enumerate(img):
                            if dom.is_zero(coeff) or r not in rpos:
                                continue
                            cur = M[roff + rpos[r]][coff + col]
                            M[roff + rpos[r]][coff + col] = dom.add(cur, dom.mul(c, coeff))
        mats[n] = M
    return dims, mats


def _rank_of(M, dom):
    if not M or not M[0]:
        return 0
    from quiverlab.fields.linalg import rank as _r
    return _r(M, dom)


def ghms_homology_dims(A, top, window=4):
    """``[dim HH_0, ..., dim HH_top]`` via the GHMS resolution -- Domain-general and fast.

    ``dim HH_n = dim C_n - rank d_n - rank d_{n+1}``.
    """
    res = GHMSResolution(A, top=top + 1, window=window)
    dom = A.domain
    dims, mats = _collapse_complex(res, top, side="hom")
    out = []
    for n in range(top + 1):
        out.append(dims[n] - _rank_of(mats.get(n, []), dom) - _rank_of(mats.get(n + 1, []), dom))
    return out


def ghms_cohomology_dims(A, top, window=4):
    """``[dim HH^0, ..., dim HH^top]`` via ``Hom_{A^e}(P_., A)``.

    ``delta^n : C^n -> C^{n+1}`` is the transpose-side of the same terms with the ``a.w.b``
    order, so ``dim HH^n = dim C^n - rank delta^n - rank delta^{n-1}``.
    """
    res = GHMSResolution(A, top=top + 1, window=window)
    dom = A.domain
    dims, mats = _collapse_complex(res, top, side="coh")
    # delta^n : C^n -> C^{n+1} is carried by d_{n+1}; its matrix from _collapse_complex is
    # stored under key n+1 with rows C^n and cols C^{n+1} -- transpose the roles.
    def r(n):
        M = mats.get(n + 1, [])
        return _rank_of(M, dom)
    out = []
    for n in range(top + 1):
        out.append(dims[n] - r(n) - (r(n - 1) if n else 0))
    return out
