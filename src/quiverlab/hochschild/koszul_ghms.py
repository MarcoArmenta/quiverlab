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
from quiverlab.fields.linalg import nullspace, rref

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
