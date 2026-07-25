"""Hom and End spaces of right A-modules over any Domain (spec §3.6, §5 c.7).

phi: M -> N is a right-module map iff N.action[b] @ phi = phi @ M.action[b] for every
generator b (arrows and idempotents). Column-stacking vec: the constraint per b is
(I_{dimM} (x) N.action[b] - M.action[b]^T (x) I_{dimN}) vec(phi) = 0. dim Hom = dim ker."""
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm


def _generators(M):
    labels = [f"e_{v}" for v in M.algebra.quiver.vertices]
    labels += list(M.algebra.quiver.arrows)
    return labels


def hom_space(M, N):
    dom = M.domain
    dm, dn = M.dim, N.dim
    Im = lm.identity(dm, dom)
    In = lm.identity(dn, dom)
    blocks = []
    for b in _generators(M):
        Nb, Mb = N.action[b], M.action[b]
        left = lm.kron(Im, Nb, dom)                     # I_dm (x) N.action[b]
        right = lm.kron(lm.transpose(Mb), In, dom)      # M.action[b]^T (x) I_dn
        blocks.append([[dom.sub(left[i][j], right[i][j]) for j in range(dm * dn)]
                       for i in range(dm * dn)])
    if not blocks:
        # _generators(M) always yields the vertex idempotents (every quiver has at
        # least one vertex) plus the arrows, so `blocks` is never empty. If it is,
        # the constraint matrix is missing entirely and the "no constraints" answer
        # would be the FULL space (dim dm*dn), not dim 0 -- returning zeros(0, ...)
        # here would silently understate Hom. Raise loudly instead of guessing.
        raise QuiverlabError(
            "hom_space: generator list unexpectedly empty (no idempotents or arrows)",
            hint="every quiver has at least one vertex, so this branch is unreachable; "
                 "if it fires, the module's algebra/quiver is malformed")
    stacked = lm.vstack(blocks)
    ker = lm.kernel_columns(stacked, dom)
    # reshape each kernel vector (length dm*dn, column-stacked) into a dn x dm matrix
    homs = []
    for z in ker:
        phi = [[z[j * dn + i] for j in range(dm)] for i in range(dn)]
        homs.append(phi)
    return homs


def hom_dim(M, N):
    return len(hom_space(M, N))


def end_dim(M):
    return hom_dim(M, M)


_ISO_SEARCH_CAP = 50000


def _prime_field_size(dom):
    """p if dom is a prime field GF(p) whose elements are ints 0..p-1; else None
    (char-0 field or GF(p^n), where a plain int ladder does not enumerate)."""
    from quiverlab.fields.primefield import PrimeField
    if isinstance(dom, PrimeField):
        return dom.p
    return None


def _combine(H, coeffs, dom):
    r, n = len(H), len(H[0])
    m = len(H[0][0])
    out = [[dom.zero() for _ in range(m)] for _ in range(n)]
    for c, mat in zip(coeffs, H):
        if dom.is_zero(c):
            continue
        for i in range(n):
            oi, mi = out[i], mat[i]
            for j in range(m):
                oi[j] = dom.add(oi[j], dom.mul(c, mi[j]))
    return out


def _search_invertible(H, dom, n, vals):
    """Try every nonzero coefficient vector over `vals`; return True on the first
    invertible combination. Caller decides whether `vals` is exhaustive."""
    import itertools
    r = len(H)
    for idx in itertools.product(range(len(vals)), repeat=r):
        if all(t == 0 for t in idx):
            continue
        C = _combine(H, [vals[t] for t in idx], dom)
        if lm.mat_rank(C, dom) == n:
            return True
    return False


def _ladder_size(r):
    k = 2
    while (k + 1) ** r <= _ISO_SEARCH_CAP:
        k += 1
    return k


def _generic_full_rank_char0(H, dom, n):
    """Exact generic rank of span{H_k} over a char-0 field: rank of sum_k t_k H_k
    over the rational function field (sympy, exact). By Noether-Deuring this equals
    the base-field iso decision (iso over the closure <=> iso over the base field
    for f.d. modules). Returns True iff the generic rank is n (invertible exists)."""
    import sympy
    r = len(H)
    ts = sympy.symbols(f"t0:{r}")
    S = sympy.zeros(n, n)
    for k in range(r):
        for i in range(n):
            row = H[k][i]
            for j in range(n):
                if not dom.is_zero(row[j]):
                    S[i, j] += ts[k] * sympy.sympify(row[j])
    return S.rank() == n


def is_isomorphic(M, N):
    """True iff M and N are isomorphic right A-modules, certified by an EXACT
    invertible module map (positive certificate -- never a false positive). Decides
    over the BASE field (not the algebraic closure).

    A negative answer is returned only when it is certain: distinct dimension /
    dimension vectors, or an exhaustive search over a small prime field found no
    isomorphism. When neither a witness nor an exhaustive refutation is reachable
    within the search budget (large GF(p^n), or an unlucky char-0 search), it raises
    loudly rather than guess."""
    if M.dim != N.dim:
        return False
    if M.dimension_vector() != N.dimension_vector():
        return False
    n = M.dim
    if n == 0:
        return True
    H = hom_space(M, N)                      # n x n matrices (dim M = dim N)
    if not H:
        return False
    dom = M.domain
    for phi in H:                            # a basis element may already be an iso
        if lm.mat_rank(phi, dom) == n:
            return True
    r = len(H)
    p = _prime_field_size(dom)
    if p is not None and p ** r <= _ISO_SEARCH_CAP:   # exhaustive over GF(p): decisive
        vals = [dom.coerce(i) for i in range(p)]
        return _search_invertible(H, dom, n, vals)
    if dom.characteristic == 0:              # char 0: generic rank decides (Noether-Deuring)
        return _generic_full_rank_char0(H, dom, n)
    # large GF(p) or GF(p^n): bounded, positive-only search over an int ladder
    hi = _ladder_size(r) + 1
    vals = [dom.coerce(i) for i in range(min(p, hi) if p is not None else hi)]
    if _search_invertible(H, dom, n, vals):  # positive-only decisive
        return True
    raise QuiverlabError(
        "is_isomorphic: could not certify over this field within the search budget",
        hint="the modules share dim and dimension vector but no base-field "
             "isomorphism was found; crosscheck via QPA or enlarge the search")
