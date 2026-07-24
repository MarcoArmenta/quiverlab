"""Frobenius / Nakayama / symmetry over ANY exact Domain (Plan 19).

DECISION (exact, any field): a basic split algebra — path-type basis — is
self-injective, equivalently Frobenius, iff every soc(e_v A) is simple and
v |-> vertex(soc(e_v A)) is a permutation (Nakayama; Skowronski–Yamagata,
Frobenius Algebras I). Both directions are Domain linear algebra:
soc(e_v A) = {x in e_v A : x r = 0}, and sufficiency is the dimension count
P_v -> I(soc P_v) with sum dim P_v = dim A = sum dim I(S_v).

FORM: the socle-dual covector is tried first, then deterministic fallbacks;
nondegeneracy is VERIFIED (Gram rank = dim A), so a returned form is never
wrong. nu = G^{-1} G^T (columns = images); the defining identity
lambda(ab) = lambda(b nu(a)) and multiplicativity are asserted by the test
battery.

SYMMETRY: A symmetric iff Frobenius and nu is INNER. (a) nontrivial
Nakayama vertex permutation -> False (inner automorphisms fix
primitive-idempotent classes); (b) else search span(U),
U = {u : nu(a) u = u a}, for an invertible element on a grid of sample
values per coordinate: det of the left-multiplication is a polynomial of
degree <= dim A on U, so (Schwartz–Zippel) vanishing on a grid with
> dim A distinct values per variable forces vanishing on U — the sweep is
CONCLUSIVE. If the Domain cannot supply enough distinct samples, or the
sweep exceeds `budget` combinations, raise loudly (never a silent wrong
answer)."""
import itertools

from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import nullspace, rank, solve
from quiverlab.invariants.pathbasis import path_type_basis


def _right_socles(A, idem, rad, src):
    """{v: basis of soc(e_v A) as full-coordinate vectors}."""
    dom = A.domain
    out = {}
    for v in idem:
        mine = [v] + [i for i in rad if src[i] == v]
        rows = []
        for b in rad:
            for t in range(A.dim):
                rows.append([A.T[i][b][t] for i in mine])
        if rows:
            ker = nullspace(rows, dom)
        else:                                   # no radical at all: soc = e_v A
            ker = [[dom.one() if k == pos else dom.zero()
                    for k in range(len(mine))] for pos in range(len(mine))]
        socs = []
        for vec in ker:
            full = [dom.zero()] * A.dim
            for pos, i in enumerate(mine):
                full[i] = vec[pos]
            socs.append(full)
        out[v] = socs
    return out


def nakayama_data(A):
    """(ok, perm, gens): ok = the socle criterion holds; perm maps each
    vertex-idempotent index v to the vertex of soc(e_v A); gens[v] = the
    1-dim socle generator (full coordinates)."""
    dom = A.domain
    idem, rad, src, tgt = path_type_basis(A, "the Frobenius test")
    socs = _right_socles(A, idem, rad, src)
    perm, gens = {}, {}
    for v in idem:
        if len(socs[v]) != 1:
            return False, None, None
        s = socs[v][0]
        nz = [i for i, c in enumerate(s) if not dom.is_zero(c)]
        targets = {tgt[i] if i in tgt else i for i in nz}
        if len(targets) != 1:
            raise QuiverlabError(
                "1-dimensional socle generator is not corner-homogeneous — "
                "impossible for a path-type basis (bug)")
        perm[v] = targets.pop()
        gens[v] = s
    if len(set(perm.values())) != len(idem):
        return False, None, None
    return True, perm, gens


def is_frobenius_generic(A):
    """Exact, conclusive, any Domain (socle criterion on a path-type basis)."""
    ok, _perm, _gens = nakayama_data(A)
    return ok


def _lam_of(lam, vec, dom):
    acc = dom.zero()
    for t, c in enumerate(vec):
        if not dom.is_zero(c):
            acc = dom.add(acc, dom.mul(lam[t], c))
    return acc


def frobenius_form_generic(A):
    """(lam, G): a VERIFIED-nondegenerate Frobenius covector (lam(x) =
    sum_t lam[t] x_t) and its Gram matrix G[i][j] = lam(f_i f_j).
    Raises QuiverlabError if A is not Frobenius, and loudly if no candidate
    passes the exact rank check (never returns an unverified form)."""
    dom = A.domain
    m = A.dim
    ok, _perm, gens = nakayama_data(A)
    if not ok:
        raise QuiverlabError(
            "algebra is not Frobenius (socle criterion), so it has no "
            "Frobenius form / Nakayama automorphism")

    def gram(lam):
        return [[_lam_of(lam, A.T[i][j], dom) for j in range(m)] for i in range(m)]

    socdual = [dom.zero()] * m
    for s in gens.values():
        for t, c in enumerate(s):
            socdual[t] = dom.add(socdual[t], c)
    cands = [socdual]
    cands.extend([[dom.one() if t == c else dom.zero() for t in range(m)]
                  for c in range(m)])
    cands.append([dom.one()] * m)
    for lam in cands:
        G = gram(lam)
        if rank(G, dom) == m:
            return lam, G
    raise QuiverlabError(
        "algebra IS Frobenius (socle criterion) but no deterministic "
        "candidate covector was nondegenerate — please report this algebra",
        hint="the socle-dual functional is expected to work for every "
             "path-type basis")


def nakayama_automorphism_generic(A):
    """nu as a matrix of Domain elements (columns = images): N = G^{-1} G^T,
    i.e. the solution of G N = G^T, from the defining identity
    lam(ab) = lam(b nu(a))."""
    lam, G = frobenius_form_generic(A)
    dom = A.domain
    m = A.dim
    N = [[dom.zero()] * m for _ in range(m)]
    for j in range(m):
        col = solve(G, [G[j][i] for i in range(m)], dom)   # (G^T) column j
        for i in range(m):
            N[i][j] = col[i]
    return N


def _sample_values(dom, need):
    """Up to `need` pairwise-distinct Domain elements, deterministically:
    integer coercions first (0, 1, 2, ...). Returns fewer than `need` only
    when the Domain's prime subfield is exhausted (small finite fields)."""
    vals = []
    k = 0
    while len(vals) < need and k < 4 * need + 8:
        cand = dom.coerce(k)
        if all(not dom.eq(cand, v) for v in vals):
            vals.append(cand)
        k += 1
    return vals


def is_symmetric_generic(A, budget=4096):
    """A symmetric iff Frobenius and nu is inner — exact; loud when the
    sample sweep cannot be made conclusive within `budget`."""
    ok, perm, _gens = nakayama_data(A)
    if not ok:
        return False
    if any(v != w for v, w in perm.items()):
        return False                    # inner autos fix vertex classes
    dom = A.domain
    m = A.dim
    N = nakayama_automorphism_generic(A)
    # U = {u : nu(f_a) u = u f_a for all a}
    rows = []
    for a in range(m):
        nua = [N[i][a] for i in range(m)]
        for t in range(m):
            row = []
            for k in range(m):
                left = dom.zero()       # coeff of u_k in (nu(f_a) u)_t
                for i in range(m):
                    if not dom.is_zero(nua[i]):
                        left = dom.add(left, dom.mul(nua[i], A.T[i][k][t]))
                row.append(dom.sub(left, A.T[k][a][t]))
            rows.append(row)
    U = nullspace(rows, dom)
    if not U:
        return False

    def invertible(u):
        L = [[dom.zero()] * m for _ in range(m)]   # left multiplication by u
        for i, c in enumerate(u):
            if dom.is_zero(c):
                continue
            for j in range(m):
                for t, x in enumerate(A.T[i][j]):
                    if not dom.is_zero(x):
                        L[t][j] = dom.add(L[t][j], dom.mul(c, x))
        return rank(L, dom) == m

    from quiverlab.fields.primefield import PrimeField
    samples = _sample_values(dom, m + 1)
    # Conclusive iff > dim A distinct samples (Schwartz–Zippel), OR the sweep
    # enumerates the WHOLE coefficient space — which integer coercions give
    # only over the prime field itself. Over GF(p^n), n > 1, with p <= dim A,
    # integer samples span just the prime subfield of each coefficient: that
    # is NOT the whole space, so it must refuse loudly, never guess.
    whole_field = isinstance(dom, PrimeField) and len(samples) == dom.characteristic
    if len(samples) < m + 1 and not whole_field:
        raise QuiverlabError(
            "symmetry test: the Domain supplies too few distinct sample "
            "values for a conclusive Schwartz–Zippel sweep",
            hint="decide symmetry over GF(p) via the engine, or over a "
                 "larger field")
    if len(samples) ** len(U) > budget:
        raise QuiverlabError(
            f"symmetry test: the conclusive sweep needs "
            f"{len(samples)}^{len(U)} > budget={budget} combinations",
            hint="raise the budget")
    for coeffs in itertools.product(samples, repeat=len(U)):
        u = [dom.zero()] * m
        for cf, base in zip(coeffs, U):
            if dom.is_zero(cf):
                continue
            for t, x in enumerate(base):
                u[t] = dom.add(u[t], dom.mul(cf, x))
        if invertible(u):
            return True
    return False
