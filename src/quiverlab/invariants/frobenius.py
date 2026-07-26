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

SYMMETRY (Plan 29): A is symmetric iff it admits a NONDEGENERATE symmetric
(trace) form lambda, lambda(ab) = lambda(ba) — equivalently A ~= DA as
bimodules (Skowronski–Yamagata, Frobenius Algebras I). Decision: (a) not
Frobenius or a nontrivial Nakayama vertex permutation -> False (weakly-
symmetric fails); (b) else CERTIFY symmetry by an EXACT witness — some
lambda in the trace-form space whose Gram ``G[i][j] = lambda(f_i f_j)`` has full
rank (no false positives, no sweep); (c) if no cheap witness is found decide
the negative by a CONCLUSIVE Schwartz–Zippel sweep of the trace-form space
(det Gram has degree <= dim A, so > dim A distinct samples per coordinate is
conclusive), raising loudly when the Domain is too small or the sweep exceeds
`budget` (never a silent wrong answer). This replaced a nu-is-INNER search
whose exhaustive grid raised on dim >= 4 and whose GF(p) delegator wrongly
required the engine's Nakayama automorphism to be the IDENTITY MATRIX (a
sufficient-not-necessary test that returned a SILENT WRONG False on
multi-vertex Brauer stars). WEAK SYMMETRY: Frobenius with the identity
Nakayama permutation (soc P_v = top P_v for every projective)."""
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
    sum_t lam_t x_t) and its Gram matrix G_{ij} = lam(f_i f_j).
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


def _trace_form_space(A):
    """Basis of the SYMMETRIC linear forms lambda with lambda(ab) = lambda(ba)
    for all a, b — equivalently Hom_{A-A}(A, DA) = (A/[A, A])^*. Presentation-
    free (structure constants only): each constraint is ``lambda . (T[i][j] -
    T[j][i]) = 0``, and i < j suffices (the (j, i) pair is its negation)."""
    dom = A.domain
    m = A.dim
    rows = []
    for i in range(m):
        for j in range(i + 1, m):
            row = [dom.sub(A.T[i][j][t], A.T[j][i][t]) for t in range(m)]
            if any(not dom.is_zero(c) for c in row):
                rows.append(row)
    if not rows:                                  # commutative: every form is a trace form
        return [[dom.one() if k == t else dom.zero() for k in range(m)]
                for t in range(m)]
    return nullspace(rows, dom)


def _gram(A, lam):
    """Gram matrix ``G[i][j] = lambda(f_i f_j)`` of a linear form lam (full coords)."""
    dom = A.domain
    m = A.dim
    return [[_lam_of(lam, A.T[i][j], dom) for j in range(m)] for i in range(m)]


def _combine(S, coeffs, dom, m):
    """The form sum_k coeffs[k] * S[k] in full coordinates."""
    out = [dom.zero()] * m
    for cf, base in zip(coeffs, S):
        if dom.is_zero(cf):
            continue
        for t, x in enumerate(base):
            if not dom.is_zero(x):
                out[t] = dom.add(out[t], dom.mul(cf, x))
    return out


def _symmetric_witness_candidates(S, dom):
    """Deterministic, cheap coefficient vectors over the trace-form basis S; if
    ANY yields a full-rank Gram, A is symmetric (an exact certificate — never a
    false positive). Each basis form, their sum, small moment curves, and
    (1, 2, ..., r) suffice for every symmetric algebra met in practice; the
    conclusive sweep is the backstop."""
    r = len(S)
    one, zero = dom.one(), dom.zero()
    for k in range(r):                                    # each basis form alone
        yield [one if t == k else zero for t in range(r)]
    if r > 1:
        yield [one] * r                                   # the sum of all basis forms
        for b in (2, 3, 5):                               # moment curves (1, b, b^2, ...)
            bb = dom.coerce(b)
            cur, vec = one, []
            for _ in range(r):
                vec.append(cur)
                cur = dom.mul(cur, bb)
            yield vec
        yield [dom.coerce(k + 1) for k in range(r)]       # (1, 2, ..., r)


def is_symmetric_generic(A, budget=4096):
    """A is symmetric iff it admits a NONDEGENERATE symmetric (trace) form —
    exact over any Domain; loud when the negative sweep cannot be made
    conclusive within `budget`.

    ROOT CAUSE fixed (Plan 29): the former GF(p) delegator declared A symmetric
    only when the engine's chosen Nakayama automorphism was the IDENTITY MATRIX.
    nu is defined only up to inner automorphism, so that is SUFFICIENT-but-not-
    NECESSARY: multi-vertex symmetric Nakayama kZ_n/J^L with n | (L-1) (Brauer
    stars kZ_2/J^3, kZ_3/J^4, kZ_4/J^5 — QPA IsSymmetricAlgebra = true) have a
    non-identity engine nu yet ARE symmetric, so the old test returned a SILENT
    WRONG False. The old generic route was correct but its exhaustive sweep of
    the nu-twisted centralizer raised on dim >= 4. Both are cured: symmetry is
    now certified by an EXACT nondegenerate-trace-form witness (cheap, no sweep),
    and only the negative goes to the conclusive sweep (or a loud raise). The
    criterion is classical (Skowronski–Yamagata, Frobenius Algebras I):
    symmetric <=> A ~= DA as bimodules <=> some trace form is nondegenerate."""
    dom = A.domain
    m = A.dim
    # nakayama_data preserves the honest refusal on presentation-less algebras
    # (path-basis-needing, consistent with is_frobenius) AND supplies the cheap
    # conclusive negatives — not Frobenius, or a nontrivial Nakayama vertex
    # permutation — that stay correct even over fields too small for the sweep.
    ok, perm, _gens = nakayama_data(A)
    if not ok:
        return False                                 # not Frobenius => not symmetric
    if any(v != w for v, w in perm.items()):
        return False                                 # weakly-symmetric fails => not symmetric
    S = _trace_form_space(A)
    if not S:
        return False
    # EXACT positive certificate: an explicit full-rank Gram of a trace form.
    for coeffs in _symmetric_witness_candidates(S, dom):
        if rank(_gram(A, _combine(S, coeffs, dom, m)), dom) == m:
            return True
    # No cheap witness: decide the negative conclusively (Schwartz–Zippel over
    # the trace-form space) or raise loudly if the Domain is too small.
    return _symmetric_sweep(A, S, budget)


def _symmetric_sweep(A, S, budget):
    """Conclusive Schwartz–Zippel decision: does some lambda in the trace-form
    space S have a full-rank Gram? Returns True/False only when conclusive
    (> dim A distinct samples, or the whole prime field), else raises loudly —
    det Gram has degree <= dim A on S, so a grid with > dim A distinct values per
    coordinate that finds no witness proves none exists (never a silent guess)."""
    from quiverlab.fields.primefield import PrimeField
    dom = A.domain
    m = A.dim
    samples = _sample_values(dom, m + 1)
    # Conclusive iff > dim A distinct samples (Schwartz–Zippel), OR the sweep
    # enumerates the WHOLE prime field (integer coercions exhaust GF(p) exactly).
    # Over GF(p^n), n > 1, with p <= dim A, integer samples span only the prime
    # subfield — NOT the whole space — so it must refuse loudly, never guess.
    whole_field = isinstance(dom, PrimeField) and len(samples) == dom.characteristic
    if len(samples) < m + 1 and not whole_field:
        raise QuiverlabError(
            "symmetry test: the Domain supplies too few distinct sample values "
            "for a conclusive Schwartz–Zippel sweep of the trace-form space",
            hint="decide symmetry over a larger field (e.g. QQ or a bigger GF(p))")
    if len(samples) ** len(S) > budget:
        raise QuiverlabError(
            f"symmetry test: the conclusive sweep needs {len(samples)}^{len(S)} "
            f"> budget={budget} combinations",
            hint="raise the budget (a symmetric algebra is normally caught by the "
                 "nondegenerate-form witness before this sweep)")
    for coeffs in itertools.product(samples, repeat=len(S)):
        if rank(_gram(A, _combine(S, coeffs, dom, m)), dom) == m:
            return True
    return False


def is_weakly_symmetric_generic(A):
    """A is WEAKLY symmetric iff Frobenius (self-injective) with the IDENTITY
    Nakayama permutation — soc(P_v) ~= top(P_v) = S_v for every indecomposable
    projective (the permutation data is in ``nakayama_data``). Exact over any
    Domain (socle criterion on a path-type basis); loud refusal on a
    presentation-less algebra.

    For self-injective Nakayama kZ_n/J^L this is n | (L - 1) (Skowronski–
    Yamagata, Frobenius Algebras I, EMS 2011); for Nakayama algebras weakly
    symmetric <=> symmetric. Every symmetric algebra is weakly symmetric; the
    converse fails (exterior / quantum complete intersections have the identity
    permutation but a non-inner Nakayama automorphism, so are NOT symmetric)."""
    ok, perm, _gens = nakayama_data(A)
    if not ok:
        return False
    return all(v == w for v, w in perm.items())
