"""Krull--Schmidt decomposition of finite-dimensional A-modules (Plan 30 Part A).

``M = (+) M_i^{m_i}`` into indecomposables, unique up to isomorphism and permutation
(Krull--Schmidt). Exact, Domain-generic, and HOUSE-HONEST: every decomposition we
return is CERTIFIED -- each summand is provably indecomposable and the summands
reassemble ``M`` (the splits are exact kernel direct sums) -- and when neither a split
nor an indecomposability certificate is reachable within the search budget we RAISE
loudly (mirroring ``modules.hom.is_isomorphic``'s refusal style) rather than emit an
unproven verdict. There is never a silent wrong "indecomposable".

ALGORITHM (Fitting splitting via minimal polynomials + a local-endomorphism certificate).
``End_A(M) = hom_space(M, M)`` is a finite-dimensional k-algebra of ``n x n`` matrices
(``n = dim M``) acting faithfully on ``M``. For an endomorphism ``phi`` with minimal
polynomial ``m(x)`` over the ground field k:

  * if ``m = f * g`` with ``f, g`` COPRIME in ``k[x]`` then, since ``m(phi) = 0``,
    the Chinese Remainder / Fitting decomposition gives
    ``M = ker f(phi) (+) ker g(phi)`` and BOTH kernels are A-SUBMODULES because
    ``f(phi)`` and ``g(phi)`` commute with the A-action (``phi`` is an A-map, so any
    polynomial in ``phi`` is too, hence its kernel is A-stable). Equivalently, with
    Bezout ``u*f + v*g = 1`` the element ``e = (v*g)(phi)`` is an idempotent in
    ``End_A(M)`` with ``im(e) = ker f(phi)``; we split via the two kernels directly, so
    no explicit idempotent matrix is needed.
  * ``m`` factors as ``prod p_i^{e_i}``; a nontrivial coprime split exists iff ``m`` is
    NOT a prime power (>= 2 distinct irreducible factors). We scan a bounded family of
    ``phi`` (the ``End`` basis and small combinations); every split strictly lowers the
    dimension of each part, so the recursion terminates.

No split found within budget ==> CERTIFY that ``End_A(M)`` is local (Fitting's lemma:
``M`` indecomposable <=> ``End_A(M)`` local <=> ``End_A(M)/rad`` a division algebra):

  * ``dim End_A(M) = 1`` ==> ``End = k * id``, a field ==> local. Rigorous in EVERY
    characteristic; this is the certificate every simple module (and most indecomposable
    projectives/injectives here) lands on.
  * else ``dim(End/rad)`` is read off the trace form ``T[i][j] = tr(H_i H_j)`` on the
    natural module ``M`` (``rad(End) = End^{perp}`` w.r.t. this form when ``char k = 0``
    or ``char k > n`` -- the classical Dickson/Cohen--Ivanyos--Wales bound; the strict
    ``char > n`` is necessary, e.g. ``k[x]/(x^p)`` in ``M_p`` has a degenerate trace
    form). When the bound holds and ``rank T = 1`` then ``End/rad = k`` (a field) ==>
    local ==> indecomposable.
  * otherwise (``rank T > 1``, or ``char p <= n`` so the trace form is unreliable) we
    RAISE the loud undecided error naming the situation -- crosscheck via QPA or work
    over a larger characteristic. NEVER a silent verdict.

tau-ADDITIVITY (docstring note carried onto ``Module.decompose``): the Auslander--Reiten
translate is additive, ``tau((+) M_i) = (+) tau M_i`` (and likewise ``tau^-``); the
Plan 30 reporting layer uses this to certify each translate summand-wise -- decompose
first, translate the indecomposable summands.
"""
from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import hom_space, is_isomorphic
from quiverlab.modules.radtopsoc import submodule

# Default cap on the number of candidate endomorphisms the split search evaluates
# before it gives up and falls through to the locality certificate. Small modules
# (the whole tested zoo) split -- or certify -- far inside this; the cap only bounds
# the pathological large-End case, where we refuse loudly rather than churn.
_DEFAULT_BUDGET = 512


# ---------------------------------------------------------------------------
# Exact polynomials over the Domain (ascending coefficient lists a0 + a1 x + ...)
# ---------------------------------------------------------------------------
def _min_poly_coeffs(phi, dom):
    """The exact minimal polynomial of the ``n x n`` matrix ``phi`` over ``dom``, as a
    monic ascending coefficient list ``[a0, ..., a_{d-1}, 1]`` (degree ``d = 1..n``).

    Krylov / first-dependence: find the least ``d`` for which ``phi^d`` lies in the span
    of ``I, phi, ..., phi^{d-1}`` (as vectors in ``k^{n^2}``); the dependence coefficients
    give ``phi^d = sum a_i phi^i``, i.e. ``m(x) = x^d - sum a_i x^i``. Cayley--Hamilton
    guarantees termination by ``d = n``. Pure Domain linear algebra -- no sympy, exact
    over every field."""
    n = len(phi)
    def _vec(mat):
        return [mat[i][j] for i in range(n) for j in range(n)]
    cur = lm.identity(n, dom)               # phi^0 = I
    basis_vecs = [_vec(cur)]                 # vec(phi^0), vec(phi^1), ...
    for _ in range(1, n + 1):
        cur = lm.matmul(phi, cur, dom)       # next power phi^k
        target = _vec(cur)
        # solve (columns = vec(phi^0..phi^{k-1})) * x = vec(phi^k)
        B = lm.cols_to_matrix(basis_vecs)
        x = linalg.solve(B, target, dom)
        if x is not None:                    # phi^k = sum x_i phi^i  -> minimal poly found
            return [dom.neg(c) for c in x] + [dom.one()]
        basis_vecs.append(target)
    # Unreachable: phi^n is dependent on lower powers by Cayley--Hamilton. Guard loudly
    # so a linear-algebra regression can never silently corrupt a decomposition.
    raise QuiverlabError(
        "minimal polynomial exceeded degree dim M (Cayley--Hamilton violated)",
        hint="this indicates a bug in the exact linear-algebra kernel, not user input")


def _poly_mul(a, b, dom):
    """Product of two ascending coefficient lists over ``dom``."""
    if not a or not b:
        return []
    out = [dom.zero() for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        if dom.is_zero(ai):
            continue
        for j, bj in enumerate(b):
            out[i + j] = dom.add(out[i + j], dom.mul(ai, bj))
    return out


def _poly_pow(a, e, dom):
    """``a`` raised to the non-negative integer power ``e`` (binary exponentiation)."""
    result = [dom.one()]
    base = list(a)
    while e:
        if e & 1:
            result = _poly_mul(result, base, dom)
        e >>= 1
        if e:
            base = _poly_mul(base, base, dom)
    return result


def _poly_eval_matrix(coeffs, phi, dom):
    """``m(phi)`` for the ascending polynomial ``coeffs``, as an ``n x n`` matrix
    (Horner from the top coefficient down)."""
    n = len(phi)
    acc = lm.zeros(n, n, dom)
    for c in reversed(coeffs):               # Horner: acc <- phi*acc + c*I
        acc = lm.matmul(phi, acc, dom)
        for i in range(n):
            acc[i][i] = dom.add(acc[i][i], c)
    return acc


# ---------------------------------------------------------------------------
# Factoring the minimal polynomial over the ground field (sympy bridge)
# ---------------------------------------------------------------------------
def _factoring_supported(dom):
    """True iff we can factor a univariate polynomial over ``dom``'s ground field
    exactly (prime field GF(p), the rationals QQ, or a CC number field QQ<alpha>).
    GF(p^n) and other domains are out of the certified factoring budget -- the split
    search is skipped there and the answer rests on the trace-form / dim-End certificate
    or a loud refusal."""
    from quiverlab.fields.complexfield import SympyExactDomain
    from quiverlab.fields.primefield import PrimeField
    from quiverlab.fields.rationals import RationalField
    return isinstance(dom, (PrimeField, RationalField, SympyExactDomain))


def _factor_min_poly(coeffs, dom):
    """Factor the monic minimal polynomial ``coeffs`` (ascending) over ``dom``'s ground
    field, returning ``[(factor_coeffs, multiplicity), ...]`` with each factor a monic
    ascending Domain-coefficient list. Requires ``_factoring_supported(dom)``.

    The three exact ground fields map to the matching sympy factoring domain:
    ``GF(p)`` via ``modulus=p``, ``QQ`` via ``domain='QQ'``, and a CC number field via
    the algebra's own sympy algebraic domain ``dom.sdom`` (so ``QQ<alpha>`` factors
    over ``QQ<alpha>``). Coefficients round-trip through the Domain, so the returned
    factors are exact Domain data."""
    import sympy
    from quiverlab.fields.complexfield import SympyExactDomain
    from quiverlab.fields.primefield import PrimeField
    from quiverlab.fields.rationals import RationalField
    x = sympy.Symbol("x")

    if isinstance(dom, PrimeField):
        expr = sum(int(c) * x ** i for i, c in enumerate(coeffs))
        poly = sympy.Poly(expr, x, modulus=dom.p)
        _, facs = poly.factor_list()
        return [([dom.coerce(int(a)) for a in reversed(f.all_coeffs())], mult)
                for f, mult in facs]

    if isinstance(dom, RationalField):
        from fractions import Fraction
        expr = sum(sympy.Rational(c.numerator, c.denominator) * x ** i
                   for i, c in enumerate(coeffs))
        poly = sympy.Poly(expr, x, domain="QQ")
        _, facs = poly.factor_list()
        out = []
        for f, mult in facs:
            fc = [Fraction(int(a.p), int(a.q)) if hasattr(a, "q") else Fraction(int(a))
                  for a in reversed(f.all_coeffs())]
            out.append(([dom.coerce(a) for a in fc], mult))
        return out

    if isinstance(dom, SympyExactDomain):
        sdom = dom.sdom
        expr = sum(sdom.to_sympy(c) * x ** i for i, c in enumerate(coeffs))
        poly = sympy.Poly(expr, x, domain=sdom)
        _, facs = poly.factor_list()
        # Poly.all_coeffs() yields sympy EXPRESSIONS (1, -sqrt(2), -1, ...); dom.coerce
        # re-embeds each in the number field (never sdom.to_sympy, which expects a domain
        # element and trips on a plain sympy Integer such as NegativeOne).
        return [([dom.coerce(a) for a in reversed(f.all_coeffs())], mult)
                for f, mult in facs]

    # Guarded by _factoring_supported at the call site; reaching here is a bug.
    raise QuiverlabError(
        f"cannot factor over {dom} within the certified budget",
        hint="decomposition factoring supports GF(p), QQ, and CC number fields")


# ---------------------------------------------------------------------------
# The split search
# ---------------------------------------------------------------------------
def _lin_combo(H, coeffs, dom):
    """The endomorphism ``sum coeffs[i] * H[i]`` as a matrix."""
    n = len(H[0])
    out = lm.zeros(n, n, dom)
    for c, mat in zip(coeffs, H):
        if dom.is_zero(c):
            continue
        for i in range(n):
            oi, mi = out[i], mat[i]
            for j in range(n):
                oi[j] = dom.add(oi[j], dom.mul(c, mi[j]))
    return out


def _candidate_endomorphisms(H, dom, budget):
    """Yield up to ``budget`` deterministic candidate endomorphisms from the ``End``
    basis ``H``: every basis element first (a single basis element already splits every
    module in the tested zoo), then small pairwise combinations ``H_i +/- H_j`` and a
    coefficient ladder ``sum (i+1) H_i`` to reach idempotents that are not themselves
    basis elements. Deterministic order == reproducible decompositions."""
    r = len(H)
    count = 0
    for i in range(r):                       # 1) singletons
        if count >= budget:
            return
        yield H[i]
        count += 1
    one = dom.one()
    for i in range(r):                       # 2) pairwise sums / differences
        for j in range(i + 1, r):
            if count >= budget:
                return
            yield _lin_combo([H[i], H[j]], [one, one], dom)
            count += 1
            if count >= budget:
                return
            yield _lin_combo([H[i], H[j]], [one, dom.neg(one)], dom)
            count += 1
    if count < budget and r > 1:             # 3) a coefficient ladder over all of H
        yield _lin_combo(H, [dom.coerce(i + 1) for i in range(r)], dom)


def _split(M, phi, factors):
    """Split ``M = ker f(phi) (+) ker g(phi)`` for a nontrivial coprime factorization,
    with ``f`` = the first irreducible factor to its full multiplicity and ``g`` = the
    product of the rest. Returns two A-submodules. The dimension identity
    ``dim X + dim Y == dim M`` is the exact Fitting/CRT certificate (asserted)."""
    dom = M.domain
    f = _poly_pow(factors[0][0], factors[0][1], dom)
    g = [dom.one()]
    for fc, mult in factors[1:]:
        g = _poly_mul(g, _poly_pow(fc, mult, dom), dom)
    Xcols = lm.kernel_columns(_poly_eval_matrix(f, phi, dom), dom)   # f(phi)-primary part
    Ycols = lm.kernel_columns(_poly_eval_matrix(g, phi, dom), dom)   # complementary part
    assert Xcols and Ycols and len(Xcols) + len(Ycols) == M.dim, (
        "Fitting split failed to give a direct-sum decomposition "
        f"(dim ker f={len(Xcols)}, dim ker g={len(Ycols)}, dim M={M.dim}); "
        "the coprime factorization or the minimal polynomial is inconsistent")
    X = submodule(M, Xcols, name=f"{M.name}(1)")
    Y = submodule(M, Ycols, name=f"{M.name}(2)")
    return X, Y


def _try_split(M, H, budget):
    """Search for a Fitting split of ``M`` within budget. Returns ``(X, Y)`` submodules
    or ``None`` (no split found -- either indecomposable, or beyond the search budget;
    the caller certifies which). Skipped entirely when factoring is unsupported over the
    domain (then only the dim-End/trace certificate or a loud refusal decides)."""
    dom = M.domain
    if not _factoring_supported(dom):
        return None
    for phi in _candidate_endomorphisms(H, dom, budget):
        m = _min_poly_coeffs(phi, dom)
        if len(m) <= 2:                      # degree <= 1: scalar, no coprime split
            continue
        factors = _factor_min_poly(m, dom)
        if len(factors) >= 2:                # >= 2 distinct irreducibles => coprime split
            return _split(M, phi, factors)
    return None


# ---------------------------------------------------------------------------
# The locality certificate (M indecomposable <=> End_A(M) local)
# ---------------------------------------------------------------------------
def _trace_form_rank(H, dom):
    """Rank of the natural trace form ``T[i][j] = tr(H_i H_j)`` on ``M`` (an ``r x r``
    Gram matrix, ``r = dim End``). Equals ``dim(End/rad End)`` when ``char = 0`` or
    ``char > dim M`` (Dickson / Cohen--Ivanyos--Wales)."""
    r = len(H)
    n = len(H[0]) if H else 0
    T = lm.zeros(r, r, dom)
    for i in range(r):
        for j in range(r):
            prod = lm.matmul(H[i], H[j], dom)
            s = dom.zero()
            for d in range(n):
                s = dom.add(s, prod[d][d])
            T[i][j] = s
    return lm.mat_rank(T, dom)


def _certify_local(M, H):
    """Certify that ``End_A(M)`` is local (so ``M`` is indecomposable), or RAISE loudly
    when no certificate is reachable. Called only after the split search found nothing.

    Certificates, in order:
      * ``dim End = 1`` ==> ``End = k * id`` is a field ==> local (any characteristic).
      * trace-form rank ``= dim(End/rad) = 1`` ==> ``End/rad = k`` ==> local, valid when
        the trace form computes the radical: ``char 0`` or ``char > dim M``.
    Anything else is undecided within budget and refuses loudly."""
    r = len(H)
    if r <= 1:                               # End = k*id (r == 1); r == 0 only for dim 0
        return True
    dom = M.domain
    char = dom.characteristic
    n = M.dim
    trace_rigorous = (char == 0) or (char > n)
    if not trace_rigorous:
        raise QuiverlabError(
            f"decompose: cannot certify indecomposability of {M.name} in characteristic "
            f"{char} <= dim M = {n} within budget (the trace-form radical is unreliable "
            "when char <= dim, e.g. it over-counts on k[x]/(x^p)), and no Fitting split "
            "was found",
            hint="crosscheck via QPA (DecomposeModuleWithMultiplicities / "
                 "IsIndecomposableModule) or recompute over a characteristic > dim M")
    rank_T = _trace_form_rank(H, dom)
    if rank_T == 1:                          # End/rad = k, a division algebra => local
        return True
    raise QuiverlabError(
        f"decompose: could not certify {M.name} within budget -- End/rad has dimension "
        f"{rank_T} > 1 (End is not obviously local) yet the bounded Fitting search found "
        "no split; the module may be decomposable with a splitting element outside the "
        "searched combinations, or genuinely indecomposable over a division algebra",
        hint="enlarge budget, or crosscheck via QPA "
             "(DecomposeModuleWithMultiplicities / IsIndecomposableModule)")


# ---------------------------------------------------------------------------
# Public surface
# ---------------------------------------------------------------------------
def _indecomposable_summands(M, budget):
    """Flat list of certified-indecomposable submodules whose direct sum is ``M``.
    Recurses on Fitting splits; certifies every leaf (raising loudly if a leaf can be
    neither split nor certified). Dimension strictly drops on each split, so this
    terminates."""
    if M.dim == 0:
        return []
    H = hom_space(M, M)
    split = _try_split(M, H, budget)
    if split is not None:
        X, Y = split
        return _indecomposable_summands(X, budget) + _indecomposable_summands(Y, budget)
    _certify_local(M, H)                     # raises unless a locality certificate holds
    return [M]


def decompose(M, budget=_DEFAULT_BUDGET):
    """Krull--Schmidt decomposition of ``M`` into indecomposables with multiplicities:
    a list ``[(M_i, m_i), ...]`` with each ``M_i`` a certified-indecomposable summand,
    ``m_i`` its multiplicity, and ``(+) M_i^{m_i} ~ M``. Summands are grouped up to
    isomorphism by the exact ``is_isomorphic`` certificate; distinct entries are
    pairwise non-isomorphic. Order is deterministic (first appearance in the recursion).

    Raises loudly (never a silent wrong answer) when a summand can be neither split nor
    certified indecomposable within budget -- see the module docstring for the honest
    scope (``char <= dim`` over GF(p), or ``End/rad`` of dimension ``> 1``).

    tau-additivity: ``tau(M) = (+) tau(M_i)^{m_i}`` (the AR translate is additive), so a
    translate of a decomposable module is certified summand-wise."""
    summands = _indecomposable_summands(M, budget)
    groups = []                              # list of [representative, multiplicity]
    for s in summands:
        for g in groups:
            if is_isomorphic(g[0], s):
                g[1] += 1
                break
        else:
            groups.append([s, 1])
    return [(g[0], g[1]) for g in groups]


def is_indecomposable(M, budget=_DEFAULT_BUDGET):
    """True iff ``M`` is indecomposable, certified (End local); False iff a Fitting split
    exists. The zero module is NOT indecomposable (``False``). Raises loudly when neither
    a split nor a locality certificate is reachable within budget -- same honest scope as
    :func:`decompose`."""
    if M.dim == 0:
        return False                         # the zero module is the empty direct sum
    H = hom_space(M, M)
    if _try_split(M, H, budget) is not None:
        return False
    return _certify_local(M, H)
