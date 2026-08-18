"""Plan 54 Task A -- the decidable Nakayama-automorphism hypothesis gate.

The BV transport ``(double-dagger)`` is a theorem only when ``nu`` is "nice". The
record (R2) names three computable conditions; Plan 54 settles them exactly and
picks deliberately (spec s2.4):

    Fact (perfect field): ``ord(nu) < inf and char k does not divide ord(nu)``
    <=> ``nu is semisimple of finite order``. So Volkov's condition (1405.5155,
    ``char k does not divide ord nu``) is CONTAINED in LZZ's (1405.5325,
    ``nu semisimple``); the difference {semisimple, infinite order} is char-0 only.

The primary decidable gate is therefore LZZ semisimplicity, decided by ONE exact
test: ``nu`` is semisimple <=> its minimal polynomial is squarefree. Over a
PERFECT field (GF(p), GF(p^n), QQ -- all perfect) every irreducible factor is
separable, so ``squarefree <=> gcd(minpoly, minpoly') is a unit`` is conclusive
and exact (no float). Over GF(p) this coincides with Volkov's ``p does not
divide ord(nu)``; both phrasings are reported in the provenance.

Routing (per instance, in order; ``BVHypothesis.label`` records which fired):

1. ``is_symmetric()``            -> "symmetric" (Tradler AIF 2008; nu inner, no twist).
2. ``is_frobenius()`` + ss nu    -> "semisimple" (Lambre-Zhou-Zimmermann 2016; Volkov).
3. Frobenius, self-inj. Nakayama, NON-semisimple nu (necessarily char | ord nu)
   -> loud refusal: the general Bian-Itagaki-Kou-Lyu-Zhou (arXiv:2603.04834)
   construction is out of v1 scope (the e=1 case k[x]/(x^N) is symmetric, route 1).
4. else                          -> loud refusal (not Frobenius, or non-semisimple nu).

nu is defined only up to inner automorphism; this gate tests the concrete
representative ``Algebra.nakayama_automorphism()`` returns. Semisimplicity of a
WRONG representative would fail the downstream bracket arbiter (a loud refusal),
never return a silent wrong Delta -- the arbiter is the per-instance CORRECTNESS
GATE, but note precisely what it certifies: it pins Delta MODULO cup-derivations
(the data the BV relation ``[a,b] = eps(Delta(aUb) - Delta a U b - (-1)^p a U Delta b)``
constrains), not every last coordinate; a Delta and a cup-derivation-shifted Delta
are indistinguishable to it (spec s2.4 up-to-inner caveat).
"""
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError

_ORDER_CAP = 512


# ---------------------------------------------------------------------------
# exact univariate-polynomial / matrix helpers over a Domain (no float, no sympy)
# ---------------------------------------------------------------------------
def _matmul(dom, X, Y):
    n = len(X)
    k = len(Y)
    w = len(Y[0]) if Y else 0
    out = [[dom.zero()] * w for _ in range(n)]
    for i in range(n):
        Xi = X[i]
        for t in range(k):
            xit = Xi[t]
            if dom.is_zero(xit):
                continue
            Yt = Y[t]
            oi = out[i]
            for j in range(w):
                ytj = Yt[j]
                if not dom.is_zero(ytj):
                    oi[j] = dom.add(oi[j], dom.mul(xit, ytj))
    return out


def _identity(dom, m):
    return [[dom.one() if i == j else dom.zero() for j in range(m)] for i in range(m)]


def _poly_trim(p, dom):
    d = -1
    for i, c in enumerate(p):
        if not dom.is_zero(c):
            d = i
    return list(p[:d + 1])


def _poly_deg(p, dom):
    d = -1
    for i, c in enumerate(p):
        if not dom.is_zero(c):
            d = i
    return d


def _poly_gcd(a, b, dom):
    """Monic gcd of two coefficient lists (low->high) over the field ``dom``."""
    a, b = _poly_trim(a, dom), _poly_trim(b, dom)
    while b:
        db = _poly_deg(b, dom)
        r = list(a)
        while _poly_deg(r, dom) >= db:
            dr = _poly_deg(r, dom)
            lead = dom.mul(r[dr], dom.inv(b[db]))
            for i in range(db + 1):
                r[dr - db + i] = dom.sub(r[dr - db + i], dom.mul(lead, b[i]))
            r = _poly_trim(r, dom)
        a, b = b, r
    return _poly_trim(a, dom)


def minimal_polynomial(M, dom):
    """Minimal polynomial of the square matrix ``M`` (Domain entries) as a
    coefficient list ``[c_0, ..., c_d]`` (low->high, monic ``c_d = 1``). Found by
    the first linear dependency among ``I, M, M^2, ...`` -- exact linear algebra
    over ``dom`` (fields.linalg.solve), no sympy, any exact Domain."""
    from quiverlab.fields.linalg import solve
    m = len(M)
    M = [[dom.coerce(M[i][j]) for j in range(m)] for i in range(m)]

    def flat(X):
        return [X[i][j] for i in range(m) for j in range(m)]

    I = _identity(dom, m)
    cols = [flat(I)]          # powers M^0, M^1, ... flattened as columns
    cur = I
    for d in range(1, m + 1):
        cur = _matmul(dom, cur, M)
        v = flat(cur)
        Amat = [[cols[c][r] for c in range(len(cols))] for r in range(m * m)]
        x = solve(Amat, list(v), dom)
        if x is not None:                       # M^d = sum x_i M^i
            coeffs = [dom.neg(x[i]) for i in range(len(cols))]
            coeffs.append(dom.one())
            return coeffs
        cols.append(v)
    # I..M^m are always dependent (Cayley-Hamilton); we never fall through
    raise QuiverlabError("minimal polynomial not found within degree dim "
                         "(internal invariant) -- please report this matrix")


def nu_is_semisimple(nu_matrix, domain):
    """``True`` iff ``nu`` (a square matrix of Domain elements) is semisimple,
    decided EXACTLY over a perfect field by ``squarefree(minpoly)`` via
    ``gcd(minpoly, minpoly') is a unit``. Conclusive, no float."""
    dom = domain
    mp = minimal_polynomial(nu_matrix, dom)
    deriv = [dom.mul(dom.coerce(i), mp[i]) for i in range(1, len(mp))]
    deriv = _poly_trim(deriv, dom)
    if not deriv:                                # minpoly' == 0
        # In char p a non-constant minpoly with zero derivative is a poly in x^p,
        # hence (over a perfect field) a p-th power => NOT squarefree, unless the
        # minpoly is a unit (impossible: degree >= 1). So non-semisimple.
        return _poly_deg(mp, dom) <= 0
    g = _poly_gcd(mp, deriv, dom)
    return _poly_deg(g, dom) <= 0                # squarefree <=> gcd is a constant


def _nu_order(nu_matrix, dom):
    """The multiplicative order of ``nu`` (bounded power search, cap ``_ORDER_CAP``);
    ``None`` past the cap. Provenance only -- never load-bearing."""
    m = len(nu_matrix)
    M = [[dom.coerce(nu_matrix[i][j]) for j in range(m)] for i in range(m)]
    I = _identity(dom, m)
    cur = M
    for k in range(1, _ORDER_CAP + 1):
        if all(dom.eq(cur[i][j], I[i][j]) for i in range(m) for j in range(m)):
            return k
        cur = _matmul(dom, cur, M)
    return None


# ---------------------------------------------------------------------------
# the gate
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class BVHypothesis:
    applies: bool
    route: str | None            # "symmetric" | "semisimple" | None
    label: str                   # human provenance (the s2.4 strings)
    nu_semisimple: bool
    nu_order: int | None
    nu_inner: bool               # symmetric => True
    refusal: str | None          # the loud message when applies is False


_NOT_FROBENIUS = (
    "BV operator needs a Frobenius algebra whose Nakayama automorphism is "
    "semisimple (or a symmetric algebra); this algebra is not Frobenius -- no "
    "implemented BV hypothesis applies")

_BIKLZ_BLOCKED = (
    "BV operator: self-injective Nakayama with NON-semisimple Nakayama "
    "automorphism (char | ord nu) -- the general Bian-Itagaki-Kou-Lyu-Zhou "
    "construction (arXiv:2603.04834) is NOT implemented in v1 (out of scope); "
    "note the e = 1 case k[x]/(x^N) is symmetric and served by route (1)")

_FROB_NONSS = (
    "BV operator needs a Frobenius algebra whose Nakayama automorphism is "
    "semisimple (or a symmetric algebra); this algebra is Frobenius with a "
    "non-semisimple Nakayama automorphism (char | ord nu) -- no implemented BV "
    "hypothesis applies")


def _is_nakayama_quiver(A):
    """True iff the quiver is a Nakayama quiver: every vertex has at most one
    incoming and at most one outgoing arrow (linear A_n or a single cycle Z_e)."""
    q = A.quiver
    if q is None:
        return False
    outdeg = {v: 0 for v in q.vertices}
    indeg = {v: 0 for v in q.vertices}
    for (s, t) in q.arrows.values():
        outdeg[s] = outdeg.get(s, 0) + 1
        indeg[t] = indeg.get(t, 0) + 1
    return all(outdeg[v] <= 1 for v in q.vertices) and \
        all(indeg[v] <= 1 for v in q.vertices)


def classify_bv(A):
    """Classify ``A`` for the BV operator per spec s2.4, returning a
    :class:`BVHypothesis` (never raises for a presented algebra: an unsupported
    algebra yields ``applies=False`` with a loud ``refusal`` string)."""
    # 1. symmetric (nu inner) -- the P52-free Tradler anchor.
    if A.is_symmetric():
        return BVHypothesis(
            applies=True, route="symmetric",
            label="symmetric (Tradler AIF 2008)",
            nu_semisimple=True, nu_order=1, nu_inner=True, refusal=None)
    # not symmetric: needs Frobenius + semisimple nu.
    if not A.is_frobenius():
        return BVHypothesis(
            applies=False, route=None, label="not Frobenius",
            nu_semisimple=False, nu_order=None, nu_inner=False,
            refusal=_NOT_FROBENIUS)
    nu = A.nakayama_automorphism()
    ss = nu_is_semisimple(nu, A.domain)
    order = _nu_order(nu, A.domain)
    if ss:
        ord_txt = str(order) if order is not None else "finite"
        label = ("semisimple Nakayama automorphism (Lambre-Zhou-Zimmermann 2016; "
                 f"Volkov 2016 over GF(p): char does not divide ord nu = {ord_txt})")
        return BVHypothesis(
            applies=True, route="semisimple", label=label,
            nu_semisimple=True, nu_order=order, nu_inner=False, refusal=None)
    # Frobenius with a non-semisimple nu -> BIKLZ-blocked (self-inj Nakayama) or generic.
    if _is_nakayama_quiver(A) and A.is_selfinjective():
        refusal = _BIKLZ_BLOCKED
    else:
        refusal = _FROB_NONSS
    return BVHypothesis(
        applies=False, route=None,
        label="Frobenius, non-semisimple Nakayama automorphism",
        nu_semisimple=False, nu_order=order, nu_inner=False, refusal=refusal)
