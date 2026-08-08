"""Exact spectral radius and Mahler measure of a Coxeter polynomial (spec section 5
component 8). Reimplements the deleted hanlab float layer (`_roots_abs`, `spectral_radius`,
`mahler_measure`, which used mpmath `Poly.nroots`) with EXACT sympy algebraic numbers.

No floats in this module (tests/test_no_floats.py enforces): magnitudes are `sympy.Abs`
of exact `CRootOf` roots, comparisons use `.is_positive`, and cyclotomic input short-
circuits to the exact integer 1. The bank floats survive only as test oracles.

SOUNDNESS. We must count every root of modulus > 1, including COMPLEX ones. `real_roots`
alone is NOT enough -- it silently drops complex off-circle roots (t^4-7t^3+16t^2-7t+1 has
a complex pair of modulus 3.5465 that real_roots misses); "real roots suffice" is a theorem
only for tree/bipartite = HEREDITARY quivers (A'Campo). Full `all_roots` is exact but
pathologically slow on a high-degree irreducible factor (Lehmer, > 2 min). So we DECIDE,
with no complex-root isolation, whether real roots suffice, via the self-inversive
y = z + 1/z substitution:

  * cyclotomic short-circuit -> 1;
  * q = non-cyclotomic part (reciprocal, even degree 2s, no roots at +-1); q has a complex
    off-circle root iff Q(y) = q(z)/z^s (y = z + 1/z, degree s, via Dickson polynomials)
    has a non-real root -- decided by counting distinct real roots of Q (Sturm) against Q's
    squarefree degree. Equal -> real_roots(q) is exact and fast (Branch A: hereditary /
    Salem / Lehmer); unequal -> all_roots(q) fallback (Branch B: non-hereditary, q small).

spectral_radius = max |root| over q's off-circle roots; mahler = |lc| * prod of |root| over
roots with |root| > 1."""
import sympy as sp

from quiverlab.engine.coxeter_spectrum import is_cyclotomic_product

_T = sp.symbols("t")
_Y = sp.symbols("y_spec")


def _as_expr(poly):
    """Accept a `sympy.Poly` (what `Algebra.coxeter_polynomial()` returns, in the same
    `t = Symbol('t')` this module's `_T` uses) as well as a bare expr. A Poly has no
    `.expand`, so `sp.expand(Poly)` raises a raw `AttributeError` downstream -- convert
    it to its expr up front so `spectral_radius(A.coxeter_polynomial())` just works. None
    and plain exprs pass through unchanged.

    A `Poly` whose generator is NOT the expected Coxeter variable `_T` is out of contract
    (nothing in the repo produces one -- `coxeter_polynomial` always builds in `t`). Refuse
    it LOUDLY here rather than `.as_expr()`-ing a foreign variable and drifting into a
    silent `None` from the degree/cyclotomic machinery downstream."""
    if isinstance(poly, sp.Poly):
        if poly.gens != (_T,):
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(
                f"spectral_radius/mahler_measure expect a polynomial in the Coxeter "
                f"variable {_T!r}; got a sympy.Poly in {poly.gens}",
                hint="pass Algebra.coxeter_polynomial() (a Poly in t) or a bare sympy "
                     "expression in t")
        return poly.as_expr()
    return poly


def _degree_below_one(poly):
    """True for a constant or zero polynomial (degree < 1). Such an input has no
    reciprocal-root structure (the y = z + 1/z machinery assumes even degree >= 2), so
    both public entry points guard on it (mirrors coxeter_spectrum's `P.degree() == 0`
    branch); sp.Poly(0, _T).degree() is -oo, also < 1."""
    return sp.Poly(sp.expand(poly), _T).degree() < 1


def _noncyclotomic_part(poly):
    """Product (with multiplicity) of the non-cyclotomic irreducible factors of poly, or
    None if poly is a product of cyclotomics."""
    q = sp.Integer(1)
    for fac, mult in sp.factor_list(sp.expand(poly), _T)[1]:
        if not is_cyclotomic_product(fac):
            q = q * sp.Poly(fac, _T).as_expr() ** mult
    return sp.Poly(q, _T) if q != sp.Integer(1) else None


def _is_reciprocal(qp):
    c = qp.all_coeffs()
    return c == c[::-1] or c == [-x for x in c[::-1]]


def _reciprocal_to_y(qp):
    """For a reciprocal q of even degree 2s: Q(y) with q(z) = z^s * Q(z + 1/z), built from
    Dickson polynomials D_k(y) = z^k + z^{-k} (D_0=2, D_1=y, D_k = y*D_{k-1} - D_{k-2});
    q(z)/z^s = a_s + sum_{k=1..s} a_{s+k} D_k(y)."""
    a = list(reversed(qp.all_coeffs()))              # a[i] = coeff of z^i
    s = qp.degree() // 2
    D = [sp.Integer(2), _Y]
    for k in range(2, s + 1):
        D.append(sp.expand(_Y * D[k - 1] - D[k - 2]))
    Q = sp.Integer(a[s])
    for k in range(1, s + 1):
        Q = Q + a[s + k] * D[k]
    return sp.Poly(sp.expand(Q), _Y)


def _real_roots_suffice(qp):
    """True iff every off-circle root of the reciprocal q is REAL -- decided by counting
    DISTINCT real roots of Q(y) (Sturm on the squarefree part) against its squarefree
    degree, with no complex-root isolation. Non-reciprocal / odd-degree q (non-Coxeter
    input) -> False (be safe: fall back to all_roots)."""
    if qp.degree() % 2 or not _is_reciprocal(qp):
        return False
    Qsf = _reciprocal_to_y(qp).sqf_part()            # squarefree: count_roots == distinct
    return Qsf.count_roots() == Qsf.degree()


def _off_circle_roots(qp):
    """Exact roots of q with |z| > 1: real_roots when they provably suffice, else all_roots
    on the (small) non-cyclotomic factor q."""
    roots = sp.real_roots(qp) if _real_roots_suffice(qp) else qp.all_roots()
    return [r for r in roots if (sp.Abs(r) - 1).is_positive]


def spectral_radius(poly):
    """max_i |alpha_i| over the roots of poly, EXACT. 1 on cyclotomic input; None on a
    constant/zero polynomial (degree < 1), matching the poly-is-None convention."""
    if poly is None:
        return None
    poly = _as_expr(poly)
    if _degree_below_one(poly):
        return None
    if is_cyclotomic_product(poly):
        return sp.Integer(1)
    q = _noncyclotomic_part(poly)
    if q is None:
        return sp.Integer(1)
    best = sp.Integer(1)
    for r in _off_circle_roots(q):
        if (sp.Abs(r) - best).is_positive:
            best = sp.Abs(r)
    return best


def mahler_measure(poly):
    """|lc| * prod over roots with |alpha| > 1 of |alpha|, EXACT. 1 on cyclotomic input;
    None on a constant/zero polynomial (degree < 1), matching the poly-is-None convention."""
    if poly is None:
        return None
    poly = _as_expr(poly)
    if _degree_below_one(poly):
        return None
    if is_cyclotomic_product(poly):
        return sp.Integer(1)
    q = _noncyclotomic_part(poly)
    if q is None:
        return sp.Integer(1)
    m = sp.Abs(sp.Poly(poly, _T).LC())               # rational-LC-safe (never Integer(abs()))
    for r in _off_circle_roots(q):
        m = m * sp.Abs(r)
    return sp.simplify(m)


# =====================================================================
# Plan 58: certification + labelling + count layer (NO new root-finding;
# wraps/reuses spectral_radius, mahler_measure, is_cyclotomic_product and the
# private _noncyclotomic_part / _off_circle_roots helpers above). Every value is
# EXACT -- minimal polynomial + rational isolating interval, never a float.
# =====================================================================

def certify_real_algebraic(alpha) -> dict:
    """Certify an exact non-negative REAL sympy value ``alpha`` (what
    ``spectral_radius`` / ``mahler_measure`` return on a real-dominant spectrum:
    ``Integer`` / ``Rational`` / a radical ``Add``/``Pow`` / a real ``CRootOf``) as
    an algebraic number: its minimal polynomial over ``QQ`` (primitive ``ZZ[x]``,
    high->low degree), a RATIONAL isolating interval, a Sturm-certified unique root
    index, and the exact ``latex``/``value`` display forms.  Never a float; the
    rational interval IS the numeric localisation (§D2).

    CALLER CONTRACT (§D2/§D3): ``alpha`` must be a real radical/CRootOf/rational.
    ``sympy.minimal_polynomial`` is intractable (~121 s, measured) on the
    complex-modulus form ``sqrt(CRootOf*CRootOf)`` a complex-dominant spectrum
    yields, and it does NOT raise there -- it hangs -- so the assembler gates on the
    shipped ``_real_roots_suffice`` predicate BEFORE calling this and never hands a
    complex-modulus value in.  This function's only loud refusal is the
    non-unique-location guard below.

    Raises ``QuiverlabError`` if ``alpha`` is ``None``, if ``minimal_polynomial``
    raises, or if the root cannot be uniquely located among the isolating
    intervals (never a guessed root)."""
    from quiverlab.errors import QuiverlabError
    if alpha is None:
        raise QuiverlabError(
            "certify_real_algebraic: no value to certify (alpha is None)",
            hint="the Coxeter polynomial has degree < 1")
    x = sp.Symbol("x")
    if alpha.is_Rational:
        # minimal polynomial of p/q over QQ, primitive ZZ[x]: q*x - p
        p, qd = int(alpha.p), int(alpha.q)
        return {
            "minpoly": [qd, -p],
            "degree": 1,
            "interval": [str(alpha), str(alpha)],
            "root_index": 0,
            "is_rational": True,
            "latex": sp.latex(alpha),
            "value": str(alpha),
        }
    try:
        m = sp.minimal_polynomial(alpha, x, polys=True)
    except Exception as exc:                          # pragma: no cover - defensive
        raise QuiverlabError(
            "certify_real_algebraic: minimal_polynomial failed on %r" % (alpha,),
            hint="alpha must be a real algebraic number (radical / real CRootOf / "
                 "rational)") from exc
    ivs = m.intervals()                               # ascending, disjoint, rational
    hits = []
    for k, (endpoints, _mult) in enumerate(ivs):
        a_k, b_k = endpoints
        if (alpha - a_k).is_nonnegative and (b_k - alpha).is_nonnegative:
            hits.append((k, a_k, b_k))
    if len(hits) != 1:
        raise QuiverlabError(
            "certify_real_algebraic: could not uniquely locate %r among the %d "
            "isolating intervals (%d candidate hits)" % (alpha, len(ivs), len(hits)),
            hint="alpha must be a real root of its minimal polynomial")
    k, a_k, b_k = hits[0]
    return {
        "minpoly": [int(c) for c in m.all_coeffs()],
        "degree": m.degree(),
        "interval": [str(a_k), str(b_k)],
        "root_index": k,
        "is_rational": False,
        "latex": sp.latex(alpha),
        "value": str(alpha),
    }


def cyclotomic_factorization(poly) -> list:
    """Factor ``poly`` over ``ZZ`` and label each monic irreducible factor with its
    cyclotomic index (§D1).  Returns ``[{factor, latex, multiplicity,
    cyclotomic_index: int|None}, ...]`` -- one entry per ``(factor, multiplicity)``
    pair from ``sympy.factor_list`` (repeated cyclotomic factors keep their
    multiplicity).

    For a factor of degree ``d`` the search runs ``n in {1, ..., 2 d^2}`` for
    ``totient(n) == d`` and ``Poly(factor) == cyclotomic_poly(n)``: since
    ``phi(n) >= sqrt(n/2)``, ``phi(n) = d`` forces ``n <= 2 d^2``, so the bound is
    exact -- consistent with (and slightly tighter than) the ``range(1, 2 d^2 + 3)``
    ``is_cyclotomic_product`` already searches.  The search starts at ``n = 1`` so
    ``Phi_1 = t - 1`` and ``Phi_2 = t + 1`` are labelled like any other index; a bare
    factor ``t`` (root 0) matches no cyclotomic polynomial and is labelled ``None``."""
    poly = _as_expr(poly)
    _const, factors = sp.factor_list(sp.expand(poly), _T)
    out = []
    for fac, mult in factors:
        fp = sp.Poly(fac, _T)
        if fp.LC() < 0:                               # normalise sign (monic)
            fp = sp.Poly(-fp.as_expr(), _T)
        d = fp.degree()
        index = None
        for nn in range(1, 2 * d * d + 1):
            if sp.totient(nn) == d and fp == sp.Poly(sp.cyclotomic_poly(nn, _T), _T):
                index = nn
                break
        out.append({
            "factor": str(fp.as_expr()),
            "latex": sp.latex(fp.as_expr()),
            "multiplicity": int(mult),
            "cyclotomic_index": index,
        })
    return out


def off_circle_root_count(poly) -> int:
    """Exact count of roots of ``poly`` with ``|z| > 1`` (§D4): ``0`` on a
    constant / cyclotomic polynomial (all roots on the unit circle), else
    ``len(_off_circle_roots(_noncyclotomic_part(poly)))`` -- the shipped, sound
    helper (it counts complex off-circle roots too, via its ``all_roots``
    fallback)."""
    poly = _as_expr(poly)
    if _degree_below_one(poly):
        return 0
    q = _noncyclotomic_part(poly)
    if q is None:                                     # product of cyclotomics
        return 0
    return len(_off_circle_roots(q))
