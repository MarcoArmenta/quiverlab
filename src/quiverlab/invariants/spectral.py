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
    """max_i |alpha_i| over the roots of poly, EXACT. 1 on cyclotomic input."""
    if poly is None:
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
    """|lc| * prod over roots with |alpha| > 1 of |alpha|, EXACT. 1 on cyclotomic input."""
    if poly is None:
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
