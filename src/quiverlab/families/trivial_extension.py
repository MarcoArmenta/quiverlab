"""Trivial extension T(A) = A |x D(A), D(A) = Hom_k(A, k) (spec §3.4). Basis:
a_0..a_{n-1} then a_0^*..a_{n-1}^*; (a,f)(b,g) = (ab, a.g + f.b), symmetric algebra."""
from quiverlab.core.algebra import Algebra


def TrivialExtension(A):
    dom = A.domain
    n = A.dim
    m = 2 * n
    zero, one = dom.zero(), dom.one()

    def basis(i):
        v = [zero] * n
        v[i] = one
        return v

    # left/right action of A on itself in coordinates: prod[i][k] = a_i * a_k
    prod = [[A.multiply(basis(i), basis(k)) for k in range(n)] for i in range(n)]
    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    # (a_i, 0)(a_k, 0) = (a_i a_k, 0)
    for i in range(n):
        for k in range(n):
            for t in range(n):
                T[i][k][t] = prod[i][k][t]
    # (a_i,0)(0, a_k^*) = (0, a_i . a_k^*) ; (a_i . g)(c) = g(c a_i) => coeff on a_l^* is [a_l a_i]_k
    for i in range(n):
        for k in range(n):
            for l in range(n):
                T[i][n + k][n + l] = prod[l][i][k]
    # (0, a_i^*)(a_k, 0) = (0, a_i^* . a_k) ; (f . b)(c) = f(b c) => coeff on a_l^* is [a_k a_l]_i
    for i in range(n):
        for k in range(n):
            for l in range(n):
                T[n + i][k][n + l] = prod[k][l][i]
    # (0, *)(0, *) = 0  (D(A) is a square-zero ideal) -- already zero.
    unit = [zero] * m
    for t in range(n):
        unit[t] = A.unit[t]
    la = A.basis_labels or [f"a{i}" for i in range(n)]
    labels = list(la) + [f"{lbl}*" for lbl in la]
    Text = Algebra(dom, T, unit, basis_labels=labels)
    Text._family_citations = ("assem_book",)
    return Text
