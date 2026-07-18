"""Lower a general (non-monomial) kQ/I to a Plan-01 structure-constant Algebra
(spec §3.3, §5 components 3-4). The basis is the certified irreducible-path basis
from the reduction system; structure constants are products reduced to normal form.

The basis is ordered EXACTLY as core.monomial.build_monomial_algebra: trivial
vertex idempotents first (in quiver.vertices order), then irreducible words sorted
(len, word). On a monomial input the reduction system's forbidden words equal the
relation words and reduction sends every reducible product to 0, so the resulting
Algebra is elementwise-equal to the monomial route (a required cross-check)."""
from quiverlab.core.algebra import Algebra
from quiverlab.groebner.system import build_reduction_system


def groebner_algebra(quiver, relations, field, degree_bound=None, trace=None):
    rs = build_reduction_system(quiver, relations, field,
                                degree_bound=degree_bound, trace=trace)
    dom = rs.domain
    zero, one = dom.zero(), dom.one()

    basis = [("e", v) for v in quiver.vertices] + [("p", w) for w in rs.irreducibles]
    index = {b: i for i, b in enumerate(basis)}
    m = len(basis)

    def src(b):
        return b[1] if b[0] == "e" else quiver.word_source(b[1])

    def tgt(b):
        return b[1] if b[0] == "e" else quiver.word_target(b[1])

    def product_vector(bi, bj):
        """Coordinate vector of b_i * b_j (concatenate, then reduce to normal form)."""
        vec = [zero] * m
        if tgt(bi) != src(bj):
            return vec
        if bi[0] == "e":                       # e * y = y
            vec[index[bj]] = one
            return vec
        if bj[0] == "e":                       # x * e = x
            vec[index[bi]] = one
            return vec
        nf = rs.reduce({bi[1] + bj[1]: one})   # normal form of the concatenated word
        for word, coeff in nf.items():
            vec[index[("p", word)]] = coeff
        return vec

    T = [[product_vector(bi, bj) for bj in basis] for bi in basis]
    unit = [zero] * m
    for v in quiver.vertices:
        unit[index[("e", v)]] = one
    labels = [f"e_{b[1]}" if b[0] == "e" else "*".join(b[1]) for b in basis]
    return Algebra(dom, T, unit, basis_labels=labels, _quiver=quiver,
                   _relations=list(relations))
