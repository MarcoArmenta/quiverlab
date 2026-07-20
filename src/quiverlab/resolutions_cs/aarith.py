"""A-arithmetic for the CS resolution over any Domain. A path word maps to its normal-form
A-vector by FOLDING core.Algebra.multiply over the arrow generators (A's structure constants
already encode reduction mod I, so this IS the map β/π). The vertex/irreducible basis order
matches groebner_algebra: vertices first, then rs.irreducibles."""


class AArith:
    def __init__(self, A, rs):
        self.A = A
        self.rs = rs
        self.dom = A.domain
        Q = rs.quiver
        self._vertices = list(Q.vertices)
        self._basis_words = [("v", v) for v in self._vertices] + [("p", w) for w in rs.irreducibles]
        self._word_index = {w: i for i, w in enumerate(self._basis_words)}
        self._arrow_vec = {name: self._lookup_vec(("p", (name,))) for name in Q.arrows}
        self._vertex_idem = {v: self._lookup_vec(("v", v)) for v in self._vertices}

    def _lookup_vec(self, tagged_word):
        v = [self.dom.zero()] * self.A.dim
        v[self._word_index[tagged_word]] = self.dom.one()
        return v

    def basis_word(self, idx):
        return self._basis_words[idx]

    def mul(self, u, v):
        return self.A.multiply(u, v)

    def vertex_vec(self, vtx):
        return list(self._vertex_idem[vtx])

    def path_vec(self, word):
        """Normal-form A-vector of a NONEMPTY path (arrow-name tuple), left-to-right.
        Callers use vertex_vec for the empty path e_v."""
        if len(word) == 0:
            raise ValueError("path_vec needs a nonempty path; use vertex_vec for e_v")
        acc = list(self._arrow_vec[word[0]])
        for name in word[1:]:
            acc = self.A.multiply(acc, self._arrow_vec[name])
        return acc

    def corner(self, o, t, side):
        """Basis indices j with b_j in the corner: side="hom" -> e_t A e_o (paths t->o);
        side="coh" -> e_o A e_t (paths o->t)."""
        left, right = (t, o) if side == "hom" else (o, t)
        el, er = self._vertex_idem[left], self._vertex_idem[right]
        return [j for j in range(self.A.dim)
                if self.A.multiply(el, self.A.multiply(self.A._basis_vec(j), er)) == self.A._basis_vec(j)]
