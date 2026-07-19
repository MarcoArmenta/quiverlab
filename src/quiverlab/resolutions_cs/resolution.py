"""The Chouhy-Solotar differential (arXiv:1406.2300 §4, §6), Domain-generic. A TERM is
(coeff, a_word, target_word, c_word) meaning coeff·(a ⊗ target ⊗ c) in P_{n-1}, where
a_word,c_word are PATHS (arrow-name tuples; () = the unit) reduced to normal form only at
collapse time, and target_word identifies a chain in S_{n-1}. Composition is LEFT TO RIGHT."""
from quiverlab.resolutions_cs.ambiguities import SSequence
from quiverlab.resolutions_cs.aarith import AArith


class ChouhySolotarResolution:
    def __init__(self, A, rs, max_degree, max_cells=4_000_000):
        self.A = A
        self.rs = rs
        self.dom = A.domain
        self.ss = SSequence(rs, max_degree, max_cells)
        self.ar = AArith(A, rs)
        self._chain_index = {}
        self._d_cache = {}

    def to_int(self, c):
        return self.dom.to_int(c) if hasattr(self.dom, "to_int") else int(c)

    def _one(self):
        return self.dom.one()

    def _neg(self, c):
        return self.dom.neg(c)

    def _chain(self, degree, word):
        idx = self._chain_index.get(degree)
        if idx is None:
            idx = {c.word: c for c in self.ss.S(degree)}
            self._chain_index[degree] = idx
        return idx.get(tuple(word))

    # -- Fox derivative φ_0 (CS §6), left-to-right --------------------------
    def _fox(self, word, coeff):
        return [(coeff, word[:k], (word[k],), word[k + 1:]) for k in range(len(word))]

    # -- leading map δ_n (CS f_{n-1}), quiverlab index ----------------------
    def delta_terms(self, n, chain):
        if n == 1:
            return self._d1_terms(chain)
        one, none = self._one(), self._neg(self._one())
        blocks = chain.blocks
        if n % 2 == 1:                                            # CS f_{n-1} EVEN: 2-term map
            u0, ulast = blocks[0], blocks[-1]
            P = tuple(x for blk in blocks[1:] for x in blk)      # u_1..u_{n-1} in S_{n-1}
            Q = tuple(x for blk in blocks[:-1] for x in blk)     # u_0..u_{n-2} in S_{n-1}
            return [(one, u0, P, ()), (none, (), Q, ulast)]
        prev = {c.word for c in self.ss.S(n - 1)}                # n even: CS f_{n-1} ODD, big sum
        w, out = chain.word, []
        for i in range(len(w)):
            for j in range(i + 1, len(w) + 1):
                if w[i:j] in prev:
                    out.append((one, w[:i], w[i:j], w[j:]))
        return out

    # -- d_1 (arrows → vertices), CS §6 -------------------------------------
    def _d1_terms(self, chain):
        one, none = self._one(), self._neg(self._one())
        name = chain.word[0]
        return [(one, (name,), ("__v__", chain.t), ()),         # α ⊗ e_t ⊗ 1
                (none, (), ("__v__", chain.o), (name,))]         # 1 ⊗ e_o ⊗ α

    # -- d_2 = φ_0(s) − φ_0(β(s)), CS §6 (β = FULL normal form, Pillar-4 fix) --
    def _d2_terms(self, chain):
        s = chain.word
        terms = list(self._fox(s, self._one()))                  # φ_0(s)
        for word, coeff in self.rs.normal_form(s).items():       # β(s) = Σ λ_i b_i (fully reduced dict)
            for (c, a, tw, cc) in self._fox(word, coeff):
                terms.append((self._neg(c), a, tw, cc))          # − φ_0(β(s))
        return terms

    def d_terms(self, n, chain):
        if n == 1:
            return self._d1_terms(chain)
        if n == 2:
            return self._d2_terms(chain)
        return self._d_general(n, chain)                         # Task 6

    def _d_general(self, n, chain):
        raise NotImplementedError("filled in Task 6 (order-condition-pinned correction)")
