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
        self._require_in_scope()      # public + ungated (edit #1): a direct odd-n≥3 call on a
        # non-quadratic non-monomial presentation would return silently-wrong even-δ terms
        # (left-vs-right decomposition divergence); refuse at the same boundary as _d_general.
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

    def _idx_word(self, idx):
        tag, val = self.ar.basis_word(idx)
        return () if tag == "v" else val

    def _lower_generators(self, n, chain):
        """Basis of ⟨\\overline{𝓛}_{n-1}^{≺}(σ)⟩: single-term lists (1, b, p.word, b') with
        p∈S_{n-1}, b,b'∈B, and the loop path b·p·b' strictly ≺ σ (parallel to σ). Finite by DCC."""
        order, skey, one = self.rs.order, self.rs.order.key(chain.word), self._one()
        gens = []
        for p in self.ss.S(n - 1):
            for bi in self.ar.corner(chain.o, p.o, "coh"):        # b : o(σ)->o(p)
                bw = self._idx_word(bi)
                for ci in self.ar.corner(p.t, chain.t, "coh"):    # b' : t(p)->t(σ)
                    cw = self._idx_word(ci)
                    if order.key(bw + p.word + cw) < skey:
                        gens.append((one, bw, p.word, cw))
        return gens

    def _require_in_scope(self):
        """RESTRICT (edit #1): the uncollapsed leading map is provably CS's δ only for
        quadratic tips (v_n = u_0; CS Prop, TeX 772) or monomial presentations (no
        correction; collapsed map = validated Bardzell). A non-quadratic non-monomial
        presentation raises NotImplementedError at this exact boundary (spec §6)."""
        quadratic = all(len(w) <= 2 for w in self.rs.leading_words())
        monomial = all(len(r.tail) == 0 for r in self.rs.rules)
        if not quadratic and not monomial:
            raise NotImplementedError(
                "CS is certified for quadratic tips (S ⊆ Q_2; CS Prop, arXiv:1406.2300 TeX "
                "line 772) or monomial presentations; this presentation has a tip of length "
                ">= 3 AND a nonzero tail (non-quadratic non-monomial). There v_n ≠ u_0 and the "
                "uncollapsed leading map feeding the correction solve is not provably δ_CS. "
                "Lifting this needs the right_decomposition upgrade (a Plan-04 stretch item); "
                "boundary per spec §6 risk register.")

    def _d_general(self, n, chain):
        from quiverlab.resolutions_cs.pelt import terms_to_pelt, apply_lower
        self._require_in_scope()                         # NotImplementedError at the exact boundary
        key = (n, chain.word)
        if key in self._d_cache:
            return self._d_cache[key]
        delta = self.delta_terms(n, chain)
        gens = self._lower_generators(n, chain)
        if not gens:
            self._d_cache[key] = delta
            return delta
        dom = self.dom
        rhs_pe = apply_lower(self, n, terms_to_pelt(self, delta))
        cols = [apply_lower(self, n, terms_to_pelt(self, [g])) for g in gens]
        keys = sorted(set(rhs_pe) | {k for col in cols for k in col})
        M = [[col.get(k, dom.zero()) for col in cols] for k in keys]
        rhs = [dom.neg(rhs_pe.get(k, dom.zero())) for k in keys]
        gamma = self._solve(M, rhs, len(gens))
        if gamma is None:
            raise NotImplementedError(
                f"CS correction linear system is inconsistent at degree {n}, chain "
                f"{chain.word}: this admissible algebra needs the higher CS homotopy "
                f"correction, outside quiverlab v1's construction (spec §6 risk register)")
        terms = list(delta)
        for coeff, (c1, a, tw, cc) in zip(gamma, gens):
            if not dom.is_zero(coeff):
                terms.append((dom.mul(coeff, c1), a, tw, cc))
        self._d_cache[key] = terms
        return terms

    def _solve(self, M, rhs, ncols):
        """Return γ (any solution of M·γ = rhs over the domain) or None if inconsistent.
        Uses fields.linalg.solve; if that requires square/consistent input, fall back to a
        particular solution (augmented RREF) + nullspace freedom (exact, over the domain)."""
        from quiverlab.fields.linalg import solve
        if not M:                                        # no equations -> any γ; choose zeros
            return [self.dom.zero()] * ncols
        return solve(M, rhs, self.dom)                   # returns a solution or None (inconsistent)

    def _cochains(self, n):
        # cochains beyond the certified resolution range are empty by construction; the
        # differential into them is the zero map. The coh differential deliberately reads
        # one degree past max_degree (matrix(max_degree,"coh") touches S(max_degree+1)); make
        # that empty cochain space explicit here so SSequence.S can raise loudly out of range.
        if n > self.ss.max_degree:
            return []
        return self.ss.S(n)

    def _basis(self, n, side):
        return [(ch, j) for ch in self._cochains(n) for j in self.ar.corner(ch.o, ch.t, side)]

    def dim_C(self, n, side):
        return len(self._basis(n, side))

    def matrix(self, n, side):
        from quiverlab.resolutions_cs.pelt import _resolve_chain, _vecs
        dom = self.dom
        if side == "hom":
            rows, cols = self._basis(n - 1, "hom"), self._basis(n, "hom")
        else:
            rows, cols = self._basis(n + 1, "coh"), self._basis(n, "coh")
        ridx = {(ch.word, j): i for i, (ch, j) in enumerate(rows)}
        M = [[dom.zero()] * len(cols) for _ in range(len(rows))]
        if side == "hom":
            for cj, (sigma, j) in enumerate(cols):
                ej = self.ar.A._basis_vec(j)
                for (coeff, a_word, tw, c_word) in self.d_terms(n, sigma):
                    a_vec, c_vec = _vecs(self, sigma, a_word, c_word)
                    val = self.ar.mul(c_vec, self.ar.mul(ej, a_vec))     # b·w·a  (homology collapse)
                    tw_word = _resolve_chain(self, tw).word
                    for p, vp in enumerate(val):
                        if not dom.is_zero(vp) and (tw_word, p) in ridx:
                            r = ridx[(tw_word, p)]
                            M[r][cj] = dom.add(M[r][cj], dom.mul(coeff, vp))
        else:
            for cj, (sigma, j) in enumerate(cols):                       # δ^n: C^n -> C^{n+1}
                ej = self.ar.A._basis_vec(j)
                for tau in self._cochains(n + 1):                        # empty beyond cap (see _cochains)
                    for (coeff, a_word, tw, c_word) in self.d_terms(n + 1, tau):
                        if _resolve_chain(self, tw).word != sigma.word:
                            continue
                        a_vec, c_vec = _vecs(self, tau, a_word, c_word)
                        val = self.ar.mul(a_vec, self.ar.mul(ej, c_vec)) # a·w·b  (cohomology collapse)
                        for p, vp in enumerate(val):
                            if not dom.is_zero(vp) and (tau.word, p) in ridx:
                                r = ridx[(tau.word, p)]
                                M[r][cj] = dom.add(M[r][cj], dom.mul(coeff, vp))
        return M

    def _matmul(self, A_, B_):
        dom = self.dom
        if not A_ or not B_:
            return []
        inner, rows_out, cols_out = len(B_), len(A_), len(B_[0])
        out = [[dom.zero()] * cols_out for _ in range(rows_out)]
        for i in range(rows_out):
            for k in range(inner):
                aik = A_[i][k]
                if dom.is_zero(aik):
                    continue
                for j in range(cols_out):
                    out[i][j] = dom.add(out[i][j], dom.mul(aik, B_[k][j]))
        return out

    def assert_dd_zero(self, upto, side):
        dom = self.dom
        for n in range(2, upto + 1):
            if side == "hom":
                prod = self._matmul(self.matrix(n - 1, "hom"), self.matrix(n, "hom"))
            else:
                prod = self._matmul(self.matrix(n, "coh"), self.matrix(n - 1, "coh"))
            if any(not dom.is_zero(x) for row in prod for x in row):
                raise AssertionError(
                    f"CS d²≠0 at degree {n} (side={side}); the correction solve failed to "
                    "close — a bug, never an approximation")

    def assert_order_condition(self, upto):
        """CS Theorem 4.1 condition (2): every correction term's loop path b·p·b' is strictly ≺ σ."""
        order = self.rs.order
        for n in range(3, upto + 1):
            for sigma in self.ss.S(n):
                lead = {(a, tw, cc): c for (c, a, tw, cc) in self.delta_terms(n, sigma)}
                for (c, a, tw, cc) in self.d_terms(n, sigma):
                    if lead.get((a, tw, cc)) == c:
                        continue                                 # a leading term
                    if not order.key(a + tw + cc) < order.key(sigma.word):
                        raise AssertionError(
                            f"order condition violated at degree {n}, chain {sigma.word}: "
                            f"correction {(a, tw, cc)} is not ≺ σ")
