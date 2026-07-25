"""Plan 20 Task 1: the double-PELT model of P_p ⊗_A P_q and its Koszul-signed
tensor differential — the ambient of the lifted diagonal Δ: P → P ⊗_A P (Task 2).

A DOUBLE-PELT element (total chain-degree n) is a dict

    (a_idx, tau_word, mid_idx, rho_word, c_idx) -> coeff

meaning  Σ coeff · ( b_{a_idx} ⊗ τ ⊗ b_{mid_idx} ⊗ ρ ⊗ b_{c_idx} )  in
P_p ⊗_A P_q  (p = deg τ, q = deg ρ, p+q = n), where b_• are A-basis elements and

    P_p ⊗_A P_q = ⊕_{τ∈S_p, ρ∈S_q}  A e_{o(τ)} ⊗ (e_{t(τ)} A e_{o(ρ)}) ⊗ e_{t(ρ)} A .

So b_{a_idx} ∈ A e_{o(τ)} (path ending at o(τ)), b_{mid_idx} ∈ e_{t(τ)} A e_{o(ρ)}
(corner path t(τ)→o(ρ)), b_{c_idx} ∈ e_{t(ρ)} A (path starting at t(ρ)).  A chain
word is a real path tuple for degree ≥ 1 and ("__v__", v) for the degree-0
vertex chain σ_v — the same convention as `resolutions_cs/pelt.py`, which this
module extends by gluing two single-PELTs along the middle A-factor.

The tensor differential is the Koszul-signed differential of a tensor product of
complexes (of A^e-projectives),

    d^{P⊗P}_n  =  d_p ⊗ 1  +  (−1)^p · 1 ⊗ d_q          on the P_p ⊗_A P_q summand,

the sign living on the SECOND summand only.  Each summand is applied by gluing a
`res.d_terms` output (coeff, a', tw, c') into a slot exactly as `pelt.apply_lower`
glues into a single-PELT:
  * slot 1 (d_p ⊗ 1):  a' multiplies a_idx from the RIGHT (b_a·a'),  c' multiplies
    mid_idx from the LEFT (c'·b_mid);
  * slot 2 (1 ⊗ d_q):  a' multiplies mid_idx from the RIGHT (b_mid·a'),  c'
    multiplies c_idx from the LEFT (c'·b_c), with the extra (−1)^p.

Domain-generic (all arithmetic through A.domain / AArith); exact only; no floats;
no engine imports.  Basis/matrix orderings are deterministic (sorted S-sequence,
ascending corner indices, nested iteration) so downstream Δ is byte-reproducible.
"""
from quiverlab.fields.linalg import reduce_mod_nullspace, solve
from quiverlab.resolutions_cs.pelt import _resolve_chain, _vecs, _accum


class TensorComplex:
    """d^{P⊗P} on double-PELTs over a fixed `ChouhySolotarResolution`.

    Public surface (Task 2 consumes all three):
      * ``basis(p, q, o, t)``      — the (o,t)-corner k-basis of P_p ⊗_A P_q;
      * ``tensor_matrix(n, o, t)`` — d^{P⊗P}_n on the (o,t)-corner of ⊕_{p+q=n};
      * ``apply_tensor_d(n, dpelt)`` — the differential of a double-PELT element.
    Plus ``generators(n, o, t)`` — the free A^e-generators of that corner.
    """

    def __init__(self, res):
        self.res = res
        self.ar = res.ar
        self.dom = res.dom
        self._chain_cache = {}
        self._basis_cache = {}
        self._full_basis_cache = {}
        self._matrix_cache = {}
        self._gen_cache = {}
        self._diag_cache = {}

    # -- chains -------------------------------------------------------------
    def _chain(self, word):
        ch = self._chain_cache.get(word)
        if ch is None:
            ch = _resolve_chain(self.res, word)
            self._chain_cache[word] = ch
        return ch

    @staticmethod
    def _chain_word(ch):
        """Storage word for a chain: real path for degree ≥ 1, ("__v__", v) for a
        vertex (matches pelt/d_terms target words)."""
        return ("__v__", ch.o) if ch.degree == 0 else ch.word

    def _vidx(self, v):
        """Basis index of the vertex idempotent e_v."""
        return self.ar._word_index[("v", v)]

    # -- k-bases ------------------------------------------------------------
    def basis(self, p, q, o, t):
        """Deterministic k-basis of the (o,t)-corner e_o (P_p ⊗_A P_q) e_t as
        (a_idx, tau_word, mid_idx, rho_word, c_idx) tuples:
        a_idx ∈ corner(o, o(τ)), mid_idx ∈ corner(t(τ), o(ρ)), c_idx ∈ corner(t(ρ), t)."""
        key = (p, q, o, t)
        cached = self._basis_cache.get(key)
        if cached is not None:
            return cached
        ar, res = self.ar, self.res
        out = []
        for tau in res.ss.S(p):
            tau_word = self._chain_word(tau)
            a_corner = ar.corner(o, tau.o, "coh")            # paths o -> o(τ)
            for rho in res.ss.S(q):
                rho_word = self._chain_word(rho)
                mid_corner = ar.corner(tau.t, rho.o, "coh")  # paths t(τ) -> o(ρ)
                c_corner = ar.corner(rho.t, t, "coh")        # paths t(ρ) -> t
                for ai in a_corner:
                    for mi in mid_corner:
                        for ci in c_corner:
                            out.append((ai, tau_word, mi, rho_word, ci))
        self._basis_cache[key] = out
        return out

    def _full_basis(self, n, o, t):
        """Concatenated k-basis of the (o,t)-corner of ⊕_{p+q=n} P_p ⊗_A P_q
        (p ascending) — the row/column index set of tensor_matrix."""
        if n < 0:
            return []
        key = (n, o, t)
        cached = self._full_basis_cache.get(key)
        if cached is not None:
            return cached
        out = []
        for p in range(n + 1):
            out.extend(self.basis(p, n - p, o, t))
        self._full_basis_cache[key] = out
        return out

    def generators(self, n, o, t):
        """Free A^e-generators of the (o,t)-corner of ⊕_{p+q=n} P_p ⊗_A P_q:
        vertex idempotents in the two outer slots, mid ranging over corner(t(τ),o(ρ)),
        with o(τ)=o and t(ρ)=t.  (d^{P⊗P})² = 0 on these ⟺ on the whole module."""
        key = (n, o, t)
        cached = self._gen_cache.get(key)
        if cached is not None:
            return cached
        ar, res = self.ar, self.res
        ev, et = self._vidx(o), self._vidx(t)
        out = []
        for p in range(n + 1):
            q = n - p
            for tau in res.ss.S(p):
                if tau.o != o:
                    continue
                tau_word = self._chain_word(tau)
                for rho in res.ss.S(q):
                    if rho.t != t:
                        continue
                    rho_word = self._chain_word(rho)
                    for mi in ar.corner(tau.t, rho.o, "coh"):
                        out.append((ev, tau_word, mi, rho_word, et))
        self._gen_cache[key] = out
        return out

    # -- the tensor differential -------------------------------------------
    def apply_tensor_d(self, n, dpelt):
        """d^{P⊗P}_n of a double-PELT element of total degree n, returned as a
        double-PELT of total degree n-1.  Koszul sign (−1)^p on the 1⊗d_q summand."""
        res, ar, dom = self.res, self.ar, self.dom
        one = dom.one()
        out = {}
        for (ai, tau_word, mi, rho_word, ci), coeff in dpelt.items():
            if dom.is_zero(coeff):
                continue
            tau, rho = self._chain(tau_word), self._chain(rho_word)
            p, q = tau.degree, rho.degree
            assert p + q == n, f"double-PELT key {(tau_word, rho_word)} has degree {p+q}, not {n}"
            a_vec = ar.A._basis_vec(ai)
            mid_vec = ar.A._basis_vec(mi)
            c_vec = ar.A._basis_vec(ci)

            # slot 1:  d_p ⊗ 1   (no d_0)
            if p >= 1:
                for (dc, a_word, tw, c_word) in res.d_terms(p, tau):
                    a2, c2 = _vecs(res, tau, a_word, c_word)
                    left = ar.mul(a_vec, a2)            # b_{a}·a'
                    newmid = ar.mul(c2, mid_vec)        # c'·b_{mid}
                    base = dom.mul(coeff, dc)
                    for aj, av in enumerate(left):
                        if dom.is_zero(av):
                            continue
                        for mj, mv in enumerate(newmid):
                            if dom.is_zero(mv):
                                continue
                            val = dom.mul(base, dom.mul(av, mv))
                            _accum(out, (aj, tw, mj, rho_word, ci), val, dom)

            # slot 2:  (−1)^p · 1 ⊗ d_q   (no d_0)
            if q >= 1:
                sign = one if p % 2 == 0 else dom.neg(one)
                signed = dom.mul(coeff, sign)
                for (dc, a_word, tw, c_word) in res.d_terms(q, rho):
                    a2, c2 = _vecs(res, rho, a_word, c_word)
                    newmid = ar.mul(mid_vec, a2)        # b_{mid}·a'
                    right = ar.mul(c2, c_vec)           # c'·b_{c}
                    base = dom.mul(signed, dc)
                    for mj, mv in enumerate(newmid):
                        if dom.is_zero(mv):
                            continue
                        for cj, cv in enumerate(right):
                            if dom.is_zero(cv):
                                continue
                            val = dom.mul(base, dom.mul(mv, cv))
                            _accum(out, (ai, tau_word, mj, tw, cj), val, dom)
        return out

    def tensor_matrix(self, n, o, t):
        """Dense matrix of d^{P⊗P}_n on the (o,t)-corner: rows indexed by
        _full_basis(n-1, o, t), columns by _full_basis(n, o, t).  Cached per
        (n, o, t).  Column j is apply_tensor_d applied to the j-th corner basis
        element; every image term lands in the same corner (a bimodule map), so it
        must index a row — an escape is a bug and raises."""
        key = (n, o, t)
        cached = self._matrix_cache.get(key)
        if cached is not None:
            return cached
        dom = self.dom
        cols = self._full_basis(n, o, t)
        rows = self._full_basis(n - 1, o, t)
        ridx = {tup: i for i, tup in enumerate(rows)}
        M = [[dom.zero()] * len(cols) for _ in range(len(rows))]
        for cj, tup in enumerate(cols):
            image = self.apply_tensor_d(n, {tup: dom.one()})
            for rtup, val in image.items():
                ri = ridx.get(rtup)
                if ri is None:
                    raise AssertionError(
                        f"d^{{P⊗P}}_{n} image term {rtup} escaped the ({o},{t})-corner "
                        f"row basis — a bug in the corner bookkeeping, never an approximation")
                M[ri][cj] = dom.add(M[ri][cj], val)
        self._matrix_cache[key] = M
        return M

    # -- the lifted diagonal Δ: P → P ⊗_A P (Plan 20 Task 2) ----------------
    def _zeta(self, n, sigma, prev):
        """RHS ζ(σ) = Δ_{n−1}(d_n σ) as a double-PELT of total degree n−1.

        `d_n σ = Σ (coeff, a', τ, c')` in P_{n−1}; Δ_{n−1} is an A^e-module (bimodule)
        map, so Δ_{n−1}(a'·τ·c') = a'·Δ_{n−1}(τ)·c'.  For each term the produced a'
        LEFT-multiplies the a-slot of Δ_{n−1}(τ) (b_a → a'·b_a) and c' RIGHT-multiplies
        the c-slot (b_c → b_c·c'); the middle slot and both chain words are untouched.
        This is the OUTER bimodule action — distinct from `apply_tensor_d`'s INTERIOR
        gluing, which sends a' into the a-slot from the right and c' into the mid-slot
        from the left.  ζ lands in the (o(σ), t(σ)) corner of ⊕_{p+q=n−1} P_p ⊗_A P_q."""
        res, ar, dom = self.res, self.ar, self.dom
        out = {}
        for (coeff, a_word, tw, c_word) in res.d_terms(n, sigma):
            if dom.is_zero(coeff):
                continue
            a_vec, c_vec = _vecs(res, sigma, a_word, c_word)   # a' : o(σ)→o(τ),  c' : t(τ)→t(σ)
            for (ai, tau_w, mi, rho_w, ci), kappa in prev[tw].items():
                base = dom.mul(coeff, kappa)
                if dom.is_zero(base):
                    continue
                left = ar.mul(a_vec, ar.A._basis_vec(ai))       # a'·b_a  (LEFT action)
                right = ar.mul(ar.A._basis_vec(ci), c_vec)      # b_c·c'  (RIGHT action)
                for aj, av in enumerate(left):
                    if dom.is_zero(av):
                        continue
                    for cj, cv in enumerate(right):
                        if dom.is_zero(cv):
                            continue
                        val = dom.mul(base, dom.mul(av, cv))
                        _accum(out, (aj, tau_w, mi, rho_w, cj), val, dom)
        return out

    def diagonal(self, n):
        """Δ_n as {chain-word: double-PELT} for σ ∈ S_n, cached and recursive.

        Base: Δ_0(σ_v) = e_v ⊗ σ_v ⊗ e_v ⊗ σ_v ⊗ e_v (corner idempotent indices).
        Step: per σ solve the lifting equation d^{P⊗P}_n · Δ_n(σ) = Δ_{n−1}(d_n σ)
        over the (o(σ),t(σ))-corner — coefficient matrix `tensor_matrix(n,o,t)`,
        RHS = ζ(σ) coordinates over `_full_basis(n−1,o,t)`; `solve` + `reduce_mod_nullspace`
        pin the canonical (free-variables-zero) representative, so Δ is byte-reproducible.
        An inconsistent solve is the identical scope edge `_d_general` flags — raise the
        same loud NotImplementedError, never a fallback."""
        cached = self._diag_cache.get(n)
        if cached is not None:
            return cached
        res, dom = self.res, self.dom
        out = {}
        if n == 0:
            one = dom.one()
            for sv in res.ss.S(0):
                w = self._chain_word(sv)
                ev = self._vidx(sv.o)
                out[w] = {(ev, w, ev, w, ev): one}
            self._diag_cache[0] = out
            return out
        prev = self.diagonal(n - 1)
        for sigma in res.ss.S(n):
            o, t = sigma.o, sigma.t
            zeta = self._zeta(n, sigma, prev)
            cols = self._full_basis(n, o, t)
            rows = self._full_basis(n - 1, o, t)
            ridx = {tup: i for i, tup in enumerate(rows)}
            rhs = [dom.zero()] * len(rows)
            for tup, val in zeta.items():
                i = ridx.get(tup)
                if i is None:
                    raise AssertionError(
                        f"ζ(σ) term {tup} escaped the ({o},{t})-corner row basis at "
                        f"degree {n}, chain {sigma.word} — a corner-bookkeeping bug, "
                        f"never an approximation")
                rhs[i] = dom.add(rhs[i], val)
            M = self.tensor_matrix(n, o, t)
            if not M:                                          # no equations: choose the zero lift
                x = [dom.zero()] * len(cols) if all(dom.is_zero(v) for v in rhs) else None
            else:
                x = solve(M, rhs, dom)
            if x is None:
                raise NotImplementedError(
                    f"diagonal lift-solve is inconsistent at degree {n}, chain "
                    f"{sigma.word}: this admissible algebra needs the higher CS homotopy "
                    f"correction for the diagonal, outside quiverlab v1's construction "
                    f"(spec §6 risk register)")
            x = reduce_mod_nullspace(x, M, dom)
            dpelt = {}
            for coeff, tup in zip(x, cols):
                if not dom.is_zero(coeff):
                    dpelt[tup] = coeff
            out[self._chain_word(sigma)] = dpelt
        self._diag_cache[n] = out
        return out


def diagonal(res, n):
    """Comparison-lifted diagonal Δ_n on the CS resolution `res`: {chain-word:
    double-PELT} for σ ∈ S_n.  Caches a `TensorComplex` on `res` so repeated calls
    (and the recursion) reuse degrees, and the result is byte-reproducible."""
    tc = getattr(res, "_tensor_complex", None)
    if tc is None:
        tc = TensorComplex(res)
        res._tensor_complex = tc
    return tc.diagonal(n)
