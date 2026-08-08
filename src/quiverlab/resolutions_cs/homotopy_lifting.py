"""Plan 51: the Negron-Witherspoon / Volkov homotopy-lifting tower ψ_η on the
Chouhy-Solotar resolution -- the past-window, any-Domain carrier of the
Gerstenhaber bracket (assembled in `bracket.py`).

A CS cochain η ∈ C^n is a bimodule cocycle η: P_n → A (coordinates over
`res._basis(n, "coh")`).  Its HOMOTOPY LIFTING ψ_η is the degree -(n-1) bimodule
map ψ_η: P_k → P_{k-n+1} solving, as an identity of Hom-complex differentials
(Oke arXiv:2103.12331 Def 3.1 / Volkov arXiv:1610.05741 Def 4),

        d_P ∘ ψ_η  −  (−1)^{n-1} ψ_η ∘ d_P   =   (η ⊗ 1_P − 1_P ⊗ η) ∘ Δ_P ,

restricted on a free generator σ ∈ S_k (ψ_η(σ) the unknown PELT of the
(o(σ), t(σ))-corner of P_{k-n+1}) to the FULL corner linear system

    D_corner(k-n+1, o, t) · ψ_η(σ) = [(η⊗1 − 1⊗η)Δ_k(σ)] + (−1)^{n-1} ψ_η(d_k σ)

    * D_corner(m, o, t): the per-corner RESOLUTION-differential matrix of
      d_P: P_m → P_{m-1} (the single-complex sibling of
      diagonal.TensorComplex.tensor_matrix; assembled off `res.d_terms` through
      pelt.apply_lower, NOT res.matrix which is the Hom-dual coboundary);
    * RHS: the diagonal half-collapse (η eats a degree-n tensor factor of the
      shipped Δ, the surviving factor lands in P_{k-n}) plus the already-solved
      lower ψ term (the outer bimodule action, the sibling of TensorComplex._zeta);
    * `fields.linalg.solve` on the FULL corner system (like diagonal.diagonal),
      then `reduce_mod_nullspace` pins the free-variables-zero canonical
      representative -> ψ is byte-reproducible.

The RHS is a d_P-boundary at every degree Δ was built (Volkov Lemma 2: for a
cocycle η, (η⊗1 − 1⊗η) is null-homotopic on P⊗_A P), so the solve is ALWAYS
consistent there -- the `solve is None` branch is a defensive AssertionError
("cannot occur once Δ is built; a bug, never an approximation"), NOT a
user-facing scope boundary (Plan 51 DD1a).  The ONE genuine CS scope edge is Δ's
own NotImplementedError, raised while BUILDING Δ (before any ψ-solve).  A
non-cocycle input raises a DISTINCT typed CocycleError (never Δ's message).

Domain-generic (all arithmetic through res.dom / AArith); exact only; no floats;
no engine imports.  Determinism follows from Δ's byte-reproducibility and the
deterministic `_basis`/`corner` orderings.
"""
from quiverlab.errors import CocycleError
from quiverlab.fields.linalg import reduce_mod_nullspace, solve
from quiverlab.resolutions_cs.cup import _cochain_evaluator
from quiverlab.resolutions_cs.diagonal import diagonal
from quiverlab.resolutions_cs.pelt import _accum, _resolve_chain, _vecs, apply_lower


# --------------------------------------------------------------------------- #
# the single-complex corner machinery (sibling of diagonal.TensorComplex)     #
# --------------------------------------------------------------------------- #
class CornerComplex:
    """d_P on single-PELTs over a fixed ChouhySolotarResolution -- the
    single-complex sibling of `diagonal.TensorComplex`.  A single-PELT element of
    chain-degree m is a dict (a_idx, chain_word, c_idx) -> coeff (pelt.py's model),
    and this caches the per-corner k-basis and the per-corner differential matrix.
    """

    def __init__(self, res):
        self.res = res
        self.ar = res.ar
        self.dom = res.dom
        self._basis_cache = {}
        self._matrix_cache = {}

    def p_basis(self, m, o, t):
        """Deterministic k-basis of the (o,t)-corner e_o P_m e_t as
        (a_idx, tau_word, c_idx) tuples: a_idx ∈ corner(o, o(τ)), c_idx ∈
        corner(t(τ), t) for τ ∈ S_m.  The chain word is `tau.word` -- the SAME
        single-PELT convention `pelt.apply_lower`/`terms_to_pelt` use (`()` for a
        degree-0 vertex chain; the vertex is fixed by the corner indices), so the
        matrix rows/cols line up with apply_lower's output keys."""
        if m < 0:
            return []
        key = (m, o, t)
        cached = self._basis_cache.get(key)
        if cached is not None:
            return cached
        ar, res = self.ar, self.res
        out = []
        for tau in res.ss.S(m):
            for ai in ar.corner(o, tau.o, "coh"):
                for ci in ar.corner(tau.t, t, "coh"):
                    out.append((ai, tau.word, ci))
        self._basis_cache[key] = out
        return out

    def apply_d(self, m, pelt):
        """d_m: P_m → P_{m-1} on a single-PELT of chain-degree m.  `apply_lower`
        applies d_{N-1} to a chain-degree-(N-1) PELT, so d_m = apply_lower(·, m+1)."""
        return apply_lower(self.res, m + 1, pelt)

    def d_matrix(self, m, o, t):
        """Dense matrix of d_m on the (o,t)-corner: columns = p_basis(m, o, t),
        rows = p_basis(m-1, o, t).  The single-complex sibling of
        TensorComplex.tensor_matrix; an image term outside the corner row basis is
        a bookkeeping bug and raises (d_P is a bimodule map -> corners preserved)."""
        key = (m, o, t)
        cached = self._matrix_cache.get(key)
        if cached is not None:
            return cached
        dom = self.dom
        cols = self.p_basis(m, o, t)
        rows = self.p_basis(m - 1, o, t)
        ridx = {tup: i for i, tup in enumerate(rows)}
        M = [[dom.zero()] * len(cols) for _ in range(len(rows))]
        for cj, tup in enumerate(cols):
            image = self.apply_d(m, {tup: dom.one()})
            for rtup, val in image.items():
                ri = ridx.get(rtup)
                if ri is None:
                    raise AssertionError(
                        f"d_P({m}) image term {rtup} escaped the ({o},{t})-corner "
                        f"row basis -- a corner-bookkeeping bug, never an approximation")
                M[ri][cj] = dom.add(M[ri][cj], val)
        self._matrix_cache[key] = M
        return M


def _corner_complex(res):
    cc = getattr(res, "_corner_complex", None)
    if cc is None:
        cc = CornerComplex(res)
        res._corner_complex = cc
    return cc


def D_corner(res, m, o, t):
    """The per-corner resolution-differential matrix of d_P: P_m → P_{m-1} on the
    (o,t)-corner (Plan 51 Task A.0): columns = the corner k-basis of P_m, rows =
    the corner k-basis of P_{m-1}, entries assembled by applying `res.d_terms` glued
    through `pelt.apply_lower`.  This is NOT `res.matrix(·, "coh")` (the Hom-dual
    coboundary δ: C^n→C^{n+1}, the wrong object and direction) -- it is the
    single-complex sibling of `diagonal.TensorComplex.tensor_matrix`.
    Cached per (res, m, o, t); certified by D_corner(m-1)·D_corner(m) == 0."""
    return _corner_complex(res).d_matrix(m, o, t)


def _tensor_complex(res):
    """The cached TensorComplex on res (shared with diagonal()); a degree-0 build
    guarantees it exists.  Provides _chain / _chain_word used by the half-collapse."""
    diagonal(res, 0)
    return res._tensor_complex


# --------------------------------------------------------------------------- #
# the homotopy-lifting tower ψ_η                                              #
# --------------------------------------------------------------------------- #
class HomotopyLifting:
    """The canonical homotopy lifting ψ_η of a CS cocycle η ∈ C^n, bound to `res`
    and cached degreewise.  `apply(k, σ)` is ψ_η(σ) ∈ P_{k-n+1} as a single-PELT;
    `tower(k)` is the whole degree-k layer {chain-word: PELT}.  `image_cochain(k)`
    is the raw ψ layer keyed for a caller wanting the whole degree at once."""

    def __init__(self, res, eta_vec, n):
        if n < 1:
            raise ValueError(
                f"homotopy_lifting needs a positive cochain degree n >= 1; got {n}")
        self.res = res
        self.ar = res.ar
        self.dom = res.dom
        self.n = n
        self.eta_vec = list(eta_vec)
        self.cc = _corner_complex(res)
        self.tc = _tensor_complex(res)
        self._eta_at = _cochain_evaluator(res, eta_vec, n)
        self._tower_cache = {}
        self._check_cocycle()

    def _check_cocycle(self):
        """Raise a DISTINCT CocycleError (never the diagonal's scope message) if
        δη != 0.  δ^n = res.matrix(n, "coh") (rows C^{n+1}, cols C^n); δη = M·η."""
        res, dom = self.res, self.dom
        M = res.matrix(self.n, "coh")
        for row in M:
            acc = dom.zero()
            for j, mij in enumerate(row):
                if not dom.is_zero(mij):
                    acc = dom.add(acc, dom.mul(mij, dom.coerce(self.eta_vec[j])))
            if not dom.is_zero(acc):
                raise CocycleError(
                    f"homotopy_lifting needs a cocycle; delta eta != 0 at degree "
                    f"{self.n}",
                    hint="pass a representative of an HH^n class (delta eta = 0)")

    def apply(self, k, sigma):
        """ψ_η(σ) for a chain σ ∈ S_k, as a single-PELT over the (o(σ),t(σ))-corner
        of P_{k-n+1} (empty when k-n+1 < 0)."""
        return self.tower(k).get(self.tc._chain_word(sigma), {})

    def image_cochain(self, k):
        """The whole degree-k ψ layer {chain-word: PELT} (== tower(k))."""
        return self.tower(k)

    def tower(self, k):
        """ψ_η on S_k as {chain-word: single-PELT over P_{k-n+1}}, cached and
        recursive.  The base (k < n-1: target P_{<0} = 0; k = n-1: target P_0, RHS
        vanishes) falls out of the general FULL-corner lift-solve via the
        no-equations branch (empty D_corner) with the free-variables-zero lift."""
        if k < 0:
            return {}
        cached = self._tower_cache.get(k)
        if cached is not None:
            return cached
        res, dom, n = self.res, self.dom, self.n
        prev = self.tower(k - 1)
        diag_k = diagonal(res, k)
        sign_psi = dom.one() if (n - 1) % 2 == 0 else dom.neg(dom.one())
        out = {}
        for sigma in res.ss.S(k):
            o, t = sigma.o, sigma.t
            m = k - n + 1                       # target degree of ψ(σ)
            # RHS = (η⊗1 − 1⊗η)Δ_k(σ)  +  (−1)^{n-1} ψ(d_k σ)   in P_{k-n} = P_{m-1}
            rhs_pelt = self._half_collapse(sigma, diag_k)
            for key, val in self._psi_on_diff(k, sigma, prev).items():
                _accum(rhs_pelt, key, dom.mul(sign_psi, val), dom)
            rows = self.cc.p_basis(m - 1, o, t)
            cols = self.cc.p_basis(m, o, t)
            ridx = {tup: i for i, tup in enumerate(rows)}
            rhs = [dom.zero()] * len(rows)
            for tup, val in rhs_pelt.items():
                i = ridx.get(tup)
                if i is None:
                    raise AssertionError(
                        f"psi RHS term {tup} escaped the ({o},{t})-corner row basis "
                        f"at degree {k}, chain {sigma.word} -- a corner-bookkeeping "
                        f"bug, never an approximation")
                rhs[i] = dom.add(rhs[i], val)
            M = self.cc.d_matrix(m, o, t)
            if not M:                            # no equations: choose the zero lift
                x = ([dom.zero()] * len(cols)
                     if all(dom.is_zero(v) for v in rhs) else None)
            else:
                x = solve(M, rhs, dom)
            if x is None:
                raise AssertionError(
                    f"homotopy-lifting solve is inconsistent at degree {k}, chain "
                    f"{sigma.word}: this cannot occur once Delta is built (Volkov "
                    f"Lemma 2, Plan 51 DD1a) -- a bug, never an approximation")
            x = reduce_mod_nullspace(x, M, dom)
            psi = {}
            for coeff, tup in zip(x, cols):
                if not dom.is_zero(coeff):
                    psi[tup] = coeff
            out[self.tc._chain_word(sigma)] = psi
        self._tower_cache[k] = out
        return out

    def _koszul(self, e):
        """(−1)^e as a Domain element."""
        return self.dom.one() if e % 2 == 0 else self.dom.neg(self.dom.one())

    def _half_collapse(self, sigma, diag_k):
        """(η⊗1 − 1⊗η) Δ_k(σ) as a single-PELT over the (o(σ),t(σ))-corner of
        P_{k-n}.  For a double-PELT key (a, τ, mid, ρ, c): if deg τ = n then η eats
        τ (η⊗1) leaving the survivor (b_a·η(τ)·b_mid) ⊗ ρ ⊗ b_c; if deg ρ = n then
        η eats ρ (−1⊗η) leaving b_a ⊗ τ ⊗ (b_mid·η(ρ)·b_c).  The a·w·b interior
        collapse is the cup convention (`cup.native_cup`)."""
        res, ar, dom, tc, n = self.res, self.ar, self.dom, self.tc, self.n
        eta_at = self._eta_at
        dpelt = diag_k.get(tc._chain_word(sigma), {})
        out = {}
        for (ai, tau_word, mi, rho_word, ci), coeff in dpelt.items():
            if dom.is_zero(coeff):
                continue
            tau = tc._chain(tau_word)
            rho = tc._chain(rho_word)
            b_a = ar.A._basis_vec(ai)
            b_mid = ar.A._basis_vec(mi)
            b_c = ar.A._basis_vec(ci)
            if tau.degree == n:                             # η ⊗ 1
                eta_tau = eta_at(tau)
                if not all(dom.is_zero(v) for v in eta_tau):
                    lvec = ar.mul(b_a, ar.mul(eta_tau, b_mid))     # ∈ e_o A e_{o(ρ)}
                    for aj, av in enumerate(lvec):
                        if not dom.is_zero(av):
                            _accum(out, (aj, rho.word, ci),
                                   dom.mul(coeff, av), dom)
            if rho.degree == n:                             # − (1 ⊗ η)
                eta_rho = eta_at(rho)
                if not all(dom.is_zero(v) for v in eta_rho):
                    rvec = ar.mul(b_mid, ar.mul(eta_rho, b_c))     # ∈ e_{t(τ)} A e_t
                    # Koszul sign of the tensor of maps: (1⊗η)(x⊗y) = (−1)^{n·|x|}
                    # x·η(y), |x| = deg of the surviving first factor τ.
                    scoeff = dom.mul(coeff, self._koszul(n * tau.degree))
                    for cj, cv in enumerate(rvec):
                        if not dom.is_zero(cv):
                            _accum(out, (ai, tau.word, cj),
                                   dom.neg(dom.mul(scoeff, cv)), dom)
        return out

    def _psi_on_diff(self, k, sigma, prev):
        """ψ_η(d_k σ) as a single-PELT over the (o(σ),t(σ))-corner of P_{k-n}.
        d_k σ = Σ (coeff, a', tw, c') in P_{k-1}; ψ is a bimodule map, so
        ψ(a'·tw·c') = a'·ψ(tw)·c' (the OUTER action, the sibling of
        TensorComplex._zeta): a' left-multiplies the a-slot, c' right-multiplies
        the c-slot of the cached ψ(tw)."""
        res, ar, dom, tc = self.res, self.ar, self.dom, self.tc
        out = {}
        for (coeff, a_word, tw, c_word) in res.d_terms(k, sigma):
            if dom.is_zero(coeff):
                continue
            a_vec, c_vec = _vecs(res, sigma, a_word, c_word)
            target = _resolve_chain(res, tw)
            psi_tw = prev.get(tc._chain_word(target), {})
            if not psi_tw:
                continue
            for (ai, chw, ci), kappa in psi_tw.items():
                base = dom.mul(coeff, kappa)
                if dom.is_zero(base):
                    continue
                left = ar.mul(a_vec, ar.A._basis_vec(ai))
                right = ar.mul(ar.A._basis_vec(ci), c_vec)
                for aj, av in enumerate(left):
                    if dom.is_zero(av):
                        continue
                    for cj, cv in enumerate(right):
                        if dom.is_zero(cv):
                            continue
                        _accum(out, (aj, chw, cj),
                               dom.mul(base, dom.mul(av, cv)), dom)
        return out


def homotopy_lifting(res, eta_vec, n):
    """The canonical homotopy lifting ψ_η of the CS cocycle η ∈ C^n (coordinates
    over res._basis(n, "coh")) as a `HomotopyLifting` bound to `res`.  Raises
    `CocycleError` if δη != 0; the ψ-solve is consistent at every degree Δ was
    built (Volkov Lemma 2, Plan 51 DD1a)."""
    return HomotopyLifting(res, eta_vec, n)
