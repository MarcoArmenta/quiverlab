"""Plan 21: the native Hochschild cap product on the Chouhy–Solotar resolution.

The cap reuses the SAME lifted diagonal Δ: P → P ⊗_A P (`diagonal.py`, Plan 20) the
native cup is read off — no new construction.  Where the cup is the COHOMOLOGY-side
`a·w·b` collapse of Δ (two interior cochain values), the cap is the HOMOLOGY-side
`b·w·a` collapse: it caps a cochain f ∈ C^p with a chain z ∈ C_n and lands in C_{n-p}.

The covariance flip is the heart of it, exactly as in the minimal engine's Plan-16
note (`engine/resolutions_minimal.py`): cohomology collapses `a·w·b` while HOMOLOGY
collapses `b·w·a` on SWAPPED corner tags.  For a double-PELT key
(a_idx, τ.word, mid_idx, ρ.word, c_idx) of the (p, n-p)-component of Δ_n(σ)
(τ.degree = p, ρ.degree = n-p), with

    b_a = e_{o(σ)}Ae_{o(τ)},  b_mid = e_{t(τ)}Ae_{o(ρ)},  b_c = e_{t(ρ)}Ae_{t(σ)}

the cochain f eats the degree-p FIRST factor τ (f(τ) ∈ e_{o(τ)}Ae_{t(τ)}), the
degree-(n-p) factor ρ survives, and — collapsing `x ⊗_{A^e} P_{n-p}` to C_{n-p} the
same op-twisted way `resolution.py`'s homology differential does (`x ⊗ (L⊗[·]⊗R) →
R·x·L`) — the value on ρ is

    (f ∩ z)(ρ)  =  Σ_Δ  coeff · b_c · x · b_a · f(τ) · b_mid       (SIGN-FREE)

where x = z(σ) ∈ e_{t(σ)}Ae_{o(σ)} is the chain's value (its `hom`-corner
coordinates).  The result lands in e_{t(ρ)}Ae_{o(ρ)} = ρ's `hom` corner and is read
against `corner(o(ρ), t(ρ), "hom")`.  No bar object is ever built, so the cap works at
ANY degree (past the bar-comparison window).

Why sign-free.  All Koszul sign already lives inside Δ's tensor differential (the ONE
`(−1)^p` on `1⊗d_q`); the collapse adds none — the bar cap is likewise sign-free
(`engine/tt_calculus.py:29-33`).  The arbiters: the in-window native ≡ transported
`Comparison.cap_of_cs_classes` (the anchor, non-commutative QCI distinguishes the
op-twist from `a·w·b`), the unit cap `1 ∩ z = z`, the module identity
`(z∩f)∩g ~ z∩(f∪g)`, and the exact cap-Leibniz `b(f∩z) = (−1)^{p+1}(δf∩z) +
(−1)^p(f∩bz)`.

Domain-generic (arithmetic through `res.dom` / AArith); exact; no floats; no engine
imports.  Determinism follows from Δ's byte-reproducibility and the deterministic
`_basis`/`corner` orderings.
"""
from quiverlab.resolutions_cs.cup import _cochain_evaluator
from quiverlab.resolutions_cs.diagonal import diagonal


def native_cap(res, f_vec, p, z_vec, n):
    """The native CS cap f ∩ z of a cochain f ∈ C^p (coordinates over
    `res._basis(p, "coh")`) with a chain z ∈ C_n (coordinates over
    `res._basis(n, "hom")`), returned as a coordinate vector over
    `res._basis(n-p, "hom")`.

    Iterates the lifted diagonal Δ_n(σ) for every σ ∈ S_n, keeps its
    (p, n-p)-component (τ.degree = p, ρ.degree = n-p), evaluates f on τ, and collapses
    per the sign-free op-twisted `b_c · x · b_a · f(τ) · b_mid` formula (x = z's value
    on σ), reading the resulting A-vector off against corner(o(ρ), t(ρ), "hom").

    Raises ``ValueError`` for p > n (a cap into negative degree — bar convention,
    `engine.tt_calculus.cap_cochain`).  Δ is cached on `res`, so a cap reuses it."""
    if p > n:
        raise ValueError(
            f"cap product needs p <= n (C^{p} ∩ C_{n} lands in C_{n - p} < 0)")
    ar, dom = res.ar, res.dom
    diag = diagonal(res, n)                       # {chain-word: double-PELT}, cached
    tc = res._tensor_complex                      # set by diagonal(); its chain cache
    f_at = _cochain_evaluator(res, f_vec, p)

    # z as {σ.word: A-vector value on the hom corner of σ}
    z_val = {}
    for coord, (ch, j) in zip(z_vec, res._basis(n, "hom")):
        if dom.is_zero(coord):
            continue
        vec = z_val.get(ch.word)
        if vec is None:
            vec = [dom.zero()] * res.A.dim
            z_val[ch.word] = vec
        vec[j] = dom.add(vec[j], coord)

    out_basis = res._basis(n - p, "hom")
    out_pos = {(ch.word, k): i for i, (ch, k) in enumerate(out_basis)}
    out = [dom.zero()] * len(out_basis)

    for sigma_word, dpelt in diag.items():
        sigma = tc._chain(sigma_word)
        x = z_val.get(sigma.word)
        if x is None:                             # z has no support on this σ
            continue
        for (ai, tau_word, mi, rho_word, ci), coeff in dpelt.items():
            if dom.is_zero(coeff):
                continue
            tau = tc._chain(tau_word)
            if tau.degree != p:
                continue
            rho = tc._chain(rho_word)
            if rho.degree != n - p:
                continue
            fval = f_at(tau)
            if all(dom.is_zero(v) for v in fval):
                continue
            b_a = ar.A._basis_vec(ai)
            b_mid = ar.A._basis_vec(mi)
            b_c = ar.A._basis_vec(ci)
            # op-twisted b·w·a homology collapse (R·x·L, L = b_a·f(τ)·b_mid, R = b_c);
            # sign-free — all Koszul sign is inside Δ (matches resolution.py:187).
            prod = ar.mul(b_c, ar.mul(x, ar.mul(b_a, ar.mul(fval, b_mid))))
            for k in ar.corner(rho.o, rho.t, "hom"):
                v = prod[k]
                if not dom.is_zero(v):
                    pos = out_pos[(rho.word, k)]
                    out[pos] = dom.add(out[pos], dom.mul(coeff, v))
    return out
