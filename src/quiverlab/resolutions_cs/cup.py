"""Plan 20 Task 3: the native Hochschild cup product on the Chouhy–Solotar resolution.

Given the lifted diagonal Δ: P → P ⊗_A P (`diagonal.py`, Task 2) as a chain map
against the Koszul-signed tensor differential, the cup of two CS cochains
f ∈ C^p, g ∈ C^q (coordinate vectors over `ChouhySolotarResolution._basis(·,"coh")`)
is read off degreewise, with NO bar object ever built:

    (f∪g)(σ) = Σ_Δ  coeff · ar.mul(b_a, ar.mul(f(τ), ar.mul(b_mid, ar.mul(g(ρ), b_c))))

for σ ∈ S_{p+q}, where Δ_{p+q}(σ) is the double-PELT dict with keys
(a_idx, τ.word, mid_idx, ρ.word, c_idx); ONLY the keys whose τ has degree p and ρ
has degree q — the (p,q)-component — contribute.  b_a = e_{o(σ)}Ae_{o(τ)},
b_mid = e_{t(τ)}Ae_{o(ρ)}, b_c = e_{t(ρ)}Ae_{t(σ)} are A-basis vectors (from the
double-PELT indices); f(τ) = Σ_j f[(τ,j)]·basis_vec(j) ∈ e_{o(τ)}Ae_{t(τ)} and
g(ρ) ∈ e_{o(ρ)}Ae_{t(ρ)} are the cochain values.  The resulting A-vector lands in
e_{o(σ)}Ae_{t(σ)}, and its coordinates over corner(o(σ), t(σ), "coh") pair with σ
to give the output cochain over `_basis(p+q, "coh")`.

The evaluation is SIGN-FREE — the a·w·b cohomology collapse (`resolution.py:201`)
generalized to two interior values, matching the engine's bar cup convention.  All
Koszul sign lives in the tensor differential inside Δ; the Leibniz identity
δ(f∪g) = (δf)∪g + (−1)^p f∪(δg) is the arbiter (see test_native_cup.py).

Domain-generic (arithmetic through `res.dom` / AArith); exact; no floats; no engine
imports.  Determinism follows from Δ's byte-reproducibility and the deterministic
`_basis`/`corner` orderings.
"""
from quiverlab.resolutions_cs.diagonal import diagonal
from quiverlab.resolutions_cs.pelt import _resolve_chain


def _cochain_evaluator(res, vec, deg):
    """Return `evaluate(chain)` -> the A-vector value f(chain) ∈ e_o A e_t of the
    cochain `vec` (coordinates over `_basis(deg,"coh")`) on a degree-`deg` chain.

    f(chain) = Σ_{j ∈ corner(o,t,"coh")} vec[(chain,j)] · basis_vec(j); since
    basis_vec(j) is the j-th unit A-vector, the value simply carries vec's coordinate
    at each corner index j.  Keyed by (chain.word, j): for degree ≥ 1 the word is
    unique, and for the degree-0 vertex chains (word = ()) the corner index j alone
    fixes the vertex (e_v-supported), so there is no collision."""
    ar, dom = res.ar, res.dom
    basis = res._basis(deg, "coh")
    coeff_of = {(ch.word, j): vec[i] for i, (ch, j) in enumerate(basis)}

    def evaluate(chain):
        val = [dom.zero()] * res.A.dim
        for j in ar.corner(chain.o, chain.t, "coh"):
            c = coeff_of.get((chain.word, j))
            if c is not None and not dom.is_zero(c):
                val[j] = c
        return val

    return evaluate


def native_cup(res, f_vec, p, g_vec, q):
    """The native CS cup f∪g of two cochains f ∈ C^p, g ∈ C^q, returned as a
    coordinate vector over `res._basis(p+q, "coh")`.

    Iterates the lifted diagonal Δ_{p+q}(σ) for every σ ∈ S_{p+q}, keeps its
    (p,q)-component (τ.degree = p, ρ.degree = q), evaluates f on τ and g on ρ,
    collapses per the sign-free a·w·b formula, and reads the resulting A-vector off
    against corner(o(σ), t(σ), "coh").  Δ is cached on `res` (module `diagonal`), so
    repeated cups reuse it."""
    ar, dom = res.ar, res.dom
    n = p + q
    diag = diagonal(res, n)                       # {chain-word: double-PELT}, cached
    f_at = _cochain_evaluator(res, f_vec, p)
    g_at = _cochain_evaluator(res, g_vec, q)

    out_basis = res._basis(n, "coh")
    out_pos = {(ch.word, j): i for i, (ch, j) in enumerate(out_basis)}
    out = [dom.zero()] * len(out_basis)

    for sigma_word, dpelt in diag.items():
        sigma = _resolve_chain(res, sigma_word)
        acc = [dom.zero()] * res.A.dim            # (f∪g)(σ) ∈ e_{o(σ)}Ae_{t(σ)}
        for (ai, tau_word, mi, rho_word, ci), coeff in dpelt.items():
            if dom.is_zero(coeff):
                continue
            tau = _resolve_chain(res, tau_word)
            if tau.degree != p:
                continue
            rho = _resolve_chain(res, rho_word)
            if rho.degree != q:
                continue
            fval = f_at(tau)
            if all(dom.is_zero(v) for v in fval):
                continue
            gval = g_at(rho)
            if all(dom.is_zero(v) for v in gval):
                continue
            b_a = ar.A._basis_vec(ai)
            b_mid = ar.A._basis_vec(mi)
            b_c = ar.A._basis_vec(ci)
            # a·w·b collapse with two interior values (sign-free, bar convention)
            prod = ar.mul(b_a, ar.mul(fval, ar.mul(b_mid, ar.mul(gval, b_c))))
            for k, pv in enumerate(prod):
                if not dom.is_zero(pv):
                    acc[k] = dom.add(acc[k], dom.mul(coeff, pv))
        for j in ar.corner(sigma.o, sigma.t, "coh"):
            v = acc[j]
            if not dom.is_zero(v):
                pos = out_pos[(sigma.word, j)]
                out[pos] = dom.add(out[pos], v)
    return out
