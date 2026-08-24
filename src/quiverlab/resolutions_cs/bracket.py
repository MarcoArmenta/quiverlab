"""Plan 51: the native Hochschild Gerstenhaber bracket on the Chouhy-Solotar
resolution, past the bar window, over any exact Domain -- the sibling of
`cup.native_cup` / `cap.native_cap`, assembled from two homotopy liftings
(`homotopy_lifting.py`) with NO bar object and NO brace/circle machinery.

Given CS cochains f ∈ C^p, g ∈ C^q (coordinate vectors over
`res._basis(·, "coh")`), the bracket cochain is read off on generators
σ ∈ S_{p+q-1} with no further solve (Negron-Witherspoon / Oke Thm 3.5):

    [f, g](σ)  =  f( ψ_g(σ) )  −  (−1)^{(p-1)(q-1)} · g( ψ_f(σ) )

where ψ_g: P_{p+q-1} → P_p and ψ_f: P_{p+q-1} → P_q are the homotopy liftings, and
f/g are evaluated on the resulting single-PELTs by the bimodule-map collapse
Σ κ · b_a · f(chain) · b_c (the a·w·b cup convention).  The result is a
(p+q-1)-cocycle whose HH class is choice-independent (the NW/Volkov theorem); the
canonical ψ makes the returned representative byte-reproducible (the native_cup
precedent -- only the class is invariant).

DD4 -- the sign is ARBITRATED, not assumed.  The `(−1)^{(p-1)(q-1)}` factor is the
Oke/NW convention (matching the homotopy lifting's `(−1)^{n-1}`); it is DECIDED by
the in-window anchor `native_bracket ≡ Comparison.bracket_of_cs_classes` over GF(p)
(the QuantumCI (2,2) odd-exponent case is the load-bearing discriminator).  The
anchor holds with this sign on the tested fixtures, so it stands; any future
deviation is recorded here and folded back into the research doc.

Degree-0 insertion is out of scope (p, q >= 1), same as Plan 35.  Domain-generic;
exact; no floats; no engine imports.
"""
from quiverlab.resolutions_cs.cup import _cochain_evaluator
from quiverlab.resolutions_cs.homotopy_lifting import homotopy_lifting
from quiverlab.resolutions_cs.pelt import _resolve_chain


def _eval_cochain_on_pelt(res, eval_fn, pelt):
    """Evaluate a bimodule-map cochain (via `eval_fn` = `cup._cochain_evaluator`)
    on a single-PELT {(a_idx, chain_word, c_idx): κ}: Σ κ · b_a · f(chain) · b_c,
    the A-vector value of the cochain on the PELT (the a·w·b collapse)."""
    ar, dom = res.ar, res.dom
    acc = [dom.zero()] * res.A.dim
    for (ai, chw, ci), kappa in pelt.items():
        if dom.is_zero(kappa):
            continue
        fval = eval_fn(_resolve_chain(res, chw))
        if all(dom.is_zero(v) for v in fval):
            continue
        prod = ar.mul(ar.A._basis_vec(ai), ar.mul(fval, ar.A._basis_vec(ci)))
        for kk, pv in enumerate(prod):
            if not dom.is_zero(pv):
                acc[kk] = dom.add(acc[kk], dom.mul(kappa, pv))
    return acc


def native_bracket(res, f_vec, p, g_vec, q, psi_f=None, psi_g=None):
    """The native CS Gerstenhaber bracket [f, g] of cochains f ∈ C^p, g ∈ C^q, as a
    coordinate vector over `res._basis(p+q-1, "coh")`:
    ``[f,g](σ) = f(ψ_g(σ)) − (−1)^{(p-1)(q-1)} g(ψ_f(σ))`` for σ ∈ S_{p+q-1}.

    Builds the two homotopy liftings ψ_f (of f, degree p) and ψ_g (of g, degree q)
    unless PREBUILT ones are supplied via `psi_f`/`psi_g` — a table builder reuses
    one lifting per class across a whole row/column instead of rebuilding it per pair
    (DD1's "two ψ towers per class-pair, cached per class"; the canonical solve makes
    the cache a pure byte-identical win).  No bar object -- any degree, any Domain.
    Requires p, q >= 1 (ValueError otherwise -- degree-0 insertion out of scope)."""
    if p < 1 or q < 1:
        raise ValueError(
            f"the Gerstenhaber bracket needs p, q >= 1 (degree-0 insertion is out "
            f"of scope); got p={p}, q={q}")
    ar, dom = res.ar, res.dom
    n = p + q - 1
    if psi_f is None:
        psi_f = homotopy_lifting(res, f_vec, p)     # ψ_f: P_n → P_{n-p+1} = P_q
    if psi_g is None:
        psi_g = homotopy_lifting(res, g_vec, q)     # ψ_g: P_n → P_{n-q+1} = P_p
    f_at = _cochain_evaluator(res, f_vec, p)
    g_at = _cochain_evaluator(res, g_vec, q)
    sign = dom.one() if ((p - 1) * (q - 1)) % 2 == 0 else dom.neg(dom.one())

    out_basis = res._basis(n, "coh")
    out_pos = {(ch.word, j): i for i, (ch, j) in enumerate(out_basis)}
    out = [dom.zero()] * len(out_basis)
    for sigma in res.ss.S(n):
        val_fg = _eval_cochain_on_pelt(res, f_at, psi_g.apply(n, sigma))
        val_gf = _eval_cochain_on_pelt(res, g_at, psi_f.apply(n, sigma))
        for j in ar.corner(sigma.o, sigma.t, "coh"):
            v = dom.add(val_fg[j], dom.neg(dom.mul(sign, val_gf[j])))
            if not dom.is_zero(v):
                out[out_pos[(sigma.word, j)]] = v
    return out
