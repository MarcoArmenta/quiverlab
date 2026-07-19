"""Uncollapsed P-element arithmetic for the correction linear-solve and the gates.
A P-ELEMENT (chain-degree n) is a dict (a_idx, chain_word, c_idx) -> coeff, where a_idx,c_idx
are A-basis indices (b,b' ∈ B), i.e. Σ coeff·(b ⊗ chain ⊗ b'). Built from a term list
(coeff, a_word, chain_word, c_word) by reducing a_word,c_word to normal-form A-vectors (π=β)
and expanding over B."""


def _resolve_chain(res, chain_word):
    if isinstance(chain_word, tuple) and len(chain_word) == 2 and chain_word[0] == "__v__":
        for c in res.ss.S(0):
            if c.o == chain_word[1]:
                return c
        raise KeyError(("vertex chain", chain_word))
    for n in range(res.ss.max_degree + 1):
        c = res._chain(n, chain_word)
        if c is not None:
            return c
    raise KeyError(("chain not found", chain_word))


def _accum(out, key, val, dom):
    cur = out.get(key)
    tot = val if cur is None else dom.add(cur, val)
    if dom.is_zero(tot):
        out.pop(key, None)
    else:
        out[key] = tot


def _vecs(res, ch, a_word, c_word):
    a = res.ar.vertex_vec(ch.o) if len(a_word) == 0 else res.ar.path_vec(a_word)
    c = res.ar.vertex_vec(ch.t) if len(c_word) == 0 else res.ar.path_vec(c_word)
    return a, c


def terms_to_pelt(res, term_list):
    dom, out = res.dom, {}
    for (coeff, a_word, chain_word, c_word) in term_list:
        ch = _resolve_chain(res, chain_word)
        a_vec, c_vec = _vecs(res, ch, a_word, c_word)
        for ai, av in enumerate(a_vec):
            if dom.is_zero(av):
                continue
            for ci, cv in enumerate(c_vec):
                if dom.is_zero(cv):
                    continue
                _accum(out, (ai, ch.word, ci), dom.mul(coeff, dom.mul(av, cv)), dom)
    return out


def apply_lower(res, n, pelt):
    """d_{n-1} applied to a P-element at chain-degree n-1, returning chain-degree n-2. For a
    key (b, σ, b') with d_{n-1}(1⊗σ⊗1) = Σ (c, a, τ, c'): image += c·(b·a ⊗ τ ⊗ c'·b')."""
    ar, dom, out = res.ar, res.dom, {}
    for (bi, chain_word, ci), coeff in pelt.items():
        ch = _resolve_chain(res, chain_word)
        bvec, cvec = ar.A._basis_vec(bi), ar.A._basis_vec(ci)
        for (c, a_word, tw, c_word) in res.d_terms(n - 1, ch):
            a_vec, c_vec = _vecs(res, ch, a_word, c_word)
            left, right = ar.mul(bvec, a_vec), ar.mul(c_vec, cvec)
            tch, base = _resolve_chain(res, tw), dom.mul(coeff, c)
            for ai, av in enumerate(left):
                if dom.is_zero(av):
                    continue
                for cj, cv in enumerate(right):
                    if dom.is_zero(cv):
                        continue
                    _accum(out, (ai, tch.word, cj), dom.mul(base, dom.mul(av, cv)), dom)
    return out
