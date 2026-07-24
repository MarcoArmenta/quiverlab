"""Betti numbers of the minimal A^e resolution over ANY exact Domain
(Plan 19), via the E-relative (Cibils) complex.

For A = E (+) r split basic (path-type basis, E separable), the relative
bar resolution A (x)_E r^{(x)_E n} (x)_E A is an A^e-projective resolution
of A; applying (x)_{A^e} (E (x) E) kills the two outer faces (their factor
acts through A/r on E), leaving the SMALL complex

    T_n = r^{(x)_E n},   d(r_1 (x) ... (x) r_n)
        = sum_{i=1}^{n-1} (-1)^i r_1 (x)..(x) r_i r_{i+1} (x)..(x) r_n

whose H_n = Tor_n^{A^e}(A, E (x) E). Minimality of the true minimal
resolution kills ITS induced differential on E (x) E coefficients, so
dim H_n = the number of A^e-generators of P_n — the engine's rks[n]
(engine/resolutions_minimal.py) — over every field. GF(p) parity is gated
in tests/invariants/test_betti_generic.py."""
from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.fields.linalg import rank
from quiverlab.invariants.pathbasis import path_type_basis

_GUARD_HINT = ("the relative-Tor chain count grows with the algebra's own "
               "complexity; raise max_cells only if you know what you are doing")


def _relative_complex(A, top, max_cells):
    """(chains, mats): chains[n] = composable radical index tuples (n >= 1;
    chains[0] is the list of idempotent indices standing for E), mats[n] =
    the matrix of d_n : T_n -> T_{n-1} for n >= 2 (d_1 = 0)."""
    dom = A.domain
    idem, rad, src, tgt = path_type_basis(A, "complexity")
    idemset = set(idem)
    chains = {0: [(v,) for v in idem], 1: [(i,) for i in rad]}
    for n in range(2, top + 2):
        chains[n] = [ch + (j,) for ch in chains[n - 1]
                     for j in rad if src[j] == tgt[ch[-1]]]
    mats = {}
    for n in range(2, top + 2):
        rows, cols = len(chains[n - 1]), len(chains[n])
        if rows * cols > max_cells:
            raise DepthLimitError(
                f"relative Tor d_{n}: {rows} x {cols} entries "
                f"(> max_cells = {max_cells})", hint=_GUARD_HINT)
        ridx = {ch: k for k, ch in enumerate(chains[n - 1])}
        M = [[dom.zero()] * cols for _ in range(rows)]
        for ci, ch in enumerate(chains[n]):
            for i in range(1, n):
                prod = A.T[ch[i - 1]][ch[i]]
                for t, cf in enumerate(prod):
                    if dom.is_zero(cf):
                        continue
                    if t in idemset:
                        raise QuiverlabError(
                            "product of radical basis elements has an "
                            "idempotent coordinate — basis is not path-type")
                    r = ridx[ch[:i - 1] + (t,) + ch[i + 1:]]
                    val = cf if i % 2 == 0 else dom.neg(cf)
                    M[r][ci] = dom.add(M[r][ci], val)
        mats[n] = M
    return chains, mats


def relative_betti_numbers(A, top, max_cells=4_000_000):
    """[rks_0, ..., rks_top]: A^e-generator counts of the minimal resolution
    of A as a bimodule, over A.domain (any exact field)."""
    dom = A.domain
    chains, mats = _relative_complex(A, top, max_cells)
    ranks = {0: 0, 1: 0}
    for n in range(2, top + 2):
        M = mats[n]
        ranks[n] = rank(M, dom) if M and M[0] else 0
    return [len(chains[n]) - ranks[n] - ranks[n + 1] for n in range(top + 1)]
