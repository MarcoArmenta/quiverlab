# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""General minimal projective A^e-bimodule resolution via iterated syzygies.

This is the *resolution-independent* engine: it computes
    dim HH_n(A; F_p) = dim Tor_n^{A^e}(A, A) over F_p
for ANY finite-dimensional algebra A given by its structure constants, by building
the **minimal** free A^e-resolution

    ... -> (A^e)^{r_n} --d_n--> (A^e)^{r_{n-1}} -> ... -> A^e --mu--> A -> 0

one degree at a time (kernel -> minimal generators modulo rad(A^e) -> next term),
then applying A (x)_{A^e} (-) to get the Hochschild chain complex.  Since Tor does
not depend on the chosen projective resolution, this computes exactly the same HH_*
that the Chouhy-Solotar resolution targets -- but it is valid for *every* algebra,
not only the two closed-form CS families, and its terms are the SMALLEST possible
(the minimal number r_n of generators), so it reaches far past the normalized bar
complex's `dim C_n = m*(m-1)^n` blow-up.

WHY THIS IS THE RIGHT TOOL FOR THE OPEN ZONE.  The genuinely-open cases of Han's
conjecture (non-monomial, non-complete-intersection, non-Koszul algebras) have no
closed-form resolution a human can write down; the general CS differential is an
intricate recursive comparison map.  Computing the minimal resolution by syzygies
sidesteps that: it is mechanical, exact, and -- crucially -- *validatable* against
the normalized bar complex on small cases and against the CS closed forms.

CHARACTERISTIC.  A minimal resolution is characteristic-specific (the rank r_n and
the generators depend on F_p), so the whole resolution is rebuilt per prime.  The
large prime PRIME = 32003 is the faithful char-0 proxy; small primes expose torsion.

MEMORY.  The local host is small, so `minimal_homology_dims` takes a `max_term_dim`
budget: if a term (A^e)^{r_n} would exceed `m^2 * r_n > max_term_dim` k-dimensions,
the build stops and the returned dims are flagged truncated (the degrees already
computed are exact).

Validated against `hh_engine.hochschild_homology_dims` (the bar oracle) and against
the closed-form `ChouhySolotarResolution`; see tests/test_minimal_resolution.py.
"""

import numpy as np

from quiverlab.engine import _kernels
from quiverlab.engine.hh_engine import rank_mod_p


# ---------------------------------------------------------------------------
# A^e arithmetic (left multiplication), mod a given prime p
# ---------------------------------------------------------------------------
class AeEngine:
    """A^e = A (x) A^op arithmetic over F_p, precomputing the left-multiplication
    matrices of the A^e basis elements e_a (x) e_b on A^e (dimension m^2)."""

    def __init__(self, A, p):
        self.A = A
        self.p = p
        self.m = A.m
        self.m2 = self.m * self.m
        m = self.m
        self.T = A.T % p
        T = self.T
        # Lb[a, b] : the (m2 x m2) matrix of left-mult by (e_a (x) e_b) on A^e.
        # (e_a (x) e_b) . (e_p (x) e_q) = (e_a e_p) (x) (e_q e_b).
        self.Lb = np.zeros((m, m, self.m2, self.m2), dtype=np.int64)
        for a in range(m):
            for b in range(m):
                for p_ in range(m):
                    for q_ in range(m):
                        col = p_ * m + q_
                        aep = T[a, p_, :]      # e_a * e_p
                        eqb = T[q_, b, :]      # e_q * e_b
                        for ia in np.nonzero(aep % p)[0]:
                            for ib in np.nonzero(eqb % p)[0]:
                                self.Lb[a, b, ia * m + ib, col] = (
                                    self.Lb[a, b, ia * m + ib, col]
                                    + aep[ia] * eqb[ib]) % p

    def apply_block(self, a, b, vec, r):
        """Left-multiply by (e_a (x) e_b) on a vector in (A^e)^r, block by block."""
        m2 = self.m2
        out = np.zeros(m2 * r, dtype=np.int64)
        Lab = self.Lb[a, b]
        for blk in range(r):
            out[blk * m2:(blk + 1) * m2] = (Lab @ vec[blk * m2:(blk + 1) * m2]) % self.p
        return out


# ---------------------------------------------------------------------------
# Exact F_p linear algebra (nullspace, span tests)
# ---------------------------------------------------------------------------
def nullspace_mod_p(M, p):
    """Basis of the right nullspace { x : M x = 0 } over F_p, as a list of vectors."""
    if _kernels.USE_KERNELS:
        K = _kernels.nullspace_kernel(np.ascontiguousarray(M.astype(np.int64)), p)
        return [K[i] for i in range(K.shape[0])]
    return _nullspace_mod_p_py(M, p)


def _nullspace_mod_p_py(M, p):
    """Basis of the right nullspace { x : M x = 0 } over F_p, as a list of vectors."""
    M = M % p
    rows, cols = M.shape
    A = M.copy().astype(np.int64)
    where = {}
    pivots = []
    r = 0
    for c in range(cols):
        piv = -1
        for i in range(r, rows):
            if A[i, c] % p:
                piv = i
                break
        if piv < 0:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(rows):
            if i != r and A[i, c] % p:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        where[c] = r
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in where]
    basis = []
    for f in free:
        v = np.zeros(cols, dtype=np.int64)
        v[f] = 1
        for c in pivots:
            v[c] = (-A[where[c], f]) % p
        basis.append(v % p)
    return basis


# ---------------------------------------------------------------------------
# rad(A) and rad(A^e) over F_p (for picking minimal generators)
# ---------------------------------------------------------------------------
def radical_basis(A, p):
    """A basis of rad(A) for a LOCAL algebra, characteristic-independently.

    A finite-dimensional local algebra splits as A = k.1 (+) rad(A), and in the
    unit-adapted basis (where 1_A is the basis vector f_t) the radical is exactly the
    span of the non-unit basis vectors {f_i : i != t}.  This is correct in every
    characteristic, unlike the trace form (which degenerates when p divides the
    relevant dimensions, e.g. k[x]/(x^2) at p=2).

    All algebras in this engine's intended use -- the open-zone local non-monomial
    zoo, plus the validation algebras k[x]/(x^a) and the quantum CIs -- are local with
    unit = f_t, so this is exact.  (A genuinely multi-vertex algebra would need its
    idempotents; that is out of scope here and would be flagged by the bar-oracle
    cross-check failing.)
    """
    m = A.m
    out = []
    for i in range(m):
        if i != A.t:
            v = np.zeros(m, dtype=np.int64)
            v[i] = 1
            out.append(v)
    return out


def _rad_ae_columns(eng, radA):
    """Spanning set of rad(A^e) = rad A (x) A + A (x) rad A^op as vectors in A^e."""
    m = eng.m
    p = eng.p
    cols = []
    for r in radA:                                  # rad A (x) A
        for b in range(m):
            v = np.zeros(eng.m2, dtype=np.int64)
            for a in np.nonzero(r % p)[0]:
                v[a * m + b] = (v[a * m + b] + r[a]) % p
            cols.append(v)
    for r in radA:                                  # A (x) rad A
        for a in range(m):
            v = np.zeros(eng.m2, dtype=np.int64)
            for b in np.nonzero(r % p)[0]:
                v[a * m + b] = (v[a * m + b] + r[b]) % p
            cols.append(v)
    return cols


def _independent_modulo(span_cols, candidates, p):
    """Candidates independent modulo span(span_cols). Returns the chosen vectors."""
    if _kernels.USE_KERNELS:
        L = len(candidates[0])
        if isinstance(span_cols, np.ndarray):
            # radK arrives straight from _build_radK as a contiguous (nspan, L) int64 array,
            # already reduced mod p -> pass it through with NO copy (the kernel re-reduces per
            # entry anyway).  This is what removes the duplicate radK from the peak.
            span_arr = span_cols
        else:
            nspan = len(span_cols)
            span_arr = (np.array([v % p for v in span_cols], dtype=np.int64).reshape(nspan, L)
                        if nspan else np.zeros((0, L), dtype=np.int64))
        cand_arr = np.array([v % p for v in candidates], dtype=np.int64).reshape(
            len(candidates), L)
        idx = _kernels.independent_modulo_kernel(span_arr, cand_arr, p)
        return [candidates[int(i)] for i in idx]
    return _independent_modulo_py(span_cols, candidates, p)


def _independent_modulo_py(span_cols, candidates, p):
    """Greedily pick the candidates that are independent modulo the span of
    `span_cols`, using one incremental RREF.  Returns the chosen candidate vectors.

    This replaces the O(#cand) repeated full-rank tests with a single elimination."""
    if span_cols:
        basis = [v.copy() % p for v in span_cols]
    else:
        basis = []
    # reduce `basis` to echelon (list of (pivot_col, row))
    pivot_of = {}
    ech = []
    for v in basis:
        w = v.copy() % p
        for (pc, row) in ech:
            if w[pc]:
                w = (w - w[pc] * row) % p
        nz = np.nonzero(w)[0]
        if nz.size:
            pc = int(nz[0])
            inv = pow(int(w[pc]), p - 2, p)
            w = (w * inv) % p
            ech.append((pc, w))
            pivot_of[pc] = w
    chosen = []
    for cand in candidates:
        w = cand.copy() % p
        for (pc, row) in ech:
            if w[pc]:
                w = (w - w[pc] * row) % p
        nz = np.nonzero(w)[0]
        if nz.size:
            pc = int(nz[0])
            inv = pow(int(w[pc]), p - 2, p)
            w = (w * inv) % p
            ech.append((pc, w))
            chosen.append(cand)
    return chosen


# ---------------------------------------------------------------------------
# rad(A^e) . ker syzygy-image build (extracted hot loop)
# ---------------------------------------------------------------------------
def _rad_ab_pairs(radAe, m):
    """Each rad(A^e) column has exactly one nonzero, at flat index a*m+b; return the
    (a, b) pairs as an (nrad, 2) int64 array (so the kernel can apply_block directly)."""
    pairs = np.empty((len(radAe), 2), dtype=np.int64)
    for ri, v in enumerate(radAe):
        idx = int(np.nonzero(v)[0][0])
        pairs[ri, 0] = idx // m
        pairs[ri, 1] = idx % m
    return pairs


def _build_radK_py(eng, ker, cur_r, p):
    """rad(A^e).ker as a list of vectors (pure-Python reference)."""
    m = eng.m
    m2 = eng.m2
    radA = radical_basis(eng.A, p)
    radAe = _rad_ae_columns(eng, radA)
    radK = []
    for k in ker:
        for rv in radAe:
            E = rv.reshape(m, m)
            img = np.zeros(m2 * cur_r, dtype=np.int64)
            for a in range(m):
                for b in range(m):
                    if E[a, b] % p:
                        img = (img + E[a, b] * eng.apply_block(a, b, k, cur_r)) % p
            radK.append(img)
    return radK


def _build_radK(eng, rad_ab, ker, cur_r, p):
    if _kernels.USE_KERNELS:
        ker_arr = np.array(ker, dtype=np.int64).reshape(len(ker), eng.m2 * cur_r)
        # Return the contiguous (nnull*nrad, L) array directly -- it is already reduced mod p.
        # The old `[out[i] for ...]` list was then rebuilt by _independent_modulo into a full
        # duplicate `span_arr`, doubling the radK footprint (the dominant ~3x peak).  Handing
        # the array straight through lets _independent_modulo consume it with NO copy.
        return _kernels.radK_kernel(eng.Lb, rad_ab, ker_arr, cur_r, eng.m2, p)
    return _build_radK_py(eng, ker, cur_r, p)


# ---------------------------------------------------------------------------
# d_n differential build (extracted hot loop)
# ---------------------------------------------------------------------------
def _build_Dn_py(eng, gens, cur_r, p):
    """Matrix of d_n : (A^e)^{r_n} -> (A^e)^{cur_r} as an int64 array of shape
    (m2*cur_r, m2*r_n) (pure-Python reference). Column j*m2 + (a*m+b) is
    apply_block(a, b, gens[j], cur_r)."""
    m = eng.m
    m2 = eng.m2
    r_n = len(gens)
    Dn = np.zeros((m2 * cur_r, m2 * r_n), dtype=np.int64)
    for j, g in enumerate(gens):
        for a in range(m):
            for b in range(m):
                Dn[:, j * m2 + (a * m + b)] = eng.apply_block(a, b, g, cur_r)
    return Dn


def _build_Dn(eng, gens, cur_r, p):
    if _kernels.USE_KERNELS:
        gens_arr = np.array(gens, dtype=np.int64).reshape(len(gens), eng.m2 * cur_r)
        return _kernels.Dn_kernel(eng.Lb, gens_arr, cur_r, eng.m, eng.m2, p)
    return _build_Dn_py(eng, gens, cur_r, p)


# ---------------------------------------------------------------------------
# The minimal resolution
# ---------------------------------------------------------------------------
def _init_resolution(A, p):
    """Set up the minimal-resolution state: engine, rad(A^e) spanning pairs, d_0
    augmentation.  Returns a mutable state dict advanced by _advance_resolution."""
    eng = AeEngine(A, p)
    m = A.m
    m2 = eng.m2
    T = eng.T
    radA = radical_basis(A, p)
    radAe = _rad_ae_columns(eng, radA)
    rad_ab_pairs = _rad_ab_pairs(radAe, m)
    aug = np.zeros((m, m2), dtype=np.int64)
    for a in range(m):
        for b in range(m):
            prod = T[a, b, :]
            for c in np.nonzero(prod % p)[0]:
                aug[c, a * m + b] = (aug[c, a * m + b] + prod[c]) % p
    return {"eng": eng, "rad_ab_pairs": rad_ab_pairs, "cur": aug, "cur_r": 1,
            "rks": {0: 1}, "cols": {0: None}, "n": 0}


def _advance_resolution(state, p, max_term_dim, max_transient_bytes):
    """Compute the next degree (n = state['n']+1) of the minimal resolution in place.
    Returns {'status', 'radK_bytes', 'r_n'} -- see the plan's interface block for the
    exact state mutations per status."""
    eng = state["eng"]
    rad_ab_pairs = state["rad_ab_pairs"]
    m2 = eng.m2
    n = state["n"] + 1
    cur = state["cur"]
    cur_r = state["cur_r"]
    ker = nullspace_mod_p(cur, p)
    if not ker:
        state["rks"][n] = 0
        state["cols"][n] = []
        state["n"] = n
        return {"status": "terminated", "radK_bytes": None, "r_n": 0}
    radK_bytes = None
    if max_transient_bytes is not None:
        nrad = rad_ab_pairs.shape[0]                  # 2m(m-1)
        radK_bytes = len(ker) * nrad * (m2 * cur_r) * 8
        if radK_bytes > max_transient_bytes:
            # break BEFORE building d_n -> last known differential is d_{n-1}; leave
            # state['n'] = n-1 so truncated_at = state['n'].
            return {"status": "memory", "radK_bytes": radK_bytes, "r_n": None}
    radK = _build_radK(eng, rad_ab_pairs, ker, cur_r, p)
    gens = _independent_modulo(radK, ker, p)
    del radK                                          # free the big transient before _build_Dn
    r_n = len(gens)
    state["rks"][n] = r_n
    state["cols"][n] = gens
    if r_n == 0:
        state["n"] = n
        return {"status": "terminated", "radK_bytes": radK_bytes, "r_n": 0}
    if m2 * r_n > max_term_dim:
        state["n"] = n                                # cols[n] set -> truncated_at = n
        return {"status": "term", "radK_bytes": radK_bytes, "r_n": r_n}
    state["cur"] = _build_Dn(eng, gens, cur_r, p)
    state["cur_r"] = r_n
    state["n"] = n
    return {"status": "ok", "radK_bytes": radK_bytes, "r_n": r_n}


def minimal_resolution(A, N, p, max_term_dim=20000, max_transient_bytes=None):
    """Build the minimal free A^e-resolution of A over F_p up to degree N+1.

    Returns (rks, cols, eng, truncated_at):
      rks[n]          : r_n, the number of free A^e-generators of P_n;
      cols[n]         : list of r_n column vectors of d_n, each in (A^e)^{r_{n-1}};
      eng             : the AeEngine (mod p);
      truncated_at    : None if complete to N+1, else the highest degree whose
                        differential cols are known when a budget forced a stop -- the
                        term-dim cap m^2 * r_n > max_term_dim, OR (if max_transient_bytes
                        is set) the predicted radK transient. HH_* is exact through
                        truncated_at - 1 either way (the shared rule in
                        `minimal_homology_dims`).

    `max_transient_bytes` (bytes; None = off) bounds the PEAK allocation, not the final
    term. The minimal-generator selection at each degree first materializes the dense
    int64 `radK = rad(A^e) . ker` of shape (len(ker)*2m(m-1), m^2*cur_r) -- typically the
    single largest array in the whole build, one m^2-blocking factor LARGER than the term
    m^2*r_n that max_term_dim bounds. We predict its size BEFORE `_build_radK` allocates
    it and truncate gracefully if it would exceed the budget; otherwise a high-dim build
    OOM-kills the process (an uncatchable SIGKILL on the cluster -> the shard records
    nothing and times out) instead of recording an honest partial result.
    """
    state = _init_resolution(A, p)
    truncated_at = None
    for _ in range(1, N + 2):
        res = _advance_resolution(state, p, max_term_dim, max_transient_bytes)
        st = res["status"]
        if st == "terminated":
            for nn in range(state["n"] + 1, N + 2):
                state["rks"][nn] = 0
                state["cols"][nn] = []
            break
        if st == "memory":
            truncated_at = state["n"]                 # = n-1
            break
        if st == "term":
            truncated_at = state["n"]                 # = n
            break
    return state["rks"], state["cols"], state["eng"], truncated_at


def _contracted_degree(eng, gens_n, r_nm1, n):
    """A (x)_{A^e} d_n : one degree of the contracted complex.  `gens_n` = cols[n] (list of
    r_n generator blocks), `r_nm1` = rks[n-1].  Returns dbar_n of shape (m*r_nm1, m*r_n)."""
    m = eng.m
    m2 = eng.m2
    p = eng.p
    T = eng.T
    r_n = len(gens_n)
    M = np.zeros((m * r_nm1, m * r_n), dtype=np.int64)
    for j, g in enumerate(gens_n):
        for a in range(m):
            col = j * m + a
            for blk in range(r_nm1):
                w = g[blk * m2:(blk + 1) * m2]
                acc = np.zeros(m, dtype=np.int64)
                for uu in range(m):
                    for vv in range(m):
                        cf = w[uu * m + vv]
                        if cf % p == 0:
                            continue
                        mid = T[vv, a, :]          # e_v * e_a
                        for s in np.nonzero(mid % p)[0]:
                            acc = (acc + cf * mid[s] * T[s, uu, :]) % p   # (e_v e_a) e_u
                M[blk * m:(blk + 1) * m, col] = (M[blk * m:(blk + 1) * m, col] + acc) % p
    return M


def _contracted_complex(A, rks, cols, eng, N):
    """Apply A (x)_{A^e} (-) to every degree 1..N+1; returns {n: dbar_n}."""
    Dbar = {}
    for n in range(1, N + 2):
        Dbar[n] = _contracted_degree(eng, cols.get(n, []) or [], rks.get(n - 1, 0), n)
    return Dbar


def minimal_homology_dims(A, N, primes=(32003,), max_term_dim=20000,
                          max_transient_bytes=None):
    """dim HH_n(A; F_p) for n=0..N, each prime in `primes`, via the minimal A^e
    resolution (rebuilt per prime).  Returns {p: [dim HH_0, ..., dim HH_M]} where
    M = N unless a budget (term-dim cap or `max_transient_bytes` peak-memory guard)
    truncated the build earlier (then the list is shorter and exact up to its last
    entry; HH at the truncation degree is omitted because d_{n+1} is unknown there)."""
    out = {}
    for p in primes:
        rks, cols, eng, trunc = minimal_resolution(
            A, N, p, max_term_dim=max_term_dim, max_transient_bytes=max_transient_bytes)
        Dbar = _contracted_complex(A, rks, cols, eng, N)
        m = eng.m
        # if truncated at degree t, we know d_1..d_t but not d_{t+1}, so HH_n is exact
        # only for n <= t-1.
        last = (trunc - 1) if trunc is not None else N
        dims = []
        for n in range(0, last + 1):
            dimn = m * rks.get(n, 0)
            rn = rank_mod_p(Dbar[n], p) if (n >= 1 and rks.get(n, 0) > 0) else 0
            rnp1 = rank_mod_p(Dbar[n + 1], p) if rks.get(n + 1, 0) > 0 else 0
            dims.append(int(dimn - rn - rnp1))
        out[p] = dims
    return out


def hochschild_dimension(A, N, p=32003, max_term_dim=20000, max_transient_bytes=None):
    """The largest n <= N with r_n > 0 in the minimal A^e resolution over F_p, or
    'N+ (truncated)' info.  proj.dim_{A^e}(A) (the Hochschild dimension) is finite iff
    the minimal resolution terminates (r_n = 0 from some point on); this returns
    (pd, terminated) where `terminated` is True iff r_n hit 0 at pd (so pd is the
    actual Hochschild dimension), False iff the resolution was still nonzero at N
    (evidence of infinite Hochschild dimension / global dimension)."""
    rks, cols, eng, trunc = minimal_resolution(
        A, N, p, max_term_dim=max_term_dim, max_transient_bytes=max_transient_bytes)
    nonzero = [n for n in range(0, N + 2) if rks.get(n, 0) > 0]
    pd = max(nonzero) if nonzero else 0
    terminated = any(rks.get(n, 0) == 0 for n in range(1, N + 2)) and trunc is None
    return pd, terminated
