# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""
Hochschild homology engine for finite-dimensional algebras over a field.

We compute dim_k HH_n(A) = dim Tor_n^{A^e}(A,A) via the NORMALIZED bar complex
    C_n = A (x) (A/k.1)^{(x) n},   differential b (cyclic Hochschild differential),
and take ranks of the differentials over F_p for several primes p.  Rank over a
large prime is used as a faithful proxy for the characteristic-0 dimension; small
primes detect characteristic-specific behaviour.

dim HH_n = dim ker(b_n) - rank(b_{n+1}) = (dim C_n - rank b_n) - rank b_{n+1}.

An algebra is specified by:
  m            : dimension over k (basis indices 0..m-1)
  T            : integer 3-array of structure constants, T[i,j,:] = (e_i * e_j) in the basis
  unit         : integer vector giving 1_A in the basis (must have a coordinate == 1)

Everything is exact integer arithmetic reduced mod p.
"""

import numpy as np
import itertools


# ----------------------------------------------------------------------
# Exact rank over F_p (dense Gaussian elimination, vectorised per pivot)
# ----------------------------------------------------------------------
def rank_mod_p(M, p):
    """Rank of integer matrix M over F_p."""
    if M.size == 0:
        return 0
    A = (M.astype(np.int64) % p)
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        # find pivot in column c at or below row r
        piv = -1
        for i in range(r, rows):
            if A[i, c] % p != 0:
                piv = i
                break
        if piv == -1:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)  # modular inverse (p prime)
        A[r] = (A[r] * inv) % p
        # eliminate this column from all other rows
        col = A[:, c].copy()
        col[r] = 0
        nz = np.nonzero(col)[0]
        if nz.size:
            A[nz] = (A[nz] - np.outer(col[nz], A[r])) % p
        r += 1
        if r == rows:
            break
    return r


# ----------------------------------------------------------------------
# Algebra wrapper: change to a unit-adapted basis so 1_A is a basis vector,
# letting the reduced space A/k.1 be "drop coordinate t".
# ----------------------------------------------------------------------
class Algebra:
    """A finite-dimensional algebra in a unit-adapted basis.

    Constructed from `(m, T, unit)`:
      * `m`    -- dimension over k (basis indices 0..m-1);
      * `T`    -- integer structure-constant 3-tensor, `T[i,j,:] = e_i * e_j`;
      * `unit` -- integer coordinate vector of `1_A` (must have a coordinate == 1).

    The constructor changes to a basis `f_i = B[:,i]` in which `1_A` is the single
    basis vector `f_t` (where `unit[t] == 1`). This is what lets multi-vertex
    algebras (`1_A = sum e_i`) work transparently: the reduced space `A/k.1` is then
    just "drop coordinate t". The structure constants are recomputed exactly in the
    f-basis (`B` is required unimodular, so the inverse is integral).

    Attributes:
        m: Dimension over k.
        t: Index of the unit basis vector.
        T: Structure constants in the f-basis, shape `(m, m, m)`.
        R: Reduced basis indices (all but `t`).
        Rpos: Map from index to reduced position in the reduced basis.
        mr: Reduced dimension `m - 1`.
    """

    def __init__(self, m, T, unit, name="A"):
        self.name = name
        self.m = m
        unit = np.array(unit, dtype=np.int64)
        # find t with unit[t] == 1
        ts = np.nonzero(unit == 1)[0]
        assert ts.size > 0, "unit must have a coordinate equal to 1"
        t = int(ts[0])
        self.t = t
        # Original-basis inputs, kept for vertex-idempotent detection (Plan 13:
        # resolutions_minimal reads the orthogonal idempotents off the unit's
        # 1-coordinates; the f-basis change below only touches column t, so the
        # non-t vertex vectors stay standard basis vectors).
        self.unit_input = unit.copy()
        self.T_input = np.array(T, dtype=np.int64, copy=True)
        # change-of-basis matrix B: identity with column t replaced by unit
        B = np.eye(m, dtype=np.int64)
        B[:, t] = unit
        Binv = np.round(np.linalg.inv(B.astype(float))).astype(np.int64)
        assert np.array_equal(B @ Binv, np.eye(m, dtype=np.int64)), "B not unimodular"
        # recompute structure constants in new basis f_i = B[:,i]
        # f_i * f_j = sum_{p,q} B[p,i] B[q,j] T[p,q,:]  (old coords) -> Binv @ (...)
        Tnew = np.zeros((m, m, m), dtype=np.int64)
        for i in range(m):
            for j in range(m):
                v = np.zeros(m, dtype=np.int64)
                bi = B[:, i]
                bj = B[:, j]
                pi = np.nonzero(bi)[0]
                pj = np.nonzero(bj)[0]
                for a in pi:
                    for b in pj:
                        v += bi[a] * bj[b] * T[a, b, :]
                Tnew[i, j, :] = Binv @ v
        self.T = Tnew
        # reduced index set R = all basis indices except t
        self.R = [i for i in range(m) if i != t]
        self.Rpos = {idx: k for k, idx in enumerate(self.R)}  # idx -> position in reduced basis
        self.mr = m - 1
        # vertex idempotent indices (original basis; None if the basis is not
        # recognizably path-type -- resolutions_minimal then guards loudly)
        from quiverlab.engine.resolutions_minimal import _vertex_indices
        self.vertices = _vertex_indices(self)

    def mult_full(self, i, j):
        """e_i * e_j as a full length-m vector (f-basis coords)."""
        return self.T[i, j, :]


# ----------------------------------------------------------------------
# Normalized bar complex: enumerate basis of C_n and build b_n as a matrix
# basis element of C_n: (a0, r1,...,rn) with a0 in 0..m-1, ri in R
# ----------------------------------------------------------------------
def cn_basis(alg, n):
    if n == 0:
        return [(a0,) for a0 in range(alg.m)]
    Rt = alg.R
    out = []
    for a0 in range(alg.m):
        for tail in itertools.product(Rt, repeat=n):
            out.append((a0,) + tail)
    return out


def differential_matrix(alg, n, basis_n, index_nm1):
    """Matrix of b_n : C_n -> C_{n-1}.  shape (dim C_{n-1}, dim C_n)."""
    m = alg.m
    t = alg.t
    rows = len(index_nm1)
    cols = len(basis_n)
    M = np.zeros((rows, cols), dtype=np.int64)
    if n == 0:
        return M  # b_0 = 0
    for cidx, gen in enumerate(basis_n):
        a0 = gen[0]
        rs = gen[1:]  # length n, entries in R
        # term i=0:  (a0 * r1 ; r2,...,rn) in A-slot, sign +
        prod = alg.mult_full(a0, rs[0])  # full vector
        nz = np.nonzero(prod)[0]
        for a0p in nz:
            key = (int(a0p),) + tuple(rs[1:])
            M[index_nm1[key], cidx] += prod[a0p]
        # terms i=1..n-1: interior multiplications, projected to reduced slot
        for i in range(1, n):
            prod = alg.mult_full(rs[i - 1], rs[i])  # full vector; project: drop coord t
            sign = -1 if (i % 2 == 1) else 1
            nz = np.nonzero(prod)[0]
            for idx in nz:
                if idx == t:
                    continue  # killed in A/k.1
                newtail = list(rs)
                # replace positions i-1,i (0-based within rs) by single 'idx'
                merged = rs[:i - 1] + (int(idx),) + rs[i + 1:]
                key = (a0,) + merged
                M[index_nm1[key], cidx] += sign * prod[idx]
        # term i=n: (rn * a0 ; r1,...,r_{n-1}) in A-slot, sign (-1)^n
        prod = alg.mult_full(rs[n - 1], a0)
        sign = -1 if (n % 2 == 1) else 1
        nz = np.nonzero(prod)[0]
        for a0p in nz:
            key = (int(a0p),) + tuple(rs[:n - 1])
            M[index_nm1[key], cidx] += sign * prod[a0p]
    return M


def hochschild_homology_dims(alg, N, primes=(32003, 2, 3, 5), resolution=None):
    """Return dict p -> list [dim HH_0,...,dim HH_N].

    resolution: a Resolution backend supplying term_basis / differential_matrix
    (default: the normalized bar complex). Faster/smaller backends are cross-checked
    against this default -- see notes/B1_chouhy_solotar_plan.md.
    """
    from quiverlab.engine.resolutions import _default
    resolution = _default(resolution)
    # precompute bases and indices up to N+1
    bases = {}
    indices = {}
    for n in range(0, N + 2):
        bases[n] = resolution.term_basis(alg, n)
        indices[n] = {g: i for i, g in enumerate(bases[n])}
    # differentials b_n for n=1..N+1
    ranks = {p: {} for p in primes}
    Bmats = {}
    for n in range(1, N + 2):
        Bmats[n] = resolution.differential_matrix(alg, n, bases[n], indices[n - 1])
    for p in primes:
        for n in range(1, N + 2):
            ranks[p][n] = rank_mod_p(Bmats[n], p)
    out = {}
    for p in primes:
        dims = []
        for n in range(0, N + 1):
            dimCn = len(bases[n])
            rn = ranks[p].get(n, 0) if n >= 1 else 0  # rank b_n  (b_0 = 0)
            rnp1 = ranks[p].get(n + 1, 0)
            dims.append(dimCn - rn - rnp1)
        out[p] = dims
    return out


# ----------------------------------------------------------------------
# Builders for test algebras
# ----------------------------------------------------------------------
def truncated_polynomial(a, coeff_field_note=""):
    """k[x]/(x^a). basis 1,x,...,x^{a-1}. unit = e_0."""
    m = a
    T = np.zeros((m, m, m), dtype=np.int64)
    for i in range(m):
        for j in range(m):
            if i + j < m:
                T[i, j, i + j] = 1
    unit = np.zeros(m, dtype=np.int64)
    unit[0] = 1
    return Algebra(m, T, unit, name=f"k[x]/(x^{a})")


def two_gen_local(rel_x2, rel_y2, rel_yx, name):
    """
    Local algebra k<x,y> with basis {1, x, y, xy}, dimension 4, where the
    products x*x, y*y, y*x are prescribed as linear combinations of the basis
    (each given as a length-4 vector in basis order [1,x,y,xy]).  We always set
    x*y = xy (the fourth basis vector).  Associativity is NOT verified here -- callers
    must run check_associative (the scan drivers do).
    basis index: 0=1, 1=x, 2=y, 3=xy.
    """
    m = 4
    T = np.zeros((m, m, m), dtype=np.int64)
    # unit multiplication
    for j in range(m):
        T[0, j, j] = 1
        T[j, 0, j] = 1
    # x*y = xy
    T[1, 2, :] = np.array([0, 0, 0, 1], dtype=np.int64)
    # prescribed relations
    T[1, 1, :] = np.array(rel_x2, dtype=np.int64)   # x*x
    T[2, 2, :] = np.array(rel_y2, dtype=np.int64)   # y*y
    T[2, 1, :] = np.array(rel_yx, dtype=np.int64)   # y*x
    # products involving xy (=basis 3): x*xy, xy*x, y*xy, xy*y, xy*xy
    # Compute them associatively from the above where forced; default 0 (radical^3 = 0).
    # We take the convention rad^3 = 0 so xy * anything in rad = 0 and (xy)^2 = 0,
    # x*xy = x*(x*y)=(x*x)*y, etc. To stay consistent we DERIVE these:
    def mult_vec(u, v):
        out = np.zeros(m, dtype=np.int64)
        for i in np.nonzero(u)[0]:
            for j in np.nonzero(v)[0]:
                out += u[i] * v[j] * T[i, j, :]
        return out
    ex = np.array([0, 1, 0, 0]); ey = np.array([0, 0, 1, 0]); exy = np.array([0, 0, 0, 1])
    # x*xy = (x*x)*y
    T[1, 3, :] = mult_vec(mult_vec(ex, ex), ey)
    # xy*x = x*(y*x)
    T[3, 1, :] = mult_vec(ex, mult_vec(ey, ex))
    # y*xy = (y*x)*y
    T[2, 3, :] = mult_vec(mult_vec(ey, ex), ey)
    # xy*y = x*(y*y)
    T[3, 2, :] = mult_vec(ex, mult_vec(ey, ey))
    # xy*xy = x*(y*x)*y
    T[3, 3, :] = mult_vec(mult_vec(ex, mult_vec(ey, ex)), ey)
    unit = np.array([1, 0, 0, 0], dtype=np.int64)
    return Algebra(m, T, unit, name=name)


def check_associative(alg, p=32003):
    """Verify (e_i e_j) e_k = e_i (e_j e_k) mod p on basis."""
    m = alg.m
    T = alg.T % p
    for i in range(m):
        for j in range(m):
            for k in range(m):
                left = np.zeros(m, dtype=np.int64)
                ij = T[i, j, :]
                for a in np.nonzero(ij)[0]:
                    left += ij[a] * T[a, k, :]
                right = np.zeros(m, dtype=np.int64)
                jk = T[j, k, :]
                for b in np.nonzero(jk)[0]:
                    right += jk[b] * T[i, b, :]
                if np.any((left - right) % p != 0):
                    return False, (i, j, k)
    return True, None
