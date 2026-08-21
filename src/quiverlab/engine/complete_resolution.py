"""Complete resolution of A over A^e, and Tate-Hochschild dimensions (Plan 76 / R3).

A complete resolution of the diagonal bimodule A over B = A^e is an ACYCLIC complex of
finitely generated projective B-modules

    T :  ... -> T_1 -> T_0 -> T_-1 -> T_-2 -> ...

that coincides with a projective resolution of A in high degrees.  Tate-Hochschild
(= singular Hochschild) cohomology is its cohomology in EVERY integer degree,

    HHhat^m(A) = H^m( Hom_{A^e}(T, A) ),      HHhat_m(A) = H_m( A (x)_{A^e} T ),

so it extends ordinary HH below degree 0.  Above the Gorenstein dimension d of A^e it
agrees with the ordinary theory (Bergh-Jorgensen): HHhat^m = HH^m for m >= d+1, and
d = 0 for self-injective A, so agreement starts at m = 1.  Degree 0 does NOT agree:
HHhat^0 is the STABLE centre, a proper quotient of Z(A) = HH^0.

SCOPE (v1): SELF-INJECTIVE A only.  Then B = A^e is self-injective, so the B-dual of a
projective is projective and the dual-splice below is well defined.  A Gorenstein but
non-self-injective algebra needs a periodicity-extension complete resolution instead --
DEFERRED, refused loudly (see `_require_scope`).

THE CONSTRUCTION (derived, then arbitrated by exactness -- see below).  Write
eps_{v,w} = e_v (x) e_w^op for the corner idempotents of B, so the indecomposable
projectives are B eps_{v,w} = A e_v (x) e_w A.  Let lambda be a Frobenius form of A,
nu its Nakayama automorphism (lambda(ab) = lambda(b nu(a))) and pi the VERTEX Nakayama
permutation (soc(e_v A) = S_{pi(v)}).  Two functors:

  * Psi(M) = D(M) = Hom_k(M, k) made a LEFT B-module through the swap anti-automorphism
    a (x) b^op |-> b (x) a^op of B.  Psi is an exact contravariant self-equivalence with
    Psi^2 = id, so it turns projective resolutions into injective coresolutions; since B
    is self-injective those are again projective.
  * theta = nu^{-1} (x) id, the twist with theta^*(A) = D(A).  (The identification is
    phi(a)(m) = lambda(am), which satisfies phi(nu^{-1}(x) a y) = x phi(a) y -- hence
    nu^{-1}, NOT nu.  This is the whole content of the twist and the ONLY place the
    non-symmetric case differs from the symmetric one, where nu = id kills it.)

Because Psi(theta^*(A)) = Psi(D(A)) = A, applying Psi to the minimal projective
resolution P_* -> theta^*(A) = D(A) yields an injective coresolution of A.  So

    T_n = P_n  (n >= 0),        T_{-n-1} = Psi(theta^*(P_n))  (n >= 0),

with, on corner tags and on differential entries (c in B, coordinates c[a*m+b] for
f_a (x) f_b^op),

    tag   (v, w)  |->  (pi(w), v)                          # `_dual_tag`
    entry f_a (x) f_b^op  |->  nu(f_b) (x) f_a^op          # `_Phi`, and the matrix
                                                           # of the differential is
                                                           # transposed (Psi is
                                                           # contravariant)

and the splice joint d_0 : T_0 = P_0 -> T_-1 is the composite P_0 ->>  A  >-> T_-1 of
the augmentation with the injective envelope.  Its matrix is an exact solve: the
u-block of the image of the generator eps_{v,v} is Z^{(v)} . eps_{dual tag of u} with

    Z^{(v)} = G^{-1} R^{(v)} G^{-1},   ``R^{(v)}[i][j]`` = lambda( nu^{-1}(f_j) f_i e_v ),

``G[i][j]`` = lambda(f_i f_j) the (nondegenerate) Gram matrix of the form.

ARBITRATION -- WHAT ACTUALLY PINS THE CONVENTIONS.  d.d = 0 is NOT sufficient here and
must not be trusted as the arbiter: the negative half is carried as B-module maps in
AMBIENT (A^e)^r coordinates, and the composite is taken there WITHOUT reference to the
tags, so a wrong pi leaves d.d = 0 intact (measured; pinned by
tests/engine/test_complete_resolution_p76.py::test_dd_zero_alone_does_not_discriminate_pi).
What a wrong pi does break is (a) the CORNER TYPING of the dual half -- entries stop
lying in eps_t . B . eps_s for their declared tags, gated at build time by
`_assert_corner_typed` -- and (b) exactness, gated by `assert_acyclic`.  Together with
the positive-degree agreement against the shipped HH engine those are the certificates;
the twist direction (nu vs nu^{-1}) was settled by them, not assumed.

Over GF(p), int64; the sibling of `engine/resolutions_minimal.py`, whose corner
machinery (tags, corner bases, the Hom / (x) collapses) is reused verbatim so the
negative half costs no new linear algebra.
"""
import numpy as np

from quiverlab.engine.adapter import to_engine
from quiverlab.engine.hh_engine import rank_mod_p
from quiverlab.engine.resolutions_minimal import (_cohomology_degree,
                                                  _contracted_degree,
                                                  _corner_cohomology_degree,
                                                  _corner_contracted_degree,
                                                  minimal_resolution)
from quiverlab.errors import QuiverlabError

__all__ = ["CompleteResolution", "complete_resolution", "nakayama_permutation",
           "tate_cohomology_dims", "tate_homology_dims"]


# ---------------------------------------------------------------------------
# exact GF(p) helpers
# ---------------------------------------------------------------------------
def _inv_mod_p(M, p, what="matrix"):
    """Exact inverse mod p by Gauss-Jordan; loud when singular."""
    n = M.shape[0]
    aug = np.concatenate([M % p, np.eye(n, dtype=np.int64)], axis=1)
    row = 0
    for col in range(n):
        piv = -1
        for i in range(row, n):
            if aug[i, col] % p:
                piv = i
                break
        if piv < 0:
            raise QuiverlabError(
                f"{what} is singular mod {p}, so the complete resolution cannot be built",
                hint="the Frobenius form must stay nondegenerate over the chosen prime")
        aug[[row, piv]] = aug[[piv, row]]
        aug[row] = (aug[row] * pow(int(aug[row, col]), p - 2, p)) % p
        for i in range(n):
            if i != row and aug[i, col] % p:
                aug[i] = (aug[i] - aug[i, col] * aug[row]) % p
        row += 1
    return aug[:, n:] % p


def _left_mult_matrix(T, vec, p):
    """(m, m) matrix of x |-> vec . x on A (columns = images of the basis)."""
    m = T.shape[0]
    out = np.zeros((m, m), dtype=np.int64)
    for c in np.nonzero(vec % p)[0]:
        out = (out + int(vec[c]) * T[c, :, :].T) % p
    return out


def _right_mult_matrix(T, vec, p):
    """(m, m) matrix of x |-> x . vec on A (columns = images of the basis)."""
    m = T.shape[0]
    out = np.zeros((m, m), dtype=np.int64)
    for c in np.nonzero(vec % p)[0]:
        out = (out + int(vec[c]) * T[:, c, :].T) % p
    return out


# ---------------------------------------------------------------------------
# the vertex Nakayama permutation (public helper)
# ---------------------------------------------------------------------------
def nakayama_permutation(A):
    """{vertex: vertex} with soc(e_v A) = S_{pi(v)}, keyed by QUIVER VERTEX LABEL.

    This is the VERTEX-level permutation, NOT the d x d matrix that
    `Algebra.nakayama_automorphism()` returns.  It is the identity exactly when A is
    weakly symmetric (in particular for every symmetric A), and it is what the dual
    half's corner tags are permuted by.  Raises when A is not self-injective (there is
    no such permutation) or has no path-type basis.
    """
    from quiverlab.invariants.frobenius import nakayama_data
    ok, perm, _gens = nakayama_data(A)
    if not ok:
        raise QuiverlabError(
            "algebra is not self-injective (socle criterion), so it has no Nakayama "
            "vertex permutation")
    labels = _vertex_labels(A)
    return {labels[v]: labels[w] for v, w in perm.items()}


def _vertex_labels(A):
    """{basis index of e_v: quiver vertex label}, matched on the basis labels."""
    if A.quiver is None or not A.basis_labels:
        raise QuiverlabError(
            "the complete resolution needs a quiver presentation (path-type basis); "
            "this algebra has none",
            hint="present the algebra via Quiver(...).algebra(...)")
    out = {}
    for v in A.quiver.vertices:
        want = f"e_{v}"
        hits = [i for i, lab in enumerate(A.basis_labels) if lab == want]
        if len(hits) != 1:
            raise QuiverlabError(
                f"vertex idempotent '{want}' is not a single basis vector -- the basis "
                f"is not path-type, so the corner tags cannot be named")
        out[hits[0]] = v
    return out


# ---------------------------------------------------------------------------
# scope gate
# ---------------------------------------------------------------------------
def _require_scope(A):
    """Self-injective, path-presented, over GF(p).  Every boundary is loud."""
    _vertex_labels(A)                                  # presentation gate, loud
    if not A.is_selfinjective():
        from quiverlab.modules.homdims import gorenstein_dimension
        try:
            gor = gorenstein_dimension(A)
        except Exception:                              # pragma: no cover - advisory only
            gor = None
        if gor is not None and any(d is not None for d in _as_pair(gor)):
            raise QuiverlabError(
                "Tate-Hochschild cohomology of a Gorenstein but NOT self-injective "
                "algebra is DEFERRED: the D(P_n)-dual splice needs A^e self-injective "
                "for the dual half to be projective, which fails here",
                hint="the periodicity-extension complete resolution (Usui) is the named "
                     "follow-up task; the eventual-periodicity certificate is still "
                     "computed for this algebra")
        raise QuiverlabError(
            "Tate-Hochschild cohomology needs the complete resolution, which requires "
            "A^e to be Gorenstein; this algebra is not self-injective and no Gorenstein "
            "certificate was produced")


def _as_pair(gor):
    if isinstance(gor, (tuple, list)):
        return tuple(gor)
    return (gor,)


# ---------------------------------------------------------------------------
# the complete resolution
# ---------------------------------------------------------------------------
class CompleteResolution:
    """T_n for -N-1 <= n <= N+1 with its differentials, over GF(p).

    Attributes mirror the minimal engine's so the shipped corner collapses apply
    verbatim:
      ranks[n]  : number of projective generators of T_n;
      tags[n]   : one corner tag (v, w) per generator (basis indices), or None on the
                  single-vertex free path, where every term is B-free;
      diffs[n]  : the r_n images of the generators of T_n in the ambient (A^e)^{r_{n-1}}
                  -- exactly the format `minimal_resolution` returns for `cols`.
    """

    def __init__(self, algebra, eng, ctx, p, top, ranks, tags, diffs):
        self.algebra = algebra
        self.eng = eng
        self.ctx = ctx
        self.p = p
        self.top = top
        self.ranks = ranks
        self.tags = tags
        self.diffs = diffs

    # -- shape ------------------------------------------------------------
    def term_rank(self, n):
        """Number of projective generators of T_n."""
        return int(self.ranks.get(n, 0))

    def term_dim(self, n):
        """dim_k T_n."""
        if self.ctx is None:
            return self.eng.m2 * self.term_rank(n)
        return int(sum(self.ctx.corner_basis[tg].shape[1]
                       for tg in self.tags.get(n, [])))

    def degrees(self):
        return sorted(self.diffs)

    # -- self-certs -------------------------------------------------------
    def assert_dd_zero(self):
        """d_{n-1} . d_n = 0 for every consecutive pair, INCLUDING both splice joints.

        Computed in ambient A^e coordinates, so it certifies that T is a complex --
        and NOTHING about the corner tags (see the module docstring; `assert_acyclic`
        is what sees those)."""
        eng, p, m, m2 = self.eng, self.p, self.eng.m, self.eng.m2
        for n in sorted(self.diffs):
            if (n - 1) not in self.diffs:
                continue
            r_prev = self.term_rank(n - 2)
            for j, g in enumerate(self.diffs[n]):
                acc = np.zeros(m2 * max(r_prev, 1), dtype=np.int64)
                for i in range(self.term_rank(n - 1)):
                    c = g[i * m2:(i + 1) * m2]
                    for a in range(m):
                        for b in range(m):
                            cf = int(c[a * m + b]) % p
                            if cf:
                                acc = (acc + cf * eng.apply_block(
                                    a, b, self.diffs[n - 1][i], max(r_prev, 1))) % p
                if np.any(acc % p):
                    raise AssertionError(
                        f"d_{n - 1} . d_{n} != 0 on generator {j} -- the splice is not a "
                        f"complex")

    def _term_matrix(self, n):
        """Full k-linear matrix of d_n : T_n -> ambient(T_{n-1})."""
        eng, p, m, m2 = self.eng, self.p, self.eng.m, self.eng.m2
        r_prev = max(self.term_rank(n - 1), 1)
        cols = []
        for j, g in enumerate(self.diffs.get(n, [])):
            basis = (self.ctx.corner_basis[self.tags[n][j]] if self.ctx is not None
                     else np.eye(m2, dtype=np.int64))
            for c in range(basis.shape[1]):
                x = basis[:, c]
                col = np.zeros(m2 * r_prev, dtype=np.int64)
                for a in range(m):
                    for b in range(m):
                        cf = int(x[a * m + b]) % p
                        if cf:
                            col = (col + cf * eng.apply_block(a, b, g, r_prev)) % p
                cols.append(col)
        if not cols:
            return np.zeros((m2 * r_prev, 0), dtype=np.int64)
        return np.stack(cols, axis=1) % p

    def assert_acyclic(self, window=None):
        """T is exact at every degree in `window` (rank-nullity, tag-sensitive).

        This is the arbiter for the Nakayama twist: with a wrong vertex permutation the
        dual half's entries fall outside their declared corners, the expanded columns
        collapse and exactness fails loudly."""
        lo, hi = window if window is not None else (-self.top, self.top)
        ranks = {}
        for n in (set(range(lo, hi + 2)) & set(self.diffs)):
            M = self._term_matrix(n)
            ranks[n] = int(rank_mod_p(M, self.p)) if M.size else 0
        for n in range(lo, hi + 1):
            if n not in ranks or (n + 1) not in ranks:
                continue
            defect = self.term_dim(n) - ranks[n] - ranks[n + 1]
            if defect != 0:
                raise AssertionError(
                    f"the complete resolution is NOT exact at degree {n}: "
                    f"dim T_{n} = {self.term_dim(n)}, rank d_{n} = {ranks[n]}, "
                    f"rank d_{n + 1} = {ranks[n + 1]} (homology {defect})")

    # -- collapses --------------------------------------------------------
    def cochain_matrix(self, n):
        """delta^{n-1} : Hom(T_{n-1}, A) -> Hom(T_n, A) (the Plan-16 a.w.b side)."""
        if n not in self.diffs:
            return None
        if self.ctx is None:
            return _cohomology_degree(self.eng, self.diffs[n],
                                      self.term_rank(n - 1), n)
        return _corner_cohomology_degree(self.eng, self.ctx, self.diffs[n],
                                         self.tags.get(n, []),
                                         self.tags.get(n - 1, []), self.p)

    def chain_matrix(self, n):
        """A (x) d_n : A (x) T_n -> A (x) T_{n-1} (the b.w.a side, swapped tags)."""
        if n not in self.diffs:
            return None
        if self.ctx is None:
            return _contracted_degree(self.eng, self.diffs[n],
                                      self.term_rank(n - 1), n)
        return _corner_contracted_degree(self.eng, self.ctx, self.diffs[n],
                                         self.tags.get(n, []),
                                         self.tags.get(n - 1, []), self.p)

    def cochain_dim(self, n):
        """dim_k Hom_{A^e}(T_n, A) = sum over generators of dim e_v A e_w."""
        if self.ctx is None:
            return self.eng.m * self.term_rank(n)
        return int(sum(self.ctx.corner_dim_A((tg[1], tg[0]))
                       for tg in self.tags.get(n, [])))

    def chain_dim(self, n):
        """dim_k A (x)_{A^e} T_n = sum over generators of dim e_w A e_v."""
        if self.ctx is None:
            return self.eng.m * self.term_rank(n)
        return int(sum(self.ctx.corner_dim_A(tg) for tg in self.tags.get(n, [])))


def complete_resolution(A, N, p, max_term_dim=20000, _skip_corner_gate=False):
    """Build the complete resolution of A over A^e in the window [-N-1, N+1], GF(p).

    `A` is a quiver-presented, SELF-INJECTIVE `quiverlab.Algebra` over GF(p).  The
    positive half is the shipped minimal A^e-resolution; the negative half is its
    Nakayama-twisted A^e-dual (module docstring).  `_skip_corner_gate` is for the
    self-cert battery only -- it disables the build-time corner-typing assertion so a
    test can exhibit what d.d = 0 alone fails to catch.
    """
    _require_scope(A)
    labels = _vertex_labels(A)
    multi = len(labels) > 1
    Aq = A if multi else A.unit_adapted()          # Plan-18: never unit-adapt multi-vertex
    eng_alg = to_engine(Aq)
    p = int(p)

    rks, cols, eng, truncated = minimal_resolution(eng_alg, N + 1, p,
                                                   max_term_dim=max_term_dim)
    if truncated is not None:
        raise QuiverlabError(
            f"the minimal A^e-resolution of this algebra was truncated at degree "
            f"{truncated} by the term-dimension budget, so the complete resolution "
            f"cannot be spliced to degree {N}",
            hint="raise max_term_dim, or ask for a smaller top degree")
    ctx = eng.corner_ctx
    m, m2, T = eng.m, eng.m2, eng.T % p

    lam, G, Ginv, nu, nu_inv = _form_data(Aq, eng_alg, eng, p)
    pi = _pi_indices(A, Aq, labels, ctx, eng_alg)

    def dual_tag(tg):
        return (pi[tg[1]], tg[0])                  # (v, w) |-> (pi(w), v)

    def Phi(c):                                    # f_a (x) f_b |-> nu(f_b) (x) f_a
        return ((nu @ (c % p).reshape(m, m).T) % p).reshape(m2)

    tags = {n: list(eng.corner_tags.get(n, [])) for n in range(N + 2)} if ctx else {}
    ranks = {n: int(rks.get(n, 0)) for n in range(N + 2)}
    diffs = {n: list(cols.get(n) or []) for n in range(1, N + 2)}

    # negative terms T_{-n-1} = Psi(theta^* P_n)
    for n in range(N + 1):
        ranks[-n - 1] = int(rks.get(n, 0))
        if ctx is not None:
            tags[-n - 1] = [dual_tag(tg) for tg in eng.corner_tags.get(n, [])]

    # negative differentials d_{-n} : T_{-n} -> T_{-n-1}, the transpose of d_n
    for n in range(1, N + 1):
        r_gen, r_blk = ranks.get(n, 0), ranks.get(n - 1, 0)
        out = []
        for i in range(r_blk):
            vec = np.zeros(m2 * max(r_gen, 1), dtype=np.int64)
            for j in range(r_gen):
                vec[j * m2:(j + 1) * m2] = Phi(cols[n][j][i * m2:(i + 1) * m2])
            out.append(vec)
        diffs[-n] = out

    diffs[0] = _splice_joint(eng, ctx, T, m, m2, p, G, Ginv, lam, nu_inv,
                             tags if ctx is not None else None, ranks, eng_alg)

    res = CompleteResolution(A, eng, ctx, p, N, ranks, tags, diffs)
    if ctx is not None and not _skip_corner_gate:
        _assert_corner_typed(res)
    return res


def _form_data(Aq, eng_alg, eng, p):
    """(lambda, G, G^-1, nu, nu^-1) in the ENGINE basis.

    `frobenius_form_generic` works in Aq's basis; the engine's f-basis differs only in
    the unit column (f_t = 1_A, f_i = e_i otherwise), so the covector transports by
    lambda(f_t) = lambda(1_A) and is unchanged elsewhere.  G and nu are then rebuilt
    from the engine's own structure constants -- no matrix conjugation to get wrong."""
    from quiverlab.invariants.frobenius import frobenius_form_generic
    lam_pub, _G_pub = frobenius_form_generic(Aq)
    m, t = eng.m, eng_alg.t
    unit = np.array([int(c) for c in Aq.unit], dtype=np.int64)
    lam = np.array([int(c) % p for c in lam_pub], dtype=np.int64)
    lam[t] = int((unit * np.array([int(c) % p for c in lam_pub])).sum() % p)
    T = eng.T % p
    G = np.array([[int(lam @ T[i, j, :]) % p for j in range(m)] for i in range(m)],
                 dtype=np.int64)
    if rank_mod_p(G, p) != m:
        raise QuiverlabError(
            f"the Frobenius form of this algebra degenerates mod {p}, so the "
            f"A^e-duality the complete resolution is built from does not exist here",
            hint="choose a prime that does not divide the form's discriminant")
    Ginv = _inv_mod_p(G, p, "the Frobenius Gram matrix")
    nu = (Ginv @ G.T) % p                          # lambda(ab) = lambda(b nu(a))
    nu_inv = _inv_mod_p(nu, p, "the Nakayama automorphism")
    return lam, G, Ginv, nu, nu_inv


def _pi_indices(A, Aq, labels, ctx, eng_alg):
    """The Nakayama vertex permutation as a map on ENGINE vertex indices.

    Read from `nakayama_permutation` (the socle computation) so that a monkeypatched /
    wrong permutation really reaches the tags -- the corner-typing gate is then a
    genuine cross-check against the form-derived nu, not a tautology."""
    if ctx is None:
        return {v: v for v in (eng_alg.vertices or [0])}
    perm = nakayama_permutation(A)
    idx_of = {lab: i for i, lab in labels.items()}
    out = {}
    for i, lab in labels.items():
        if lab not in perm:
            raise QuiverlabError(
                f"the Nakayama permutation does not name vertex {lab}")
        out[i] = idx_of[perm[lab]]
    if sorted(out.values()) != sorted(out):
        raise QuiverlabError("the Nakayama vertex map is not a permutation")
    return out


def _splice_joint(eng, ctx, T, m, m2, p, G, Ginv, lam, nu_inv, tags, ranks, eng_alg):
    """d_0 : T_0 = P_0 ->> A >-> T_-1, the norm map of the splice.

    The image of the generator eps_{v,v} is determined by
    Lambda(swap(y) . z) = lambda( mu(theta(y)) . e_v ) for all y, which in the basis
    reads (G Z G)_{i,j} = lambda(nu^{-1}(f_j) f_i e_v), i.e. Z = G^-1 R G^-1; the
    u-block is then Z cut down by the target block's corner idempotent."""
    r_m1 = ranks.get(-1, 0)
    vertices = list(ctx.vertices) if ctx is not None else [None]
    src_tags = tags[0] if tags is not None else [None] * ranks.get(0, 0)
    tgt_tags = tags[-1] if tags is not None else [None] * r_m1
    out = []
    for blk, tg0 in enumerate(src_tags):
        ev = (ctx.idem[tg0[0]] % p if ctx is not None
              else _unit_vector(m, eng_alg.t))
        R = np.zeros((m, m), dtype=np.int64)
        Rev = _right_mult_matrix(T, ev, p)
        for j in range(m):
            nuj = nu_inv[:, j]                     # nu^{-1}(f_j)
            left = _left_mult_matrix(T, nuj, p)    # x |-> nu^{-1}(f_j) . x
            for i in range(m):
                col = (Rev @ left[:, i]) % p       # nu^{-1}(f_j) f_i e_v
                R[i, j] = int(lam @ col) % p
        Z = (Ginv @ ((R @ Ginv) % p)) % p
        vec = np.zeros(m2 * max(r_m1, 1), dtype=np.int64)
        for u, tgu in enumerate(tgt_tags):
            z = Z.reshape(m2)
            if ctx is not None:
                z = _sandwich(ctx, T, m, p, z, tg0, tgu)
            vec[u * m2:(u + 1) * m2] = z % p
        out.append(vec)
    return out


def _unit_vector(m, t):
    v = np.zeros(m, dtype=np.int64)
    v[t] = 1
    return v


def _sandwich(ctx, T, m, p, z, left_tag, right_tag):
    """eps_left . z . eps_right in B, with B's convention
    (e_a (x) e_b)(e_c (x) e_d) = (e_a e_c) (x) (e_d e_b)."""
    el, er = ctx.idem[left_tag[0]] % p, ctx.idem[left_tag[1]] % p
    fl, fr = ctx.idem[right_tag[0]] % p, ctx.idem[right_tag[1]] % p
    Z = (z % p).reshape(m, m)
    first = (_right_mult_matrix(T, fl, p) @ _left_mult_matrix(T, el, p)) % p
    second = (_left_mult_matrix(T, fr, p) @ _right_mult_matrix(T, er, p)) % p
    return ((first @ ((Z @ second.T) % p)) % p).reshape(m * m)


def _assert_corner_typed(res):
    """Every differential entry lies in eps_{target tag} . B . eps_{source tag}.

    This is where a WRONG Nakayama permutation dies: the dual half's entries are built
    from the form-derived nu, its tags from the socle-derived pi, and the two must
    agree.  d.d = 0 cannot see this (module docstring)."""
    ctx, eng, p, m = res.ctx, res.eng, res.p, res.eng.m
    m2 = eng.m2
    T = eng.T % p
    for n, gens in res.diffs.items():
        src = res.tags.get(n, [])
        tgt = res.tags.get(n - 1, [])
        for j, g in enumerate(gens):
            if j >= len(src):
                continue
            for i, tg in enumerate(tgt):
                c = g[i * m2:(i + 1) * m2] % p
                if not np.any(c):
                    continue
                back = _sandwich(ctx, T, m, p, c, src[j], tg)
                if not np.array_equal(back % p, c % p):
                    raise QuiverlabError(
                        f"the complete resolution is not corner-typed at degree {n} "
                        f"(generator {j}, block {i}): an entry declared to live in "
                        f"eps_{src[j]} . A^e . eps_{tg} does not",
                        hint="the Nakayama vertex permutation used for the dual tags "
                             "disagrees with the Nakayama automorphism the entries "
                             "were twisted by")


# ---------------------------------------------------------------------------
# Tate dimensions
# ---------------------------------------------------------------------------
def _homology_dims(res, N, side):
    """dim H^m / H_m of the collapsed complete complex for -N <= m <= N."""
    p = res.p
    mats, rks = {}, {}
    for n in res.degrees():
        M = res.cochain_matrix(n) if side == "coh" else res.chain_matrix(n)
        mats[n] = M
        rks[n] = int(rank_mod_p(M, p)) if (M is not None and M.size) else 0
    dim_of = res.cochain_dim if side == "coh" else res.chain_dim
    out = {}
    for mdeg in range(-N, N + 1):
        # cohomology: delta^{m-1} comes from d_m, delta^m from d_{m+1}; the collapsed
        # complex is the same shape on both sides, so one formula serves both.
        out[mdeg] = int(dim_of(mdeg) - rks.get(mdeg, 0) - rks.get(mdeg + 1, 0))
    return out


def tate_cohomology_dims(A, N, p, max_term_dim=20000, verify_positive=True):
    """{m: dim HHhat^m(A)} for -N <= m <= N, over GF(p), via the complete resolution.

    `verify_positive` runs the bookkeeping drift gate described in
    `_assert_positive_agreement` -- read its docstring for what that check does and
    does NOT establish."""
    res = complete_resolution(A, N + 1, p, max_term_dim=max_term_dim)
    dims = _homology_dims(res, N, "coh")
    if verify_positive:
        _assert_positive_agreement(A, dims, N, p, max_term_dim)
    return dims


def tate_homology_dims(A, N, p, max_term_dim=20000):
    """{m: dim HHhat_m(A)} for -N <= m <= N, over GF(p), via the complete resolution."""
    res = complete_resolution(A, N + 1, p, max_term_dim=max_term_dim)
    return _homology_dims(res, N, "hom")


def _assert_positive_agreement(A, dims, N, p, max_term_dim):
    """Drift gate: HHhat^m == HH^m for m >= 1 against the minimal A^e engine.

    HONEST SCOPE -- this is a BOOKKEEPING gate, not an independent oracle.  In degrees
    m >= 1 the Tate cochain complex is literally the minimal resolution's cochain
    complex, so agreement is expected by construction; what the gate catches is a
    window/offset error in the negative indexing leaking into the positive half.  The
    genuine cross-engine oracle is the comparison against the shipped bar/fast HH
    engine, which lives in the test battery."""
    from quiverlab.engine.resolutions_minimal import minimal_cohomology_dims
    labels = _vertex_labels(A)
    Aq = A if len(labels) > 1 else A.unit_adapted()
    ref = minimal_cohomology_dims(to_engine(Aq), N, primes=(p,),
                                  max_term_dim=max_term_dim)[p]
    for mdeg in range(1, min(N, len(ref) - 1) + 1):
        if dims[mdeg] != ref[mdeg]:
            raise QuiverlabError(
                f"Tate-Hochschild drift: HHhat^{mdeg} = {dims[mdeg]} but the minimal "
                f"A^e resolution gives HH^{mdeg} = {ref[mdeg]}, while the "
                f"Bergh-Jorgensen threshold makes them equal for degree >= 1")
