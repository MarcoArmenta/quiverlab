"""Module Ext^n via a minimal projective resolution + Hom, and global dimension
(spec section 3.5-3.6). Ext^n_A(M,N) = H^n(Hom_A(Q_*, N)); gl.dim A = sup_v pd(S_v)."""
from dataclasses import dataclass

from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import hom_space
from quiverlab.modules.resolution import minimal_resolution


def _delta_matrix(Hn, Hn1, dn1, dom):
    """Matrix of delta^n: Hom(Q_n, N) -> Hom(Q_{n+1}, N), phi |-> phi @ d_{n+1}, in the
    given Hom-basis coordinates. Column j = coords of (Hn[j] @ dn1) in the Hn1 basis."""
    if not Hn or not Hn1:
        return lm.zeros(len(Hn1), len(Hn), dom)
    # flatten each Hom basis matrix (dn x dq) into a vector; build the coordinate solve
    def flat(mat):
        return [x for row in mat for x in row]
    basisH1 = lm.cols_to_matrix([flat(h) for h in Hn1])
    cols = []
    for phi in Hn:
        comp = lm.matmul(phi, dn1, dom)                # phi @ d_{n+1}: (dn x dq_{n+1})
        x = lm.solve_columns(basisH1, lm.cols_to_matrix([flat(comp)]), dom)[0]
        cols.append(x)
    return lm.cols_to_matrix(cols)


def ext_dims(A, M, N, top):
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1)
    Qs = [t.module for t in terms]
    Homs = [hom_space(Q, N) if (Q is not None and Q.dim) else [] for Q in Qs]
    # delta^n uses d_{n+1} = dmats[n+1]
    deltas = []
    for n in range(len(Qs) - 1):
        dn1 = dmats[n + 1]
        deltas.append(_delta_matrix(Homs[n], Homs[n + 1], dn1, dom)
                      if (dn1 and dn1[0]) else lm.zeros(len(Homs[n + 1]), len(Homs[n]), dom))
    out = []
    for n in range(top + 1):
        cn = len(Homs[n]) if n < len(Homs) else 0
        r_n = lm.mat_rank(deltas[n], dom) if n < len(deltas) and deltas[n] and deltas[n][0] else 0
        r_nm1 = lm.mat_rank(deltas[n - 1], dom) if (n - 1) >= 0 and (n - 1) < len(deltas) \
            and deltas[n - 1] and deltas[n - 1][0] else 0
        out.append(cn - r_n - r_nm1)
    return out


def ext(A, M, N, n):
    return ext_dims(A, M, N, n)[n]


@dataclass
class GlobalDimension:
    value: int
    exact: bool

    def __int__(self):
        return self.value

    def __eq__(self, other):
        if isinstance(other, GlobalDimension):
            return (self.value, self.exact) == (other.value, other.exact)
        if isinstance(other, int):
            return self.exact and self.value == other
        return NotImplemented

    def __repr__(self):
        if self.exact:
            return f"gl.dim = {self.value} (exact)"
        return (f">= {self.value} (certified lower bound; not resolved within depth "
                f"{self.value})")


def global_dimension(A, bound=32):
    """sup over simples of projective dimension. Exact if every simple resolves within
    `bound`; otherwise a certified lower bound (some simple has pd >= bound)."""
    best, exact = 0, True
    for v in A.quiver.vertices:
        res = A.simple(v).projective_resolution(bound)
        pd = res.pd()
        if pd is None:                     # not resolved within bound -> lower bound
            exact = False
            best = max(best, bound)
        else:
            best = max(best, pd)
    return GlobalDimension(best, exact)


def is_selfinjective(A):
    """A is self-injective iff every indecomposable projective P_v is injective, iff
    soc(P_v) is SIMPLE for every v and v -> (socle vertex) is a BIJECTION (the Nakayama
    permutation). Exact over any field (uses the socle, not the GF(p) engine); for a
    finite-dimensional algebra self-injective == Frobenius (spec section 3.5)."""
    if A.quiver is None:
        from quiverlab.errors import QuiverlabError
        raise QuiverlabError("is_selfinjective needs the quiver presentation",
                             hint="construct via Quiver.algebra(...)")
    socle_vertex = {}
    for v in A.quiver.vertices:
        soc = A.projective(v).socle()
        dv = soc.dimension_vector()
        support = [w for w, d in dv.items() if d > 0]
        if soc.dim != 1 or len(support) != 1:
            return False                     # socle not simple -> P_v not injective
        socle_vertex[v] = support[0]
    return len(set(socle_vertex.values())) == len(list(A.quiver.vertices))
