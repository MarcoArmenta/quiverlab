"""HH^1 as a Lie algebra = Der(A)/Inn(A) with the commutator bracket (Plan 70 / R11).

The first Hochschild cohomology of a finite-dimensional algebra is the
outer-derivation Lie algebra ``HH^1(A) = Der(A)/Inn(A)`` with the bracket
``[D, E] = D o E - E o D`` (Gerstenhaber 1963, degree 1). This module computes it
**field-generally over any exact Domain** from the algebra's own structure-constant
tensor ``A.T`` -- the Leibniz null space + the matrix commutator -- entirely
independent of the window-bounded Gerstenhaber-bracket engine (Plan 35 / P51), which
enters only as a cross-engine oracle.

The **field-general** invariants (dimension, the bracket structure constants in the
Der/Inn representative basis, the derived and lower-central series, and the
solvable / nilpotent / abelian / perfect verdicts) are exact linear algebra in ANY
characteristic. The **characteristic-0** classification (solvable radical, Levi
decomposition, sl2-count, toral rank) rests on the Killing form and Cartan's /
Levi's / Weyl's theorems -- theorems of characteristic 0 with NO positive-char
analogue -- and is a loud ``QuiverlabError`` (or ``None`` + a ``char0_note``) over
``char > 0``. See ``docs/plans/2026-08-08-plan-70-hh1-lie.md``.

References: RSS 1903.12145 (the Ext-quiver solvability criterion), Strametz 2006
(monomial HH^1-Lie), Eisele-Raedschelders 1903.07380 (the sl2-count from Kronecker
subquivers), Gerstenhaber 1963 (degree-1 bracket = commutator), de Graaf 2000
(radical / Levi / simple-ideal algorithms)."""
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import nullspace, rank, rref, solve

#: The hard cap on ``A.dim`` for the Der solve. The Leibniz system is ``d^2``
#: unknowns / ``d^3`` equations; the exact QQ solve is empirically ~ d^5.4
#: (measured: dim 16 = 4.3 s, dim 20 = 13 s, dim 24 = 36 s; GF(p) markedly faster).
DEFAULT_MAXDIM = 48

_REFERENCES = ["rss_hh1_lie", "eisele_raedschelders", "strametz_hh1_lie",
               "gerstenhaber1963", "assem_book"]


# ---------------------------------------------------------------------------
# helpers over an exact Domain
# ---------------------------------------------------------------------------
def _span_basis(vecs, dom):
    """A reduced-row-echelon basis of the span of ``vecs`` (rows)."""
    if not vecs:
        return []
    R, piv = rref(vecs, dom)
    return [R[r] for r in range(len(piv))]


def _rank(vecs, dom):
    return rank(vecs, dom) if vecs else 0


def _matmul_commutator_flat(u, v, d, dom):
    """The flattened commutator ``[U, V] = U V - V U`` of two ``d x d`` matrices
    given in row-major flattened form (``x[a*d + b] = X[a][b]``)."""
    W = [dom.zero()] * (d * d)
    for a in range(d):
        ad = a * d
        for c in range(d):
            s = dom.zero()
            for b in range(d):
                uv = dom.mul(u[ad + b], v[b * d + c])
                vu = dom.mul(v[ad + b], u[b * d + c])
                s = dom.add(s, dom.sub(uv, vu))
            W[ad + c] = s
    return W


def _guard_budget(A, budget, what):
    d = A.dim
    if d > budget:
        raise QuiverlabError(
            f"{what}: dim A = {d} exceeds the budget {budget} -- the derivation "
            f"solve is {d * d} unknowns in {d * d * d} equations (~ d^5.4 over QQ)",
            hint=f"raise the budget above {d} only if you can afford the exact solve, "
                 f"or restrict to a subalgebra")


# ---------------------------------------------------------------------------
# Der(A), Inn(A)
# ---------------------------------------------------------------------------
def derivations(A, *, budget=DEFAULT_MAXDIM):
    """A basis of ``Der(A)`` as flattened ``d x d`` matrices over ``A.domain``
    (row-major, column ``j`` = coordinates of ``D(e_j)``): the exact null space of
    the Leibniz system, for all ``i, j, k``,

        sum_m T[i][j][m] D[k][m] - sum_p D[p][i] T[p][j][k] - sum_q D[q][j] T[i][q][k] = 0.

    The unit forces ``D(1) = 0`` automatically."""
    _guard_budget(A, budget, "derivations")
    dom, d, T = A.domain, A.dim, A.T
    rows = []
    for i in range(d):
        Ti = T[i]
        for j in range(d):
            for k in range(d):
                row = [dom.zero()] * (d * d)
                for m in range(d):                       # +T[i][j][m] * D[k][m]
                    row[k * d + m] = dom.add(row[k * d + m], Ti[j][m])
                for p in range(d):                       # -T[p][j][k] * D[p][i]
                    row[p * d + i] = dom.sub(row[p * d + i], T[p][j][k])
                for q in range(d):                       # -T[i][q][k] * D[q][j]
                    row[q * d + j] = dom.sub(row[q * d + j], Ti[q][k])
                if any(not dom.is_zero(c) for c in row):
                    rows.append(row)
    if not rows:                                          # no constraints: every matrix
        return [[dom.one() if u == v else dom.zero() for u in range(d * d)]
                for v in range(d * d)]
    return nullspace(rows, dom)


def inner_derivations(A, *, budget=DEFAULT_MAXDIM):
    """A basis of ``Inn(A) = {ad_x : x in A}`` as flattened ``d x d`` matrices.
    ``ad_{e_i}`` has column ``j = (T[i][j][a] - T[j][i][a])_a``; ``Inn`` is an ideal
    of ``Der`` (``[D, ad_x] = ad_{D(x)}``)."""
    _guard_budget(A, budget, "inner_derivations")
    dom, d, T = A.domain, A.dim, A.T
    vecs = []
    for i in range(d):
        v = [dom.zero()] * (d * d)
        for b in range(d):
            for a in range(d):
                v[a * d + b] = dom.sub(T[i][b][a], T[b][i][a])
        if any(not dom.is_zero(c) for c in v):
            vecs.append(v)
    return _span_basis(vecs, dom)


# ---------------------------------------------------------------------------
# HH^1 representatives + bracket structure constants (der_inn basis)
# ---------------------------------------------------------------------------
def _hh1_reps(Der, Inn, dom):
    """Extend the ``Inn`` basis to a basis of ``Der``; the added ``Der`` vectors
    represent ``HH^1 = Der/Inn`` (their images form a basis of the quotient)."""
    cur = _span_basis(Inn, dom)
    reps = []
    for v in Der:
        test = _span_basis(cur + [v], dom)
        if len(test) > len(cur):
            reps.append(v)
            cur = test
    return reps


def _coords_in(basis_vecs, w, dom):
    """coords ``c`` with ``sum_t c[t] basis_vecs[t] = w`` (basis independent, w in span)."""
    if not basis_vecs:
        return []
    n = len(w)
    G = [[basis_vecs[t][r] for t in range(len(basis_vecs))] for r in range(n)]
    c = solve(G, list(w), dom)
    if c is None:
        raise QuiverlabError(
            "HH^1-Lie: a bracket left the span of Der -- Inn is not an ideal "
            "(bug or an inconsistent structure-constant tensor)")
    return c


def _bracket_constants(reps, Inn, dom, d):
    """The bracket structure constants ``c[i][j][k]`` of ``HH^1`` in the ``reps``
    (Der/Inn) basis: ``[reps_i, reps_j] = sum_k c[i][j][k] reps_k`` mod ``Inn``."""
    m = len(reps)
    inn = _span_basis(Inn, dom)
    basis = inn + reps                        # Inn first, reps last
    c = [[[dom.zero()] * m for _ in range(m)] for _ in range(m)]
    for i in range(m):
        for j in range(m):
            w = _matmul_commutator_flat(reps[i], reps[j], d, dom)
            coo = _coords_in(basis, w, dom)
            rep_coords = coo[len(inn):]
            for k in range(m):
                c[i][j][k] = rep_coords[k]
    return c


# ---------------------------------------------------------------------------
# abstract Lie analysis on structure constants c[i][j][k] (m basis elements)
# ---------------------------------------------------------------------------
def _bracket_coords(x, y, c, m, dom):
    """``[x, y]`` in coordinates for ``x, y`` m-vectors, given constants ``c``."""
    out = [dom.zero()] * m
    for p in range(m):
        if dom.is_zero(x[p]):
            continue
        for q in range(m):
            if dom.is_zero(y[q]):
                continue
            f = dom.mul(x[p], y[q])
            cpq = c[p][q]
            for k in range(m):
                if not dom.is_zero(cpq[k]):
                    out[k] = dom.add(out[k], dom.mul(f, cpq[k]))
    return out


def _series(c, m, dom, central):
    """The derived series (``central=False``: L^{(k+1)} = [L^{(k)}, L^{(k)}]) or the
    lower-central series (``central=True``: L^{[k+1]} = [L, L^{[k]}]) dimension
    profile, terminating at 0 (reached the trivial subalgebra) or on stabilization."""
    ident = [[dom.one() if u == v else dom.zero() for u in range(m)] for v in range(m)]
    Vk = ident
    dims = [m]
    while True:
        left = ident if central else Vk
        brs = [_bracket_coords(a, b, c, m, dom) for a in left for b in Vk]
        Vn = _span_basis(brs, dom)
        dn = len(Vn)
        dims.append(dn)
        if dn == 0 or dn == dims[-2]:
            break
        Vk = Vn
    return dims


def _jacobi_holds(c, m, dom):
    """[[x,y],z] + [[y,z],x] + [[z,x],y] = 0 on the basis (self-cert)."""
    e = [[dom.one() if u == v else dom.zero() for u in range(m)] for v in range(m)]
    for i in range(m):
        for j in range(m):
            for k in range(m):
                t1 = _bracket_coords(_bracket_coords(e[i], e[j], c, m, dom), e[k], c, m, dom)
                t2 = _bracket_coords(_bracket_coords(e[j], e[k], c, m, dom), e[i], c, m, dom)
                t3 = _bracket_coords(_bracket_coords(e[k], e[i], c, m, dom), e[j], c, m, dom)
                for a in range(m):
                    if not dom.is_zero(dom.add(dom.add(t1[a], t2[a]), t3[a])):
                        return False
    return True


# ---------------------------------------------------------------------------
# the report dataclass
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HH1Lie:
    """A frozen report of the first Hochschild cohomology as a Lie algebra."""
    dim: int
    dim_der: int
    dim_inn: int
    basis: str
    constants: list
    abelian: bool
    perfect: bool
    solvable: bool
    nilpotent: bool
    derived_series_dims: list
    lower_central_dims: list
    characteristic: int
    algebraically_closed: bool
    radical_dim: object = None
    levi_dim: object = None
    sl2_count: object = None
    toral_rank: object = None
    semisimple: object = None
    simple: object = None
    levi_type: object = None
    base_change_note: object = None
    char0_note: object = None
    status: str = "complete"
    note: str = ""
    references: tuple = ()


def _field_general(A, budget):
    """(dim, dim_der, dim_inn, constants c, series, verdicts) over any field."""
    dom, d = A.domain, A.dim
    Der = derivations(A, budget=budget)
    Inn = inner_derivations(A, budget=budget)
    dim_der, dim_inn = len(Der), len(Inn)
    # self-cert: Inn subset Der (adding any ad to Der must not raise the rank)
    if _rank(Der + Inn, dom) != dim_der:
        raise QuiverlabError("HH^1-Lie self-cert failed: Inn is not contained in Der")
    reps = _hh1_reps(Der, Inn, dom)
    m = len(reps)
    c = _bracket_constants(reps, Inn, dom, d)
    ds = _series(c, m, dom, central=False)
    lcs = _series(c, m, dom, central=True)
    solvable = ds[-1] == 0
    nilpotent = lcs[-1] == 0
    abelian = len(ds) >= 2 and ds[1] == 0
    perfect = len(ds) >= 2 and ds[1] == ds[0]
    return dict(dim=m, dim_der=dim_der, dim_inn=dim_inn, c=c, reps=reps,
                derived_series_dims=ds, lower_central_dims=lcs,
                solvable=solvable, nilpotent=nilpotent, abelian=abelian, perfect=perfect)


def hh1_lie_structure(A, *, budget=DEFAULT_MAXDIM, require_char0=False):
    """The full HH^1-Lie report (:class:`HH1Lie`).

    The field-general block (dim, bracket structure constants in the ``der_inn``
    basis, derived / lower-central series, solvable / nilpotent / abelian / perfect)
    is computed over any exact Domain. The char-0 classification (radical, Levi,
    sl2-count, toral rank) is filled only when ``A.domain.characteristic == 0``;
    over ``char > 0`` those fields are ``None`` with a ``char0_note`` (and
    ``require_char0=True`` raises loudly)."""
    dom = A.domain
    fg = _field_general(A, budget)
    m, c = fg["dim"], fg["c"]
    consts = [[[dom.to_str(c[i][j][k]) for k in range(m)] for j in range(m)]
              for i in range(m)]
    common = dict(
        dim=m, dim_der=fg["dim_der"], dim_inn=fg["dim_inn"], basis="der_inn",
        constants=consts, abelian=fg["abelian"], perfect=fg["perfect"],
        solvable=fg["solvable"], nilpotent=fg["nilpotent"],
        derived_series_dims=fg["derived_series_dims"],
        lower_central_dims=fg["lower_central_dims"],
        characteristic=dom.characteristic,
        algebraically_closed=bool(dom.is_algebraically_closed),
        references=tuple(_REFERENCES))
    if dom.characteristic != 0:
        if require_char0:
            raise QuiverlabError(
                "Levi/radical/sl2-count/toral-rank need characteristic 0: "
                "Cartan/Levi/Weyl have NO positive-characteristic analogue, so the "
                "classification is unjustified for any char > 0 (regardless of dim)",
                hint="solvable/nilpotent are returned in every characteristic")
        note = (
            "radical / Levi / sl2-count / toral-rank need characteristic 0 "
            "(Cartan/Levi/Weyl have no positive-characteristic analogue); "
            "solvable / nilpotent are still computed. In char p, "
            "HH^1(k[x]/(x^p)) = W_1 (Jacobson-Witt), on which the Killing form is "
            "degenerate and Levi's theorem has no content.")
        return HH1Lie(char0_note=note, status="complete", **common)
    # char 0: the classification (radical / Levi / sl2-count) is filled by Task 3;
    # until then the fields default to None (green intermediate slice).
    return HH1Lie(status="complete", **common)


# ---------------------------------------------------------------------------
# thin field-general verdicts
# ---------------------------------------------------------------------------
def is_solvable_hh1(A, *, budget=DEFAULT_MAXDIM):
    """Is ``HH^1(A)`` a solvable Lie algebra? Any exact field."""
    return _field_general(A, budget)["solvable"]


def is_nilpotent_hh1(A, *, budget=DEFAULT_MAXDIM):
    """Is ``HH^1(A)`` a nilpotent Lie algebra? Any exact field."""
    return _field_general(A, budget)["nilpotent"]
