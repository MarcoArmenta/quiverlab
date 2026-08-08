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
    cls = _classify_char0(c, m, dom)
    return HH1Lie(status=cls.pop("status"), note=cls.pop("note"),
                  base_change_note=_base_change_note(dom), **cls, **common)


# ---------------------------------------------------------------------------
# characteristic-0 classification: radical, Levi, sl2-count, toral rank
# ---------------------------------------------------------------------------
def _identity_vecs(m, dom):
    return [[dom.one() if u == v else dom.zero() for u in range(m)] for v in range(m)]


def _killing(c, m, dom):
    """The Killing form ``kappa(i, j) = tr(ad_i ad_j) = sum_{a,b} c[i][b][a] c[j][a][b]``."""
    K = [[dom.zero()] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            s = dom.zero()
            for a in range(m):
                for b in range(m):
                    cib = c[i][b][a]
                    if not dom.is_zero(cib):
                        s = dom.add(s, dom.mul(cib, c[j][a][b]))
            K[i][j] = s
    return K


def _quotient_by_ideal(c, m, ideal_basis, dom):
    """Structure constants ``cS`` (and dim ``mS``) of ``L / I`` for the ideal ``I``
    spanned by ``ideal_basis``: a vector-space complement of ``I`` carries the
    induced bracket ``[s_i, s_j] mod I``."""
    ideal = _span_basis(ideal_basis, dom)
    nrad = len(ideal)
    cur = list(ideal)
    reps = []
    for v in _identity_vecs(m, dom):
        test = _span_basis(cur + [v], dom)
        if len(test) > len(cur):
            reps.append(v)
            cur = test
    mS = len(reps)
    basis = ideal + reps
    cS = [[[dom.zero()] * mS for _ in range(mS)] for _ in range(mS)]
    for i in range(mS):
        for j in range(mS):
            w = _bracket_coords(reps[i], reps[j], c, m, dom)
            coo = _coords_in(basis, w, dom)
            for k in range(mS):
                cS[i][j][k] = coo[nrad + k]
    return cS, mS


def _num_simple_ideals(cS, mS, dom):
    """The number of simple ideals of a SEMISIMPLE Lie algebra ``S`` = the dimension
    of the space of invariant symmetric bilinear forms on ``S`` (each simple ideal
    contributes exactly the multiples of its own Killing form; cross terms vanish by
    invariance + perfectness). A clean exact linear-algebra count (de Graaf)."""
    if mS == 0:
        return 0
    idx, cnt = {}, 0
    for p in range(mS):
        for q in range(p, mS):
            idx[(p, q)] = cnt
            cnt += 1

    def vi(p, q):
        return idx[(p, q)] if p <= q else idx[(q, p)]

    rows = []
    for a in range(mS):
        for b in range(mS):
            for cc in range(mS):
                row = [dom.zero()] * cnt
                for k in range(mS):
                    # invariance: sum_k cS[a][b][k] B[k][cc] + sum_k cS[a][cc][k] B[b][k] = 0
                    row[vi(k, cc)] = dom.add(row[vi(k, cc)], cS[a][b][k])
                    row[vi(b, k)] = dom.add(row[vi(b, k)], cS[a][cc][k])
                if any(not dom.is_zero(x) for x in row):
                    rows.append(row)
    return len(nullspace(rows, dom)) if rows else cnt


def _classify_char0(c, m, dom):
    """The char-0 classification dict (radical / Levi / sl2-count / toral rank /
    type), on the Killing form and de Graaf's algorithms."""
    K = _killing(c, m, dom)
    LL = _span_basis([c[i][j][:] for i in range(m) for j in range(m)], dom)
    # rad(L) = [L, L]^perp = { x : kappa(x, y) = 0 for all y in [L, L] } (de Graaf)
    rows = []
    for y in LL:
        row = [dom.zero()] * m
        for i in range(m):
            s = dom.zero()
            for jj in range(m):
                if not dom.is_zero(y[jj]):
                    s = dom.add(s, dom.mul(y[jj], K[i][jj]))
            row[i] = s
        rows.append(row)
    rad_basis = nullspace(rows, dom) if rows else _identity_vecs(m, dom)
    radical_dim = len(rad_basis)
    levi_dim = m - radical_dim
    semisimple = (radical_dim == 0)
    # the semisimple Levi factor S = L / rad(L)
    cS, mS = _quotient_by_ideal(c, m, rad_basis, dom)
    if mS and rank(_killing(cS, mS, dom), dom) != mS:
        raise QuiverlabError(
            "HH^1-Lie char-0 classification: L/rad(L) is not semisimple (its Killing "
            "form is degenerate) -- the radical computation is inconsistent")
    r = _num_simple_ideals(cS, mS, dom)
    note = ""
    if levi_dim == 0:
        sl2_count, toral_rank, levi_type, simple, status = 0, 0, "0", False, "complete"
    elif levi_dim == 3 * r:                      # every simple ideal is a 3-dim form of sl2
        sl2_count, toral_rank = r, r
        levi_type = "+".join(["A1"] * r)
        simple = (radical_dim == 0 and r == 1)
        status = "complete"
    else:                                        # a simple ideal of dim > 3: type not identified
        sl2_count = toral_rank = levi_type = None
        simple = (radical_dim == 0 and r == 1)
        status = "levi_incomplete"
        note = (f"the Levi factor S (dim {levi_dim}, {r} simple ideal(s)) is not a sum "
                "of sl2's; its type / sl2-count / toral rank are not identified "
                "(honest-scope, Plan 70 Task 3 -- only dim S in {0, 3} is pinned)")
    return dict(radical_dim=radical_dim, levi_dim=levi_dim, sl2_count=sl2_count,
                toral_rank=toral_rank, semisimple=semisimple, simple=simple,
                levi_type=levi_type, char0_note=None, status=status, note=note)


def _base_change_note(dom):
    """The base-change provenance, attached over EVERY char-0 quiverlab domain
    (gated on the exact ARITHMETIC field, not the formal ``is_algebraically_closed``
    flag). Every quiverlab char-0 domain -- including the DEFAULT ``CC`` -- computes
    in exact QQ, so an anisotropic-over-QQ form of sl2 can be undercounted."""
    field = dom.name
    if bool(dom.is_algebraically_closed):        # the DEFAULT CC: formally closed, exact QQ
        return (f"sl2-count / toral rank count SPLIT sl2 forms over the exact arithmetic "
                f"field ({field} is formally algebraically closed but computes in exact "
                f"QQ arithmetic); base change to the algebraic closure may increase the "
                f"count -- an anisotropic 3-dim simple factor splits only over k-bar "
                f"(no fabricated 'algebraically closed' claim on the count)")
    return (f"sl2-count / toral rank count SPLIT sl2 forms over {field}; base change to "
            f"the algebraic closure may increase the count -- an anisotropic 3-dim simple "
            f"factor (a form of sl2) splits only over k-bar")


# ---------------------------------------------------------------------------
# thin field-general verdicts
# ---------------------------------------------------------------------------
def is_solvable_hh1(A, *, budget=DEFAULT_MAXDIM):
    """Is ``HH^1(A)`` a solvable Lie algebra? Any exact field."""
    return _field_general(A, budget)["solvable"]


def is_nilpotent_hh1(A, *, budget=DEFAULT_MAXDIM):
    """Is ``HH^1(A)`` a nilpotent Lie algebra? Any exact field."""
    return _field_general(A, budget)["nilpotent"]


# ---------------------------------------------------------------------------
# Task 4: the RSS Ext-quiver criterion + the Gerstenhaber cross-engine oracle
# ---------------------------------------------------------------------------
def rss_solvable_certificate(A):
    """The RSS 1903.12145 Ext-quiver solvability criterion (arbitrary characteristic):
    if the Gabriel quiver has NO loops (``Ext^1(S,S) = 0``) and NO parallel arrows
    (``dim Ext^1(S,T) <= 1``) then ``HH^1(A)`` is SOLVABLE. A pure quiver predicate.

    Returns ``{solvable_by_criterion, has_loops, has_parallel, note}``. When the
    criterion is silent (loops or parallel arrows present) ``solvable_by_criterion``
    is ``None`` (the criterion is one-directional -- its failure says nothing).
    Loud refusal on a structure-constant-only algebra (no Gabriel quiver)."""
    if A.quiver is None:
        raise QuiverlabError(
            "rss_solvable_certificate needs the Gabriel quiver -- a structure-constant-only "
            "algebra has no Ext-quiver to test",
            hint="build the algebra from a quiver (kQ/I)")
    # A.quiver IS the Gabriel quiver for an admissible presentation (its arrows are the
    # Ext^1(S_i, S_j) basis); read loops / parallel arrows off it directly -- a pure
    # quiver predicate over ANY characteristic (unlike gabriel_quiver()'s char-sensitive
    # basic-ization).
    gq = A.quiver
    pair_count, has_loops, has_parallel = {}, False, False
    for name in gq.arrows:
        s, t = gq.source(name), gq.target(name)
        if s == t:
            has_loops = True
        else:
            pair_count[(s, t)] = pair_count.get((s, t), 0) + 1
            if pair_count[(s, t)] >= 2:
                has_parallel = True
    applies = (not has_loops) and (not has_parallel)
    if applies:
        note = ("RSS 1903.12145: no loops and no parallel arrows in the Gabriel quiver "
                "(Ext^1(S,S) = 0, dim Ext^1(S,T) <= 1) => HH^1 is solvable, in any "
                "characteristic")
    else:
        reasons = ([("loops (Ext^1(S,S) != 0)")] if has_loops else []) + \
                  ([("parallel arrows (dim Ext^1(S,T) > 1)")] if has_parallel else [])
        note = ("the RSS no-loops/no-parallel criterion does not apply (" +
                " and ".join(reasons) + "); it is silent -- HH^1 may or may not be solvable")
    return {"solvable_by_criterion": True if applies else None,
            "has_loops": has_loops, "has_parallel": has_parallel, "note": note}


def _lie_invariants_from_products(g):
    """Read the HH^1 Lie algebra off a Plan-35 ``HHProducts`` degree-(1,1) bracket
    table and derive the BASIS-INDEPENDENT invariants (dim, solvable, nilpotent, the
    derived-series dims). The Gerstenhaber bracket route is GF(p) in-window (Plan 35 /
    P51); the field is recovered from the products' basis provenance (``bar/GF(p)``).
    Structure constants are NOT compared across engines (basis-dependent, Plan-35 rule)."""
    import re

    from quiverlab.fields import GF

    t = g.tables.get((1, 1))
    if t is None:
        raise QuiverlabError(
            "HH^1-Lie cross-engine: the products object has no degree-(1,1) bracket table")
    match = re.search(r"GF\((\d+)\)", str(g.basis) or "")
    if match is None:
        raise QuiverlabError(
            f"HH^1-Lie cross-engine: cannot read the field from basis {g.basis!r} "
            f"(the Gerstenhaber bracket route is GF(p))")
    dom = GF(int(match.group(1)))
    consts, m = t.constants, t.dims[0]
    # c[i][j][k] = coeff of basis element k in [b_i, b_j] (ProductTable: constants[k][i][j])
    c = [[[dom.coerce(int(consts[k][i][j])) for k in range(m)] for j in range(m)]
         for i in range(m)]
    ds = _series(c, m, dom, central=False)
    lcs = _series(c, m, dom, central=True)
    return dict(dim=m, solvable=ds[-1] == 0, nilpotent=lcs[-1] == 0,
                derived_series_dims=ds, lower_central_dims=lcs)
