"""HH^*(A) as a graded Lie module over HH^1(A) (Plan 71 / R12).

The whole Hochschild cohomology ``HH^*(A) = (+)_n HH^n(A)`` is a graded module over
the Lie algebra ``HH^1(A) = Der(A)/Inn(A)`` (Plan 70), with the action the degree-1
component of the Gerstenhaber bracket ``[.,.]: HH^1 (x) HH^n -> HH^n``. On cochains
this is the **Lie derivative** of a derivation ``D`` (Gerstenhaber 1963, degree 1):

    (L_D f)(a_1, ..., a_n) = D(f(a_1, ..., a_n)) - sum_i f(a_1, ..., D(a_i), ..., a_n).

Because ``p = 1`` collapses every Koszul sign in the shipped circle product
(``engine/tt_calculus.py``), ``L_D`` is byte-identical to the shipped
``gerstenhaber_bracket_cochain(alg, 1, n, D, f)`` -- so the field-general route needs
only ``D`` as a ``k``-linear map ``A -> A`` (P70's ``derivations(A)``) and the
normalized bar cochain complex over ``A.domain`` (``hochschild/bar.py``). No
resolution engine, no bracket window.

- **Any exact field**: the action matrices ``rho_n(D) in End(HH^n)`` in a Der/Inn
  basis of ``HH^1``, the graded module ``(+)_n HH^n``, and the self-certified module
  axiom ``rho_n([D,E]) = [rho_n(D), rho_n(E)]`` + inner-acts-zero
  ``rho_n(ad_x) = 0`` (the action factors through ``HH^1 = Der/Inn``).
- **char 0 (a hard, loud gate, P70's)**: the weight/torus decomposition of each
  ``HH^n`` w.r.t. a maximal torus of ``HH^1`` (Task 3) and the indecomposable
  Lie-module summand decomposition via ``modules/decompose.py`` (Task 4, governed
  INDEPENDENTLY by ``decompose``'s own ``char 0 or char > d_n`` guard).

See ``docs/plans/2026-08-08-plan-71-hh-lie-module.md``.
"""
from dataclasses import dataclass

from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.fields.linalg import nullspace, rank, solve
from quiverlab.hochschild.bar import coboundary_matrix, _cochain_basis
from quiverlab.modules import linalg_mod as lm
from quiverlab.invariants.hh1_lie import (
    DEFAULT_MAXDIM,
    _hh1_reps,
    _matmul_commutator_flat,
    derivations,
    inner_derivations,
)

#: Registry citation keys (reuse the existing ``alsolotar_toupie``; P71 adds
#: ``mnprs_special_biserial`` + ``csss_gentle_tt``; P70 keys reused for the rest).
_REFERENCES = ["alsolotar_toupie", "mnprs_special_biserial", "csss_gentle_tt",
               "gerstenhaber1963", "rss_hh1_lie", "degraaf_lie"]


# ---------------------------------------------------------------------------
# the cohomology quotient of the normalized bar complex (reps + coboundary image)
# ---------------------------------------------------------------------------
def _cols_to_matrix(cols):
    """A list of column vectors -> the row-major matrix whose columns they are."""
    if not cols:
        return []
    L = len(cols[0])
    return [[cols[c][r] for c in range(len(cols))] for r in range(L)]


def _cohomology_quotient(B, n, max_cells):
    """``(reps, image_cols)`` for ``HH^n(B)`` on the normalized bar cochain complex
    over ``B.domain``: ``reps`` are cocycles kept independent modulo ``image_cols`` =
    the columns of the coboundary ``delta^{n-1}`` (the coboundaries), so reps and image
    share one orientation (the Plan-35 ``_generic_homology_quotient`` cohomology twin).
    ``dim HH^n = len(reps)``."""
    dom = B.domain
    Dn, ncols, nrows = coboundary_matrix(B, n, max_cells)
    if nrows and ncols:
        cycles = nullspace(Dn, dom)
    else:                                            # no target rows: every cochain a cocycle
        cycles = [[dom.one() if i == j else dom.zero() for j in range(ncols)]
                  for i in range(ncols)]
    if n == 0:
        image = []
    else:
        Dm, ncols_m, nrows_m = coboundary_matrix(B, n - 1, max_cells)
        image = [[Dm[r][c] for r in range(nrows_m)] for c in range(ncols_m)] \
            if (nrows_m and ncols_m) else []
    reps, base = [], list(image)
    r0 = rank(_cols_to_matrix(base), dom) if base else 0
    for v in cycles:
        rr = rank(_cols_to_matrix(base + [v]), dom)
        if rr > r0:
            reps.append(v)
            base.append(v)
            r0 = rr
    return reps, image


def _class_coords(vec, reps, image_cols, dom):
    """Coordinates of ``vec`` in the class basis ``reps`` modulo ``image_cols``: solve
    ``[image | reps] x = vec`` and read off the reps segment. Loud when ``vec`` is not a
    (co)cycle representative (the descent check). Domain elements (not stringified)."""
    cols = list(image_cols) + list(reps)
    if not cols:
        if any(not dom.is_zero(dom.coerce(v)) for v in vec):
            raise QuiverlabError("HH-Lie-module: L_D landed outside the zero space")
        return []
    Mat = [[dom.coerce(cols[c][r]) for c in range(len(cols))] for r in range(len(vec))]
    x = solve(Mat, [dom.coerce(v) for v in vec], dom)
    if x is None:
        raise QuiverlabError(
            "HH-Lie-module: L_D r did not descend to HH^n (not in the cocycle span) -- "
            "the Lie derivative is not a chain map here (a bug -- report it)")
    return list(x[len(image_cols):])


# ---------------------------------------------------------------------------
# the Lie derivative L_D on the normalized degree-n cochain space
# ---------------------------------------------------------------------------
def _apply_lie_derivative(B, Dflat, n, gvec):
    """``L_D g`` as a coordinate vector over ``_cochain_basis(m, n)`` for the derivation
    ``Dflat`` (flattened ``d x d``, ``Dflat[a*m + b] = <e_a, D(e_b)>``) on ``B`` (unit at
    index 0). The two terms of the Lie derivative:

      term 1 (D o g):  <e_t, D(g(e_J))>              = sum_s Dflat[t*m+s] g[(s, J)]
      term 2 (g o D): -sum_i g(..., proj D(e_{J[i]}), ...)
                     = -sum_i sum_{k=1}^{m-1} Dflat[k*m+J[i]] g[(t, J with pos i -> k)]

    Inserted derivation values are reduced mod ``k.1`` (drop index 0 -- ``k`` starts at
    1) exactly as the shipped normalized circle product does; the output value lives in
    the full ``A`` (``t`` ranges over ``0..m-1``)."""
    dom = B.domain
    m = B.dim
    basis = _cochain_basis(m, n)
    idx = {b: i for i, b in enumerate(basis)}
    out = [dom.zero()] * len(basis)
    for o, (t, J) in enumerate(basis):
        val = dom.zero()
        tm = t * m
        for s in range(m):                            # term 1: D o g
            c = Dflat[tm + s]
            if not dom.is_zero(c):
                val = dom.add(val, dom.mul(c, gvec[idx[(s, J)]]))
        for i in range(n):                            # term 2: - g o D
            ji = J[i]
            for k in range(1, m):
                c = Dflat[k * m + ji]
                if dom.is_zero(c):
                    continue
                Jk = J[:i] + (k,) + J[i + 1:]
                val = dom.sub(val, dom.mul(c, gvec[idx[(t, Jk)]]))
        out[o] = val
    return out


def _rho(B, Dflat, n, reps, image):
    """``rho_n(D)`` as a ``d_n x d_n`` matrix (column ``j`` = class coords of
    ``L_D reps[j]``) over ``B.domain``."""
    dom = B.domain
    d = len(reps)
    cols = [_class_coords(_apply_lie_derivative(B, Dflat, n, r), reps, image, dom)
            for r in reps]
    return [[cols[j][a] for j in range(d)] for a in range(d)]


def _is_zero_matrix(M, dom):
    return all(dom.is_zero(M[i][j]) for i in range(len(M)) for j in range(len(M[0]) if M else 0))


def _mat_commutator(X, Y, dom):
    d = len(X)
    out = [[dom.zero()] * d for _ in range(d)]
    for i in range(d):
        for j in range(d):
            s = dom.zero()
            for k in range(d):
                s = dom.add(s, dom.sub(dom.mul(X[i][k], Y[k][j]), dom.mul(Y[i][k], X[k][j])))
            out[i][j] = s
    return out


def _mat_equal(X, Y, dom):
    return all(dom.is_zero(dom.sub(X[i][j], Y[i][j]))
               for i in range(len(X)) for j in range(len(X[0]) if X else 0))


# ---------------------------------------------------------------------------
# Task 3: char-0 maximal torus of HH^1 + weight decomposition of HH^n
#   (NET-NEW code -- P70 provides NONE of this. Cartan subalgebra (de Graaf
#    Engel/Fitting) -> ad-semisimple part -> radical-toral extension ->
#    rational simultaneous diagonalization; each step self-certified.)
# ---------------------------------------------------------------------------
def _identity_vecs(m, dom):
    return [[dom.one() if i == j else dom.zero() for j in range(m)] for i in range(m)]


def _ad_matrix(c, m, dom, x):
    """``ad_x`` on ``HH^1`` as an ``m x m`` matrix (col ``j`` = ``[x, e_j]``):
    ``(ad_x)[k][j] = sum_i x[i] c[i][j][k]``."""
    M = [[dom.zero()] * m for _ in range(m)]
    for i in range(m):
        xi = x[i]
        if dom.is_zero(xi):
            continue
        ci = c[i]
        for j in range(m):
            cij = ci[j]
            for k in range(m):
                if not dom.is_zero(cij[k]):
                    M[k][j] = dom.add(M[k][j], dom.mul(xi, cij[k]))
    return M


def _combine(basis, coords, dom, d):
    out = [dom.zero()] * d
    for j, cj in enumerate(coords):
        if dom.is_zero(cj):
            continue
        bj = basis[j]
        for t in range(d):
            out[t] = dom.add(out[t], dom.mul(cj, bj[t]))
    return out


def _restrict(M, basis, dom, d):
    """Matrix of ``M`` (``d x d``) restricted to the invariant subspace ``span(basis)``
    (columns = ambient vecs), in ``basis`` coordinates; loud if not invariant."""
    s = len(basis)
    colsT = [[basis[j][r] for j in range(s)] for r in range(d)]
    out_cols = []
    for j in range(s):
        Mb = lm.matvec(M, basis[j], dom)
        coeff = solve(colsT, Mb, dom)
        if coeff is None:
            raise QuiverlabError("HH-Lie-module torus: subspace is not invariant "
                                 "(a bug in the simultaneous diagonalization)")
        out_cols.append(coeff)
    return [[out_cols[j][i] for j in range(s)] for i in range(s)]


def _is_nilpotent_mat(M, dom):
    P = M
    for _ in range(len(M) - 1):
        P = lm.matmul(P, M, dom)
    return all(dom.is_zero(P[i][j]) for i in range(len(P)) for j in range(len(P)))


def _min_poly_factors(M, dom):
    """Factor the minimal polynomial of ``M`` over ``dom`` -> ``[(factor, mult)]``
    (ascending-coeff monic factors). Reuses the decompose char-0/GF(p) factoring."""
    from quiverlab.modules.decompose import _factor_min_poly, _min_poly_coeffs
    return _factor_min_poly(_min_poly_coeffs(M, dom), dom)


def _is_ad_semisimple(c, m, dom, x):
    """Is ``ad_x`` semisimple over ``k`` (diagonalizable) -- min poly squarefree AND
    split into linear factors? (``ad_x = 0`` counts: min poly ``x``, one linear factor.)"""
    ad = _ad_matrix(c, m, dom, x)
    if all(dom.is_zero(ad[i][j]) for i in range(m) for j in range(m)):
        return True
    for fac, mult in _min_poly_factors(ad, dom):
        if mult > 1 or len(fac) != 2:            # repeated (nilpotent part) or non-linear
            return False
    return True


# --- de Graaf Cartan subalgebra via the Engel / Fitting-null recursion ---
def _subalgebra_ad(c, m, dom, x, sub):
    """``ad_x`` restricted to the subalgebra ``span(sub)`` (``sub`` = ambient coord vecs),
    as a matrix in ``sub`` coordinates; loud if ``sub`` is not ``ad_x``-stable."""
    s = len(sub)
    subT = [[sub[j][r] for j in range(s)] for r in range(m)]
    cols = []
    for j in range(s):
        w = [dom.zero()] * m                     # [x, sub[j]] in ambient coords
        bj = sub[j]
        for a in range(m):
            xa = x[a]
            if dom.is_zero(xa):
                continue
            ca = c[a]
            for b in range(m):
                bb = bj[b]
                if dom.is_zero(bb):
                    continue
                f = dom.mul(xa, bb)
                cab = ca[b]
                for k in range(m):
                    if not dom.is_zero(cab[k]):
                        w[k] = dom.add(w[k], dom.mul(f, cab[k]))
        coeff = solve(subT, w, dom)
        if coeff is None:
            raise QuiverlabError("HH-Lie-module Cartan: bracket left the subalgebra "
                                 "(a bug in the Engel recursion)")
        cols.append(coeff)
    return [[cols[j][i] for j in range(s)] for i in range(s)]


def _find_nonnilpotent(c, m, dom, sub):
    """An element ``x in span(sub)`` (ambient coords) with ``ad_x|_sub`` NOT nilpotent,
    or ``None`` (``span(sub)`` is nilpotent -> a Cartan). Deterministic ladder: basis
    elements, then small integer combinations (an exact field is infinite / large)."""
    s = len(sub)
    ladders = list(_identity_vecs(s, dom))
    for a in range(s):                           # pairwise sums
        for b in range(a + 1, s):
            ladders.append([dom.add(dom.one() if t == a else dom.zero(),
                                    dom.one() if t == b else dom.zero()) for t in range(s)])
    ladders.append([dom.coerce(t + 1) for t in range(s)])   # coefficient ladder
    for coeff in ladders:
        x = _combine(sub, coeff, dom, m)
        if not _is_nilpotent_mat(_subalgebra_ad(c, m, dom, x, sub), dom):
            return x
    return None


def _fitting_null(c, m, dom, x, sub):
    """Generalized 0-eigenspace of ``ad_x`` within ``span(sub)`` (a subalgebra
    containing a Cartan), returned as ambient coord vectors."""
    s = len(sub)
    adsub = _subalgebra_ad(c, m, dom, x, sub)
    P = adsub
    for _ in range(s - 1):
        P = lm.matmul(P, adsub, dom)
    return [_combine(sub, v, dom, m) for v in nullspace(P, dom)]


def _cartan_subalgebra(c, m, dom):
    """A Cartan subalgebra of ``HH^1`` (ambient coord basis) by the Engel/Fitting
    recursion. Self-cert: the result is nilpotent (no non-nilpotent element) -- the
    loop exits exactly then."""
    if m == 0:
        return []
    sub = _identity_vecs(m, dom)
    while True:
        x = _find_nonnilpotent(c, m, dom, sub)
        if x is None:
            return sub
        sub = _fitting_null(c, m, dom, x, sub)


def _commutes(c, m, dom, x, y):
    for k in range(m):
        s = dom.zero()
        for i in range(m):
            if dom.is_zero(x[i]):
                continue
            for j in range(m):
                if not dom.is_zero(y[j]) and not dom.is_zero(c[i][j][k]):
                    s = dom.add(s, dom.mul(dom.mul(x[i], y[j]), c[i][j][k]))
        if not dom.is_zero(s):
            return False
    return True


def _independent(vecs, dom):
    return rank(vecs, dom) == len(vecs)


def _maximal_torus(c, m, dom):
    """A maximal ad-diagonalizable abelian subalgebra ``t`` of ``HH^1`` (ambient coord
    basis) and a provenance string. Route: the ad-semisimple part of a Cartan
    subalgebra, then a greedy radical-toral extension by ad-semisimple ``HH^1`` basis
    elements commuting with the current torus (the ``k[x]/x^n`` grading ``x d`` -- a
    radical toral element P70 scoped OUT -- is captured here). Self-cert: ``t`` is
    abelian and every generator is ad-semisimple."""
    if m == 0:
        return [], "trivial (HH^1 = 0)"
    H = _cartan_subalgebra(c, m, dom)
    torus = []
    for h in H:                                  # ad-semisimple part of the Cartan
        if _is_ad_semisimple(c, m, dom, h) and all(_commutes(c, m, dom, h, t) for t in torus) \
                and _independent(torus + [h], dom):
            torus.append(h)
    for e in _identity_vecs(m, dom):             # radical-toral greedy extension
        if _is_ad_semisimple(c, m, dom, e) and all(_commutes(c, m, dom, e, t) for t in torus) \
                and _independent(torus + [e], dom):
            torus.append(e)
    # self-cert: abelian + each ad-semisimple
    for a in range(len(torus)):
        for b in range(len(torus)):
            if not _commutes(c, m, dom, torus[a], torus[b]):
                raise QuiverlabError("HH-Lie-module torus self-cert failed: not abelian")
    prov = (f"maximal torus (dim {len(torus)}) = the ad-semisimple part of a Cartan "
            "subalgebra of HH^1 + radical-toral extension; basis-dependent (the "
            "normalization of each generator is a choice -- weight LABELS are provenance, "
            "the P70 sl2-triple NON-NORMATIVE precedent)")
    return torus, prov


def _rational_eigenspaces(M, dom, s):
    """Full eigenspaces of the ``s x s`` matrix ``M`` over ``k``: ``[(lam, [ker vecs])]``
    for each rational eigenvalue, or ``("anisotropic", None)`` if a factor is non-linear
    (eigenvalue outside ``k``), or ``("nonsemisimple", None)`` if not diagonalizable
    (a repeated min-poly factor -- eigenspace dims would not fill ``s``)."""
    if s == 0:
        return []
    facs = _min_poly_factors(M, dom)
    out = []
    for fac, mult in facs:
        if len(fac) != 2:
            return ("anisotropic", None)
        if mult > 1:
            return ("nonsemisimple", None)
        # root of a1*x + a0 = -a0/a1 (the factor may be non-monic over QQ, e.g. 2x-1)
        lam = dom.neg(dom.mul(fac[0], dom.inv(fac[1])))
        E = [[dom.sub(M[i][j], lam if i == j else dom.zero()) for j in range(s)]
             for i in range(s)]
        out.append((lam, nullspace(E, dom)))
    return out


def _simultaneous_weights(mats, dom, d):
    """Simultaneous eigenspace decomposition of the commuting family ``mats`` on
    ``k^d``: ``{weight-tuple (Domain elts): dim}``, or ``("anisotropic"|"nonsemisimple",
    None)`` on the honest fallbacks. ``d = 0`` -> ``{}``."""
    if d == 0:
        return {}
    spaces = [([], _identity_vecs(d, dom))]      # (weight prefix, ambient basis of subspace)
    for M in mats:
        new = []
        for wt, basis in spaces:
            s = len(basis)
            if s == 0:
                continue
            eig = _rational_eigenspaces(_restrict(M, basis, dom, d), dom, s)
            if isinstance(eig, tuple) and eig[0] in ("anisotropic", "nonsemisimple"):
                return eig
            for lam, kervecs in eig:
                amb = [_combine(basis, v, dom, d) for v in kervecs]
                if amb:
                    new.append((wt + [lam], amb))
        spaces = new
    result = {}
    for wt, basis in spaces:
        key = tuple(wt)
        result[key] = result.get(key, 0) + len(basis)
    return result


# ---------------------------------------------------------------------------
# Task 4: indecomposable Lie-module summands (via decompose_representation) + iso grouping
# ---------------------------------------------------------------------------
def _rep_hom_space(g1, g2, dom, d):
    """Basis of ``{X in M_d : X g1[i] = g2[i] X for all i}`` (intertwiners g1 -> g2)."""
    rows = []
    for a, b in zip(g1, g2):
        for i in range(d):
            for l in range(d):
                row = [dom.zero()] * (d * d)
                for k in range(d):
                    row[i * d + k] = dom.add(row[i * d + k], a[k][l])
                    row[k * d + l] = dom.sub(row[k * d + l], b[i][k])
                if any(not dom.is_zero(x) for x in row):
                    rows.append(row)
    flats = nullspace(rows, dom) if rows else \
        [[dom.one() if t == u else dom.zero() for t in range(d * d)] for u in range(d * d)]
    return [[[v[i * d + j] for j in range(d)] for i in range(d)] for v in flats]


def _rep_isomorphic(s1, s2, dom):
    """Are the two :class:`RepSummand`-like ``(dim, gens)`` isomorphic? Certified by an
    INVERTIBLE intertwiner (simultaneous similarity). Non-iso => no invertible element
    (never a false positive); an undecidable case returns False (conservative -- may
    over-count classes, never merges non-iso ones)."""
    if s1["dim"] != s2["dim"]:
        return False
    d = s1["dim"]
    if d == 0:
        return True
    homs = _rep_hom_space(s1["gens"], s2["gens"], dom, d)
    if not homs:
        return False
    cands = list(homs)                                # singletons
    for i in range(len(homs)):                        # + pairwise sums
        for j in range(i + 1, len(homs)):
            cands.append([[dom.add(homs[i][a][b], homs[j][a][b]) for b in range(d)]
                          for a in range(d)])
    return any(rank(X, dom) == d for X in cands)


def _group_summands(flat, dom):
    """Group a flat list of RepSummand into ``[{"dim", "label", "mult"}]`` by the exact
    representation-iso certificate; distinct entries are pairwise non-isomorphic.
    ``label`` is always ``None`` (never guessed)."""
    groups = []                                       # [ {"dim","gens","mult"} ]
    for s in flat:
        entry = {"dim": s.dim, "gens": [list(map(list, g)) for g in s.gens]}
        for g in groups:
            if _rep_isomorphic(g, entry, dom):
                g["mult"] += 1
                break
        else:
            entry["mult"] = 1
            groups.append(entry)
    return [{"dim": g["dim"], "label": None, "mult": g["mult"]} for g in groups]


def _summand_tables(per_degree, dom, budget):
    """Per-degree indecomposable-summand decomposition, governed INDEPENDENTLY by
    decompose's own char guard (``char 0 or char > d_n``). ``summands[n]`` is a list of
    ``{"dim","label","mult"}`` parts, or ``{"error": note}`` where decompose refuses
    (``char p <= d_n`` with no ``dim End = 1`` certificate) -- surfaced, never a crash
    and never a guessed decomposition."""
    from quiverlab.modules.decompose import decompose_representation
    out = []
    for pd in per_degree:
        dn = len(pd["reps"])
        if dn == 0:
            out.append([])
            continue
        if not pd["rhos"]:                            # HH^1 = 0: trivial action, d_n lines
            out.append([{"dim": 1, "label": None, "mult": dn}])
            continue
        try:
            flat = decompose_representation(pd["rhos"], dom, budget=budget)
            out.append(_group_summands(flat, dom))
        except QuiverlabError as exc:
            out.append({"error": str(exc)})
    return out


def _weight_tables(B, per_degree, torus, reps_hh1, dom):
    """Per-degree weight tables over ``k`` (char 0). Returns ``(tables, base_change_note)``:
    ``tables[n]`` is ``{"n", "torus_rank", "weights": [[ [str(lam)...], dim ]...] }`` where
    weights split rationally, else a per-degree note; ``base_change_note`` is set when an
    anisotropic torus is met (P70 base-change precedent -- no fabricated split)."""
    m = len(reps_hh1)
    tables = []
    base_change_note = None
    for n, pd in enumerate(per_degree):
        reps, image, rhos = pd["reps"], pd["image"], pd["rhos"]
        dn = len(reps)
        entry = {"n": n, "torus_rank": len(torus)}
        if dn == 0 or not torus:
            entry["weights"] = ([[[], dn]] if dn else [])
            tables.append(entry)
            continue
        mats = []                                # rho_n(h) for each torus generator h
        for h in torus:
            R = [[dom.zero()] * dn for _ in range(dn)]
            for i in range(m):
                hi = h[i]
                if dom.is_zero(hi):
                    continue
                Ri = rhos[i]
                for a in range(dn):
                    for b in range(dn):
                        if not dom.is_zero(Ri[a][b]):
                            R[a][b] = dom.add(R[a][b], dom.mul(hi, Ri[a][b]))
            mats.append(R)
        w = _simultaneous_weights(mats, dom, dn)
        if isinstance(w, tuple):                 # ("anisotropic"|"nonsemisimple", None)
            entry["weights"] = None
            if w[0] == "anisotropic":
                entry["note"] = ("a torus eigenvalue lies outside the base field -- "
                                 "weights are unavailable over k (base change to k-bar "
                                 "needed); no fabricated split")
                base_change_note = ("some HH^n carries an anisotropic torus action; its "
                                    "weights split only over the algebraic closure "
                                    "(P70 base-change precedent -- no fabricated split)")
            else:
                entry["note"] = ("a torus generator does not act semisimply on this HH^n "
                                 "(a repeated eigenvalue) -- weights are not read")
        else:
            entry["weights"] = [[[dom.to_str(l) for l in wt], dim] for wt, dim in w.items()]
        tables.append(entry)
    return tables, base_change_note


# ---------------------------------------------------------------------------
# the public primitive + the report
# ---------------------------------------------------------------------------
def lie_derivative_on_hh(A, D, n, *, max_cells=4_000_000):
    """``rho_n(D)``: the action of the derivation ``D`` on ``HH^n(A)`` as a ``d_n x d_n``
    matrix over ``A.domain`` (``d_n = dim HH^n``). ``D`` is a flattened ``d x d``
    derivation matrix in the ``A.unit_adapted()`` basis (``D[a*m + b] = <e_a, D(e_b)>``,
    the ``derivations(A.unit_adapted())`` convention). Field-general: no resolution."""
    B = A.unit_adapted()
    reps, image = _cohomology_quotient(B, n, max_cells)
    if not reps:
        return []
    return _rho(B, D, n, reps, image)


@dataclass(frozen=True)
class HHLieModule:
    """A frozen report of ``HH^*`` as a graded Lie module over ``HH^1`` (Plan 71).
    ``__bool__`` is NOT defined (a data report, like Plan-35's ``HHProducts``)."""
    top: int
    hh_dims: list                 # [dim HH^0, ..., dim HH^top]  (== A.hochschild_cohomology)
    hh1_dim: int                  # dim HH^1 (the acting Lie algebra; P70)
    basis: str                    # provenance: "der_inn/bar"
    action: list                  # per n: {"n", "dim", "gens": [flattened rho_n(D) as str]}
    module_axiom_ok: bool         # rho_n([D,E]) == [rho_n(D), rho_n(E)] on all pairs
    inner_acts_zero: bool         # rho_n(ad_x) == 0 (factors through Der/Inn)
    characteristic: int
    weights: object = None        # per n weight table (char 0 only), else None + char0_note
    torus_provenance: object = None
    summands: object = None       # per n indecomposable summands (decompose's own guard)
    weight_base_change_note: object = None
    char0_note: object = None
    status: str = "complete"
    note: str = ""
    references: tuple = ()


def _field_general_action(B, top, max_cells):
    """The field-general core: (hh_dims, hh1_dim, per-degree quotient + rho matrices,
    module_axiom_ok, inner_acts_zero). Reps of HH^1 (Der/Inn) act; the commutator of
    derivations is P70's ``_matmul_commutator_flat``."""
    dom = B.domain
    d = B.dim
    Der = derivations(B)
    Inn = inner_derivations(B)
    reps_hh1 = _hh1_reps(Der, Inn, dom)
    hh1_dim = len(reps_hh1)
    hh_dims, per_degree, action = [], [], []
    axiom_ok, inner_ok = True, True
    for n in range(top + 1):
        reps, image = _cohomology_quotient(B, n, max_cells)
        dn = len(reps)
        hh_dims.append(dn)
        rhos = [_rho(B, Dm, n, reps, image) for Dm in reps_hh1] if dn else []
        per_degree.append({"reps": reps, "image": image, "rhos": rhos})
        action.append({"n": n, "dim": dn,
                       "gens": [[dom.to_str(R[i][j]) for i in range(dn) for j in range(dn)]
                                for R in rhos]})
        if not dn:
            continue
        for Dx in Inn:                                # inner derivations act as zero
            if not _is_zero_matrix(_rho(B, Dx, n, reps, image), dom):
                inner_ok = False
        for a in range(hh1_dim):                      # module axiom on all rep pairs
            for b in range(hh1_dim):
                comm = _matmul_commutator_flat(reps_hh1[a], reps_hh1[b], d, dom)
                lhs = _rho(B, comm, n, reps, image)
                rhs = _mat_commutator(rhos[a], rhos[b], dom)
                if not _mat_equal(lhs, rhs, dom):
                    axiom_ok = False
    return dict(hh_dims=hh_dims, hh1_dim=hh1_dim, reps_hh1=reps_hh1, Der=Der, Inn=Inn,
                per_degree=per_degree, action=action,
                module_axiom_ok=axiom_ok, inner_acts_zero=inner_ok)


def lie_module_action(A, top, *, budget=DEFAULT_MAXDIM, max_cells=4_000_000,
                      require_char0=False):
    """The full ``HH^*``-as-a-Lie-module-over-``HH^1`` report (:class:`HHLieModule`).

    The field-general block (``hh_dims``, ``hh1_dim``, the per-degree action matrices
    ``rho_n(D)`` tagged ``basis="der_inn/bar"``, the module axiom + inner-acts-zero
    self-certs) is computed over any exact field. The char-0 weight/torus decomposition
    (Task 3) and the indecomposable-summand decomposition (Task 4) are filled only over
    ``char 0`` / where ``decompose``'s own guard holds; over ``char > 0`` the weights are
    ``None`` with a ``char0_note`` (``require_char0=True`` raises loudly for the
    weight/torus accessor). ``budget`` guards ``A.dim`` (the Der solve); ``max_cells``
    guards ``top`` (the exponential bar cochain complex)."""
    dom = A.domain
    if require_char0 and dom.characteristic != 0:
        raise QuiverlabError(
            "the HH-Lie-module weight/torus decomposition needs characteristic 0: "
            "Cartan/Weyl theory has NO positive-characteristic analogue (P70's gate)",
            hint="the field-general action rho_n / module axiom are returned in every "
                 "characteristic; only the weight/torus block is char-0 gated")
    B = A.unit_adapted()
    # guard A.dim up front (the Der solve), mirroring P70's budget guard
    from quiverlab.invariants.hh1_lie import _guard_budget
    _guard_budget(B, budget, "lie_module_action")
    fg = _field_general_action(B, top, max_cells)

    char0_note = None
    if dom.characteristic != 0:
        char0_note = (
            "the weight / torus decomposition needs characteristic 0 (Cartan/Weyl have "
            "no positive-characteristic analogue, P70's gate); the field-general action "
            "rho_n, the module axiom and inner-acts-zero are computed in every "
            "characteristic, and the indecomposable-summand decomposition follows "
            "decompose's own char guard (char 0 or char > dim HH^n) INDEPENDENTLY.")

    weights = None
    torus_provenance = None
    weight_base_change_note = None
    summands = None

    # --- Task 3: char-0 weight / torus decomposition (loud char-p gate) ---
    if dom.characteristic == 0 and fg["hh1_dim"] > 0:
        from quiverlab.invariants.hh1_lie import _bracket_constants
        c = _bracket_constants(fg["reps_hh1"], fg["Inn"], dom, B.dim)
        torus, torus_provenance = _maximal_torus(c, fg["hh1_dim"], dom)
        weights, weight_base_change_note = _weight_tables(
            B, fg["per_degree"], torus, fg["reps_hh1"], dom)

    # --- Task 4: indecomposable-summand decomposition (INDEPENDENT of the weight gate;
    #     governed by decompose's own char guard char 0 or char > d_n) ---
    summands = _summand_tables(fg["per_degree"], dom, budget)

    return HHLieModule(
        top=top, hh_dims=fg["hh_dims"], hh1_dim=fg["hh1_dim"], basis="der_inn/bar",
        action=fg["action"], module_axiom_ok=fg["module_axiom_ok"],
        inner_acts_zero=fg["inner_acts_zero"], characteristic=dom.characteristic,
        weights=weights, torus_provenance=torus_provenance, summands=summands,
        weight_base_change_note=weight_base_change_note, char0_note=char0_note,
        status="complete", note="", references=tuple(_REFERENCES))


def hh_lie_module_block(A, top, *, budget=DEFAULT_MAXDIM, max_cells=4_000_000):
    """The JSON block for the ``hh_lie_module`` compute kind (Plan 71), shared by all
    three tiers (byte-identical). Ships ``hh_dims``, the verdicts, the char-0 WEIGHT
    table and the indecomposable-SUMMAND decomposition table -- the action structure
    constants are NOT shipped (they explode and are basis-dependent). An oversized
    (``A.dim > budget``) or over-``max_cells`` refusal is a clean ``{"error": ...}``
    block, never a 500."""
    refs = list(_REFERENCES)
    try:
        L = lie_module_action(A, top, budget=budget, max_cells=max_cells)
    except (QuiverlabError, DepthLimitError) as exc:
        return {"kind": "hh_lie_module", "status": "budget", "error": str(exc),
                "references": refs}
    summ = None
    if L.summands is not None:
        summ = [{"n": n, "parts": L.summands[n]} for n in range(len(L.summands))]
    return {
        "kind": "hh_lie_module", "top": L.top, "hh_dims": list(L.hh_dims),
        "hh1_dim": L.hh1_dim, "module_axiom_ok": L.module_axiom_ok,
        "inner_acts_zero": L.inner_acts_zero, "characteristic": L.characteristic,
        "weights": L.weights, "torus_provenance": L.torus_provenance,
        "summands": summ, "weight_base_change_note": L.weight_base_change_note,
        "char0_note": L.char0_note, "status": L.status, "note": L.note or None,
        "references": refs,
    }
