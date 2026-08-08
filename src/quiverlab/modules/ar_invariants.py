"""Auslander-Reiten component invariants (Plan 57 / R21): Liu's left/right degrees
of irreducible maps, sectional paths, the postprojective/preinjective/regular
partition, directing modules, the representation-directed recognizer, and the
generalized-standard flag -- all finite sweeps on the knitted AR quiver, exact.

Liu's left degree (arXiv:1704.03933 / ART 2019, corrected `in` direction). For an
irreducible ``f: X_i -> X_j`` (a class of ``rad(X_i, X_j) / rad^2``), the left
degree ``d_l(f)`` is the least ``m`` for which some indecomposable ``Z`` carries a
``g in rad^m(Z, X_i) \\ rad^{m+1}`` with ``g.then(f) in rad^{m+2}(Z, X_j)`` -- a
radical-layer DROP.  In the house left-to-right convention ``g.then(f)`` has matrix
``f.matrix @ g.matrix``.  The right degree is the formal dual with ``Z`` on the
target side (``h in rad^m(X_j, Z)``, ``f.then(h) in rad^{m+2}(X_i, Z)``); it is one
code path parameterized by which argument is fixed (``side``), not a guess.

Everything is decidable because the knit is finite and ``rad^N = 0`` (Plan 57 /
R37).  ``inf`` degree is ``None`` with a ``finite`` flag (no ``float('inf')``).
Certified iff the knit closes (rep-finite); off-scope input refuses honestly.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.ar import _unvec, _vec
from quiverlab.modules.radical import RadicalFiltration, radical_filtration


# --------------------------------------------------------------------------- #
# Liu left / right degrees.
# --------------------------------------------------------------------------- #
def _lincomb(cols, coeffs, dom):
    """``sum coeffs[k] * cols[k]`` (columns are equal-length vectors)."""
    if not cols:
        return []
    length = len(cols[0])
    out = [dom.zero()] * length
    for c, colv in zip(coeffs, cols):
        if dom.is_zero(c):
            continue
        for t in range(length):
            out[t] = dom.add(out[t], dom.mul(c, colv[t]))
    return out


def _preimage(C, R, dom):
    """A basis (in ``C``-coordinates) of ``{ c : sum_k c_k C[k] in colspan(R) }`` --
    the ``c``-parts of the nullspace of the block ``[C | -R]``."""
    p = len(C)
    if p == 0:
        return []
    negR = [[dom.neg(x) for x in rc] for rc in R]
    M = lm.cols_to_matrix(C + negR)                 # L x (p + q)
    null = lm.kernel_columns(M, dom)
    cparts = [nv[:p] for nv in null]
    if not cparts:
        return []
    piv = lm.column_space_pivots(lm.cols_to_matrix(cparts), dom)
    return [cparts[i] for i in piv]


def _default_f_vec(rf: RadicalFiltration, i, j):
    """A representative of a generator of ``rad(X_i, X_j) / rad^2(X_i, X_j)`` (the
    irreducible-map class); refuses loudly when the pair carries no irreducible map."""
    dom = rf._dom
    rad1 = rf._layer_basis(i, j, 1)
    rad2 = rf._layer_basis(i, j, 2)
    idxs = lm.independent_modulo(rad1, rad2, dom)
    if not idxs:
        raise QuiverlabError(
            "left_degree: X_%d -> X_%d carries no irreducible map "
            "(rad/rad^2 = 0 for this pair)" % (i, j),
            hint="pass a pair that is an arrow of the AR quiver")
    return rad1[idxs[0]]


def _degree_sweep(rf: RadicalFiltration, i, j, fmat, side):
    """The Liu degree sweep for an irreducible ``f`` (matrix ``fmat``, ``X_i -> X_j``)
    on ``side in {"left", "right"}``.  Returns ``(d, True, witness)`` or
    ``(None, False, None)``.  ``witness`` for the left sweep is
    ``(z, g_vec, comp_vec)`` with ``g in rad^d(X_z, X_i)`` and ``comp = g.then(f)`` in
    ``Hom(X_z, X_j)`` so the caller can assert all three memberships externally (M3)."""
    dom = rf._dom
    indecs = rf.indecs
    r = len(indecs)
    N = rf.nilpotency_index
    if N is None:
        raise QuiverlabError(
            "degree: the knit did not close (status=%s) -- degrees are defined only on "
            "the certified representation-finite domain" % rf.status)
    di, dj = indecs[i].dim, indecs[j].dim
    for d in range(1, N):                           # d = 1 .. N-1 (rad^N = 0 => termination)
        for z in range(r):
            dz = indecs[z].dim
            if side == "left":
                dom_pair, comp_pair, esc_pair = (z, i), (z, j), (z, i)
                g_rows, g_cols = di, dz             # g: X_z -> X_i  (di x dz)
            else:
                dom_pair, comp_pair, esc_pair = (j, z), (i, z), (j, z)
                g_rows, g_cols = dz, dj             # h: X_j -> X_z  (dz x dj)
            B_d = rf._layer_basis(*dom_pair, d)
            if not B_d:
                continue
            C = []
            for gk in B_d:
                gmat = _unvec(gk, g_rows, g_cols)
                comp = lm.matmul(fmat, gmat, dom) if side == "left" \
                    else lm.matmul(gmat, fmat, dom)
                C.append(_vec(comp))
            R = rf._layer_basis(*comp_pair, d + 2)
            Kc = _preimage(C, R, dom)
            if not Kc:
                continue
            gvecs = [_lincomb(B_d, c, dom) for c in Kc]
            esc = rf._layer_basis(*esc_pair, d + 1)
            idxs = lm.independent_modulo(gvecs, esc, dom)
            if idxs:
                t = idxs[0]
                witness = (z, gvecs[t], _lincomb(C, Kc[t], dom))
                return d, True, witness
    return None, False, None


def left_degree(rf: RadicalFiltration, i, j, f_vec=None):
    """Liu's left degree ``d_l`` of the irreducible-map class ``X_i -> X_j``.
    Returns ``(d, True, (z, g_vec, comp_vec))`` for a finite degree ``d`` -- the
    witness certifies the layer drop -- or ``(None, False, None)`` for infinite
    degree.  ``f_vec`` defaults to the first ``rad/rad^2`` class representative."""
    if f_vec is None:
        f_vec = _default_f_vec(rf, i, j)
    fmat = _unvec(f_vec, rf.indecs[j].dim, rf.indecs[i].dim)
    return _degree_sweep(rf, i, j, fmat, "left")


def right_degree(rf: RadicalFiltration, i, j, f_vec=None):
    """Liu's right degree ``d_r`` -- the formal dual of :func:`left_degree` with the
    target on the free side (M2, one parameterized code path, no 'flip and see').
    Returns ``(d, True, (z, h_vec, comp_vec))`` or ``(None, False, None)``."""
    if f_vec is None:
        f_vec = _default_f_vec(rf, i, j)
    fmat = _unvec(f_vec, rf.indecs[j].dim, rf.indecs[i].dim)
    return _degree_sweep(rf, i, j, fmat, "right")


def _class_reps(rf: RadicalFiltration, i, j):
    """The class representatives of ``rad(X_i, X_j) / rad^2`` (one per arrow
    multiplicity): the ``rad^1`` basis vectors independent modulo ``rad^2``."""
    dom = rf._dom
    rad1 = rf._layer_basis(i, j, 1)
    rad2 = rf._layer_basis(i, j, 2)
    idxs = lm.independent_modulo(rad1, rad2, dom)
    return [rad1[t] for t in idxs]


def _min_finite(vals):
    """The least FINITE degree in ``vals`` (``None`` = infinite), or ``None`` if every
    entry is infinite."""
    finite = [v for v in vals if v is not None]
    return min(finite) if finite else None


def degree_table(rf: RadicalFiltration, ar):
    """``{(i, j): {...}}`` over the irreducible-map arrows of the knit. A
    multiplicity-1 arrow carries ``{"d_l", "d_r", "finite_l", "finite_r"}`` (``d`` is
    ``None`` for an infinite degree). A multiplicity > 1 arrow additionally carries
    ``{"mult", "d_l_classes", "d_r_classes"}`` -- the per-class list, with ``d_l``/``d_r``
    the min (finite-preferring) over the classes (Plan-57 honest per-class report).

    # PIN (missing algebra-level fixture): every knittable representation-finite
    # quiver algebra in scope has ONLY multiplicity-1 AR arrows (`dim rad(X,Y)/rad^2 =
    # 1`) -- verified over kA_4/kA_5, D_4/D_5 and the Nakayama zoo [3,2,2]/[4,3,2,2].
    # Multiplicity > 1 needs a non-trivial AR valuation (a species / non-algebraically-
    # closed base field), which `knit_ar_quiver` (projective-seeded quiver-algebra BFS)
    # does not produce. The per-class code path is therefore covered by a synthetic
    # `_class_reps` unit test (`tests/modules/test_liu_degrees.py`); a genuine mult>1
    # algebra-level fixture is DEFERRED until the AR engine grows a species surface."""
    out = {}
    for (i, j) in ar.arrows:
        reps = _class_reps(rf, i, j)
        if len(reps) <= 1:
            f = reps[0] if reps else None
            dl, fin_l, _ = left_degree(rf, i, j, f)
            dr, fin_r, _ = right_degree(rf, i, j, f)
            out[(i, j)] = {"d_l": dl, "d_r": dr,
                           "finite_l": bool(fin_l), "finite_r": bool(fin_r)}
        else:
            dls = [left_degree(rf, i, j, f) for f in reps]
            drs = [right_degree(rf, i, j, f) for f in reps]
            dl_cls = [d for (d, _f, _w) in dls]
            dr_cls = [d for (d, _f, _w) in drs]
            out[(i, j)] = {
                "d_l": _min_finite(dl_cls), "d_r": _min_finite(dr_cls),
                "finite_l": any(fin for (_d, fin, _w) in dls),
                "finite_r": any(fin for (_d, fin, _w) in drs),
                "mult": len(reps),
                "d_l_classes": dl_cls, "d_r_classes": dr_cls}
    return out


# --------------------------------------------------------------------------- #
# Sectional paths (Bautista-Smalo: a sectional path has a nonzero composite).
# --------------------------------------------------------------------------- #
def _tau_index_map(ar):
    """``{k: m}`` with ``X_m ~ tau(X_k)`` (omitted when ``X_k`` is projective,
    ``tau = 0``), matched by the exact iso certificate, dimension-vector prefiltered."""
    from quiverlab.modules.duality import tau
    from quiverlab.modules.hom import is_isomorphic
    mods = [v["module"] for v in ar.vertices]
    out = {}
    for k, X in enumerate(mods):
        t = tau(X)
        if t.dim == 0:
            continue
        dv = t.dimension_vector()
        for m, Y in enumerate(mods):
            if Y.dim == t.dim and Y.dimension_vector() == dv and is_isomorphic(Y, t):
                out[k] = m
                break
    return out


def is_sectional(ar, path):
    """``True`` iff ``path`` (a list of vertex indices, consecutive pairs arrows of
    ``Gamma_A``) has no mesh reversal ``X_{k-1} = tau X_{k+1}``.  Sectional composites
    are nonzero (Bautista-Smalo)."""
    for k in range(len(path) - 1):
        if (path[k], path[k + 1]) not in ar.arrows:
            raise QuiverlabError(
                "is_sectional: (%s, %s) is not an arrow of Gamma_A"
                % (path[k], path[k + 1]))
    tau_idx = _tau_index_map(ar)
    for k in range(1, len(path) - 1):
        if tau_idx.get(path[k + 1]) == path[k - 1]:
            return False
    return True


def max_sectional_length(ar):
    """The maximal length (number of arrows) of a sectional path in ``Gamma_A``.
    Exact on a directed (acyclic) component; on a component with an oriented cycle the
    DFS is depth-capped at the vertex count (an honest informative bound)."""
    if not getattr(ar, "is_complete", False):
        return 0
    tau_idx = _tau_index_map(ar)
    adj = {}
    for (a, b) in ar.arrows:
        adj.setdefault(a, []).append(b)
    n = len(ar.vertices)
    best = 0

    def dfs(path):
        nonlocal best
        best = max(best, len(path) - 1)
        if len(path) > n:                           # cap: terminate on cyclic components
            return
        last = path[-1]
        for nxt in adj.get(last, []):
            if len(path) >= 2 and tau_idx.get(nxt) == path[-2]:
                continue                            # mesh reversal -> not sectional
            dfs(path + [nxt])

    for v in range(n):
        dfs([v])
    return best


# --------------------------------------------------------------------------- #
# Partition + directing + rep-directed.
# --------------------------------------------------------------------------- #
def _partition(ar):
    """``{k: "postprojective"|"preinjective"|"regular"}`` by walking each tau-orbit:
    the orbit contains a projective (``tau = 0``) => postprojective; else an injective
    (``tau^- = 0``) => preinjective; else regular.  On a representation-directed
    component every orbit reaches a projective (all postprojective, no regular)."""
    from quiverlab.modules.duality import tau, tau_minus
    mods = [v["module"] for v in ar.vertices]
    is_proj = [tau(X).dim == 0 for X in mods]
    is_inj = [tau_minus(X).dim == 0 for X in mods]
    part = {}
    for orbit in ar.tau_orbits:
        has_p = any(is_proj[k] for k in orbit)
        has_i = any(is_inj[k] for k in orbit)
        label = "postprojective" if has_p else ("preinjective" if has_i else "regular")
        for k in orbit:
            part[k] = label
    return part


def _tarjan_sccs(n, adj):
    """Strongly-connected components of the digraph ``0..n-1`` (iterative Tarjan)."""
    index = [None] * n
    low = [0] * n
    on_stack = [False] * n
    stack = []
    sccs = []
    counter = [0]
    for start in range(n):
        if index[start] is not None:
            continue
        work = [(start, 0)]                          # (node, next-neighbour pointer)
        while work:
            v, pi = work[-1]
            if pi == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                on_stack[v] = True
            neigh = adj.get(v, [])
            if pi < len(neigh):
                work[-1] = (v, pi + 1)
                w = neigh[pi]
                if index[w] is None:
                    work.append((w, 0))
                elif on_stack[w]:
                    low[v] = min(low[v], index[w])
            else:
                if low[v] == index[v]:
                    comp = []
                    while True:
                        w = stack.pop()
                        on_stack[w] = False
                        comp.append(w)
                        if w == v:
                            break
                    sccs.append(comp)
                work.pop()
                if work:
                    u = work[-1][0]
                    low[u] = min(low[u], low[v])
    return sccs


def _directing(ar):
    """``(directing_set, is_representation_directed)``.  ``X`` is directing iff its SCC
    in ``Gamma_A`` is a singleton with no self-loop -- the arrows are irreducible maps
    (nonzero non-isomorphisms), so an oriented cycle IS a Ringel cycle of nonzero
    non-isos through ``X`` (no separate nonzero-composite test is needed, and the total
    round-trip composite is 0 in ``rad^{>=N}`` yet the module is still non-directing).
    ``A`` is representation-directed iff every vertex is directing (``Gamma_A``
    acyclic)."""
    n = len(ar.vertices)
    adj = {v: [] for v in range(n)}
    self_loops = set()
    for (a, b) in ar.arrows:
        if a == b:
            self_loops.add(a)
        else:
            adj[a].append(b)
    directing = set()
    for comp in _tarjan_sccs(n, adj):
        if len(comp) == 1 and comp[0] not in self_loops:
            directing.add(comp[0])
    return directing, (len(directing) == n)


# --------------------------------------------------------------------------- #
# The frozen result container + the entry point.
# --------------------------------------------------------------------------- #
@dataclass
class ARInvariants:
    partition: dict                       # {vertex index -> "postprojective"|"preinjective"|"regular"}
    directing: set                        # directing vertex indices
    is_representation_directed: bool
    generalized_standard: bool
    nilpotency_index: object              # int | None
    max_sectional_length: int
    degrees: dict                         # {(i, j) -> {"d_l", "d_r", "finite_l", "finite_r"}}
    is_complete: bool
    status: str
    note: str = ""
    names: list = field(default_factory=list)


def ar_invariants(A, *, budget_modules=256, budget_dim=4096):
    """The AR-component invariants of ``A`` (Plan 57 / R21).  Certified iff the knit
    closes (rep-finite); an incomplete knit (self-injective / budget / error) returns
    a loud off-scope container -- never a false 'directed' or a spurious index."""
    rf = radical_filtration(A, budget_modules=budget_modules, budget_dim=budget_dim)
    ar = rf._ar
    names = list(rf.names)
    if not ar.is_complete:
        note = (rf.note or ar.note or ("out of scope (status=%s)" % ar.status)) + \
            " -- AR-component invariants are certified only on the representation-" \
            "finite (knit-complete) domain; no directing/partition/degree verdict is " \
            "issued here"
        return ARInvariants(
            partition={}, directing=set(), is_representation_directed=False,
            generalized_standard=False, nilpotency_index=None,
            max_sectional_length=0, degrees={}, is_complete=False,
            status=ar.status, note=note, names=names)
    partition = _partition(ar)
    directing, rep_directed = _directing(ar)
    degs = degree_table(rf, ar)
    msl = max_sectional_length(ar)
    note = ("representation-finite: %d indecomposables; representation-directed=%s "
            "(Gamma_A %s), generalized standard (rad^inf = 0, witness N = %s)"
            % (len(rf.indecs), rep_directed,
               "acyclic" if rep_directed else "has an oriented cycle",
               rf.nilpotency_index))
    if not rep_directed:                            # honest labeling of the partition
        note = note + " -- " + _NONDIRECTED_CAVEAT
    return ARInvariants(
        partition=partition, directing=directing,
        is_representation_directed=rep_directed,
        generalized_standard=rf.is_generalized_standard(),
        nilpotency_index=rf.nilpotency_index, max_sectional_length=msl,
        degrees=degs, is_complete=True, status="complete", note=note, names=names)


# --------------------------------------------------------------------------- #
# Shared no-code ALGEBRA-block builder (byte-identical across both runners: both
# import THIS library function).  Every value is JSON-safe; Module objects are
# dropped.  A refusal (char scope) is caught by the caller into {"error": ...}.
# The radical_filtration block builder lives in modules/radical.py (the R37 home).
# --------------------------------------------------------------------------- #
_AR_INV_REFS = ["liu_degrees", "liu_semistable", "ringel_tame"]

# The plan-mandated honesty caveat on a NON-directed component: the τ-partition is
# the clean postprojective/preinjective/regular trichotomy ONLY on a directed
# component; on a non-directed one the buckets can overlap, so the classification is
# reported as a τ-orbit reading, never the clean trichotomy. Carried as block DATA
# (like ar_quiver's `partial_note`), so it needs no i18n key.
_NONDIRECTED_CAVEAT = ("τ-orbit classification; non-directed component — buckets "
                       "may overlap")


def ar_invariants_block(A, *, budget=512):
    """The ``ar_invariants`` algebra block (Plan 57 / R21): Liu degrees, the
    postprojective/preinjective/regular partition, directing modules, the
    representation-directed recognizer, sectional-path length and the
    generalized-standard flag.  Honest off-scope refusal like
    :func:`radical_filtration_block`."""
    inv = A.ar_invariants(budget_modules=budget)
    part_counts = {}
    for lab in inv.partition.values():
        part_counts[lab] = part_counts.get(lab, 0) + 1
    degrees = {}
    for (i, j), rec in inv.degrees.items():
        # Copy the whole record: a multiplicity-1 arrow carries exactly
        # {d_l, d_r, finite_l, finite_r} (byte-unchanged), a mult>1 arrow additionally
        # carries {mult, d_l_classes, d_r_classes} (Plan-57 per-class honesty).
        degrees["%d->%d" % (i, j)] = dict(rec)
    block = {
        "kind": "ar_invariants",
        "status": inv.status,
        "complete": bool(inv.is_complete),
        "budget": int(budget),
        "representation_directed": (inv.is_representation_directed
                                    if inv.is_complete else None),
        "generalized_standard": inv.generalized_standard if inv.is_complete else None,
        "nilpotency_index": inv.nilpotency_index,
        "num_directing": len(inv.directing),
        "num_indecomposables": len(inv.names),
        "max_sectional_length": inv.max_sectional_length,
        "partition_counts": part_counts,
        "partition": {str(k): v for k, v in sorted(inv.partition.items())},
        "degrees": degrees,
        "vertices": [{"name": nm} for nm in inv.names],
        "note": inv.note,
        "latex": (r"d_\ell,\ d_r\ \text{(Liu)},\ \text{postproj/preinj/regular},\ "
                  r"\text{directing},\ \Gamma_A\ \text{acyclic}"),
        "references": list(_AR_INV_REFS),
    }
    # Honest partition caveat: on a NON-directed component the buckets can overlap, so
    # the classification is a tau-orbit reading, never the clean trichotomy (plan-
    # mandated). Block data, like ar_quiver's partial_note -- no i18n key needed.
    if inv.is_complete and not inv.is_representation_directed:
        block["partition_note"] = _NONDIRECTED_CAVEAT
    return block
