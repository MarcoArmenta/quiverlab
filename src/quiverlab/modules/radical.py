"""The radical filtration of ``mod A`` (Plan 57 / R37, the Liu-Chaio program).

Exact linear algebra on the knitted indecomposable universe (route (ii), the
arbiter): for a representation-finite ``A`` the AR knitter (Plan 41
:func:`quiverlab.modules.ar.knit_ar_quiver`) produces the finite set
``ind A = {X_1, ..., X_r}``, and the radical of ``mod A`` is the two-sided ideal
``rad(X, Y)`` of non-isomorphisms.  On indecomposables

  * ``rad(X_i, X_j) = Hom_A(X_i, X_j)`` for ``i != j`` (a map between
    non-isomorphic indecomposables is never invertible), and
  * ``rad(X_i, X_i) = rad End_A(X_i)`` (the trace-form radical of the local ring,
    char 0 or char > dim X_i).

Higher powers, in the house left-to-right convention (``p in rad^n(X_i, X_k)``,
``q in rad(X_k, X_j)``, composite matrix ``q.matrix @ p.matrix``)::

    rad^{n+1}(X_i, X_j) = sum_{k} { compose(p, q) : p in rad^n(X_i, X_k),
                                                     q in rad(X_k, X_j) }

Summing over the INDECOMPOSABLE ``X_k`` only loses nothing (Krull-Schmidt +
Hom-additivity -- every module splits and the biproduct maps are radical when the
summands are non-isomorphic).  Each ``rad^n(X_i, X_j)`` is a coordinate subspace of
``Hom(X_i, X_j)`` (column-stacked vecs), stored as a reduced column basis; its
dimension is a ``mat_rank``.  The nilpotency index ``N`` is the least ``n`` with
``rad^n = 0`` on every pair; ``rad^inf = 0`` iff the knit closed (Auslander's
representation-finiteness certificate, with ``N`` its finite witness).

Honest semi-decision contract, mirrored from :class:`~quiverlab.modules.ar.ARQuiver`:
certified iff the knit is complete (rep-finite).  A self-injective input is refused
(``status="unsupported"``); a budget-exhausted knit yields a window-restricted view
of the discovered subcategory (a lower bound, NOT ``rad(mod A)``) with NO nilpotency
or ``rad^inf`` verdict.  Exact only -- integers/dicts, no floats; ``inf`` degree is
``None`` with a boolean flag (never ``float('inf')``).
"""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.ar import _rad_basis, _unvec, _vec, knit_ar_quiver

# Window-branch guards (rep-infinite / budget input): compute the discovered
# subcategory's layers only when it is cheap, and never claim an index.
_WINDOW_MODULE_CAP = 24
_WINDOW_DEGREE_CAP = 12
_INDEX_SAFETY = 64          # rep-finite N is tiny; a runaway loop is a bug -> raise


def _reduce_cols(cols, dom):
    """A reduced (independent) column basis spanning the same space as ``cols``."""
    if not cols:
        return []
    piv = lm.column_space_pivots(lm.cols_to_matrix(cols), dom)
    return [cols[p] for p in piv]


def _all_zero(vec, dom):
    return all(dom.is_zero(x) for x in vec)


def _rad1_layer(indecs, dom):
    """``rad^1`` on every ordered pair: reduced column bases keyed ``(i, j)``.
    Reuses the P41 ``_rad_basis`` (full Hom off the diagonal, trace-form ``rad End``
    on it -- the char guard fires there over char <= dim)."""
    r = len(indecs)
    out = {}
    for i in range(r):
        for j in range(r):
            out[(i, j)] = _reduce_cols(_rad_basis(indecs[i], indecs[j]), dom)
    return out


def _compose_layer(prev, rad1, indecs, dom):
    """``rad^{n+1}`` from ``rad^n`` (``prev``) and ``rad^1`` (``rad1``): the column
    space of every composite ``compose(p, q)``, ``p in rad^n(i, k)``,
    ``q in rad^1(k, j)`` (matrix ``q @ p``), over all intermediate ``k``."""
    r = len(indecs)
    out = {}
    for i in range(r):
        di = indecs[i].dim
        for j in range(r):
            dj = indecs[j].dim
            comps = []
            for k in range(r):
                dk = indecs[k].dim
                left = prev.get((i, k), [])          # rad^n(i, k): vecs dk x di
                right = rad1.get((k, j), [])         # rad^1(k, j): vecs dj x dk
                if not left or not right:
                    continue
                for lv in left:
                    p = _unvec(lv, dk, di)           # X_i -> X_k
                    for rv in right:
                        q = _unvec(rv, dj, dk)       # X_k -> X_j
                        comps.append(_vec(lm.matmul(q, p, dom)))   # X_i -> X_j
            out[(i, j)] = _reduce_cols(comps, dom)
    return out


class RadicalFiltration:
    """The radical filtration of ``mod A``: per-pair layer dimensions
    ``dim rad^n(X_i, X_j)``, the nilpotency index ``N``, and the
    ``rad^inf = 0`` certificate (rep-finite).  See the module docstring."""

    def __init__(self, indecs, names, layer, rad, nilpotency_index,
                 rad_infinity_zero, is_complete, status, note, dom):
        self.indecs = indecs                       # list[Module] == ARQuiver vertex modules
        self.names = names                         # list[str|None] (identify_standard names)
        self.layer = layer                         # {(i, j): [dim rad^n]_{n=1..len}}
        self._rad = rad                            # {n: {(i, j): reduced column basis}}
        self.nilpotency_index = nilpotency_index   # int N | None (None off certified scope)
        self.rad_infinity_zero = rad_infinity_zero  # True iff complete | None otherwise
        self.is_complete = is_complete
        self.status = status
        self.note = note
        self._dom = dom
        self._ar = None                            # stashed ARQuiver (set by radical_filtration)

    # -- public accessors --------------------------------------------------- #
    def layer_dim(self, i, j, n):
        """``dim rad^n(X_i, X_j)`` (``0`` for ``n`` past the last stored layer)."""
        seq = self.layer.get((i, j), [])
        return seq[n - 1] if 1 <= n <= len(seq) else 0

    def pair_layers(self, X, Y):
        """``[dim rad^n(X, Y)]_{n>=1}`` for two knit indecomposables (matched by the
        exact ``is_isomorphic`` certificate, dimension-vector prefiltered)."""
        i, j = self._match(X), self._match(Y)
        if i is None or j is None:
            raise QuiverlabError(
                "radical_filtration.pair_layers: X or Y is not one of the knitted "
                "indecomposables of this algebra",
                hint="pass modules from rf.indecs (the ind A universe)")
        return list(self.layer.get((i, j), []))

    def is_generalized_standard(self):
        """``True`` in the certified (rep-finite) scope: ``rad^inf(X, Y) = 0`` for all
        ``X, Y`` (Auslander), witnessed by the finite nilpotency index; the flag is
        informative only for representation-infinite input (out of certified scope,
        where it is not decided here)."""
        return bool(self.is_complete)

    # -- internals (Task 2's degree sweep consumes these) ------------------- #
    def _layer_basis(self, i, j, n):
        layer_n = self._rad.get(n)
        return layer_n.get((i, j), []) if layer_n is not None else []

    def _in_layer(self, i, j, n, vec):
        """``True`` iff ``vec`` (a column-stacked map in ``Hom(X_i, X_j)``) lies in
        ``rad^n(X_i, X_j)``.  ``rad^n = 0`` for ``n >= N`` -- then only the zero map
        qualifies."""
        dom = self._dom
        if _all_zero(vec, dom):
            return True
        basis = self._layer_basis(i, j, n)
        if not basis:
            return False
        return lm.solve_columns(lm.cols_to_matrix(basis),
                                lm.cols_to_matrix([vec]), dom) is not None

    def _match(self, M):
        from quiverlab.modules.hom import is_isomorphic
        dv = M.dimension_vector()
        for idx, X in enumerate(self.indecs):
            if X.dim == M.dim and X.dimension_vector() == dv and is_isomorphic(X, M):
                return idx
        return None


def _window_note(k, status):
    return ("window-restricted to the %d discovered indecomposables (add of the "
            "frontier); this is a lower bound on rad(mod A), NOT the category radical "
            "-- no nilpotency index or rad^inf verdict is issued (knit status: %s)"
            % (k, status))


def radical_filtration(A, *, budget_modules=256, budget_dim=4096, _ar=None):
    """The radical filtration of ``mod A`` (route (ii); Plan 57 / R37).

    Knits the AR quiver, then builds ``rad^n(X_i, X_j)`` exactly by iterated matrix
    composition.  Complete knit (rep-finite) => certified nilpotency index and
    ``rad^inf = 0``; self-injective / budget-exhausted input => an honest
    window/refusal with no verdict.  Pass ``_ar`` to reuse an already-knitted quiver.
    """
    ar = _ar if _ar is not None else knit_ar_quiver(
        A, budget_modules=budget_modules, budget_dim=budget_dim)
    dom = A.domain
    indecs = [v["module"] for v in ar.vertices]
    names = [v["name"] for v in ar.vertices]

    if not ar.is_complete:
        # Off the certified domain: no index, no rad^inf.  Self-injective (empty
        # vertex set) surfaces the knitter's own note; a budget/error knit ships the
        # discovered subcategory's layers (cheap only) with the loud window note.
        rad, layer = {}, {}
        if ar.status == "unsupported" or not indecs:
            note = ar.note or ("out of scope (knit status: %s)" % ar.status)
        else:
            note = _window_note(len(indecs), ar.status)
            if ar.note:
                note = note + " [" + ar.note + "]"
            if len(indecs) <= _WINDOW_MODULE_CAP:
                rad, layer, _ = _build_layers(indecs, dom, degree_cap=_WINDOW_DEGREE_CAP)
        rf = RadicalFiltration(indecs, names, layer, rad, None, None,
                               False, ar.status, note, dom)
        rf._ar = ar
        return rf

    rad, layer, N = _build_layers(indecs, dom, degree_cap=None)
    note = ("representation-finite: the knit closed on %d indecomposables, so "
            "rad^inf(mod A) = 0 (Auslander) with nilpotency index N = %d the finite "
            "witness" % (len(indecs), N))
    rf = RadicalFiltration(indecs, names, layer, rad, N, True, True, "complete",
                           note, dom)
    rf._ar = ar
    return rf


def _build_layers(indecs, dom, degree_cap=None):
    """``(rad, layer, N)`` from route (ii).  ``rad[n][(i, j)]`` = reduced column basis
    of ``rad^n(X_i, X_j)``; ``layer[(i, j)]`` = ``[dim rad^n]`` for the nonzero layers;
    ``N`` = the nilpotency index (least ``n`` with ``rad^n == 0`` everywhere) or, when a
    ``degree_cap`` is given and hit first (the window), ``None``."""
    if not indecs:
        return {}, {}, 0
    r = len(indecs)
    rad = {1: _rad1_layer(indecs, dom)}
    n = 1
    while True:
        if all(not rad[n].get((i, j)) for i in range(r) for j in range(r)):
            N = n
            break
        if degree_cap is not None and n >= degree_cap:
            N = None                                # window: capped, no index claim
            break
        if n > _INDEX_SAFETY:
            raise QuiverlabError(
                "radical_filtration: nilpotency index exceeded the safety cap %d on a "
                "complete knit -- this should be impossible for a rep-finite algebra "
                "(likely a composition bug)" % _INDEX_SAFETY)
        rad[n + 1] = _compose_layer(rad[n], rad[1], indecs, dom)
        n += 1
    # Layer sequences: keep the nonzero prefix (n = 1 .. N-1 for a complete knit).
    top = (N - 1) if N is not None else max(rad)
    layer = {}
    for i in range(r):
        for j in range(r):
            seq = [len(rad[m].get((i, j), [])) for m in range(1, top + 1)]
            while seq and seq[-1] == 0:
                seq.pop()
            layer[(i, j)] = seq
    # Drop the all-zero terminal layer (rad^N) from the stored bases.
    if N is not None and N in rad:
        del rad[N]
    return rad, layer, N


# --------------------------------------------------------------------------- #
# _thm13_index -- the Chaio-Guazzelli nilpotency-index formula (arXiv:2003.04189,
# Theorem 1.3), read off the SAME knit/route-(ii) layers.  Its agreement with
# ``nilpotency_index`` is an internal-consistency self-cert (two readings of one
# engine), NOT a cross-engine oracle.
# --------------------------------------------------------------------------- #
def _match_index(indecs, M):
    from quiverlab.modules.hom import is_isomorphic
    dv = M.dimension_vector()
    for idx, X in enumerate(indecs):
        if X.dim == M.dim and X.dimension_vector() == dv and is_isomorphic(X, M):
            return idx
    return None


def _thm13_index(A, *, budget_modules=256, budget_dim=4096):
    """Chaio-Guazzelli Theorem 1.3: ``index(rad mod A) = max_a { r_a + 1 }`` with
    ``r_a`` the length (radical degree) of the nonzero path of irreducible morphisms
    from ``P_a`` to ``I_a`` through ``S_a``.  On a serial / linear-``A_n`` algebra
    ``Hom(P_a, I_a)`` is at most one dimensional and that nonzero map factors through
    the top ``S_a`` = socle of ``I_a``, so ``r_a`` is exactly its radical degree
    ``max{ m : dim rad^m(P_a, I_a) > 0 }`` -- read straight off route (ii)."""
    rf = radical_filtration(A, budget_modules=budget_modules, budget_dim=budget_dim)
    if not rf.is_complete:
        raise QuiverlabError(
            "_thm13_index: the knit did not close (status=%s) -- Theorem 1.3 applies "
            "only to a representation-finite algebra" % rf.status)
    indecs = rf.indecs
    best = 1
    for a in A.quiver.vertices:
        pi = _match_index(indecs, A.projective(a))
        ii = _match_index(indecs, A.injective(a))
        if pi is None or ii is None:
            raise QuiverlabError(
                "_thm13_index: P_%s or I_%s is not among the knitted indecomposables "
                "(the knit is incomplete or a matching failed)" % (a, a))
        seq = rf.layer.get((pi, ii), [])
        r_a = len(seq)                              # max m with dim rad^m(P_a, I_a) > 0
        best = max(best, r_a + 1)
    return best


# --------------------------------------------------------------------------- #
# _mesh_layer_dim -- the GENUINELY INDEPENDENT route (i): the ZA_n mesh closed
# form for a standard directed type-A_n component, computed by PURE interval
# combinatorics on the dimension vectors (NO Hom, NO matmul).  Scoped to linear
# A_n interval modules (the shipped kA_n crosscheck); a non-standard or non-A_n
# component is out of scope (the documented honest boundary -- a general
# functorial mesh engine + standardness detector is the named deferral).
# --------------------------------------------------------------------------- #
def _interval(dimvec):
    """``(a, b)`` for a type-``A_n`` interval module (contiguous support, all ones);
    ``None`` if the dimension vector is not such an interval (out of route-(i) scope)."""
    supp = sorted(v for v, c in dimvec.items() if c)
    if not supp:
        return None
    if any(dimvec[v] != 1 for v in supp):
        return None
    if supp != list(range(supp[0], supp[-1] + 1)):
        return None
    return supp[0], supp[-1]


def _mesh_layer_dim(ar, i, j, n):
    """``dim rad^n(X_i, X_j)`` for a standard directed ``kA_n`` component via the ZA_n
    mesh closed form.  The indecomposables are the intervals ``[a, b]``; there is a
    nonzero map ``[a, b] -> [c, d]`` iff ``c <= a <= d <= b``, its radical degree is
    ``(a - c) + (b - d)``, and (descending filtration) it contributes ``1`` to every
    ``rad^m`` with ``1 <= m <= that degree``.  Purely combinatorial -- the independent
    crosscheck of route (ii)."""
    if i == j:
        return 0
    Ii = _interval(ar.vertices[i]["dimvec"])
    Ij = _interval(ar.vertices[j]["dimvec"])
    if Ii is None or Ij is None:
        raise QuiverlabError(
            "_mesh_layer_dim: route (i) is scoped to linear A_n interval modules; "
            "vertex %d or %d is not a contiguous-support interval" % (i, j))
    a, b = Ii
    c, d = Ij
    if not (c <= a <= d <= b):
        return 0
    degree = (a - c) + (b - d)
    return 1 if 1 <= n <= degree else 0


# --------------------------------------------------------------------------- #
# The shared no-code ALGEBRA-block builder (byte-identical across both runners:
# server hpc.spec and the Pyodide docs/gui/runner both import THIS function).
# Every value JSON-safe; Module objects dropped.  A raised QuiverlabError (char
# scope) is caught by the caller into {"error": ...} (the per-block precedent).
# --------------------------------------------------------------------------- #
_RADICAL_REFS = ["chaio_liu_radical", "liu_degrees", "cmms_radsq"]


def _dimvec_json(dv):
    return {str(v): int(c) for v, c in sorted(dv.items(), key=lambda kv: str(kv[0]))}


def _layer_totals(rf):
    """``[sum_ij dim rad^n(X_i, X_j)]_{n>=1}`` (the layer profile), or ``None`` off the
    certified scope."""
    if not rf.is_complete or rf.nilpotency_index is None:
        return None
    r = len(rf.indecs)
    top = rf.nilpotency_index - 1
    return [sum(rf.layer_dim(i, j, n) for i in range(r) for j in range(r))
            for n in range(1, top + 1)]


def radical_filtration_block(A, *, budget=512):
    """The ``radical_filtration`` algebra block (Plan 57 / R37): the radical of the
    MODULE CATEGORY ``rad^n(X, Y)`` and its nilpotency index -- NOT the Loewy radical
    series spectral sequence (``radical_filtration_ss``).  Honest semi-decision: a
    complete (rep-finite) knit ships the index + ``rad^inf`` verdict; a self-injective
    or budget knit ships ``complete=False`` with the loud note and no verdict."""
    rf = A.radical_filtration(budget_modules=budget)
    return {
        "kind": "radical_filtration",
        "status": rf.status,
        "complete": bool(rf.is_complete),
        "budget": int(budget),
        "num_indecomposables": len(rf.indecs),
        "nilpotency_index": rf.nilpotency_index,
        "rad_infinity_zero": rf.rad_infinity_zero,
        "generalized_standard": (rf.is_generalized_standard()
                                 if rf.is_complete else None),
        "layer_profile": _layer_totals(rf),
        "vertices": [{"name": nm, "dimvec": _dimvec_json(X.dimension_vector())}
                     for nm, X in zip(rf.names, rf.indecs)],
        "note": rf.note,
        "latex": (r"\operatorname{rad}^{n}(X,Y)\ \text{and}\ "
                  r"N=\min\{n:\operatorname{rad}^n=0\}"),
        "references": list(_RADICAL_REFS),
    }
