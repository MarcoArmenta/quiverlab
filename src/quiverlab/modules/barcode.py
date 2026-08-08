"""Barcodes = interval decompositions of A_n / zigzag persistence modules, and
AR-quiver-indexed generalized persistence diagrams of commutative ladders CL(n<=4).

Plan 69 / R33 -- the TDA bridge, representation theory FIRST:

* A *persistence module* over ``A_n`` is a representation of the type-A quiver; by
  **Gabriel** its indecomposables are the **interval (thin) modules** ``I[b,d]``, each
  a **brick** (``dim End = 1``), so ``decompose`` (Plan 30) certifies them in ANY
  characteristic. The **barcode** is the multiset of support intervals of the
  Krull-Schmidt summands (Gabriel / Botnan-Crawley-Boevey). Forward orientation =>
  ordinary persistence; an alternating orientation => a **zigzag** module. This route
  is **field-robust over every exact domain including GF(2)** (bricks).

* A **commutative ladder** ``CL(n) = A_n [] A_2`` is representation-finite iff
  ``n <= 4`` (Escolar-Hiraoka). Its **generalized persistence diagram** is
  ``decompose(M)`` with every summand matched to its vertex in the knitted AR quiver
  (Plan 41). CL indecomposables need not be bricks, so the CL route is **char-scoped**
  (char 0 or char > dim); off scope the knit refuses (``is_complete=False``) and
  ``barcode`` surfaces it loudly.

Exact only: filtration parameter = the discrete vertex index (``int``); no floats, no
``inf`` (a top-reaching FORWARD bar is ``death = n`` with ``essential = True``).
"""
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.decompose import decompose

_REFS = ("escolar_hiraoka", "botnan_crawley_boevey", "gabriel", "assem_book")


@dataclass(frozen=True)
class Bar:
    """One bar ``[birth, death]`` (an interval summand ``I[birth, death]``).

    ``essential`` is the FORWARD-persistence flag ``death == n`` (M2): it is
    meaningful only for ``kind="persistence"`` and is **always ``False`` for
    ``kind="zigzag"``** (a zigzag has no monotone 'top', so 'alive at the top' has no
    meaning). ``dimvec`` is the summand's dimension vector (vertex-keyed, nonzero
    entries)."""
    birth: int
    death: int
    multiplicity: int
    essential: bool
    dimvec: dict


@dataclass(frozen=True)
class DiagramEntry:
    """One indecomposable class of a CL generalized persistence diagram, matched to a
    knitted AR-quiver vertex (``ar_index`` into ``ARQuiver.vertices``; ``ar_name`` its
    standard name or ``X{i}``). ``is_interval`` iff the summand's support is thin
    (contiguous 0/1) on the CL base."""
    ar_index: int
    ar_name: str
    dimvec: dict
    multiplicity: int
    is_interval: bool


@dataclass(frozen=True)
class Barcode:
    kind: str                 # "persistence" | "zigzag" | "commutative_ladder"
    line_order: tuple         # the A_n vertex order used to read intervals (or CL base)
    bars: tuple               # tuple[Bar] -- interval summands
    diagram: tuple            # CL only: tuple[DiagramEntry]; () for A_n/zigzag
    n: int                    # ladder/line length
    field: str
    field_robust: bool        # True for A_n/zigzag (bricks); False for CL (char-scoped)
    ar_status: str            # CL only: knit status ("complete" for n<=4); "" otherwise
    references: tuple


# --------------------------------------------------------------------------- #
# Line-shape classification
# --------------------------------------------------------------------------- #
def _undirected(Q):
    adj = {v: set() for v in Q.vertices}
    for (s, t) in Q.arrows.values():
        if s != t:
            adj[s].add(t)
            adj[t].add(s)
    return adj


def _line_order(Q):
    """The path order of an ``A_n`` (n>=2) underlying graph as a list, walked from the
    lower-labelled degree-1 endpoint, or ``None`` if the graph is not a simple path."""
    V = list(Q.vertices)
    adj = _undirected(Q)
    if any(len(adj[v]) > 2 for v in V):
        return None
    endpoints = [v for v in V if len(adj[v]) == 1]
    if len(endpoints) != 2:
        return None
    try:
        start = min(endpoints)
    except TypeError:
        start = min(endpoints, key=repr)
    order, prev, cur = [start], None, start
    while True:
        nxts = [w for w in adj[cur] if w != prev]
        if not nxts:
            break
        prev, cur = cur, nxts[0]
        order.append(cur)
    if len(order) != len(V):            # not connected / not a simple path
        return None
    return order


def _is_forward_line(Q, order):
    """True iff every consecutive line edge is an arrow ``order[k] -> order[k+1]``."""
    fwd = {(s, t) for (s, t) in Q.arrows.values()}
    return all((order[k], order[k + 1]) in fwd for k in range(len(order) - 1))


# --------------------------------------------------------------------------- #
# The A_n / zigzag interval route (field-robust)
# --------------------------------------------------------------------------- #
def _interval_bars(M, order, kind, budget):
    n = len(order)
    pos = {v: i + 1 for i, v in enumerate(order)}
    bars = []
    for (I, mult) in decompose(M, budget=budget):
        dv = I.dimension_vector()
        supp = sorted(pos[v] for v, d in dv.items() if d)
        if not supp or supp != list(range(supp[0], supp[-1] + 1)):
            raise QuiverlabError(
                f"barcode: summand {I.name} of {M.name} has non-contiguous support "
                f"{supp} in the line order {tuple(order)} -- the input is not "
                "A_n-shaped",
                hint="barcode consumes a representation of a type-A line A_n")
        birth, death = supp[0], supp[-1]
        essential = (kind == "persistence" and death == n)
        bars.append(Bar(birth, death, mult, essential,
                        {v: int(d) for v, d in dv.items() if d}))
    bars.sort(key=lambda b: (b.birth, b.death))
    return tuple(bars)


def _barcode_by_ranks(M):
    """The classical forward-persistence rank formula (a genuinely independent second
    route -- only ``mat_rank`` of composite arrow matrices, no ``decompose``/``Hom``):

        mult[b,d] = (r_b^d - r_{b-1}^d) - (r_b^{d+1} - r_{b-1}^{d+1}),

    with ``r_i^j = rank(M_i -> M_j)`` the composite arrow-action rank (``r_i^i = dim
    M_i``; out of range = 0). Forward ``A_n`` ONLY -- refuses on a non-forward/zigzag
    line (the naive formula does not apply there). Returns ``[(birth, death,
    multiplicity), ...]`` for the positive multiplicities."""
    A = M.algebra
    Q = A.quiver
    if Q is None:
        raise QuiverlabError("_barcode_by_ranks needs a quiver-presented algebra")
    order = _line_order(Q)
    if order is None or not _is_forward_line(Q, order):
        raise QuiverlabError(
            "_barcode_by_ranks applies only to a FORWARD A_n line (1->2->...->n); the "
            "given module is not forward-persistence",
            hint="use barcode(M) for zigzag / commutative-ladder inputs")
    n = len(order)
    dom = M.domain
    arrow_of = {(s, t): lab for lab, (s, t) in Q.arrows.items()}
    dimvec = M.dimension_vector()
    # composite action matrices r_i^j (1-based i<=j), as the product of consecutive
    # forward arrow actions (right-module: (m*a)*b = action[b] @ action[a] @ m).
    comp = {}                                        # (i,j) -> matrix M_i -> M_j
    ranks = {}
    for i in range(1, n + 1):
        cur = None
        for j in range(i, n + 1):
            if j == i:
                ranks[(i, j)] = int(dimvec[order[i - 1]])
                continue
            Aa = M.action[arrow_of[(order[j - 2], order[j - 1])]]
            cur = Aa if cur is None else lm.matmul(Aa, cur, dom)
            ranks[(i, j)] = lm.mat_rank(cur, dom)

    def r(i, j):
        if i < 1 or j > n or i > j:
            return 0
        return ranks[(i, j)]

    out = []
    for b in range(1, n + 1):
        for d in range(b, n + 1):
            mult = (r(b, d) - r(b - 1, d)) - (r(b, d + 1) - r(b - 1, d + 1))
            if mult:
                out.append((b, d, int(mult)))
    return out


# --------------------------------------------------------------------------- #
# The commutative-ladder generalized-diagram route (char-scoped)
# --------------------------------------------------------------------------- #
def _generalized_diagram(A, M, n, budget, budget_modules):
    from quiverlab.modules.ar import knit_ar_quiver
    from quiverlab.modules.hom import is_isomorphic

    if n >= 5:                                       # Escolar-Hiraoka: rep-infinite
        raise QuiverlabError(
            f"barcode: CL({n}) is representation-infinite for n >= 5 (Escolar-Hiraoka); "
            "there is no complete AR-indexed generalized persistence diagram",
            hint="use a commutative ladder CL(n) with n <= 4")
    ar = knit_ar_quiver(A, budget_modules=max(budget_modules, 400), budget_dim=8000)
    if not ar.is_complete:
        # M6: over char <= dim the knit's internal decompose of a NON-brick raises;
        # the knit CATCHES it and returns is_complete=False / status="error" (NOT a
        # silent truncated "complete"). barcode refuses HERE, keyed off is_complete.
        raise QuiverlabError(
            f"barcode: the AR knit of this commutative ladder did not complete "
            f"(status={ar.status!r}) -- CL indecomposables need not be bricks, so the "
            f"AR-indexed diagram is certified only over char 0 or char > dim "
            f"(field {M.domain}); no complete diagram exists here",
            hint="recompute over QQ (or a characteristic > dim A); the A_n/zigzag "
                 "barcode route is field-robust, the commutative-ladder route is not")
    entries = []
    for (S, mult) in decompose(M, budget=budget):
        dv = {v: int(d) for v, d in S.dimension_vector().items() if d}
        match = None
        for i, vtx in enumerate(ar.vertices):
            if vtx["dimvec"] != S.dimension_vector():
                continue                             # dim-vector prefilter
            if is_isomorphic(S, vtx["module"]):
                match = (i, vtx.get("name") or f"X{i}")
                break
        if match is None:
            raise QuiverlabError(
                f"barcode: a summand of {M.name} (dim-vector {dv}) matched no AR-quiver "
                "vertex, yet the knit is complete -- every indecomposable summand MUST "
                "be a knit vertex; this is a bug",
                hint="report this with the module")
        is_interval = all(d == 1 for d in dv.values())    # thin support (0/1 dim vector)
        entries.append(DiagramEntry(match[0], match[1], dv, mult, is_interval))
    entries.sort(key=lambda e: (e.ar_index, e.multiplicity))
    return tuple(entries), ar.status


def _base_index(name):
    """The base index ``i`` of a CommutativeLadder vertex name ``"i_j"`` (else None)."""
    try:
        return int(str(name).split("_")[0])
    except (ValueError, IndexError):
        return None


# --------------------------------------------------------------------------- #
# Public surface
# --------------------------------------------------------------------------- #
def barcode(M, *, budget=512, budget_modules=256):
    """The barcode / generalized persistence diagram of a persistence module ``M``.

    ``M`` is a representation of an ``A_n`` line (any orientation) OR a ``CL(n<=4)``:

    * ``A_n`` forward -> ``kind="persistence"``, alternating -> ``kind="zigzag"``:
      the interval decomposition (Gabriel / Botnan-Crawley-Boevey), **field-robust**
      over every exact domain (interval modules are bricks).
    * ``CL(n<=4)`` -> ``kind="commutative_ladder"``: the AR-quiver-indexed generalized
      persistence diagram (Escolar-Hiraoka), char-scoped (char 0 / char > dim).

    The single-vertex ``A_1`` case is HANDLED (one bar ``[1,1]``, multiplicity =
    dim M -- M1), not refused. Refuses loudly: a non-``A_n``/non-CL quiver, ``CL(n>=5)``,
    a presentation-less algebra, or a char-undecidable ``decompose``/``knit``."""
    A = M.algebra
    Q = getattr(A, "quiver", None)
    if Q is None:
        raise QuiverlabError(
            "barcode needs a quiver-presented algebra (kQ/I); the given algebra is "
            "presentation-less (structure constants)",
            hint="barcode is defined for representations of A_n lines and commutative "
                 "ladders CL(n)")
    field = str(M.domain)
    V = list(Q.vertices)

    # (a0) n == 1 special case (M1): a single vertex, no arrows.
    if len(V) == 1 and not Q.arrows:
        dim_m = M.dim
        bars = ()
        if dim_m:
            bars = (Bar(1, 1, dim_m, True, {V[0]: 1}),)
        return Barcode("persistence", (V[0],), bars, (), 1, field, True, "", _REFS)

    # (a) underlying graph is a path A_n (n >= 2)?
    order = _line_order(Q)
    if order is not None:
        kind = "persistence" if _is_forward_line(Q, order) else "zigzag"
        bars = _interval_bars(M, order, kind, budget)
        return Barcode(kind, tuple(order), bars, (), len(order), field, True, "", _REFS)

    # (b) a commutative ladder CL(n)?
    from quiverlab.families.commutative_ladder import is_commutative_ladder
    ok, n = is_commutative_ladder(A)
    if ok:
        diagram, ar_status = _generalized_diagram(A, M, n, budget, budget_modules)
        base_bars = tuple(b for b in _cl_interval_bars(diagram))
        return Barcode("commutative_ladder", tuple(V), base_bars, diagram, n, field,
                       False, ar_status, _REFS)

    # (c) neither a line nor a ladder -> loud refusal.
    raise QuiverlabError(
        "barcode: the algebra is neither a type-A persistence line A_n nor a "
        "commutative ladder CL(n) -- a barcode/persistence diagram is undefined",
        hint="draw a line quiver 1-2-...-n (any orientation) or use "
             "CommutativeLadder(n) for n <= 4")


def barcode_block(A, M):
    """The core data of the ``barcode`` module compute kind as a plain JSON-able dict
    (Task 5), SHARED by both runners (``hpc/spec.py`` and ``docs/gui/runner.py``) so the
    two dispatches cannot drift. ``references``/``citations`` are added by each runner.
    A refusal (non-A_n/non-CL, CL(n>=5), presentation-less, char-undecidable) is an
    ``{"error": msg}`` block -- never a raise, never a 500. Dimvec keys are stringified
    for byte-stable JSON."""
    try:
        bc = barcode(M)
    except QuiverlabError as exc:
        return {"error": str(exc)}
    return {
        "kind": bc.kind, "n": bc.n, "field": bc.field, "field_robust": bc.field_robust,
        "bars": [{"birth": b.birth, "death": b.death, "multiplicity": b.multiplicity,
                  "essential": b.essential,
                  "dimvec": {str(k): int(v) for k, v in b.dimvec.items()}}
                 for b in bc.bars],
        "diagram": ([{"ar_name": e.ar_name,
                      "dimvec": {str(k): int(v) for k, v in e.dimvec.items()},
                      "multiplicity": e.multiplicity, "is_interval": e.is_interval}
                     for e in bc.diagram] if bc.kind == "commutative_ladder" else None),
        "ar_status": bc.ar_status or None,
    }


def _cl_interval_bars(diagram):
    """The interval sub-multiset of a CL diagram, exposed as Bars over the base index
    ``i`` parsed from the ``"i_j"`` vertex names (so a CL barcode still surfaces its
    interval bars). Non-interval summands carry no Bar -- they live only in the
    diagram. The diagram is the primary CL object; ``essential`` is always ``False``
    here (the CL primary object is the diagram, not a monotone-line barcode)."""
    bars = []
    for e in diagram:
        if not e.is_interval:
            continue
        idxs = [i for i in (_base_index(v) for v in e.dimvec) if i is not None]
        if len(idxs) == len(e.dimvec) and idxs:      # all names parsed as "i_j"
            birth, death = min(idxs), max(idxs)
        else:
            birth, death = 1, len(e.dimvec)
        bars.append(Bar(birth, death, e.multiplicity, False, dict(e.dimvec)))
    bars.sort(key=lambda b: (b.birth, b.death))
    return tuple(bars)
