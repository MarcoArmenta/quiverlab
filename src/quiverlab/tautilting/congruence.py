"""The lattice theory of torsion classes (Plan 64 / R26): the finite lattice ``tors A`` as
an abstract lattice, its congruence lattice ``Con(tors A)``, the forcing order on bricks, the
canonical join representations, and the wide-subcategory poset via Enomoto's core label order.

A thin COMBINATORIAL layer over the merged Plan-45 tau-tilting engine -- NO new
representation theory. Plan 45 ships the ingredients read-only: the oriented, brick-labelled
exchange graph (= the Hasse quiver of ``tors A``, AIR Thm 2.30 / DIRRT), ``bricks`` /
``semibricks`` / ``_edge_brick`` (the iso-class-disambiguated DIRRT brick label, non-thin
safe). P64 assembles them into a finite lattice and computes:

- ``TorsionLattice`` -- order (transitive closure of the oriented Hasse), join/meet
  (``is_lattice`` self-cert), join-/meet-irreducibles, the canonical join representations
  (down-covers + their semibrick labels), and the SD / modular / distributive flags.
- ``CongruenceLattice`` -- ``Con(tors A)`` via the standard principal-cover-congruence
  algorithm (union-find closed under the join/meet compatibility to a fixed point), its
  distinct join-irreducible congruences (= the bricks, DIRRT), the forcing order on bricks,
  and ``|Con(tors A)|`` = the forcing poset's order-ideal count. Distributive by
  Funayama-Nakayama (structural).
- ``WideSubcategoryPoset`` -- Enomoto's kappa order AND core label order on ``L``, asserted
  to coincide (Enomoto's theorem, a real self-cert), isomorphic to ``(wide A, subseteq)``.

**Honest semi-decision (STRICTER than a bounded region).** Every P64 invariant is a GLOBAL
function of the whole finite lattice; on a budget-truncated exchange graph the lattice is a
prefix, not a lattice, so a partial ``Con`` / forcing order / ``#wide`` would be a LIE. P64 is
certified **iff ``A`` is tau-tilting-finite** (``eg.is_complete``, decided by the P45 BFS --
DIJ); on an incomplete graph it returns ``is_complete=False``, ``status="budget"``, a ``note``,
and OMITS every lattice invariant (never a partial-lattice value). The 2-Kronecker
(tau-tilting-infinite) is the honest-refusal oracle.

**Scope (inherited from Plan 45).** The engine runs over QQ by default; ``bricks`` /
``semibricks`` / the exchange-graph BFS are rigorous over char 0 or char > dim (Dickson/CIW),
and a brick decides ``End(B) = k`` over the algebraically-closed / char-0 base (the GF(p^n)
proper-division-ring caveat is honest-scope). Off scope the underlying engine refuses loudly.

Float-free: the lattice is finite combinatorics -- ids are ``int``, orders/partitions are
``frozenset``, counts are ``int``. All refusals are ``QuiverlabError``.

References: dirrt_lattice_torsion (DIRRT, Trans. AMS B 10 (2023)), barnard_carroll_zhu (BCZ,
Alg. Combin. 2 (2019)), enomoto_wide_ice (Enomoto, arXiv:2201.00595), marks_stovicek
(Marks-Stovicek, Bull. LMS 49 (2017)); air_tau_tilting / demonet_iyama_jasso (Plan 45)."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError

# The O(N^3) lattice-property sweep (SD / modular / distributive) and the O(2^m) order-ideal
# enumeration are exact and cheap for tau-tilting-finite small inputs; beyond these caps the
# corresponding invariant is reported as None + a note (never a hang / silent value).
_PROP_CAP = 400          # lattice-property triple sweep cap (element count)
_IDEAL_CAP = 22          # order-ideal (down-set) enumeration cap (#join-irreducible congruences)

_INCOMPLETE_NOTE = (
    "The exchange graph did not close (status={status}): A is tau-tilting-infinite or the "
    "pair budget ({budget}) was hit. tors A is not finite, so the lattice, its congruence "
    "lattice, the forcing order and the wide-subcategory poset are ALL undefined here (a "
    "partial value would be a lie, not merely incomplete). Certified only iff A is "
    "tau-tilting-finite (DIJ). Raise the budget for a genuinely finite algebra.")


# --------------------------------------------------------------------------- #
# value objects
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TorsionLattice:
    """The finite lattice ``tors A`` extracted from the P45 oriented exchange graph."""
    algebra: object
    elements: tuple              # exchange-graph vertex ids 0..N-1 (= torsion classes)
    labels: dict                 # id -> human label (P45 _label; e.g. "(P1(+)S1,{2})")
    covers: tuple                # ((upper_id, lower_id, {"brick_dimvec":{v:m},"brick_name":str|None}),...)
    order: dict                  # id -> frozenset(ids <= it)   (reflexive-transitive closure)
    top: int
    bottom: int
    join_irreducibles: tuple     # ids with a unique lower cover
    meet_irreducibles: tuple     # ids with a unique upper cover
    canonical_joins: dict        # id -> tuple of (lower_id, brick_dimvec, brick_name) (down-covers)
    is_lattice: bool
    is_semidistributive: bool
    is_distributive: bool
    is_modular: bool
    is_complete: bool
    status: str
    note: str


@dataclass(frozen=True)
class CongruenceLattice:
    """``Con(tors A)`` via principal cover-congruences + the forcing order on bricks."""
    algebra: object
    join_irreducible_congruences: tuple   # the distinct principal cover-congruences (partitions)
    forcing_order: tuple                  # ((brick_i_id, brick_j_id), ...) meaning brick_i <= brick_j
    forcing_bricks: tuple                 # {"dimvec":{v:m},"name":str|None} per J(Con) class (BCZ transport)
    size: int                             # |Con(tors A)| = #order-ideals of the forcing poset
    is_distributive: bool                 # always True (self-cert)
    is_complete: bool
    status: str
    note: str


@dataclass(frozen=True)
class WideSubcategoryPoset:
    """The wide-subcategory poset via Enomoto's core label order / kappa order."""
    algebra: object
    elements: tuple          # the distinct core-label representatives (one per wide subcategory)
    order: dict              # id -> frozenset(ids <= it) (inclusion order among them)
    labels: tuple            # per element: the brick-label set (its "simple objects"), aligned w/ elements
    core_labels: tuple       # per element: the CORE LABEL SET (frozenset of join-irreducible ids),
    #                          aligned w/ elements -- the ORDER's foundation: the poset TOP (mod A)
    #                          has the FULL brick set, the BOTTOM (0) the empty set (pins direction,
    #                          which the self-dual M3/NC(A3)/M4 label histograms alone cannot).
    size: int                # #wide
    constructions_agree: bool  # kappa order == core label order (Enomoto; self-cert)
    is_complete: bool
    status: str
    note: str


# --------------------------------------------------------------------------- #
# extraction: the oriented Hasse quiver -> lattice order + brick-labelled covers
# --------------------------------------------------------------------------- #
def _std_name(B):
    """The standard-module name ("S1"/"P1"/"I1") of a brick module, or None."""
    if B is None:
        return None
    from quiverlab.modules.hom import identify_standard
    std = identify_standard(B)
    if std is None:
        return None
    kind, v = std
    return {"simple": "S", "projective": "P", "injective": "I"}[kind] + str(v)


def _extract(A, budget):
    """Build the exchange graph and, if complete, the lattice skeleton: covers labelled by the
    ISO-CLASS brick (:func:`torsion._edge_brick`, non-thin safe), the down-set closure ``leq``
    and the up-set closure ``geq``. Returns ``(eg, None)`` when the graph did not close (drives
    the honest refusal), else ``(eg, (n, covers, leq, geq))``.

    Each cover is ``(upper, lower, brick_module, brick_dimvec, brick_name)`` -- ``upper`` has
    the larger ``Gen`` (AIR/DIRRT: mutation edges ARE the covers of ``tors A``)."""
    from quiverlab.tautilting.mutation import exchange_graph
    from quiverlab.tautilting.torsion import (_edge_brick, _torsion_universe,
                                              hasse_orientation)
    eg = exchange_graph(A, budget_pairs=budget)
    if not eg.is_complete:
        return eg, None
    orient = hasse_orientation(eg)
    n = len(eg.vertices)
    universe = _torsion_universe(A, budget=budget)
    covers = []
    for (i, j), d in orient.items():
        a, b = (i, j) if d == "down" else (j, i)      # a (upper) > b (lower)
        dv = eg.arrows[(i, j)]["brick"]
        B = _edge_brick(A, eg.vertices[i]["pair"], eg.vertices[j]["pair"], dv, universe)
        bdv = {str(k): int(v) for k, v in (B.dimension_vector().items() if B is not None
                                           else [])}
        covers.append((a, b, B, bdv, _std_name(B)))
    # Warshall-style down-closure over the covers: leq[e] = {f : f <= e}.
    leq = {e: {e} for e in range(n)}
    changed = True
    while changed:
        changed = False
        for (a, b, _B, _dv, _nm) in covers:
            if not leq[b] <= leq[a]:
                leq[a] |= leq[b]
                changed = True
    geq = {e: {f for f in range(n) if e in leq[f]} for e in range(n)}
    return eg, (n, covers, leq, geq)


def _lattice_ops(n, leq, geq):
    """Meet / join tables (``mt[x][y]`` / ``jn[x][y]``) by the greatest-lower-bound /
    least-upper-bound of the common down-/up-sets. Raises ``QuiverlabError`` if any pair lacks
    a UNIQUE meet or join -- ``tors A`` is a lattice (DIRRT), so a violation means the
    exchange graph / orientation was mis-extracted, NOT the mathematics."""
    mt = [[0] * n for _ in range(n)]
    jn = [[0] * n for _ in range(n)]
    for x in range(n):
        lx = leq[x]
        gx = geq[x]
        for y in range(n):
            common = lx & leq[y]
            cand = [z for z in common if common <= leq[z]]
            if len(cand) != 1:
                raise QuiverlabError(
                    f"tors A is not a lattice: pair ({x},{y}) has {len(cand)} maximal common "
                    "lower bounds (a DIRRT guarantee failed -- the extracted Hasse quiver / "
                    "orientation is wrong)", hint="report the algebra + this pair")
            mt[x][y] = cand[0]
            cu = gx & geq[y]
            cj = [z for z in cu if cu <= geq[z]]
            if len(cj) != 1:
                raise QuiverlabError(
                    f"tors A is not a lattice: pair ({x},{y}) has {len(cj)} minimal common "
                    "upper bounds (a DIRRT guarantee failed -- the extracted Hasse quiver / "
                    "orientation is wrong)", hint="report the algebra + this pair")
            jn[x][y] = cj[0]
    return mt, jn


def _lattice_properties(n, mt, jn):
    """The SD / modular / distributive flags via the standard identities over all triples
    (``O(N^3)``; capped at ``_PROP_CAP`` -- fine for tau-tilting-finite small ``N``; beyond
    the cap the three flags are ``None`` + the check is skipped, DIRRT guarantee SD so this is
    a self-cert, never load-bearing for correctness)."""
    if n > _PROP_CAP:
        return None, None, None
    sd_join = sd_meet = distributive = modular = True
    rng = range(n)
    for x in rng:
        for y in rng:
            for z in rng:
                # SD-join: x v y == x v z  =>  x v y == x v (y ^ z)
                if jn[x][y] == jn[x][z] and jn[x][y] != jn[x][mt[y][z]]:
                    sd_join = False
                # SD-meet: x ^ y == x ^ z  =>  x ^ y == x ^ (y v z)
                if mt[x][y] == mt[x][z] and mt[x][y] != mt[x][jn[y][z]]:
                    sd_meet = False
                # distributive: x ^ (y v z) == (x ^ y) v (x ^ z)
                if mt[x][jn[y][z]] != jn[mt[x][y]][mt[x][z]]:
                    distributive = False
                # modular: x <= z  =>  x v (y ^ z) == (x v y) ^ z
                if jn[x][z] == z and jn[x][mt[y][z]] != mt[jn[x][y]][z]:
                    modular = False
    return (sd_join and sd_meet), distributive, modular


# --------------------------------------------------------------------------- #
# congruences: principal cover-congruences + forcing order
# --------------------------------------------------------------------------- #
def _principal_congruence(a0, b0, n, mt, jn):
    """The smallest lattice congruence collapsing ``a0 == b0`` (the principal cover-congruence
    ``con(b0, a0)``): union-find closed under the compatibility rule (if ``x == y`` then
    ``x^z == y^z`` and ``x v z == y v z`` for all ``z``), iterated to a fixed point. Returns
    the partition as a ``frozenset`` of blocks."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
            return True
        return False

    union(a0, b0)
    changed = True
    while changed:
        changed = False
        reps = {}
        for e in range(n):
            reps.setdefault(find(e), []).append(e)
        for blk in reps.values():
            for i in range(len(blk)):
                for j in range(i + 1, len(blk)):
                    x, y = blk[i], blk[j]
                    for z in range(n):
                        if union(mt[x][z], mt[y][z]):
                            changed = True
                        if union(jn[x][z], jn[y][z]):
                            changed = True
    blocks = {}
    for e in range(n):
        blocks.setdefault(find(e), set()).add(e)
    return frozenset(frozenset(v) for v in blocks.values())


def _refines(ti, tj):
    """``ti <= tj`` in the partition-refinement order: every block of ``ti`` sits inside a
    block of ``tj`` (``ti`` is finer, ``con_i subseteq con_j`` as relations)."""
    return all(any(bi <= bj for bj in tj) for bi in ti)


def _order_ideals_count(pairs, m):
    """The number of order ideals (down-sets) of the poset on ``0..m-1`` given by the strict
    relations ``pairs`` (``(i, j)`` = ``i <= j``). Brute-force bitmask enumeration, capped at
    ``_IDEAL_CAP`` (returns ``None`` above the cap -- never a hang)."""
    if m > _IDEAL_CAP:
        return None
    belowmask = [0] * m
    for (i, j) in pairs:
        if i != j:
            belowmask[j] |= (1 << i)
    count = 0
    for s in range(1 << m):
        ok = True
        for e in range(m):
            if (s >> e) & 1 and (belowmask[e] & ~s):
                ok = False
                break
        count += ok
    return count


# --------------------------------------------------------------------------- #
# Enomoto: the kappa map, the extended kappa, the core label set
# --------------------------------------------------------------------------- #
def _kappa_map(n, join_irr, jstar, mt, leq):
    """The kappa map (Enomoto Def 2.8): ``kappa(j) = max{x : j ^ x == j_*}`` for each
    join-irreducible ``j`` (``j_*`` its unique lower cover). In a semidistributive lattice the
    maximum is unique (and meet-irreducible); a missing/ambiguous max raises loudly."""
    kappa = {}
    for j in join_irr:
        cand = [x for x in range(n) if mt[j][x] == jstar[j]]
        mx = [x for x in cand if all(y in leq[x] for y in cand)]
        if len(mx) != 1:
            raise QuiverlabError(
                f"kappa map: join-irreducible {j} has {len(mx)} maximal x with j^x=j_* "
                "(the lattice is not semidistributive as DIRRT requires)",
                hint="report the algebra")
        kappa[j] = mx[0]
    return kappa


def _canonical_joinands(w, lower_of_w, join_irr, jstar, mt, jn):
    """The canonical join representation of ``w`` as a set of join-irreducible ELEMENTS: for
    each lower cover ``b`` of ``w``, the unique join-irreducible ``j`` with ``j v b == w`` and
    ``j ^ b == j_*`` (the DIRRT/BCZ cover label). Raises if it is not unique."""
    js = []
    for b in lower_of_w:
        hit = [j for j in join_irr if jn[j][b] == w and mt[j][b] == jstar[j]]
        if len(hit) != 1:
            raise QuiverlabError(
                f"canonical join representation: cover ({b} <. {w}) matched {len(hit)} "
                "join-irreducibles (the lattice is not semidistributive)",
                hint="report the algebra")
        js.append(hit[0])
    return js


def _meet_reduce(xs, top, mt):
    """The meet of a set of elements (empty meet = ``top``)."""
    it = list(xs)
    if not it:
        return top
    acc = it[0]
    for z in it[1:]:
        acc = mt[acc][z]
    return acc


# --------------------------------------------------------------------------- #
# public builders
# --------------------------------------------------------------------------- #
def _lower_map(n, covers):
    """id -> list of (lower_id, brick_module, brick_dimvec, brick_name) for its down-covers."""
    low = {e: [] for e in range(n)}
    for (a, b, B, dv, nm) in covers:
        low[a].append((b, B, dv, nm))
    return low


def _upper_count(n, covers):
    """id -> number of up-covers (covers where the id is the LOWER endpoint)."""
    up = {e: 0 for e in range(n)}
    for (a, b, _B, _dv, _nm) in covers:
        up[b] += 1
    return up


def torsion_lattice(A, *, budget=512) -> TorsionLattice:
    """The finite lattice ``tors A`` (DIRRT / AIR): the order (transitive closure of the P45
    oriented Hasse quiver), join/meet (``is_lattice`` self-cert), join-/meet-irreducibles, the
    canonical join representations (down-covers + semibrick labels), and the SD / modular /
    distributive flags. Certified complete iff ``A`` is tau-tilting-finite; else an honest
    ``is_complete=False`` value with every lattice invariant omitted."""
    eg, data = _extract(A, budget)
    if data is None:
        return TorsionLattice(
            algebra=A, elements=(), labels={}, covers=(), order={}, top=-1, bottom=-1,
            join_irreducibles=(), meet_irreducibles=(), canonical_joins={},
            is_lattice=False, is_semidistributive=False, is_distributive=False,
            is_modular=False, is_complete=False, status=eg.status,
            note=_INCOMPLETE_NOTE.format(status=eg.status, budget=budget))
    n, covers, leq, geq = data
    mt, jn = _lattice_ops(n, leq, geq)
    is_sd, is_dist, is_mod = _lattice_properties(n, mt, jn)
    top = next(e for e in range(n) if len(leq[e]) == n)
    bottom = next(e for e in range(n) if len(leq[e]) == 1)
    low = _lower_map(n, covers)
    upc = _upper_count(n, covers)
    join_irr = tuple(e for e in range(n) if len(low[e]) == 1)
    meet_irr = tuple(e for e in range(n) if upc[e] == 1)
    canonical_joins = {e: tuple((b, dv, nm) for (b, _B, dv, nm) in low[e])
                       for e in range(n)}
    cover_tuples = tuple((a, b, {"brick_dimvec": dv, "brick_name": nm})
                         for (a, b, _B, dv, nm) in covers)
    return TorsionLattice(
        algebra=A, elements=tuple(range(n)),
        labels={i: eg.vertices[i]["label"] for i in range(n)},
        covers=cover_tuples, order={e: frozenset(leq[e]) for e in range(n)},
        top=top, bottom=bottom, join_irreducibles=join_irr, meet_irreducibles=meet_irr,
        canonical_joins=canonical_joins, is_lattice=True,
        is_semidistributive=bool(is_sd) if is_sd is not None else None,
        is_distributive=bool(is_dist) if is_dist is not None else None,
        is_modular=bool(is_mod) if is_mod is not None else None,
        is_complete=True, status="complete", note="")


def congruence_lattice(A, *, budget=512) -> CongruenceLattice:
    """``Con(tors A)`` via principal cover-congruences (Grätzer / Funayama-Nakayama +
    Birkhoff): the distinct principal congruences are the join-irreducible congruences (= the
    bricks, DIRRT), ordered by refinement they form the forcing order on bricks, and
    ``|Con(tors A)|`` = the number of order ideals of that poset. Distributive by construction.
    Self-cert gates: ``#J(Con) == #bricks`` and ``#join-irreducibles == #bricks`` (loud on
    violation). Certified complete iff ``A`` is tau-tilting-finite."""
    from quiverlab.tautilting.torsion import bricks
    eg, data = _extract(A, budget)
    if data is None:
        return CongruenceLattice(
            algebra=A, join_irreducible_congruences=(), forcing_order=(),
            forcing_bricks=(), size=None, is_distributive=True, is_complete=False,
            status=eg.status, note=_INCOMPLETE_NOTE.format(status=eg.status, budget=budget))
    n, covers, leq, geq = data
    mt, jn = _lattice_ops(n, leq, geq)
    low = _lower_map(n, covers)
    n_join_irr = sum(1 for e in range(n) if len(low[e]) == 1)
    # principal cover-congruence per cover, deduped by partition equality -> J(Con).
    jcon = []                                 # (theta, brick_module, brick_dimvec, brick_name)
    for (a, b, B, dv, nm) in covers:
        theta = _principal_congruence(a, b, n, mt, jn)
        if not any(t == theta for (t, _B, _dv, _nm) in jcon):
            jcon.append((theta, B, dv, nm))
    n_bricks = len(bricks(A, budget=budget))
    if len(jcon) != n_bricks:
        raise QuiverlabError(
            f"Con(tors A): {len(jcon)} distinct principal congruences but {n_bricks} bricks "
            "-- the forcing classes must be exactly the bricks (DIRRT); the extraction is "
            "wrong", hint="report the algebra")
    if n_join_irr != n_bricks:
        raise QuiverlabError(
            f"tors A: {n_join_irr} join-irreducibles but {n_bricks} bricks -- BCZ requires "
            "them in bijection; the extraction is wrong", hint="report the algebra")
    thetas = [t for (t, _B, _dv, _nm) in jcon]
    forcing = tuple((i, j) for i in range(len(jcon)) for j in range(len(jcon))
                    if i != j and _refines(thetas[i], thetas[j]))
    size = _order_ideals_count(forcing, len(jcon))
    forcing_bricks = tuple({"dimvec": dv, "name": nm} for (_t, _B, dv, nm) in jcon)
    note = ("" if size is not None else
            f"|Con| omitted: {len(jcon)} join-irreducible congruences exceed the "
            f"order-ideal enumeration cap ({_IDEAL_CAP}).")
    return CongruenceLattice(
        algebra=A, join_irreducible_congruences=tuple(thetas), forcing_order=forcing,
        forcing_bricks=forcing_bricks, size=size, is_distributive=True, is_complete=True,
        status="complete", note=note)


def _wide_data(A, budget):
    """Shared core for the wide-subcategory poset: kappa map, extended kappa, core label sets,
    and BOTH poset structures (kappa order + core label order) on the finite lattice. Returns
    ``(eg, None)`` when incomplete, else ``(eg, dict)`` with the pieces the poset needs."""
    eg, data = _extract(A, budget)
    if data is None:
        return eg, None
    n, covers, leq, geq = data
    mt, jn = _lattice_ops(n, leq, geq)
    top = next(e for e in range(n) if len(leq[e]) == n)
    low = _lower_map(n, covers)
    join_irr = [e for e in range(n) if len(low[e]) == 1]
    jstar = {e: low[e][0][0] for e in join_irr}
    kappa = _kappa_map(n, join_irr, jstar, mt, leq)
    # extended kappa: kbar(w) = /\{kappa(j) : j in canonical join rep of w}.
    kbar = {}
    for w in range(n):
        cj = _canonical_joinands(w, [b for (b, _B, _dv, _nm) in low[w]],
                                 join_irr, jstar, mt, jn)
        kbar[w] = _meet_reduce((kappa[j] for j in cj), top, mt)
    # core label set: CLS(w) = {j in JI : j <= w and kappa(j) >= w_down}, w_down = /\ lower covers.
    wdown = {w: _meet_reduce((b for (b, _B, _dv, _nm) in low[w]), top, mt) for w in range(n)}
    cls = {}
    for w in range(n):
        cls[w] = frozenset(j for j in join_irr
                           if j in leq[w] and wdown[w] in leq[kappa[j]])
    # the semibrick (simple objects of the wide subcategory) = the down-cover brick labels.
    labels = {w: tuple({"dimvec": dv, "name": nm} for (_b, _B, dv, nm) in low[w])
              for w in range(n)}

    def kle(x, y):                       # kappa order (Enomoto): x <= y and kbar(x) >= kbar(y)
        return x in leq[y] and kbar[y] in leq[kbar[x]]

    def cle(x, y):                       # core label order: CLS(x) subseteq CLS(y)
        return cls[x] <= cls[y]

    return eg, {"n": n, "cls": cls, "labels": labels, "kle": kle, "cle": cle}


def wide_subcategories(A, *, budget=512) -> WideSubcategoryPoset:
    """The wide-subcategory poset ``(wide A, subseteq)`` via Enomoto (2201.00595): the kappa
    order (extended kappa map of Barnard-Todorov-Zhu) and the core label order, computed BOTH
    and asserted to coincide (Enomoto's theorem -- a real self-cert, loud on mismatch). The
    poset elements are the DISTINCT core label sets (so ``#wide`` can be ``< #torsion`` off the
    representation-finite case -- Marks-Stovicek); ``labels`` per element are its simple-object
    bricks. Certified complete iff ``A`` is tau-tilting-finite."""
    eg, data = _wide_data(A, budget)
    if data is None:
        return WideSubcategoryPoset(
            algebra=A, elements=(), order={}, labels=(), core_labels=(), size=None,
            constructions_agree=False, is_complete=False, status=eg.status,
            note=_INCOMPLETE_NOTE.format(status=eg.status, budget=budget))
    n, cls, labels = data["n"], data["cls"], data["labels"]
    kle, cle = data["kle"], data["cle"]
    # distinct core label sets = the wide subcategories (one representative element each).
    reps = {}
    for w in range(n):
        reps.setdefault(cls[w], w)
    nodes = sorted(reps.values())
    # Enomoto's theorem: the kappa order and the core label order coincide (on representatives).
    agree = all(kle(x, y) == cle(x, y) for x in nodes for y in nodes)
    if not agree:
        raise QuiverlabError(
            "wide subcategories: the kappa order and the core label order disagree "
            "(Enomoto's theorem failed -- a construction is mis-transcribed)",
            hint="report the algebra")
    order = {x: frozenset(y for y in nodes if cls[y] <= cls[x]) for x in nodes}
    return WideSubcategoryPoset(
        algebra=A, elements=tuple(nodes), order=order,
        labels=tuple(labels[x] for x in nodes),
        core_labels=tuple(cls[x] for x in nodes), size=len(nodes),
        constructions_agree=True, is_complete=True, status="complete", note="")


# --------------------------------------------------------------------------- #
# the shared no-code block (both runners route through this)
# --------------------------------------------------------------------------- #
_CITATIONS = ["dirrt_lattice_torsion", "barnard_carroll_zhu", "enomoto_wide_ice",
              "marks_stovicek", "air_tau_tilting", "demonet_iyama_jasso"]


def congruences_block(A, *, budget=512) -> dict:
    """The ``congruences`` compute-kind payload (Plan 64): the torsion-lattice summary + the
    canonical join representations, ``Con(tors A)`` + the forcing order on bricks, and the
    wide-subcategory poset. Certified complete iff ``A`` is tau-tilting-finite; on an
    incomplete exchange graph ``lattice = congruences = wide = None`` and a ``note`` (never a
    partial-lattice lie). Consumed byte-identically by both runners (each adds ``citations``)."""
    n = len(list(A.quiver.vertices))
    L = torsion_lattice(A, budget=budget)
    block = {
        "kind": "congruences", "n": n,
        "complete": L.is_complete,
        "status": L.status,
        "references": list(_CITATIONS),
        "note": None,
    }
    if not L.is_complete:
        block["lattice"] = None
        block["congruences"] = None
        block["wide"] = None
        block["note"] = L.note
        return block
    C = congruence_lattice(A, budget=budget)
    W = wide_subcategories(A, budget=budget)
    block["lattice"] = {
        "size": len(L.elements),
        "num_covers": len(L.covers),
        "is_semidistributive": L.is_semidistributive,
        "is_modular": L.is_modular,
        "is_distributive": L.is_distributive,
        "join_irreducibles": len(L.join_irreducibles),
        "meet_irreducibles": len(L.meet_irreducibles),
        "hasse": [{"from": a, "to": b, "brick_dimvec": meta["brick_dimvec"],
                   "brick_name": meta["brick_name"]} for (a, b, meta) in L.covers],
        "canonical_joins": [
            {"element": e, "label": L.labels[e],
             "joinands": [{"brick_dimvec": dv, "brick_name": nm}
                          for (_lo, dv, nm) in L.canonical_joins[e]]}
            for e in L.elements],
    }
    block["congruences"] = {
        "size": C.size,
        "num_join_irreducibles": len(C.join_irreducible_congruences),
        "is_distributive": C.is_distributive,
        "forcing_order": {
            "bricks": [{"dimvec": fb["dimvec"], "name": fb["name"]}
                       for fb in C.forcing_bricks],
            "relations": [[i, j] for (i, j) in C.forcing_order],
        },
    }
    idx = {e: k for k, e in enumerate(W.elements)}
    block["wide"] = {
        "size": W.size,
        "constructions_agree": W.constructions_agree,
        "relations": [[idx[x], idx[y]] for x in W.elements for y in W.elements
                      if x != y and x in W.order[y]],
        "labels": [[{"dimvec": lab["dimvec"], "name": lab["name"]} for lab in labs]
                   for labs in W.labels],
    }
    if C.note:
        block["note"] = C.note
    return block
