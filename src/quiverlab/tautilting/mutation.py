"""tau-tilting mutation + the exchange-graph BFS (Plan 45 / C4, AIR 2014 Thm 2.30 /
Sec 3). Every support tau-tilting pair has exactly n neighbours (n-regular exchange
graph); the BFS from (A, 0) reaches them ALL iff the algebra is tau-tilting-finite (the
exchange graph is connected -- AIR Cor 2.38), else it hits the budget (loud,
status='budget'). Pairs dedup by g_key (AIR: g-vectors determine the pair).

The exchange at summand k is the uniform 2-term silting mutation (``_twoterm``): the cone
of the minimal left add(U)-approximation OR the cocone of the right one (U = the pair minus
summand k) -- whichever lands 2-term is THE neighbour, ARBITRATED by validation: the unique
candidate that make_pair-validates, differs from the pair, and swaps EXACTLY the k-th
g-column. Every mutation self-certifies (validated neighbour, one-column swap, involution).
Runs over QQ by default (the decompose / is_isomorphic char caveat -- rigorous over char 0
or char > dim)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.tautilting import _twoterm as tt
from quiverlab.tautilting.pairs import initial_pair, make_pair
from quiverlab.tautilting.rigid import g_columns


def mutate(pair, k):
    """The AIR mutation at the ``k``-th exchangeable summand (0..n-1: module summands in
    order, then support vertices in vertex order -- the g-matrix column order). Returns the
    UNIQUE other support tau-tilting pair sharing the almost-complete pair (``pair`` minus
    summand ``k``). Self-certified: the result validates (make_pair), its g-matrix swaps
    exactly the k-th column, and ``mutate(mutate(pair, k), k')`` recovers ``pair``."""
    A = pair.algebra
    verts = list(A.quiver.vertices)
    n = len(verts)
    if not (0 <= k < n):
        raise QuiverlabError(f"mutate: summand index {k} out of range 0..{n-1}")
    comps = tt.to_complexes(pair)
    Xk = comps[k]
    U = comps[:k] + comps[k + 1:]
    base_mods, base_supp = _pair_minus(pair, k)
    target_col = tuple(g_columns(pair)[k])
    pair_key = pair.g_key()
    for new_summand in (tt.left_mutation_summand(Xk, U),
                        tt.right_mutation_summand(Xk, U)):
        cls = tt.summand_class(new_summand)
        if cls[0] == "bad":
            continue
        if cls[0] == "module":
            mods = list(base_mods) + [cls[1]]
            supp = base_supp
        else:                                            # ("support", v)
            mods = list(base_mods)
            supp = base_supp | {cls[1]}
        try:
            cand = make_pair(A, mods, supp, check=True)
        except QuiverlabError:
            continue
        cand_key = cand.g_key()
        if (cand_key != pair_key and (pair_key - cand_key) == {target_col}
                and len(cand_key - pair_key) == 1):
            return cand
    raise QuiverlabError(
        f"mutate: no valid exchange found at summand {k} -- the 2-term silting cone/cocone "
        "produced no validating neighbour (report the pair + index)",
        hint="over char <= dim the decompose/is_isomorphic caveat can refuse; run over QQ")


def wall_normal(pair_i, pair_j):
    """The brick dim-vector labelling the exchange edge ``pair_i -- pair_j`` (King wall
    normal): the primitive non-negative integer vector orthogonal to the ``n-1`` shared
    g-columns. Exact (sympy nullspace, no floats). Returns a vertex-keyed dict, or None if
    the pairs are not adjacent."""
    import sympy
    verts = list(pair_i.algebra.quiver.vertices)
    shared = sorted(pair_i.g_key() & pair_j.g_key())
    if len(shared) != len(verts) - 1:
        return None
    if not shared:                                    # n = 1: the wall normal is e_1
        return {verts[0]: 1}
    ns = sympy.Matrix([list(c) for c in shared]).nullspace()
    if len(ns) != 1:
        return None
    vec = ns[0]
    lcm = 1
    for x in vec:
        lcm = sympy.ilcm(lcm, sympy.Rational(x).q)
    ivec = [int(sympy.Rational(x) * lcm) for x in vec]
    g = 0
    for x in ivec:
        g = sympy.igcd(g, x)
    if g:
        ivec = [x // g for x in ivec]
    if sum(ivec) < 0:                                 # brick dim-vectors are non-negative
        ivec = [-x for x in ivec]
    return {verts[i]: int(ivec[i]) for i in range(len(verts))}


def _brick_label(pair_i, pair_j):
    """The edge's brick dim-vector (wall normal) + name placeholder. The brick MODULE is
    identified later by :mod:`quiverlab.tautilting.torsion` (needs the module universe)."""
    return {"dimvec": wall_normal(pair_i, pair_j), "name": None}


def _pair_minus(pair, k):
    """The almost-complete pair (``pair`` with summand ``k`` removed): the remaining module
    summands (tuple) and support (frozenset)."""
    verts = list(pair.algebra.quiver.vertices)
    nmods = len(pair.summands)
    if k < nmods:
        mods = tuple(pair.summands[:k] + pair.summands[k + 1:])
        return mods, pair.support
    sup_sorted = sorted(pair.support, key=verts.index)
    v_removed = sup_sorted[k - nmods]
    return tuple(pair.summands), frozenset(pair.support - {v_removed})


# --------------------------------------------------------------------------- #
# the exchange graph (BFS from (A, 0)), mirroring P41's ARQuiver contract
# --------------------------------------------------------------------------- #
class ExchangeGraph:
    """The support tau-tilting exchange graph (or the loudly budget-capped prefix of a
    tau-tilting-infinite one). ``vertices`` is a list of dicts
    ``{"pair", "g_matrix", "label", "is_initial", "support", "summand_dimvecs"}``;
    ``arrows`` maps an undirected edge ``(i, j)`` (i<j) to
    ``{"brick": dimvec|None, "brick_name": str|None}`` (Task D fills the brick); ``adj[i]``
    is the neighbour set; ``is_complete`` is True iff the BFS closed (tau-tilting-finite),
    False iff budget-capped; ``status`` is ``"complete" | "budget" | "error"``; ``n_regular``
    is True iff every discovered vertex has exactly n neighbours (checked when complete)."""

    def __init__(self, vertices, arrows, adj, is_complete, status, n_regular):
        self.vertices = vertices
        self.arrows = arrows
        self.adj = adj
        self.is_complete = is_complete
        self.status = status
        self.n_regular = n_regular

    def __repr__(self):
        return ("ExchangeGraph(%d pairs, %d edges, status=%s, n_regular=%s)"
                % (len(self.vertices), len(self.arrows), self.status, self.n_regular))


def _label(pair):
    """A human label ``(P1(+)S1, {2})`` for a pair: standard names where identifiable,
    else the dim-vector. The support set follows."""
    from quiverlab.modules.hom import identify_standard
    parts = []
    for M in pair.summands:
        std = identify_standard(M)
        if std is not None:
            kind, v = std
            sym = {"simple": "S", "projective": "P", "injective": "I"}[kind]
            parts.append(f"{sym}{v}")
        else:
            dv = M.dimension_vector()
            parts.append("[" + ",".join(str(dv[w]) for w in sorted(dv)) + "]")
    mod = "(+)".join(parts) if parts else "0"
    supp = "{" + ",".join(str(v) for v in sorted(pair.support)) + "}"
    return f"({mod},{supp})"


def _record(pair, is_initial):
    return {
        "pair": pair,
        "g_matrix": pair.g_matrix(),
        "label": _label(pair),
        "is_initial": is_initial,
        "support": sorted(pair.support),
        "summand_dimvecs": [dict(M.dimension_vector()) for M in pair.summands],
    }


def exchange_graph(A, budget_pairs=512):
    """BFS the exchange graph from ``initial_pair(A)``, dedup by ``g_key``. Closes
    (``status="complete"``) iff ``A`` is tau-tilting-finite; else STOPS LOUDLY at the budget
    (``status="budget"``, ``is_complete=False``) -- never a silent partial graph. A summand
    the exchange cannot certify sets ``status="error"`` (never a silent skip)."""
    verts = list(A.quiver.vertices)
    n = len(verts)
    start = initial_pair(A)
    records = [_record(start, True)]
    index = {start.g_key(): 0}
    adj = {0: set()}
    arrows = {}
    frontier = [0]
    status, complete = "complete", True
    while frontier:
        i = frontier.pop(0)                 # breadth-first: keeps discovered modules small
        pi = records[i]["pair"]
        for k in range(n):
            try:
                pj = mutate(pi, k)
            except QuiverlabError:
                complete = False
                status = "error"
                continue
            key = pj.g_key()
            j = index.get(key)
            if j is None:
                if len(records) >= budget_pairs:
                    return ExchangeGraph(records, arrows, adj, is_complete=False,
                                         status="budget", n_regular=False)
                j = len(records)
                index[key] = j
                records.append(_record(pj, False))
                adj[j] = set()
                frontier.append(j)
            adj[i].add(j)
            adj[j].add(i)
            e = (min(i, j), max(i, j))
            if e not in arrows:
                lab = _brick_label(pi, pj)
                arrows[e] = {"brick": lab["dimvec"], "brick_name": lab["name"]}
    n_regular = complete and all(len(adj[i]) == n for i in range(len(records)))
    return ExchangeGraph(records, arrows, adj, is_complete=complete,
                         status=status, n_regular=n_regular)
