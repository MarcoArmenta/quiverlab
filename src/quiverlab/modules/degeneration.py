"""Degeneration (closure) order for representation-finite algebras (Plan 49 / C8),
computed as the equivalent HOM order (Bongartz 1996 Adv. Math. 121; Zwara 2000
Compos. Math. 121):
  M <=_deg N  <=>  O_M subset closure(O_N)  <=>  dim Hom(M, X) >= dim Hom(N, X)
                   for every indecomposable X.
The generic module (open orbit) is the MAXIMUM; the semisimple module the MINIMUM.
Convention: Hom OUT of the class ([M,X]); degenerate = lower = bigger hom-dims.
Rep-finite only -- the finite indecomposable universe is P41 knit_ar_quiver; a
rep-infinite or self-injective input refuses loudly (never a partial poset)."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.geometry import group_dim
from quiverlab.modules.hom import hom_dim, identify_standard


@dataclass
class DegenerationPoset:
    vertices: list
    covers: list
    is_complete: bool
    status: str
    note: str = ""


def _std_name(std):
    return f"{std[0][0].upper()}_{std[1]}" if std else None


def _enumerate_classes(idimvecs, target, budget):
    """All multisets (multiplicity vectors over the indecomposables) whose
    dim vectors sum to `target`. Budget-capped; returns None on overflow."""
    n = len(target)
    m = len(idimvecs)
    out = []

    def dfs(idx, remaining, chosen):
        if len(out) > budget:
            return False
        if all(x == 0 for x in remaining):
            out.append(tuple(chosen))
            return True
        if idx >= m:
            return True
        dv = idimvecs[idx]
        if any(dv[k] > remaining[k] for k in range(n)):
            return dfs(idx + 1, remaining, chosen)
        maxc = min(remaining[k] // dv[k] for k in range(n) if dv[k] > 0)
        for c in range(maxc, -1, -1):
            newrem = tuple(remaining[k] - c * dv[k] for k in range(n))
            chosen[idx] = c
            if not dfs(idx + 1, newrem, chosen):
                chosen[idx] = 0
                return False
            chosen[idx] = 0
        return True

    ok = dfs(0, target, [0] * m)
    return out if (ok and len(out) <= budget) else None


def _hasse_covers(leq):
    """Cover pairs (i, j) with i < j (i.e. leq[i][j], i != j) and NO k strictly
    between (i < k < j)."""
    n = len(leq)
    covers = []
    for i in range(n):
        for j in range(n):
            if i == j or not leq[i][j]:
                continue
            if any(k not in (i, j) and leq[i][k] and leq[k][j] for k in range(n)):
                continue
            covers.append((i, j))
    return covers


def degeneration_order(algebra, d, *, budget=256, budget_dim=4096):
    """The degeneration (= hom) order of all iso-classes of dimension vector d over
    a representation-FINITE algebra. Returns a DegenerationPoset. STOPS LOUDLY at
    the budget / on self-injective (P41 knit unsupported) / on a representation-
    infinite input -- never a silent partial poset (`is_complete=False`, `status`
    mirrors the ARQuiver contract).

    `budget` caps the AR indecomposable universe (`ar_quiver(budget_modules=...)`)
    and the iso-class multiset enumeration; `budget_dim` caps the knit by module
    dimension. A HEREDITARY algebra is representation-finite iff it is of Dynkin
    (finite) type (P38 `form_type`), so a hereditary non-Dynkin input is refused
    IMMEDIATELY -- before the (for tame algebras, pathologically slow) knit."""
    from quiverlab.modules.morphism import direct_sum
    verts = list(algebra.quiver.vertices)
    if isinstance(d, dict):
        dvec = {v: int(d[v]) for v in verts}
    else:
        dvec = {v: int(d[i]) for i, v in enumerate(verts)}
    target = tuple(dvec[v] for v in verts)
    # Fast representation-type pre-check (hereditary case): a path algebra kQ is
    # rep-finite iff Q is Dynkin. This refuses the tame/wild hereditary inputs
    # (e.g. the Kronecker) instantly rather than driving the knit to a budget cap.
    try:
        hereditary = algebra.is_hereditary()
    except QuiverlabError:
        hereditary = False
    if hereditary and algebra.form_type() != "finite":
        return DegenerationPoset(
            [], [], is_complete=False, status="unsupported",
            note=f"hereditary of {algebra.dynkin_type()} type: representation-infinite "
                 "(infinitely many indecomposables) -- no finite degeneration poset "
                 "for a fixed dimension vector; raising the budget cannot help")
    # devil's-advocate fix (2026-08-05): the guard above covers ONLY hereditary
    # input; a non-hereditary non-self-injective rep-infinite algebra (e.g.
    # k<x,y>/rad^2, dim 3) drove the knit into the default budget_dim=4096 --
    # a silent multi-minute hang. Cap the knit's DIMENSION by the target: no
    # indecomposable of a d-dimensional iso-class can exceed dim(d), so the
    # knit never needs modules past sum(d) -- any incompleteness at that cap
    # is refused loudly below, never retried at 4096.
    dim_cap = min(budget_dim, sum(target) + 1)
    ar = algebra.ar_quiver(budget_modules=budget, budget_dim=dim_cap)
    if not ar.is_complete and hereditary and algebra.form_type() == "finite":
        # Dynkin hereditary is PROVABLY rep-finite: the full-budget knit
        # terminates, and its indecomposables may legitimately exceed sum(d)
        # (they are discarded by the class enumeration below).
        ar = algebra.ar_quiver(budget_modules=budget, budget_dim=budget_dim)
    if not ar.is_complete and ar.status != "unsupported":
        return DegenerationPoset(
            [], [], is_complete=False, status="budget",
            note=f"knit incomplete at the target-derived dimension cap {dim_cap} "
                 "and the algebra is not certified representation-finite -- "
                 "refusing rather than driving an unbounded knit; pass an "
                 "explicit budget_dim to raise the cap deliberately")
    if not ar.is_complete:
        return DegenerationPoset([], [], is_complete=False, status=ar.status,
                                 note="not representation-finite / knit did not close")
    indecs = [vtx["module"] for vtx in ar.vertices]
    idimvecs = [tuple(int(vtx["dimvec"][v]) for v in verts) for vtx in ar.vertices]
    classes = _enumerate_classes(idimvecs, target, budget)
    if classes is None:
        return DegenerationPoset([], [], is_complete=False, status="budget",
                                 note="too many iso-classes of dimension d")
    # pairwise Hom matrix over the indecomposables (computed once)
    ni = len(indecs)
    H = [[hom_dim(indecs[i], indecs[j]) for j in range(ni)] for i in range(ni)]

    def homvec(a):
        return tuple(sum(a[i] * H[i][j] for i in range(ni)) for j in range(ni))

    homvecs = [homvec(a) for a in classes]
    if len(set(homvecs)) != len(homvecs):
        raise QuiverlabError("degeneration_order: two iso-classes share a hom-vector "
                             "(rep-finite => impossible; report the input)")
    n = len(classes)
    # a <=_deg b iff a has bigger-or-equal homs everywhere (a is the degeneration)
    leq = [[all(homvecs[a][j] >= homvecs[b][j] for j in range(ni))
            for b in range(n)] for a in range(n)]
    covers = _hasse_covers(leq)
    gd = group_dim(dvec)
    max_orbit = None
    vertices = []
    for a in range(n):
        end_a = sum(classes[a][i] * classes[a][j] * H[i][j]
                    for i in range(ni) for j in range(ni))
        orbit = gd - end_a
        max_orbit = orbit if max_orbit is None else max(max_orbit, orbit)
        summ = [(indecs[i], classes[a][i]) for i in range(ni) if classes[a][i]]
        mods = []
        for (mod, mult) in summ:
            mods.extend([mod] * mult)
        module = direct_sum(*mods)[0] if len(mods) > 1 else mods[0]
        vertices.append({
            "index": a,
            "dimvec": dict(dvec),
            "summands": [(_std_name(identify_standard(mod)), mult) for mod, mult in summ],
            "orbit_dim": orbit,
            "module": module,
        })
    for v in vertices:
        v["is_generic"] = (v["orbit_dim"] == max_orbit)
    return DegenerationPoset(vertices, covers, is_complete=True, status="complete")
