"""Geometry of representations (Plan 49 / C8): orbit dimensions in the
representation variety, Voigt rigidity, and the Kac canonical decomposition.
Right modules; the orbit/rigidity layer is exact over EVERY Domain
(dim End + Ext are exact); the canonical decomposition is a hereditary-Dynkin
notion with a loud refusal off scope. Float-free (int/dict)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules.hom import end_dim


def group_dim(dimvec):
    """dim GL(d) = sum_v d_v^2 -- the base-change group acting on Rep(Q, d)."""
    return sum(int(n) * int(n) for n in dimvec.values())


def representation_variety_dim(algebra, dimvec):
    """dim Rep(Q, d) = sum over arrows a: i -> j of d_i * d_j (the AMBIENT quiver
    representation variety). For kQ/I the module variety mod_A(d) is the closed
    subvariety cut by the relations; this is its ambient space -- stated honestly
    by every codim consumer."""
    total = 0
    for (s, t) in algebra.quiver.arrows.values():
        total += int(dimvec[s]) * int(dimvec[t])
    return total


def orbit_dimension(M):
    """dim of the GL(d)-orbit of M in Rep(Q, d):
        dim O_M = dim GL(d) - dim_k End_A(M) = sum_v d_v^2 - dim End_A(M).
    Aut(M) is Zariski-open in the affine space End_A(M), so
    dim Stab(M) = dim End_A(M). Holds verbatim for kQ/I."""
    return group_dim(M.dimension_vector()) - end_dim(M)


def is_rigid(M):
    """Voigt: M is rigid iff Ext^1_A(M, M) = 0, and then O_M is OPEN in the
    module variety (rigid => open orbit, in general). One Ext call."""
    return M.algebra.ext(M, M, 1) == 0


def rigidity_codim(M):
    """dim Ext^1_A(M, M). HONEST semantics (stated by every caller):
      * HEREDITARY A: EQUALS codim of the orbit closure in Rep(Q, d) -- Voigt,
        since Rep(Q, d) is smooth: codim O_M = dim End(M) - <d,d> = dim Ext^1(M,M).
      * general kQ/I: only an UPPER BOUND -- codim_{mod_A(d)} O_M <= dim Ext^1(M,M)
        (equality iff mod_A(d) is smooth at M). NEVER claim equality off hereditary.
    """
    return M.algebra.ext(M, M, 1)


def canonical_decomposition(algebra, d, *, budget=4096):
    """The Kac canonical decomposition of the dimension vector d over a HEREDITARY
    DYNKIN algebra:
        d = sum_i m_i * beta_i,   [{"root": beta_i, "multiplicity": m_i, "name": ...}, ...]
    into positive (real Schur) roots with ext(beta_i, beta_j) = 0 for all i, j in
    the support -- i.e. the generic module G = (+) M_{beta_i}^{m_i} is RIGID (Kac).
    Certified per instance (Ext^1(G, G) == 0).

    Scope (loud QuiverlabError otherwise):
      1. non-hereditary A: the canonical decomposition is a hereditary (path-algebra)
         notion -- refused (orbit_dimension / is_rigid still apply to kQ/I).
      2. hereditary but not Dynkin (Euclidean/wild): the general Schofield /
         Derksen-Weyman recursion (imaginary Schur roots, isotropic multiplicities)
         is DEFERRED to a named successor plan -- refused.
      3. Dynkin (finite type): full support, rigidity-certified."""
    from quiverlab.modules.hom import identify_standard
    from quiverlab.modules.morphism import direct_sum
    verts = list(algebra.quiver.vertices)
    if isinstance(d, dict):
        dvec = tuple(int(d[v]) for v in verts)
    else:
        dvec = tuple(int(x) for x in d)
    if not algebra.is_hereditary():
        raise QuiverlabError(
            "canonical_decomposition: A is not hereditary -- the Kac canonical "
            "decomposition is a hereditary (path-algebra) notion",
            hint="orbit_dimension and is_rigid still apply to kQ/I; the canonical "
                 "decomposition does not")
    if algebra.form_type() != "finite":
        raise QuiverlabError(
            f"canonical_decomposition: A is hereditary of {algebra.dynkin_type()} "
            "type (Euclidean/wild, not Dynkin/finite) -- the general Schofield / "
            "Derksen-Weyman recursion (imaginary Schur roots, isotropic multiplicities) "
            "is DEFERRED to a successor plan",
            hint="Dynkin (finite type) is fully supported here")
    roots = [tuple(int(x) for x in r) for r in algebra.positive_roots()]  # Dynkin-gated
    # one indecomposable per positive root (Gabriel): read them off the AR quiver
    ar = algebra.ar_quiver()
    if not ar.is_complete:
        raise QuiverlabError("canonical_decomposition: AR knitting did not close "
                             f"(status {ar.status})")
    mod_of = {}
    for vtx in ar.vertices:
        key = tuple(int(vtx["dimvec"][v]) for v in verts)
        mod_of[key] = vtx["module"]
    missing = [r for r in roots if r not in mod_of]
    if missing:
        raise QuiverlabError("canonical_decomposition: no indecomposable found for "
                             f"root(s) {missing} (AR/root mismatch)")
    # generic ext(beta, gamma) = actual Ext^1 between the unique indecomposables
    ext = {(a, b): algebra.ext(mod_of[a], mod_of[b], 1) for a in roots for b in roots}
    best = _search_canonical(roots, dvec, ext, budget)
    if best is None:
        raise QuiverlabError("canonical_decomposition: no rigid decomposition within "
                             "budget", hint="raise `budget`")
    # per-instance certificate: the direct sum of the chosen indecomposables is
    # RIGID (Kac). Build G explicitly and demand Ext^1(G, G) == 0 -- an independent
    # oracle on top of the search's pairwise-ext bookkeeping.
    mods = []
    for root, mult in best:
        mods.extend([mod_of[root]] * mult)
    G = direct_sum(*mods)[0] if len(mods) > 1 else mods[0]
    if algebra.ext(G, G, 1) != 0:
        raise QuiverlabError(
            "canonical_decomposition: rigidity certificate FAILED "
            f"(Ext^1(G, G) != 0 for the assembled generic module of {dvec})")
    out = []
    for root, mult in sorted(best):
        std = identify_standard(mod_of[root])
        name = f"{std[0][0].upper()}_{std[1]}" if std else None
        out.append({"root": root, "multiplicity": mult, "name": name})
    return out


def orbit_geometry_block(M):
    """The shared no-code GUI/report block for the ``orbit_geometry`` module kind
    (Plan 49 / C8): orbit dim, dim GL(d), dim Rep(Q,d), dim End, the Voigt rigidity
    verdict with the HONEST codimension semantics (=codim on hereditary, an upper
    bound on kQ/I), and -- on hereditary Dynkin -- the Kac canonical decomposition.
    Byte-identical across the HPC runner and the Pyodide twin (both import this ONE
    function). Orbit dim / rigidity / codim compute for ANY module over ANY kQ/I;
    the canonical decomposition is the conditional extra (else ``None`` + a note)."""
    A = M.algebra
    dv = M.dimension_vector()
    ext1 = A.ext(M, M, 1)
    hereditary = A.is_hereditary()
    block = {
        "kind": "orbit_geometry", "side": M.side,
        "dim_vector": {str(v): int(n) for v, n in dv.items()},
        "group_dim": group_dim(dv),                       # dim GL(d)
        "rep_variety_dim": representation_variety_dim(A, dv),  # dim Rep(Q,d) (ambient)
        "end_dim": end_dim(M),                            # dim End_A(M)
        "orbit_dim": orbit_dimension(M),                  # dim O_M
        "rigid": ext1 == 0,
        "ext1_self": ext1,                                # dim Ext^1(M,M) = rigidity_codim
        "hereditary": hereditary,
        "codim_semantics": "hereditary" if hereditary else "general",
        "latex": r"\dim \mathcal{O}_M = \sum_v d_v^2 - \dim_k \operatorname{End}_A(M)",
    }
    try:                                                  # Dynkin-only extra
        block["canonical_decomposition"] = canonical_decomposition(A, dv)
    except QuiverlabError as e:
        block["canonical_decomposition"] = None
        block["canonical_note"] = str(e)
    return block


def _search_canonical(roots, target, ext, budget):
    """DFS for the (unique, over Dynkin) multiset of roots summing to `target`
    whose direct sum is RIGID: for every pair beta, gamma in the chosen support
    ext(beta, gamma) = ext(gamma, beta) = 0 (real roots have ext(beta,beta)=0).
    Since dim Ext^1(G, G) = sum_{b,g} m_b m_g ext(b, g) >= 0, this is exactly
    Ext^1(G, G) = 0 (the Kac criterion)."""
    n = len(target)
    order = sorted(roots, key=lambda r: -sum(r))       # bigger roots first
    steps = {"n": 0}

    def compatible(a, b):
        return ext[(a, b)] == 0 and ext[(b, a)] == 0

    def dfs(idx, remaining, chosen):
        steps["n"] += 1
        if steps["n"] > budget:
            raise QuiverlabError("canonical_decomposition: search budget exceeded",
                                 hint="raise `budget`")
        if all(x == 0 for x in remaining):
            return list(chosen.items())
        if idx >= len(order):
            return None
        beta = order[idx]
        if any(beta[k] > remaining[k] for k in range(n)):
            return dfs(idx + 1, remaining, chosen)
        maxm = min(remaining[k] // beta[k] for k in range(n) if beta[k] > 0)
        for m in range(maxm, -1, -1):
            if m > 0:
                if ext[(beta, beta)] != 0:              # imaginary root: skip (Dynkin: never)
                    continue
                if any(not compatible(beta, g) for g in chosen):
                    continue
            newrem = tuple(remaining[k] - m * beta[k] for k in range(n))
            newchosen = dict(chosen)
            if m > 0:
                newchosen[beta] = m
            res = dfs(idx + 1, newrem, newchosen)
            if res is not None:
                return res
        return None

    return dfs(0, target, {})
