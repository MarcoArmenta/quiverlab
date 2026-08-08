"""Tilted-algebra recognizer (Plan 60 / R17). Liu-Skowronski criterion (Liu Arch. Math. 61
(1993) 12-19; Liu arXiv:1409.2054 Thm 2.6): A is tilted iff Gamma_A has a faithful cut/section
Sigma with Hom_A(X, tau Y) = 0 for all X, Y in Sigma; certified by Ringel's slice theorem
(Thm 1.9(2), LNM 1099 (4.2)): Sigma is a slice iff S = (+) Sigma is a tilting A-module with
End_A(S) hereditary, and then A = End_H(D(S)). Representation-finite exhaustive scope (the knit
refuses self-injective + rep-infinite); three theorem gates (hereditary => tilted; non-semisimple
self-injective => not; gl.dim > 2 => not) extend and speed up the verdict. Exact; char 0 / char >
dim for the tilting count + Gabriel recovery (loud otherwise).

Completeness chain (the record's "|Sigma_0| = n IMPLIED, state as consequence"): the Thm-2.6
witness is a slice (Thm 2.6); a slice is a section (Liu's Remark after Thm 1.9, Happel-Ringel
(7.1) / Ringel LNM 1099 (4.2)); a section meets each tau-orbit of its component exactly once, so
it is a one-module-per-orbit transversal with exactly n = |Q_0| modules. Enumerating transversals
of the n-orbit components is therefore COMPLETE (the witness, if it exists, is among them), and
the #orbits == n component prune is sound. A cut (Def 2.1) need NOT be a transversal -- the
transversal property is invoked only for the slice/section witness."""
from __future__ import annotations

import itertools
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.dynkin_type import dynkin_type, is_connected
from quiverlab.invariants.recognizers import is_hereditary
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import hom_dim
from quiverlab.modules.morphism import direct_sum

_TILTED_REFS = ["liu_tilted_1993", "liu_another_2014", "happel_ringel_tilted", "assem_book"]


@dataclass
class TiltedReport:
    """The verdict of :func:`tilted_check`, mirroring the P41 ``ARQuiver`` / P55
    ``LeftRightAtlas`` honest ``is_complete``/``status`` contract. ``__bool__`` is
    ``verdict == "tilted"`` so ``assert A.is_tilted()`` reads naturally."""
    algebra: object
    verdict: str                       # "tilted" | "not_tilted" | "unknown"
    reason: str                        # "hereditary"|"self_injective"|"gldim>2"|
                                       #   "faithful_section_found"|"search_exhausted"|
                                       #   "budget"|"disconnected"|"unsupported"|"error"
    slice: list = ()                   # [ {"index","name","dimvec"} , ... ]  (Sigma)
    slice_module: object = None        # the Module S = (+) Sigma
    hereditary_type: str = None        # "A_3" | "D_4" | ... | "~A_2" | "unclassified"
    hereditary_algebra: object = None  # H = End_A(S) presented (Task B)
    reconstruction: dict = None        # {"dim_A","dim_H","type","note"} (Task B)
    universe_size: int = 0
    is_complete: bool = False
    status: str = "error"              # mirrors ARQuiver / gate provenance
    note: str = ""

    def __bool__(self):
        return self.verdict == "tilted"


# --------------------------------------------------------------------------- #
# Type layer: the shipped invariants/dynkin_type.py classifier + adapters.
# --------------------------------------------------------------------------- #

def _type_str(t) -> str:
    """('A', 3) -> "A_3"; ('~A', 2) -> "~A_2"; None -> "unclassified"."""
    if t is None:
        return "unclassified"
    letter, k = t
    return f"{letter}_{k}"


def _is_finite_dynkin(t) -> bool:
    """True iff ``dynkin_type`` returned a FINITE ADE type ``('A'|'D'|'E', n)``. A finite
    Dynkin quiver is exactly the rep-FINITE hereditary case (Gabriel), so only then is the
    projective slice knit-enumerable. Euclidean (``~A``/``~D``/``~E``) / wild / unclassified
    (``None``) are rep-INFINITE -- their AR knit does not terminate, so Gate H must NEVER build
    it (the Plan-60 hang fix): the ``tilted`` verdict + type stand and the slice is omitted."""
    return t is not None and t[0] in ("A", "D", "E")


def _section_quiver(ar, section_indices):
    """The section-graph-to-quiver adapter: a :class:`Quiver` whose vertices are
    ``section_indices`` relabeled ``1..k`` and whose arrows are the AR arrows with both
    endpoints in the section (multiplicity flattened to distinct arrow names). A section is
    connected and acyclic, so the result is a tree; ``dynkin_type`` classifies it A/D/E."""
    from quiverlab.combinat.quiver import Quiver
    pos = {orig: new + 1 for new, orig in enumerate(section_indices)}
    arrows = {}
    a = 0
    for (i, j) in ar.arrows:
        if i in pos and j in pos:
            a += 1
            arrows[f"s{a}"] = (pos[i], pos[j])
    return Quiver(list(range(1, len(section_indices) + 1)), arrows)


# --------------------------------------------------------------------------- #
# The Liu-Skowronski filters (exact over every Domain -- no char caveat).
# --------------------------------------------------------------------------- #

def _action_of_label(M, label):
    """The action matrix on ``M`` of an algebra basis element ``label`` (an idempotent
    ``e_v`` or a path ``a1*...*ak``), read off the shipped ``Module`` representation --
    NOT hardcoding the per-path direction (``_action_of_word`` composes the arrow actions
    in the module's own anti-homomorphism order)."""
    if label.startswith("e_"):
        return M.action[label]
    return M._action_of_word(tuple(label.split("*")))


def _is_faithful(A, M) -> bool:
    """``M`` is faithful iff ``ann_A(M) = 0``, iff the k-linear action map
    ``A -> End_k(M)`` sending each A-basis element to its action block on ``M`` is
    injective: ``dim ann = dim A - rank``, faithful iff ``rank == dim A``. Exact over
    every Domain (no char caveat).

    REDUNDANT for correctness -- ``is_tilting_module`` already enforces faithfulness (a
    tilting module is faithful, Bongartz) -- so the airtight ``_certify_slice`` acceptance
    never relies on this filter; it exists to realize the record's Thm-2.6 "faithful cut"
    formulation and as a cheap early prune that rejects Liu's sincere-not-faithful cut
    before the expensive tilting/End checks. Never load-bearing alone."""
    dom = M.domain
    d = M.dim
    if d == 0:
        return False
    rows = []
    for label in A.basis_labels:
        blk = _action_of_label(M, label)
        rows.append([blk[i][j] for i in range(d) for j in range(d)])
    return lm.mat_rank(rows, dom) == len(A.basis_labels)


def _hom_tau_zero(A, section_modules) -> bool:
    """``Hom_A(X, tau Y) = 0`` for all ``X, Y`` in the section (Liu-Skowronski Thm 2.6).
    ``Y.tau()`` is the zero module when ``Y`` is projective, and ``hom_dim(., 0) = 0``, so
    no projective special-case is needed. Exact over every Domain."""
    taus = [Y.tau() for Y in section_modules]
    for X in section_modules:
        for tY in taus:
            if hom_dim(X, tY) != 0:
                return False
    return True


def _certify_slice(A, S):
    """Ringel Thm 1.9(2), the independent airtight acceptance: ``S = (+) Sigma`` generates
    a slice iff ``S`` is a tilting A-module AND ``H = End_A(S)`` presented is hereditary.
    Returns ``(ok, H)`` with ``H`` the presented ``End_A(S)`` (P44 ``presented_form``) or
    ``None``. ``presented_form``/``is_tilting_module`` are char 0 / char > dim (loud
    otherwise); the char-free ``_is_faithful``/``_hom_tau_zero`` filters never reach here
    for a non-tilting S because of the short-circuit."""
    from quiverlab.core.basic import presented_form
    from quiverlab.modules.endomorphism import end_algebra
    from quiverlab.modules.tilting import is_tilting_module
    if not is_tilting_module(S).is_tilting:
        return (False, None)
    H = presented_form(end_algebra(S))
    return (is_hereditary(H), H)


# --------------------------------------------------------------------------- #
# The search: components -> #orbits==n prune -> transversal enumeration ->
# faithful + Hom(X,tauY)=0 filter -> Ringel 1.9(2) certificate.
# --------------------------------------------------------------------------- #

def _components(n, arrows):
    """Connected components of ``0..n-1`` under the UNDIRECTED AR-arrow graph."""
    adj = {i: set() for i in range(n)}
    for (i, j) in arrows:
        adj[i].add(j)
        adj[j].add(i)
    seen, comps = set(), []
    for s in range(n):
        if s in seen:
            continue
        stack, comp = [s], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            seen.add(x)
            stack.extend(adj[x] - comp)
        comps.append(sorted(comp))
    return comps


def _valid_transversal(seed, ar, n, comps):
    """True iff ``seed`` is a one-module-per-tau-orbit transversal of a single n-orbit
    component (used to guard the optional P55 seed before trusting it)."""
    if seed is None or len(seed) != n:
        return False
    N = len(ar.vertices)
    if any((not isinstance(i, int)) or i < 0 or i >= N for i in seed):
        return False
    orbit_of = {}
    for oi, orbit in enumerate(ar.tau_orbits):
        for i in orbit:
            orbit_of[i] = oi
    if len({orbit_of.get(i) for i in seed}) != n:
        return False                                   # not one-per-orbit
    for comp in comps:
        cs = set(comp)
        if all(i in cs for i in seed):
            orbits_here = [o for o in ar.tau_orbits if o[0] in cs]
            return len(orbits_here) == n
    return False


def _accepts(ar, A, combo):
    """The full Liu-Skowronski + Ringel acceptance of a transversal ``combo`` (a tuple of
    indices): faithful AND Hom(X, tau Y)=0 AND Ringel 1.9(2) certificate. Returns True/False."""
    U = [v["module"] for v in ar.vertices]
    mods = [U[i] for i in combo]
    S, _, _ = direct_sum(*mods)
    if not _is_faithful(A, S):
        return False
    if not _hom_tau_zero(A, mods):
        return False
    ok, _H = _certify_slice(A, S)
    return ok


def _faithful_section_search(ar, A, budget_sections, seed=None):
    """Search the complete knit for a faithful section with ``Hom(X, tau Y)=0`` certified a
    slice by Ringel 1.9(2). Returns the section indices (a list), ``None`` (a genuine
    rep-finite refutation), or the sentinel ``"BUDGET"`` (transversal cap tripped -- an
    honest unknown, never a false ``not_tilted``)."""
    U = [v["module"] for v in ar.vertices]
    n = len(list(A.quiver.vertices))
    comps = _components(len(U), ar.arrows)

    # Optimization (P55 seed): a certified hit short-circuits the enumeration. Guarded --
    # the exhaustive knit search below stays authoritative, and an invalid seed is ignored.
    if _valid_transversal(seed, ar, n, comps) and _accepts(ar, A, tuple(seed)):
        return list(seed)

    examined = 0
    for comp in comps:
        cs = set(comp)
        orbits = [o for o in ar.tau_orbits if o[0] in cs]
        if len(orbits) != n:                           # #orbits == n prune (sound; see docstring)
            continue
        for combo in itertools.product(*orbits):
            examined += 1
            if examined > budget_sections:
                return "BUDGET"
            if _accepts(ar, A, combo):
                return list(combo)
    return None


# --------------------------------------------------------------------------- #
# Gate H helpers.
# --------------------------------------------------------------------------- #

def _locate(U, M):
    """First index ``j`` with ``U[j]`` isomorphic to ``M`` (dim-vector prefilter, then the
    exact ``is_isomorphic`` certificate), or ``None``. QQ / char-0 decisive; over large
    GF(p) ``is_isomorphic`` is positive-only and its loud refusal propagates."""
    from quiverlab.modules.hom import is_isomorphic
    if M is None or M.dim == 0:
        return None
    dv = M.dimension_vector()
    for j, X in enumerate(U):
        if X.dim == M.dim and X.dimension_vector() == dv and is_isomorphic(X, M):
            return j
    return None


def _projective_section(A, budget_modules):
    """Gate-H certified payload: for rep-finite hereditary ``A``, the indecomposable
    projectives ``{P_v}`` form the slice ``A_A`` (``End_A(A_A)`` hereditary => a genuine
    slice). Returns ``(slice_records, S)``; ``([], None)`` when the AR knit does not complete
    within ``budget_modules`` -- the ``tilted`` verdict + type still stand, the slice is
    omitted with an honest note. Gate H calls this ONLY for a finite-Dynkin (rep-finite)
    quiver; the rep-INFINITE hereditary case (Euclidean/wild, e.g. Kronecker) never reaches
    here -- its slice is omitted upstream WITHOUT building the knit (which would not
    terminate), decided from the Dynkin type of ``A``'s own quiver."""
    ar = A.ar_quiver(budget_modules=budget_modules)
    if not ar.is_complete:
        return [], None
    U = [v["module"] for v in ar.vertices]
    recs, mods = [], []
    for v in A.quiver.vertices:
        idx = _locate(U, A.projective(v))
        if idx is None:                                # a complete knit contains every projective
            return [], None
        recs.append({"index": idx, "name": ar.vertices[idx]["name"],
                     "dimvec": ar.vertices[idx]["dimvec"]})
        mods.append(U[idx])
    S, _, _ = direct_sum(*mods)
    return recs, S


# --------------------------------------------------------------------------- #
# The optional P55 seed.
# --------------------------------------------------------------------------- #

def _p55_seed(A, budget_modules):
    """The ACLV Ext-projectives-of-add-R_A section (P55 ``ext_projectives_right``) as an
    optional first transversal for the search. Best-effort and defensive -- any refusal or
    an incomplete atlas returns ``None`` and the exhaustive search proceeds unchanged."""
    try:
        from quiverlab.modules.left_right import left_right_parts
        atlas = left_right_parts(A, budget=budget_modules)
        if not atlas.is_complete:
            return None
        idx = [r["index"] for r in atlas.ext_projectives_right]
        return idx or None
    except QuiverlabError:
        return None


# --------------------------------------------------------------------------- #
# The verdict pipeline.
# --------------------------------------------------------------------------- #

def tilted_check(A, *, budget_modules=256, budget_sections=4096):
    """Decide whether ``A = kQ/I`` is tilted (``A = End_H(T)``, ``H`` hereditary, ``T``
    tilting) by the Liu-Skowronski faithful-section criterion, certified by Ringel's slice
    theorem. Theorem gates first (H: hereditary => tilted; S: non-semisimple self-injective
    => not; G: ``gl.dim >= 3`` => not), then the rep-finite exhaustive faithful-section
    search on ``Gamma_A``. Returns a :class:`TiltedReport`. Needs a quiver presentation
    (structure-constant-only ``A`` refuses loudly)."""
    if getattr(A, "quiver", None) is None:
        raise QuiverlabError(
            "tilted_check needs a quiver presentation",
            hint="build A via Quiver.algebra(...); a tilted verdict recovers the hereditary "
                 "type from the Gabriel quiver of End_A(S)")
    n = len(list(A.quiver.vertices))

    # Connectedness -- a tilted algebra is a CONNECTED End-algebra over a CONNECTED hereditary
    # algebra (ASS2006 VIII.4), so a disconnected A is not tilted by definition. Decided up
    # front for consistency: without it a disconnected hereditary A returns tilted/unclassified
    # from Gate H while the non-hereditary path would refute it. The P55 support surface feeds
    # the recognizer its CONNECTED components one at a time, so this never fires there.
    if not is_connected(A.quiver):
        return TiltedReport(A, "not_tilted", "disconnected", is_complete=True, status="complete",
                            note="disconnected quiver -- tilted algebras are connected by "
                                 "definition (a connected End-algebra over a connected "
                                 "hereditary algebra, ASS2006)")

    # Gate H -- hereditary => tilted (Happel-Ringel); type read from A's OWN quiver (instant).
    if is_hereditary(A):
        dt = dynkin_type(A.quiver)
        typ = _type_str(dt)
        # Only finite-Dynkin (rep-finite) hereditary has a knit-enumerable projective slice.
        # For Euclidean/wild (rep-INFINITE) hereditary -- e.g. the Kronecker quiver ~A_1 -- the
        # AR knit does not terminate, so we NEVER build it: the slice is honestly omitted while
        # the tilted verdict + type stand (Plan-60 hang fix).
        if _is_finite_dynkin(dt):
            slc, S = _projective_section(A, budget_modules)
        else:
            slc, S = [], None
        H = None
        recon = None
        if S is not None:
            H = A
            recon = {"dim_A": A.dim, "dim_H": A.dim, "type": typ,
                     "note": "hereditary: H = A, T = A_A is a tilting A-module, "
                             "End_A(A_A) is hereditary"}
            note = ("hereditary: slice = the indecomposable projectives (A_A), "
                    "End_A(A_A) hereditary")
        elif _is_finite_dynkin(dt):
            note = ("hereditary rep-finite but the AR knit did not complete within "
                    "budget_modules: slice omitted; verdict + type stand")
        else:
            note = ("hereditary but representation-infinite (Euclidean/wild): the AR knit is "
                    "not built (it would not terminate); the postprojective section exists but "
                    "is not knit-enumerable, so the slice is omitted -- verdict + type stand")
        return TiltedReport(A, "tilted", "hereditary", slice=slc, slice_module=S,
                            hereditary_type=typ, hereditary_algebra=H, reconstruction=recon,
                            is_complete=True, status="complete", note=note)

    # Gate S -- non-semisimple self-injective => not tilted (Gate H ate the semisimple case).
    if A.is_selfinjective():
        return TiltedReport(A, "not_tilted", "self_injective", is_complete=True,
                            status="complete",
                            note="non-semisimple self-injective => gl.dim = infinity => not "
                                 "tilted (tilted => gl.dim <= 2)")

    # Gate G -- gl.dim >= 3 => not tilted. gd.value is ALWAYS int; exact=False => a certified
    # LOWER BOUND, and a lower bound >= 3 already proves gl.dim > 2.
    gd = A.global_dimension()
    if gd.value >= 3:
        return TiltedReport(A, "not_tilted", "gldim>2", is_complete=True, status="complete",
                            note=f"gl.dim {'=' if gd.exact else '>='} {gd.value} > 2; "
                                 f"tilted => gl.dim <= 2")

    # The faithful-section search (rep-finite, non-self-injective, gl.dim <= 2).
    ar = A.ar_quiver(budget_modules=budget_modules)
    if not ar.is_complete:
        return TiltedReport(A, "unknown", ar.status, is_complete=False, status=ar.status,
                            note=(ar.note or "") + " | rep-infinite extension path: the local "
                                 "criterion arXiv:1409.2054 (Thm 2.6 on a locally-computed cut)")
    seed = _p55_seed(A, budget_modules)
    found = _faithful_section_search(ar, A, budget_sections, seed=seed)
    if found == "BUDGET":
        return TiltedReport(A, "unknown", "budget", universe_size=len(ar.vertices),
                            is_complete=False, status="budget",
                            note="section (transversal) enumeration exceeded budget_sections")
    if found is None:
        return TiltedReport(A, "not_tilted", "search_exhausted", universe_size=len(ar.vertices),
                            is_complete=True, status="complete",
                            note="no faithful section with Hom(X, tau Y)=0 (rep-finite, "
                                 "complete knit, budget-exhaustive)")
    slc = [{"index": i, "name": ar.vertices[i]["name"], "dimvec": ar.vertices[i]["dimvec"]}
           for i in found]
    S, _, _ = direct_sum(*[ar.vertices[i]["module"] for i in found])
    sec_type = _type_str(dynkin_type(_section_quiver(ar, found)))
    ok, H = _certify_slice(A, S)                        # ok is True (search already certified)
    h_type = _type_str(dynkin_type(H.quiver))
    if sec_type != h_type:                              # a knit or Gabriel-recovery bug -- loud
        raise QuiverlabError(
            f"tilted_check: section-graph type {sec_type!r} disagrees with the recovered "
            f"Gabriel type {h_type!r} of End_A(S) -- internal inconsistency")
    recon = {"dim_A": A.dim, "dim_H": H.dim, "type": sec_type,
             "note": "A = End_H(D(S)) (Ringel Thm 1.9(2)); the reported invariants are "
                     "checkable, the isomorphism is theorem-guaranteed"}
    return TiltedReport(A, "tilted", "faithful_section_found", slice=slc, slice_module=S,
                        hereditary_type=sec_type, hereditary_algebra=H, reconstruction=recon,
                        universe_size=len(ar.vertices), is_complete=True, status="complete")


# --------------------------------------------------------------------------- #
# The JSON block the two runners share (Task E).
# --------------------------------------------------------------------------- #

def tilted_check_block(A, *, budget_modules=256, budget_sections=4096):
    """The ``tilted_check`` JSON block (Task E): serialize the :class:`TiltedReport`,
    stripping the Module/Algebra objects to dim-vector / presented-quiver summaries. The
    char-scope refusal (``presented_form``/``is_isomorphic`` over char <= dim) surfaces as
    a clean ``{"kind": "tilted_check", "error": ...}`` entry (Plan-30 precedent, never a
    500). ``references`` present in every shape so ``block["citations"] = ...`` is uniform."""
    try:
        rep = tilted_check(A, budget_modules=budget_modules, budget_sections=budget_sections)
    except QuiverlabError as exc:
        return {"kind": "tilted_check", "error": str(exc), "references": list(_TILTED_REFS)}
    slice_dimvec = None
    if rep.slice_module is not None:
        slice_dimvec = {str(k): int(v)
                        for k, v in rep.slice_module.dimension_vector().items()}
    hered = None
    if rep.hereditary_algebra is not None:
        Hq = rep.hereditary_algebra.quiver
        hered = {"vertices": [str(v) for v in Hq.vertices],
                 "arrows": {name: [str(s), str(t)] for name, (s, t) in Hq.arrows.items()},
                 "relations": [str(r) for r in (rep.hereditary_algebra.relations or [])],
                 "dim": int(rep.hereditary_algebra.dim)}
    return {
        "kind": "tilted_check",
        "n": int(len(list(A.quiver.vertices))) if getattr(A, "quiver", None) else 0,
        "complete": bool(rep.is_complete),
        "status": rep.status,
        "verdict": rep.verdict,
        "reason": rep.reason,
        "slice": [{"name": r["name"],
                   "dimvec": {str(k): int(v) for k, v in r["dimvec"].items()}}
                  for r in rep.slice],
        "slice_dimvec": slice_dimvec,
        "hereditary_type": rep.hereditary_type,
        "hereditary_algebra": hered,
        "reconstruction": (None if rep.reconstruction is None else
                           {"dim_A": int(rep.reconstruction["dim_A"]),
                            "dim_H": int(rep.reconstruction["dim_H"]),
                            "type": rep.reconstruction["type"],
                            "note": rep.reconstruction["note"]}),
        "note": rep.note or None,
        "references": list(_TILTED_REFS),
    }
