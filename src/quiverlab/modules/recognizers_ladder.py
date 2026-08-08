"""The quasi-tilted / shod / weakly-shod / laura / ada recognizer ladder (Plan 61 / R18).

A thin, exact classifier over the Plan-55 substrate: every verdict is read off ONE
``left_right_parts(A)`` atlas plus the shipped ``global_dimension`` (and, for the
weakly-shod rung, ONE fresh AR knit). The five nested per-instance certificates are

    quasi-tilted  subset  shod  subset  weakly-shod  subset  laura        and        ada,

each a definite ``True``/``False`` with a witness on ``False`` and a certificate on
``True``. The characterisations computed (see the class table in the plan doc):

  * **shod** (Coelho-Lanzilotta): every indecomposable has pd <= 1 OR id <= 1.
    Decision procedure = ``atlas.complement == []`` (survey Thm 4.1, (a)<=>(b));
    cross-checked against the independent QT2 sweep over ``atlas.pd_le_1``/``id_le_1``.
  * **quasi-tilted** (Happel-Reiten-Smalo): (QT1) gl.dim <= 2 AND (QT2) shod.
    Decision = ``gl.dim <= 2 and complement == []``; cross-checked against the ACLV
    route "every indecomposable projective lies in the left part".
  * **weakly shod** (Coelho-Lanzilotta, WSA): the lengths of irreducible-morphism paths
    from an injective to a projective are bounded. On the finite AR quiver this fails iff
    some vertex on an oriented cycle (a non-trivial strongly-connected component) is both
    reachable from an injective and can reach a projective.
  * **laura** (Assem-Coelho): ``ind A \\ (L_A u R_A)`` is finite. Trivially ``True`` in
    representation-finite scope; the finite ``complement`` is the reported datum.
  * **ada** (ACLV Def 2.1): every ``P_x`` and every ``I_x`` lies in ``L_A u R_A``, i.e.
    every projective/injective placement is not ``"neither"``.

The nesting theorems are asserted as a standing MONOTONE self-cert (a ``True`` at a more
special rung under a ``False`` at a more general one is a loud ``QuiverlabError`` -- it
means the atlas or a sweep is wrong, never the mathematics), together with two cheap
gl.dim theorem gates (shod => gl.dim <= 3, quasi-tilted => gl.dim <= 2).

For **ada** algebras over an **algebraically closed** field the shipped ``HH^1`` becomes a
COMPLETE simple-connectedness oracle: ACLV Theorem B (arXiv:1102.1188) says A is simply
connected iff ``HH^1(A) = 0`` there (and then ``HH^*(A)`` reduces to the base field). The
field gate is the additive ``A.domain.is_algebraically_closed`` flag (``True`` only on the
CC working domain -- there is no ``Algebra.field``, and QQi shares CC's ``SympyExactDomain``
class, so the flag, not the class, is the sound predicate). ``dim HH^1`` is char-0
field-independent (over CC = over QQ) by flat base change -- the boundary matrices are
defined over the smaller field and rank is preserved -- so the number computed on the exact
working domain IS the C statement; the algebraically-closed hypothesis controls the theorem's
VALIDITY, not the number, hence the gate is a declaration flag rather than a recomputation.

Representation-finite, non-self-injective scope (inherited from P55/P41): the atlas is
complete iff A is rep-finite and not self-injective, else a loud ``status`` and NO rungs --
never a partial ladder. Identification (locating P_x / I_x / injectives in the universe) is
QQ / char-0 decisive; over large GF(p)/GF(p^n) ``is_isomorphic`` is positive-only and its
loud refusal propagates (never a silent wrong verdict). Exact only; no floats."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.modules.left_right import _index_in_U, left_right_parts

_LADDER_REFS = ("hrs_quasitilted", "coelho_lanzilotta_weakly_shod", "assem_coelho_laura",
                "aclv_supports", "smith_almost_laura", "bft_quasitilted_quiver",
                "organising_module_category")

_THEOREM_B = (
    "ACLV Theorem B (Assem-Castonguay-Lanzilotta-Vargas, arXiv:1102.1188): for an ada "
    "algebra over an algebraically closed field, A is simply connected if and only if "
    "HH^1(A) = 0; moreover, when this holds, the Hochschild cohomology ring HH^*(A) "
    "reduces to the base field."
)

_LADDER_ORDER = ("quasi_tilted", "shod", "weakly_shod", "laura", "ada")


@dataclass(frozen=True)
class LadderRung:
    name: str            # "quasi_tilted" | "shod" | "weakly_shod" | "laura" | "ada"
    verdict: bool
    certificate: dict    # the positive datum (e.g. {"gldim": 1, "complement_size": 0})
    witness: dict | None = None  # populated only on verdict=False


@dataclass(frozen=True)
class SimpleConnectedness:
    applicable: bool
    field_algebraically_closed: bool
    hh1_dim: int | None
    verdict: bool | None
    theorem: str
    note: str


@dataclass(frozen=True)
class RecognizerLadder:
    algebra: object
    rungs: dict                 # {name -> LadderRung}
    complement: tuple           # ({"index","name","dimvec"}, ...) -- the laura datum
    simple_connectedness: SimpleConnectedness
    gldim: int | None
    gldim_exact: bool           # False => gldim is a CERTIFIED LOWER BOUND, not a definite value
    universe_size: int
    is_complete: bool
    status: str
    note: str

    def verdict(self, name) -> bool:
        """The Boolean verdict of the named rung (KeyError if the ladder refused)."""
        return self.rungs[name].verdict


# -- reachability + SCC helpers (pure, index-based) --------------------------

def _reach_closure(n, arrows):
    """Reflexive-transitive closure of the digraph on ``n`` vertices whose edges are the
    KEYS of ``arrows`` (a ``{(i, j): _} `` map or any iterable of ``(i, j)`` pairs)."""
    reach = [[i == j for j in range(n)] for i in range(n)]
    for (i, j) in arrows:
        reach[i][j] = True
    for k in range(n):
        rk = reach[k]
        for i in range(n):
            if reach[i][k]:
                ri = reach[i]
                for j in range(n):
                    if rk[j]:
                        ri[j] = True
    return reach


def _scc_of(reach, Z, n):
    """The strongly-connected component of ``Z``: all ``W`` mutually reachable with ``Z``."""
    return [W for W in range(n) if reach[Z][W] and reach[W][Z]]


def _weakly_shod_from_digraph(n, arrows, inj, proj):
    """PURE bounded-path SCC sweep (drives the synthetic unit test and both routes).

    Returns ``(weakly_shod, witness)``. ``A`` is weakly shod iff NO vertex ``Z`` lying on an
    oriented cycle (a non-trivial SCC, size >= 2) is BOTH reachable from an injective and can
    reach a projective -- then some injective->projective path threads the cycle and its length
    is unbounded. On failure ``witness = {"cycle", "inj_source", "proj_target"}``."""
    reach = _reach_closure(n, arrows)
    for Z in range(n):
        scc = _scc_of(reach, Z, n)
        if len(scc) < 2:                                   # trivial SCC: Z on no oriented cycle
            continue
        src = next((i for i in sorted(inj) if reach[i][Z]), None)
        tgt = next((p for p in sorted(proj) if reach[Z][p]), None)
        if src is not None and tgt is not None:
            return False, {"cycle": sorted(scc), "inj_source": src, "proj_target": tgt}
    return True, None


def _inj_proj_indices(A, U):
    """The universe indices of the standard injectives / projectives of ``A`` located in
    ``U`` by :func:`_index_in_U` (QQ/char-0 decisive; a loud refusal propagates, M1)."""
    inj, proj = set(), set()
    for v in A.quiver.vertices:
        ii = _index_in_U(U, A.injective(v))
        if ii is not None:
            inj.add(ii)
        pp = _index_in_U(U, A.projective(v))
        if pp is not None:
            proj.add(pp)
    return inj, proj


# -- the five rungs ----------------------------------------------------------

def _shod(atlas) -> LadderRung:
    """shod <=> complement empty (survey Thm 4.1, (a)<=>(b)), CROSS-CHECKED against the
    independent QT2 sweep (every indec pd<=1 or id<=1 via atlas.pd_le_1/id_le_1). The two are
    equal by the theorem (the QT2-violator set is a subset of the complement, and shod <=> both
    empty); a mismatch is a loud QuiverlabError (the atlas is wrong). Witness on False = a
    complement module that violates QT2 (pd>1 AND id>1)."""
    pd_ok, id_ok = atlas.pd_le_1, atlas.id_le_1
    qt2_all = all(pd_ok[i] or id_ok[i] for i in range(len(pd_ok)))
    complement_empty = (len(atlas.complement) == 0)
    if qt2_all != complement_empty:
        raise QuiverlabError(
            "recognizer ladder: the shod two-route cross-check disagrees -- the complement is "
            f"{'empty' if complement_empty else 'non-empty'} but the QT2 sweep (every indec "
            f"pd<=1 or id<=1) is {qt2_all} (survey Thm 4.1 (a)<=>(b) violated)",
            hint="report this presentation; the P55 atlas or its pd_le_1/id_le_1 sweep is wrong")
    verdict = complement_empty
    cert = {"complement_size": len(atlas.complement), "qt2_every_indec_pd_or_id": qt2_all}
    witness = None
    if not verdict:
        chosen = None
        for r in atlas.complement:
            i = r["index"]
            if not (pd_ok[i] or id_ok[i]):
                chosen = r
                break
        if chosen is None:                                 # theorem-forbidden; defensive
            chosen = atlas.complement[0]
        i = chosen["index"]
        witness = {"module": chosen["name"], "pd_le_1": bool(pd_ok[i]), "id_le_1": bool(id_ok[i])}
    return LadderRung("shod", verdict, cert, witness)


def _quasi_tilted_aclv_route(A, atlas=None) -> bool:
    """The ACLV route: A is quasi-tilted iff every indecomposable projective lies in the left
    part (equivalently every injective in the right part). Computed from
    ``atlas.projective_placement`` -- every placement in {"L", "both"}."""
    if atlas is None:
        atlas = left_right_parts(A)
    if not atlas.is_complete or atlas.projective_placement is None:
        raise QuiverlabError(
            "quasi-tilted ACLV route: the left/right atlas is not complete, so the projective "
            "placements are unavailable",
            hint="the algebra must be representation-finite and not self-injective")
    return all(p in ("L", "both") for p in atlas.projective_placement.values())


def _quasi_tilted(atlas, gld) -> LadderRung:
    """quasi-tilted <=> (QT1) gl.dim <= 2 AND (QT2) shod (== complement empty). Cross-checked
    against the ACLV every-projective-in-left route (loud QuiverlabError if they differ)."""
    complement_empty = (len(atlas.complement) == 0)
    qt1 = (gld is not None and gld <= 2)
    verdict = qt1 and complement_empty
    aclv = _quasi_tilted_aclv_route(atlas.algebra, atlas)
    if verdict != aclv:
        raise QuiverlabError(
            "recognizer ladder: the quasi-tilted two-route cross-check disagrees -- HRS "
            f"(gl.dim<=2 and complement empty) = {verdict} but the ACLV route (every "
            f"projective in the left part) = {aclv}",
            hint="report this presentation; the P55 atlas placements or the global dimension "
                 "disagree with the HRS characterisation")
    cert = {"gldim": gld, "gldim_le_2": qt1, "complement_size": len(atlas.complement),
            "aclv_every_projective_in_left": aclv}
    witness = None
    if not verdict:
        witness = {"gldim": gld, "gldim_le_2": qt1, "complement_size": len(atlas.complement)}
    return LadderRung("quasi_tilted", verdict, cert, witness)


def _weakly_shod(A, atlas, budget) -> LadderRung:
    """PRIMARY route: the WSA bounded-path SCC sweep on a FRESH AR knit (the irreducible
    morphisms ARE the AR arrows). Complete atlas => the knit closes; a budget trip here is an
    internal inconsistency and refuses loudly."""
    ar = A.ar_quiver(budget_modules=budget)
    if not ar.is_complete:
        raise QuiverlabError(
            "recognizer ladder: the AR knit for the weakly-shod sweep did not close although "
            f"the left/right atlas is complete (status={ar.status!r})",
            hint="raise the budget; this is an internal inconsistency to report")
    U = [v["module"] for v in ar.vertices]
    n = len(U)
    inj, proj = _inj_proj_indices(A, U)
    verdict, witness = _weakly_shod_from_digraph(n, ar.arrows, inj, proj)
    cert = {"universe_size": n, "directed": witness is None and not _has_nontrivial_scc(n, ar.arrows)}
    return LadderRung("weakly_shod", verdict, cert, witness)


def _has_nontrivial_scc(n, arrows):
    reach = _reach_closure(n, arrows)
    return any(len(_scc_of(reach, Z, n)) >= 2 for Z in range(n))


def _weakly_shod_hom(atlas) -> LadderRung:
    """CROSS-CHECK route: the SAME bounded-path SCC sweep run on the Hom-closure predecessor
    matrix ``atlas._leq`` instead of the AR arrows. P55 pins Hom-closure reachability ==
    AR-arrow reachability for representation-finite algebras (rad^inf = 0), so the SCC verdict
    must agree with :func:`_weakly_shod` -- a genuine cross-engine tie of that pin.

    (NB the cyclicity criterion is membership in a non-trivial SCC of the reachability graph,
    NOT ``end_dim >= 2``: bricks lie on AR-quiver oriented cycles -- e.g. cyclic Nakayama
    algebras, whose indecomposables all have End = k yet sit on a tube -- so a non-brick test
    would wrongly report those weakly shod. See the plan's weakly-shod section.)"""
    U = list(atlas._modules)
    n = len(U)
    A = atlas.algebra
    inj, proj = _inj_proj_indices(A, U)
    # atlas._leq is already the reflexive-transitive Hom closure; feed its edges as the digraph.
    edges = [(i, j) for i in range(n) for j in range(n) if atlas._leq[i][j] and i != j]
    verdict, witness = _weakly_shod_from_digraph(n, edges, inj, proj)
    return LadderRung("weakly_shod", verdict, {"universe_size": n, "route": "hom-closure"},
                      witness)


def _laura(atlas) -> LadderRung:
    """laura <=> ind A \\ (L_A u R_A) finite. Trivially True in representation-finite scope;
    the finite complement (size + members) is the reported datum."""
    return LadderRung("laura", True, {"complement_size": len(atlas.complement)}, None)


def _ada(atlas) -> LadderRung:
    """ada (ACLV Def 2.1) <=> every projective AND every injective placement != "neither".
    Asserts the placements are populated when the atlas is complete (C2: a typed error, never
    a NoneType crash)."""
    if atlas.projective_placement is None or atlas.injective_placement is None:
        raise QuiverlabError(
            "recognizer ladder: the left/right atlas is complete but its projective/injective "
            "placements are None (P55 contract violation)",
            hint="report this presentation; a complete atlas must populate the placements")
    proj_bad = [v for v, p in atlas.projective_placement.items() if p == "neither"]
    inj_bad = [v for v, p in atlas.injective_placement.items() if p == "neither"]
    verdict = not proj_bad and not inj_bad
    cert = {"all_projectives_in_union": not proj_bad, "all_injectives_in_union": not inj_bad}
    witness = None
    if not verdict:
        witness = {"projectives_in_neither": proj_bad, "injectives_in_neither": inj_bad}
    return LadderRung("ada", verdict, cert, witness)


# -- self-certs --------------------------------------------------------------

def _assert_nesting(rungs) -> None:
    """The monotone nesting self-cert (loud on violation): quasi_tilted => shod =>
    weakly_shod => laura, quasi_tilted => ada, shod => ada, and ada => laura (rep-finite)."""
    v = {name: rungs[name].verdict for name in rungs}
    implications = [
        ("quasi_tilted", "shod"), ("shod", "weakly_shod"), ("weakly_shod", "laura"),
        ("quasi_tilted", "ada"), ("shod", "ada"), ("ada", "laura"),
    ]
    for lo, hi in implications:
        if v[lo] and not v[hi]:
            raise QuiverlabError(
                f"recognizer ladder: monotone nesting violated -- {lo} is True but {hi} is "
                "False (an internal atlas/sweep inconsistency, never the mathematics)",
                hint="report this presentation; the P55 atlas or a recognizer sweep is wrong")


def _assert_theorem_gates(rungs, gld) -> None:
    """Cheap necessary gl.dim consequences on TRUE verdicts (loud on violation): shod =>
    gl.dim <= 3 (Coelho-Lanzilotta), quasi_tilted => gl.dim <= 2 (HRS QT1)."""
    if gld is None:
        return
    if rungs["shod"].verdict and gld > 3:
        raise QuiverlabError(
            f"recognizer ladder: shod verdict True but gl.dim = {gld} > 3 (Coelho-Lanzilotta "
            "shod => gl.dim <= 3 violated)",
            hint="report this presentation; the shod verdict or the global dimension is wrong")
    if rungs["quasi_tilted"].verdict and gld > 2:
        raise QuiverlabError(
            f"recognizer ladder: quasi-tilted verdict True but gl.dim = {gld} > 2 (HRS QT1 "
            "gl.dim <= 2 violated)",
            hint="report this presentation; the quasi-tilted verdict or the gl.dim is wrong")


# -- the ada / HH^1 simple-connectedness block (ACLV Theorem B) --------------

def _theorem_b_verdict(ada, alg_closed, hh1_dim):
    """The pure Theorem-B decision: ``None`` unless A is ada over an algebraically closed
    field (the hypothesis), else ``hh1_dim == 0`` (simply connected iff HH^1 vanishes)."""
    if not ada or not alg_closed:
        return None
    if hh1_dim is None:
        return None
    return hh1_dim == 0


def _simple_connectedness(A, ada_verdict) -> SimpleConnectedness:
    """The ACLV-Theorem-B block. Emits a definite simple-connectedness verdict ONLY on the ada
    class over an algebraically closed field (``A.domain.is_algebraically_closed`` -- the
    additive flag, True only on the CC working domain; never the CC sentinel, never a class
    check). Off ada: not applicable, no HH computed. Ada over a non-alg-closed field: dim HH^1
    reported, no verdict (Theorem B's hypothesis unmet). The HH^1 call is quiet (verbose=False)
    so the ladder writes no stray Worked-steps files; a loud HH refusal is recorded honestly,
    never a crash."""
    alg_closed = bool(getattr(A.domain, "is_algebraically_closed", False))
    if not ada_verdict:
        return SimpleConnectedness(
            applicable=False, field_algebraically_closed=alg_closed, hh1_dim=None,
            verdict=None, theorem=_THEOREM_B,
            note="Theorem B applies only to ada algebras; this algebra is not ada, so HH^1 is "
                 "not a simple-connectedness certificate here.")
    try:
        hh1 = int(A.hochschild_cohomology(1, verbose=False)[1])
    except (QuiverlabError, DepthLimitError) as e:
        return SimpleConnectedness(
            applicable=(ada_verdict and alg_closed), field_algebraically_closed=alg_closed,
            hh1_dim=None, verdict=None, theorem=_THEOREM_B, note=f"HH^1 unavailable: {e}")
    verdict = _theorem_b_verdict(ada_verdict, alg_closed, hh1)
    applicable = ada_verdict and alg_closed
    if not alg_closed:
        note = (_THEOREM_B + " This field is not algebraically closed, so no "
                "simple-connectedness verdict is emitted; recompute over CC for the Theorem-B "
                "conclusion.")
    elif verdict:
        note = (_THEOREM_B + " Here HH^1(A) = 0, so A is simply connected and HH^*(A) reduces "
                "to the base field.")
    else:
        note = (_THEOREM_B + " Here HH^1(A) != 0, so A is not simply connected.")
    return SimpleConnectedness(applicable, alg_closed, hh1, verdict, _THEOREM_B, note)


# -- the primary compute -----------------------------------------------------

def recognizer_ladder(A, *, budget=256) -> RecognizerLadder:
    """Classify ``A`` against the quasi-tilted/shod/weakly-shod/laura/ada ladder off ONE
    ``left_right_parts`` atlas + ``global_dimension`` + one AR knit (the weakly-shod sweep).
    Complete iff A is representation-finite and not self-injective, else a loud ``status`` and
    NO rungs (never a partial ladder). Raises loudly (never a silent verdict) if the monotone
    nesting, a theorem gate, or a two-route cross-check is violated, or if identification
    refuses off char-scope (M1)."""
    atlas = left_right_parts(A, budget=budget)
    default_sc = SimpleConnectedness(
        applicable=False, field_algebraically_closed=False, hh1_dim=None, verdict=None,
        theorem=_THEOREM_B, note="not computed (the ladder did not complete)")
    if not atlas.is_complete:
        return RecognizerLadder(A, {}, (), default_sc, None, True, 0, False, atlas.status,
                                atlas.note or "")
    gd = A.global_dimension()
    # The verdict logic reads the (possibly lower-bound) VALUE. When gd.exact is False the
    # value is >= the resolution bound (global_dimension caps at 32, so an inexact value is
    # always >= 32 > 3 > 2): every gl.dim-thresholded rung -- quasi-tilted (needs gld <= 2)
    # and the shod theorem gate (needs gld <= 3) -- therefore resolves to the correct False by
    # its own logic, never by trusting a spurious finite number. gd.exact is threaded through
    # ONLY for honest REPORTING (a lower bound must not be printed as a definite gl.dim).
    gld = gd.value
    rungs = {
        "shod": _shod(atlas),
        "quasi_tilted": _quasi_tilted(atlas, gld),
        "weakly_shod": _weakly_shod(A, atlas, budget),
        "laura": _laura(atlas),
        "ada": _ada(atlas),
    }
    rungs = {name: rungs[name] for name in _LADDER_ORDER}   # canonical order
    _assert_nesting(rungs)
    _assert_theorem_gates(rungs, gld)
    complement = tuple({"index": r["index"], "name": r["name"], "dimvec": r["dimvec"]}
                       for r in atlas.complement)
    sc = _simple_connectedness(A, rungs["ada"].verdict)
    return RecognizerLadder(A, rungs, complement, sc, gld, gd.exact,
                            int(atlas.universe_size), True, atlas.status, atlas.note or "")


# -- the JSON block shared byte-for-byte by both runners (Task D) -------------

def _sorted_dimvec(dv):
    return {str(w): int(n) for w, n in sorted(dv.items(), key=lambda kv: str(kv[0]))}


def _rung_json(rung):
    return {"verdict": bool(rung.verdict), "certificate": rung.certificate,
            "witness": rung.witness}


def _sc_json(sc):
    return {"applicable": bool(sc.applicable),
            "field_algebraically_closed": bool(sc.field_algebraically_closed),
            "hh1_dim": (None if sc.hh1_dim is None else int(sc.hh1_dim)),
            "verdict": sc.verdict, "theorem": sc.theorem, "note": sc.note}


def _empty_ladder_block(n, status, note):
    return {"kind": "recognizer_ladder", "n": int(n), "complete": False, "status": status,
            "universe_size": None, "gldim": None, "gldim_exact": True, "ladder": {},
            "complement": [], "simple_connectedness": _sc_json(SimpleConnectedness(
                applicable=False, field_algebraically_closed=False, hh1_dim=None,
                verdict=None, theorem=_THEOREM_B, note="not computed (the ladder did not complete)")),
            "note": note, "references": list(_LADDER_REFS)}


def recognizer_ladder_block(A, *, budget=256) -> dict:
    """The JSON block shared byte-for-byte by both runners (Task D). An ALGEBRA-level kind
    (schema v1, no module block): the five-rung ladder (verdict/certificate/witness each), the
    finite laura complement (names + dim-vectors -- the `Module`s stripped), and the
    ada/HH^1 simple-connectedness block. On refusal (rep-infinite / self-injective / knit
    error) the same shape with an empty ladder + honest `status`/`note`; a char-scope
    identification refusal (large GF(p) is_isomorphic) is surfaced as a loud `error` field
    (the Plan-30 per-entry precedent -- never a 500)."""
    n = len(A.quiver.vertices) if A.quiver is not None else 0
    try:
        L = recognizer_ladder(A, budget=budget)
    except (QuiverlabError, DepthLimitError) as e:
        block = _empty_ladder_block(n, "error", None)
        block["error"] = str(e)
        return block
    if not L.is_complete:
        return _empty_ladder_block(n, L.status, L.note or None)
    ladder = {}
    for name in _LADDER_ORDER:
        rung = L.rungs[name]
        if name == "laura":
            ladder["laura"] = {"verdict": bool(rung.verdict),
                               "complement_size": int(rung.certificate["complement_size"])}
        else:
            ladder[name] = _rung_json(rung)
    return {
        "kind": "recognizer_ladder", "n": int(n), "complete": True, "status": L.status,
        "universe_size": int(L.universe_size),
        "gldim": (None if L.gldim is None else int(L.gldim)),
        "gldim_exact": bool(L.gldim_exact),
        "ladder": ladder,
        "complement": [{"name": r["name"], "dimvec": _sorted_dimvec(r["dimvec"])}
                       for r in L.complement],
        "simple_connectedness": _sc_json(L.simple_connectedness),
        "note": (L.note or None),
        "references": list(_LADDER_REFS),
    }
