"""The tau-cluster morphism category W(A) + picture group (Plan 66 / R29).

A thin COMBINATORIAL layer over the merged Plan-45 tau-tilting engine, the Plan-64
wide-subcategory poset, and the Plan-65 Jasso tau-perpendicular reduction -- NO new
representation theory (the bricks, wides and reductions are all inherited).

Objects and morphisms (Buan-Marsh, *A category of wide subcategories*, IMRN 2021):

- **objects** = the tau-perpendicular wide subcategories of ``mod A``; for a
  tau-tilting-FINITE ``A`` these are ALL the finitely many wide subcategories, so
  ``#objects == wide_subcategories(A).size`` (the Enomoto core-label-order count, P64 --
  the R29 "object count = #wide" cross-engine oracle);
- **morphisms** ``[U]: W -> W'`` = the support tau-rigid pairs ``U`` of ``W`` with
  ``W' = J_W(U)`` its tau-perpendicular category, graded by ``rank = |U|``.

Everything is read straight off ``A``'s own ``g``-fan (the fan of the support tau-rigid
pairs -- a triangulated ``(n-1)``-sphere): each face is a subset of the signed ``g``-vector
columns of a support tau-tilting facet, and the **closed star** of a face ``U`` (the faces
``V >= U``, graded by ``|V| - |U|``) IS the ``g``-fan of the reduction ``C(U)`` (the Jasso
link identity, a theorem). So

- objects = the faces deduped by the wide ``J(U)`` they cut out (``J(U)`` identified by the
  set of bricks it contains -- ``beta in J(U)`` iff ``Hom(M, beta) = 0``, ``Hom(beta, tau M)
  = 0`` and ``beta`` vanishes on the support of ``P``, for the pair ``U = (M, P)``);
- the per-object out-degree ``= |closedStar(U_W)| = #faces(C_W)``;
- the **Hanson-Igusa classifying-space cube complex** ``face_vector`` has
  ``f_0 = #objects = #wide`` (one 0-cell per wide) and
  ``f_k = #(rank-k morphisms) = sum_W #(rank-k faces of C_W's g-fan)`` (one k-cell per
  rank-k morphism), ``sum(face_vector) = morphism_count``, ``f_n = #sTt`` -- this is the
  ACTUAL topology, NOT the ``g``-fan sphere (kept separately as ``g_fan_face_vector``).

The picture group ``pi_1(|W(A)|)`` (Igusa-Todorov-Weyman arXiv:1609.02636; Hanson-Igusa):
one generator per brick, one relation per rank-2 wide (``commutation`` for ``k x k``, an
``atom``/pentagon for a connected ``kA_2`` wide -- split by ``Ext^1`` between the two simple
bricks), and the abelianization by exact Smith normal form.

**Honest gate (STRICTER than a bounded region; P65 M1).** ``W(A)`` is a FINITE category iff
``A`` is tau-tilting-finite. A partial category / face vector / picture group off a
budget-truncated (or ``status != "complete"``) exchange graph would be a LIE, so P66 is
certified iff ``eg.is_complete and eg.status == "complete"``; otherwise every invariant is
``None`` and the honest ``note`` is set.

**Scope (inherited from P45/P64/P65).** Runs over QQ by default; ``bricks`` / ``hom`` /
``tau`` / the exchange-graph BFS are rigorous over char 0 or char > dim (Dickson/CIW), a
brick decides ``End(B) = k`` over the algebraically-closed / char-0 base (the GF(p^n)
proper-division-ring caveat is honest-scope). Off scope the underlying engine refuses loudly.

Float-free: faces are ``frozenset`` of int tuples, counts are ``int``, the abelianization is
exact SNF over ``ZZ``. All refusals are ``QuiverlabError``.

References: buan_marsh_wide (IMRN 2021 -- the category), hanson_igusa (Comm. Alg. 2021 -- the
cube complex + Nakayama K(pi,1) + picture group), igusa_todorov_weyman (arXiv:1609.02636 --
the presentation), igusa_todorov_cat0 (arXiv:2203.16679 -- hereditary-Dynkin K(pi,1) /
honest-scope); enomoto_wide_ice / jasso_reduction / air_tau_tilting / demonet_iyama_jasso
(P64/P65/P45)."""
from __future__ import annotations

import weakref
from collections import Counter
from dataclasses import dataclass
from itertools import combinations

from quiverlab.errors import QuiverlabError

# A per-(algebra, budget) exchange-graph memo LOCAL to Plan 66 (P45's exchange_graph is left
# byte-unchanged -- this only avoids rebuilding the SAME graph three+ times within one
# tau_cluster_category / picture_group / tau_cluster_block call). Weak-keyed on the algebra
# object (identity-based, drops with the algebra), so there is no cross-instance bleed.
_EG_CACHE: "weakref.WeakKeyDictionary" = weakref.WeakKeyDictionary()


def _cached_eg(A, budget):
    """The support tau-tilting exchange graph of ``A`` at ``budget``, memoized per (algebra,
    budget) for the duration of the algebra's life. Byte-identical to
    ``mutation.exchange_graph(A, budget_pairs=budget)`` (it IS that call on a miss)."""
    from quiverlab.tautilting.mutation import exchange_graph
    bucket = _EG_CACHE.setdefault(A, {})
    eg = bucket.get(budget)
    if eg is None:
        eg = exchange_graph(A, budget_pairs=budget)
        bucket[budget] = eg
    return eg


def _prewarm_universe(A, eg):
    """Pre-populate ``torsion._UNIVERSE_CACHE[A]`` from an ALREADY-built complete exchange
    graph so the downstream ``bricks`` / ``wide_subcategories`` skip their own
    ``_torsion_universe`` exchange-graph rebuild (byte-identical universe -- same dedup rule,
    same source pairs)."""
    from quiverlab.modules.hom import is_isomorphic
    from quiverlab.tautilting import torsion as _t
    if A in _t._UNIVERSE_CACHE:
        return
    uni = []
    for rec in eg.vertices:
        for M in rec["pair"].summands:
            if not any(U.dim == M.dim and U.dimension_vector() == M.dimension_vector()
                       and is_isomorphic(U, M) for U in uni):
                uni.append(M)
    _t._UNIVERSE_CACHE[A] = uni

_REFERENCES = ["buan_marsh_wide", "hanson_igusa", "igusa_todorov_weyman",
               "igusa_todorov_cat0", "enomoto_wide_ice", "jasso_reduction",
               "air_tau_tilting", "demonet_iyama_jasso"]

_KPI1_UNCERTIFIED = ("not certified (K(pi,1) known only for Nakayama / hereditary Dynkin "
                     "at the cited theorems)")

_INCOMPLETE_NOTE = (
    "The exchange graph did not close (status={status}): A is tau-tilting-infinite or the "
    "pair budget ({budget}) was hit. W(A) is an INFINITE category, so its objects, "
    "morphisms, the cube-complex classifying space and the picture group are ALL undefined "
    "here (a partial category would be a lie, not merely incomplete). Certified only iff A "
    "is tau-tilting-finite (DIJ). Raise the budget for a genuinely finite algebra.")


# --------------------------------------------------------------------------- #
# value objects
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TauClusterCategory:
    """The finite category ``W(A)`` (Buan-Marsh IMRN 2021) + its Hanson-Igusa classifying
    space. See the module docstring for the field semantics; every invariant is ``None`` on a
    non-complete exchange graph (the honest refusal)."""
    algebra: object
    objects: tuple           # (id, rank, semibrick(=simple-brick dimvecs), cw_face_vector)
    object_count: object     # int == wide_subcategories(A).size  (cross-tie) | None
    morphisms: tuple         # ((source_id, target_id, rank, pair_cols), ...)
    morphism_count: object   # int | None
    out_degree: dict         # object_id -> #morphisms out (== |closedStar(U_W)| == #faces(C_W))
    face_vector: object      # HANSON-IGUSA CLASSIFYING SPACE tuple: f_0 == object_count | None
    g_fan_face_vector: object  # the g-fan/cluster SPHERE tuple: f_n == #sTt | None
    euler_characteristic: object  # int == sum((-1)^k face_vector[k]) (NOT (-1)^n) | None
    is_kpi1: object          # True (Nakayama / hered. Dynkin) | None
    kpi1_reason: str
    is_complete: bool
    status: str
    note: str


@dataclass(frozen=True)
class PictureGroup:
    """The picture-group presentation ``pi_1(|W(A)|)`` as DATA (Igusa-Todorov-Weyman;
    Hanson-Igusa). See the module docstring; every invariant is ``None`` on a non-complete
    exchange graph."""
    algebra: object
    generators: tuple        # per brick: (id, dimvec, name)
    relations: tuple         # per rank-2 wide: (type, (i,j), ext_brick_id|None, word)
    num_generators: object   # int == #bricks | None
    num_relations: object    # int == #rank-2 wides | None
    num_atom: object
    num_commutation: object
    abelianization: tuple    # SNF invariant factors > 1 (torsion; () for a free Z^r)
    abelianization_rank: object  # int == #bricks - rank(relation matrix) (H2) | None
    is_kpi1: object
    kpi1_reason: str
    is_complete: bool
    status: str
    note: str


# --------------------------------------------------------------------------- #
# the g-fan of A: faces = subsets of the signed g-vector columns of the facets
# --------------------------------------------------------------------------- #
def _require_quiver(A, what):
    if getattr(A, "quiver", None) is None:
        raise QuiverlabError(
            f"{what}: requires a quiver-presented algebra kQ/I (a structure-constant-only "
            "algebra has no tau-tilting / wide-subcategory surface)")


def _build_gfan(A, budget):
    """Build ``A``'s ``g``-fan from the exchange graph. Returns ``(eg, None)`` when the graph
    is not ``status="complete"`` (drives the honest refusal), else ``(eg, (n, faces,
    col_kind))``: ``faces`` = every subset of a facet's column set (deduped ``frozenset`` of
    column tuples); ``col_kind`` maps a column tuple to ``("mod", Module)`` (a tau-rigid
    summand) or ``("supp", vertex)`` (a killed projective ``P_v[1]``, column ``-e_v``)."""
    eg = _cached_eg(A, budget)
    if not eg.is_complete or eg.status != "complete":
        return eg, None
    _prewarm_universe(A, eg)
    verts = list(A.quiver.vertices)
    n = len(verts)
    col_kind = {}
    facets = []
    for rec in eg.vertices:
        pair = rec["pair"]
        gm = rec["g_matrix"]
        ncols = len(gm[0]) if gm else 0
        cols = [tuple(gm[r][c] for r in range(n)) for c in range(ncols)]
        nmods = len(pair.summands)
        sup_sorted = sorted(pair.support, key=verts.index)
        for ci, col in enumerate(cols):
            if col not in col_kind:
                if ci < nmods:
                    col_kind[col] = ("mod", pair.summands[ci])
                else:
                    col_kind[col] = ("supp", sup_sorted[ci - nmods])
        facets.append(frozenset(cols))
    faces = set()
    for f in facets:
        fl = list(f)
        for k in range(len(fl) + 1):
            for sub in combinations(fl, k):
                faces.add(frozenset(sub))
    return eg, (n, faces, col_kind)


def _dvt(module_or_dict):
    """A brick / module dim-vector as a tuple over the sorted vertex order."""
    dv = module_or_dict.dimension_vector() if hasattr(module_or_dict, "dimension_vector") \
        else module_or_dict
    return tuple(int(dv[v]) for v in sorted(dv, key=str))


def _perp_tables(A, col_kind, brick_list):
    """Precompute the tau-perpendicular membership predicates ONCE per INDECOMPOSABLE (never
    per face). Because ``tau`` commutes with ``(+)`` and ``Hom`` is additive, ``beta in J(M, P)``
    factors as: ``beta`` vanishes on every support vertex of ``P``, AND for every module summand
    ``M_i`` of ``M`` both ``Hom(M_i, beta) = 0`` and ``Hom(beta, tau M_i) = 0``. So for each
    distinct ``mod`` column (an indecomposable tau-rigid ``M_i``, keyed by its ``g``-column) we
    store the boolean vector ``ok_i[b] = (Hom(M_i, b) == 0 and Hom(b, tau M_i) == 0)`` over the
    bricks, and per brick the set of vertices it vanishes on. Turns O(#faces x #bricks) Hom
    calls into O(#indecomposables x #bricks)."""
    from quiverlab.modules.hom import hom_dim
    vanish = []                                   # per brick: frozenset of vertices with mult 0
    for beta in brick_list:
        dvb = beta.dimension_vector()
        vanish.append(frozenset(v for v in A.quiver.vertices if int(dvb.get(v, 0)) == 0))
    mod_ok = {}                                   # g-column tuple -> tuple(bool over bricks)
    for col, (kind, val) in col_kind.items():
        if kind != "mod":
            continue
        M = val
        tauM = M.tau()
        mod_ok[col] = tuple(hom_dim(M, beta) == 0 and hom_dim(beta, tauM) == 0
                            for beta in brick_list)
    return vanish, mod_ok


def _wide_key(U, col_kind, brick_list, tables):
    """The identity of the wide ``J(U)`` as ``frozenset`` of the indices of the bricks it
    contains (the Jasso / Buan-Marsh tau-perpendicular condition, via the precomputed
    :func:`_perp_tables`). A wide is determined by the set of bricks it contains, so this is an
    exact, iso-class-safe dedup key (never a dim-vector key -- P64 ruling 3)."""
    vanish, mod_ok = tables
    mod_cols = []
    pverts = set()
    for col in U:
        kind, val = col_kind[col]
        if kind == "mod":
            mod_cols.append(col)
        else:
            pverts.add(val)
    in_wide = []
    for bi in range(len(brick_list)):
        if pverts and not pverts <= vanish[bi]:
            continue
        if all(mod_ok[col][bi] for col in mod_cols):
            in_wide.append(bi)
    return frozenset(in_wide)


def _simple_bricks(wide_key, brick_list):
    """A best-effort semibrick label: the bricks of the wide minimal under dim-vector
    domination (``beta`` is a simple object of the wide when no other brick of the wide has a
    componentwise-smaller dim-vector). Exact for the thin / rank<=2 batteries; a label only
    (the authoritative rank is ``n - |U_W|``), never a claimed count."""
    dvs = {bi: brick_list[bi].dimension_vector() for bi in wide_key}
    simples = []
    for bi in wide_key:
        dvb = dvs[bi]
        dominated = False
        for bj in wide_key:
            if bj == bi:
                continue
            dvc = dvs[bj]
            if (all(dvc.get(v, 0) <= dvb.get(v, 0) for v in dvb)
                    and any(dvc.get(v, 0) < dvb.get(v, 0) for v in dvb)):
                dominated = True
                break
        if not dominated:
            simples.append(_dvt(dvb))
    return tuple(sorted(simples))


# --------------------------------------------------------------------------- #
# the K(pi,1) verdict -- theorem-anchored (the honest homotopy boundary)
# --------------------------------------------------------------------------- #
def _kpi1_verdict(A):
    """``(True, reason)`` iff ``|W(A)|`` is provably a ``K(pi,1)``: Nakayama (Hanson-Igusa,
    Comm. Alg. 2021) or hereditary Dynkin (Igusa-Todorov, CAT(0), arXiv:2203.16679). Else
    ``(None, ...)`` -- NEVER asserted beyond the cited theorems (the general tau-tilting-finite
    case is delicate)."""
    from quiverlab.invariants.recognizers import is_hereditary, is_nakayama
    if is_nakayama(A):
        return True, "Nakayama (Hanson-Igusa)"
    if is_hereditary(A):
        from quiverlab.invariants.dynkin_type import dynkin_type
        dt = dynkin_type(A.quiver)
        if dt and dt[0] in ("A", "D", "E"):
            return True, "hereditary Dynkin (Igusa-Todorov, CAT(0) arXiv:2203.16679)"
    return None, _KPI1_UNCERTIFIED


# --------------------------------------------------------------------------- #
# Task A: the category
# --------------------------------------------------------------------------- #
def _composition_spotcheck(objects_by_id, morphisms):
    """A structural self-cert of the category law (Buan-Marsh: W(A) is a genuine category).
    Every morphism ``W -> W'`` of rank ``k`` drops the rank by exactly ``k``
    (``rank(W) - k == rank(W')``); every object has exactly ONE rank-0 morphism (its identity
    ``[emptyset]``), a self-loop. Loud ``QuiverlabError`` on any violation (never a
    silently-wrong category)."""
    id_rank = {oid: rank for (oid, rank, _sb, _cw) in objects_by_id}
    id_identities = {oid: 0 for (oid, _r, _sb, _cw) in objects_by_id}
    for (src, tgt, k, _cols) in morphisms:
        if id_rank[src] - k != id_rank[tgt]:
            raise QuiverlabError(
                f"tau_cluster_category: morphism {src}->{tgt} rank {k} does not drop the "
                f"source rank {id_rank[src]} to the target rank {id_rank[tgt]} -- the "
                "composition/target structure is wrong", hint="report the algebra")
        if k == 0:
            if src != tgt:
                raise QuiverlabError(
                    f"tau_cluster_category: a rank-0 morphism {src}->{tgt} is not a self-loop "
                    "(the identity of an object must fix it)", hint="report the algebra")
            id_identities[src] += 1
    for oid, cnt in id_identities.items():
        if cnt != 1:
            raise QuiverlabError(
                f"tau_cluster_category: object {oid} has {cnt} identity (rank-0) morphisms, "
                "expected exactly 1", hint="report the algebra")


def tau_cluster_category(A, *, budget=512) -> TauClusterCategory:
    """The tau-cluster morphism category ``W(A)`` (Plan 66 / R29; Buan-Marsh IMRN 2021 +
    Hanson-Igusa Comm. Alg. 2021): objects = the tau-perpendicular wide subcategories
    (``object_count == wide_subcategories(A).size``, cross-tied to P64), morphisms = support
    tau-rigid pairs of the source graded by rank (out-degree = closed star = ``#faces(C_W)``),
    the Hanson-Igusa classifying-space ``face_vector`` (``f_0 = #wide``, ``f_k = #(rank-k
    morphisms)``), the separate ``g_fan_face_vector`` (the g-fan SPHERE, ``f_n = #sTt``), the
    Euler characteristic, and the theorem-anchored ``K(pi,1)`` verdict. Certified complete iff
    ``A`` is tau-tilting-finite; else an honest ``is_complete=False`` value with every
    invariant ``None`` (the P65 M1 loud refusal, never a partial-category lie)."""
    _require_quiver(A, "tau_cluster_category")
    eg, data = _build_gfan(A, budget)
    if data is None:
        # W(A) is infinite here -- the K(pi,1) verdict is moot (no finite classifying space),
        # so it is honestly None regardless of the algebra's Nakayama/Dynkin type.
        return TauClusterCategory(
            algebra=A, objects=(), object_count=None, morphisms=(), morphism_count=None,
            out_degree={}, face_vector=None, g_fan_face_vector=None,
            euler_characteristic=None, is_kpi1=None, kpi1_reason=_KPI1_UNCERTIFIED,
            is_complete=False, status=eg.status,
            note=_INCOMPLETE_NOTE.format(status=eg.status, budget=budget))
    n, faces, col_kind = data
    from quiverlab.tautilting.torsion import bricks as _bricks
    brick_list = list(_bricks(A, budget=budget))

    # 1. the wide of every face (dedup key), grouped -> objects.
    tables = _perp_tables(A, col_kind, brick_list)
    face_key = {U: _wide_key(U, col_kind, brick_list, tables) for U in faces}
    groups = {}
    for U, key in face_key.items():
        groups.setdefault(key, []).append(U)

    # 2. object ids (deterministic: by rank then key), rank = n - |U_W|, representative face.
    def _rank_of(key):
        return n - len(groups[key][0])
    keys_sorted = sorted(groups, key=lambda k: (_rank_of(k), tuple(sorted(k))))
    key_to_id = {k: i for i, k in enumerate(keys_sorted)}
    reps = {k: min(groups[k], key=lambda U: (len(U), tuple(sorted(U)))) for k in groups}

    # 3. closed star of each representative = C_W's g-fan (Jasso link identity); face_vector.
    def _graded_star(U):
        c = Counter()
        for V in faces:
            if U <= V:
                c[len(V) - len(U)] += 1
        return tuple(c[k] for k in range(max(c) + 1))

    objects = []
    out_degree = {}
    face_vec = Counter()
    star_of = {}
    for k in keys_sorted:
        oid = key_to_id[k]
        U_W = reps[k]
        rank = n - len(U_W)
        star = _graded_star(U_W)
        star_of[k] = star
        out_degree[oid] = sum(star)
        for kk, val in enumerate(star):
            face_vec[kk] += val
        objects.append((oid, rank, _simple_bricks(k, brick_list), star))
    face_vector = tuple(face_vec[kk] for kk in range(max(face_vec) + 1))
    object_count = len(groups)
    morphism_count = sum(out_degree.values())

    # 4. the explicit morphism list: each object W (rep U_W) -> J(V) for every face V >= U_W.
    morphisms = []
    for k in keys_sorted:
        src = key_to_id[k]
        U_W = reps[k]
        for V in faces:
            if U_W <= V:
                tgt = key_to_id[face_key[V]]
                extra = tuple(sorted(V - U_W))
                morphisms.append((src, tgt, len(V) - len(U_W), extra))
    morphisms = tuple(morphisms)

    # 5. loud self-certs (H1): f_0 == #wide, sum == morphism_count, f_n == #sTt; the P64 tie.
    stt = len(eg.vertices)
    g_fan = Counter(len(f) for f in faces)
    g_fan_face_vector = tuple(g_fan[kk] for kk in range(max(g_fan) + 1))
    from quiverlab.tautilting.congruence import wide_subcategories
    wide_size = wide_subcategories(A, budget=budget).size
    if object_count != wide_size:
        raise QuiverlabError(
            f"tau_cluster_category: enumerated {object_count} objects but "
            f"wide_subcategories(A).size == {wide_size} -- the tau-perpendicular enumeration "
            "or the Enomoto core-label order is wrong (they must agree: objects = #wide)",
            hint="report the algebra")
    if face_vector[0] != object_count:
        raise QuiverlabError(
            f"tau_cluster_category: face_vector[0]={face_vector[0]} != object_count="
            f"{object_count} -- f_0 must be #wide (H1); the g-fan sphere was re-emitted as "
            "the classifying space", hint="report the algebra")
    if sum(face_vector) != morphism_count:
        raise QuiverlabError(
            f"tau_cluster_category: sum(face_vector)={sum(face_vector)} != morphism_count="
            f"{morphism_count}", hint="report the algebra")
    if face_vector[-1] != stt or g_fan_face_vector[-1] != stt:
        raise QuiverlabError(
            f"tau_cluster_category: face_vector[-1]={face_vector[-1]}, "
            f"g_fan_face_vector[-1]={g_fan_face_vector[-1]} -- both top cells must equal "
            f"#sTt={stt}", hint="report the algebra")
    _composition_spotcheck(objects, morphisms)

    euler = sum((-1) ** kk * f for kk, f in enumerate(face_vector))
    kpi1, reason = _kpi1_verdict(A)
    return TauClusterCategory(
        algebra=A, objects=tuple(objects), object_count=object_count, morphisms=morphisms,
        morphism_count=morphism_count, out_degree=out_degree, face_vector=face_vector,
        g_fan_face_vector=g_fan_face_vector, euler_characteristic=euler, is_kpi1=kpi1,
        kpi1_reason=reason, is_complete=True, status="complete", note="")


# --------------------------------------------------------------------------- #
# Task B: the picture group
# --------------------------------------------------------------------------- #
def _std_name(B):
    from quiverlab.modules.hom import identify_standard
    std = identify_standard(B)
    if std is None:
        return None
    kind, v = std
    return {"simple": "S", "projective": "P", "injective": "I"}[kind] + str(v)


def _find_ext_brick(A, Bi, Bj, brick_list):
    """The extension brick of a connected rank-2 wide ``{Bi, Bj}``: the third brick, of
    dim-vector ``Bi + Bj`` (the connected case), identified by ISO-CLASS among the bricks (a
    dim-vector prefilter, then ``is_isomorphic``; never a bare dim-vector match)."""
    from quiverlab.modules.hom import is_isomorphic
    target = {}
    dvi, dvj = Bi.dimension_vector(), Bj.dimension_vector()
    for v in A.quiver.vertices:
        target[v] = int(dvi.get(v, 0)) + int(dvj.get(v, 0))
    for bi, beta in enumerate(brick_list):
        if beta.dimension_vector() != target:
            continue
        # sum brick != Bi, Bj; certify by iso-class among the same-dim-vector candidates
        try:
            if is_isomorphic(beta, Bi) or is_isomorphic(beta, Bj):
                continue
        except QuiverlabError:
            pass
        return bi
    return None


def _relation_type(A, i, j, brick_list):
    """Classify the rank-2 wide ``{B_i, B_j}`` (hom-orthogonal bricks): ``("commutation",
    None)`` when ``Ext^1(B_i, B_j) = Ext^1(B_j, B_i) = 0`` (``W ~ k x k``), else
    ``("atom", ext_brick_id)`` (``W ~ mod kA_2``, connected). At most one Ext direction is
    nonzero for a tau-tilting-finite rank-2 wide (both nonzero = Kronecker-type =
    tau-tilting-infinite) -- asserted loudly."""
    from quiverlab.modules.ext import ext_dims
    Bi, Bj = brick_list[i], brick_list[j]
    e1 = ext_dims(A, Bi, Bj, 1)[1]
    e2 = ext_dims(A, Bj, Bi, 1)[1]
    if e1 == 0 and e2 == 0:
        return "commutation", None
    if e1 != 0 and e2 != 0:
        raise QuiverlabError(
            "picture_group: a rank-2 wide has Ext^1 in BOTH directions "
            f"(Ext(B{i},B{j})={e1}, Ext(B{j},B{i})={e2}) -- that is a Kronecker-type wide, "
            "hence tau-tilting-INFINITE (W(A) would be an infinite category)",
            hint="over char <= dim the brick/Ext caveat can misreport; run over QQ")
    ext_id = _find_ext_brick(A, Bi, Bj, brick_list)
    return "atom", ext_id


def _abelianization(num_gen, relations):
    """The abelianization ``G(A)^ab = Z^{num_gen} / (relation matrix)`` by exact Smith normal
    form over ``ZZ`` (H2 correction). Each commutation relation is the zero row; each atom
    relation ``x_i + x_j = x_i + x_delta + x_j`` gives the row ``e_delta`` (killing the
    EXTENSION brick ``delta``). Several atoms may share one ``delta`` (identical rows), so
    ``rank(matrix) = #(distinct extension bricks) <= #atom`` -- NEVER ``#atom``. Returns
    ``(invariant_factors > 1, abelianization_rank = num_gen - rank)``."""
    import sympy as sp
    from sympy import ZZ
    from sympy.matrices.normalforms import invariant_factors
    rows = []
    for (typ, _pair, ext_id, _word) in relations:
        if typ == "atom" and ext_id is not None:
            row = [0] * num_gen
            row[ext_id] = 1
            rows.append(row)
        # commutation (and an atom whose ext brick could not be identified) -> no relation row
    if not rows:
        return (), num_gen
    M = sp.Matrix(rows)
    rank = M.rank()
    factors = tuple(int(d) for d in invariant_factors(M, domain=ZZ) if int(d) > 1)
    return factors, num_gen - rank


def picture_group(A, *, budget=512) -> PictureGroup:
    """The picture-group presentation ``pi_1(|W(A)|)`` as DATA (Plan 66 / R29;
    Igusa-Todorov-Weyman arXiv:1609.02636; Hanson-Igusa): generators = the bricks
    (``num_generators == #bricks``), relations = the rank-2 wides (``num_relations == #rank-2
    wides``, each typed ``commutation`` for ``k x k`` / ``atom`` for a connected ``kA_2`` via
    ``Ext^1``), and the abelianization by exact SNF (``abelianization_rank = #bricks -
    rank(relation matrix) = #bricks - #(distinct extension bricks among the atoms)``, H2 -- NOT
    ``#bricks - #atom``). Certified complete iff ``A`` is tau-tilting-finite; else an honest
    ``is_complete=False`` value with every invariant ``None``."""
    _require_quiver(A, "picture_group")
    eg = _cached_eg(A, budget)
    if not eg.is_complete or eg.status != "complete":
        return PictureGroup(
            algebra=A, generators=(), relations=(), num_generators=None, num_relations=None,
            num_atom=None, num_commutation=None, abelianization=(), abelianization_rank=None,
            is_kpi1=None, kpi1_reason=_KPI1_UNCERTIFIED, is_complete=False, status=eg.status,
            note=_INCOMPLETE_NOTE.format(status=eg.status, budget=budget))
    _prewarm_universe(A, eg)
    from quiverlab.modules.hom import hom_dim
    from quiverlab.tautilting.torsion import bricks as _bricks
    brick_list = list(_bricks(A, budget=budget))
    generators = tuple((i, _dvt(B), _std_name(B)) for i, B in enumerate(brick_list))

    def _orth(i, j):
        return hom_dim(brick_list[i], brick_list[j]) == 0 \
            and hom_dim(brick_list[j], brick_list[i]) == 0

    relations = []
    num_atom = num_commutation = 0
    for i, j in combinations(range(len(brick_list)), 2):
        if not _orth(i, j):
            continue
        typ, ext_id = _relation_type(A, i, j, brick_list)
        if typ == "atom":
            num_atom += 1
            word = _atom_word(i, j, ext_id, generators)
        else:
            num_commutation += 1
            word = _commutation_word(i, j, generators)
        relations.append((typ, (i, j), ext_id, word))
    relations = tuple(relations)
    factors, ab_rank = _abelianization(len(brick_list), relations)
    kpi1, reason = _kpi1_verdict(A)
    return PictureGroup(
        algebra=A, generators=generators, relations=relations,
        num_generators=len(brick_list), num_relations=len(relations), num_atom=num_atom,
        num_commutation=num_commutation, abelianization=factors, abelianization_rank=ab_rank,
        is_kpi1=kpi1, kpi1_reason=reason, is_complete=True, status="complete", note="")


def _gen_name(gid, generators):
    """A display token ``x_{S1}`` / ``x_{(1,1)}`` for generator ``gid``."""
    _i, dv, name = generators[gid]
    label = name if name is not None else "(" + ",".join(str(x) for x in dv) + ")"
    return "x_{%s}" % label


def _commutation_word(i, j, generators):
    """The commutation relation ``x_i x_j = x_j x_i`` (``W ~ k x k``, ITW)."""
    return "%s %s = %s %s" % (_gen_name(i, generators), _gen_name(j, generators),
                              _gen_name(j, generators), _gen_name(i, generators))


def _atom_word(i, j, ext_id, generators):
    """The atom/pentagon relation of a connected rank-2 wide ``{B_i, B_j}`` with extension
    brick ``B_delta`` (ITW slope-ordered normal form, LHS the two simple bricks, RHS all three
    bricks of the wide): ``x_i x_j = x_i x_delta x_j``. Read cyclically it is the length-5
    pentagon of the two maximal green sequences of the wide."""
    if ext_id is None:
        return "%s %s = %s %s" % (_gen_name(i, generators), _gen_name(j, generators),
                                  _gen_name(j, generators), _gen_name(i, generators))
    return "%s %s = %s %s %s" % (
        _gen_name(i, generators), _gen_name(j, generators),
        _gen_name(i, generators), _gen_name(ext_id, generators), _gen_name(j, generators))


# --------------------------------------------------------------------------- #
# the shared no-code block (both runners route through this)
# --------------------------------------------------------------------------- #
def _brick_entry(dvt_tuple, name, A):
    """A ``{"dimvec": {v: m}, "name": str|None}`` block entry from a dim-vector tuple."""
    verts = sorted(A.quiver.vertices, key=str)
    return {"dimvec": {str(v): int(m) for v, m in zip(verts, dvt_tuple) if m},
            "name": name}


def tau_cluster_block(A, *, budget=512) -> dict:
    """The shared ``tau_cluster`` compute-kind payload (Plan 66 Task D): the category summary
    (objects, morphisms, ranks), the Hanson-Igusa classifying-space cube complex (``face_vector``,
    Euler characteristic, ``K(pi,1)`` verdict) alongside the ``g``-fan/cluster SPHERE
    (``g_fan_face_vector``, ``f_n = #sTt``), and the picture-group presentation (generators,
    typed relations, abelianization). Certified complete iff ``A`` is tau-tilting-finite; on an
    incomplete exchange graph ``category = picture_group = None`` and a ``note`` (never a
    partial-category lie). Consumed byte-identically by both runners (each adds ``citations``)."""
    if getattr(A, "quiver", None) is None:
        return {"kind": "tau_cluster",
                "error": ("the tau-cluster morphism category requires a quiver-presented "
                          "algebra kQ/I (a structure-constant-only algebra has no tau-tilting "
                          "/ wide-subcategory surface)"),
                "references": list(_REFERENCES)}
    n = len(list(A.quiver.vertices))
    C = tau_cluster_category(A, budget=budget)
    block = {
        "kind": "tau_cluster", "n": n,
        "complete": C.is_complete, "status": C.status,
        "references": list(_REFERENCES), "note": None,
    }
    if not C.is_complete:
        block["category"] = None
        block["picture_group"] = None
        block["note"] = C.note
        return block
    objects_by_rank = Counter(rank for (_oid, rank, _sb, _cw) in C.objects)
    morphisms_by_rank = {k: f for k, f in enumerate(C.face_vector)}
    block["category"] = {
        "object_count": C.object_count,
        "morphism_count": C.morphism_count,
        "objects_by_rank": {int(r): int(c) for r, c in sorted(objects_by_rank.items())},
        "morphisms_by_rank": {int(r): int(c) for r, c in sorted(morphisms_by_rank.items())},
        "face_vector": list(C.face_vector),
        "g_fan_face_vector": list(C.g_fan_face_vector),
        "euler_characteristic": C.euler_characteristic,
        "is_kpi1": C.is_kpi1, "kpi1_reason": C.kpi1_reason,
    }
    G = picture_group(A, budget=budget)
    gens = [_brick_entry(dv, name, A) for (_i, dv, name) in G.generators]
    rels = []
    for (typ, (i, j), ext_id, word) in G.relations:
        gi, gj = G.generators[i], G.generators[j]
        entry = {
            "type": typ,
            "simple_bricks": [_brick_entry(gi[1], gi[2], A), _brick_entry(gj[1], gj[2], A)],
            "ext_brick": (None if ext_id is None
                          else _brick_entry(G.generators[ext_id][1],
                                            G.generators[ext_id][2], A)),
            "word": word,
        }
        rels.append(entry)
    block["picture_group"] = {
        "num_generators": G.num_generators,
        "generators": gens,
        "num_relations": G.num_relations,
        "num_atom": G.num_atom,
        "num_commutation": G.num_commutation,
        "relations": rels,
        "abelianization": list(G.abelianization),
        "abelianization_rank": G.abelianization_rank,
    }
    return block
