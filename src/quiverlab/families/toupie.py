"""Toupie algebras (Plan 59 / R35).

Artenstein-Lanzilotta-Solotar (Algebras and Representation Theory 23 (2020) 421-456,
arXiv:1803.10310): a TOUPIE algebra is a kQ/I whose quiver has a UNIQUE source, a UNIQUE
sink, and every other vertex has in-degree 1 and out-degree 1 -- a vertex-disjoint
bundle of ``a`` directed BRANCHES (parallel paths) from source to sink -- with ``I`` an
admissible ideal (monomial branch-truncations and/or non-monomial cross-branch linear
combinations of the full source->sink branch-paths). Same ordinary quiver as the
canonical algebras.

This module ships:
  * ``ToupieAlgebra(branches, relations=None, field=None, expected_dim=None)`` -- the
    presented ``kQ/I`` constructor, with a per-instance dimension certificate for the
    relation-free case (``dim = 2 + sum(l_i - 1) + sum l_i(l_i+1)/2``).
  * ``is_toupie(A)`` -- the graph-shape recognizer (CONNECTED + ACYCLIC + the degree
    pattern; both connectivity and acyclicity are load-bearing -- the path U
    oriented-cycle quiver passes the degree checks yet is neither).
  * ``toupie_branch_count(A)`` -- the number of branches ``a = out_deg(source) ==
    in_deg(sink)`` (asserts the equality).
  * ``toupie_direct_arrow_count(A)`` -- the number of arrows DIRECTLY source->sink
    (length-1 branches). THIS is the ``sl_a`` ``a`` (ALS Thm 6.5), NOT the branch count.
  * ``toupie_sl_a_lower_bound(A)`` -- ``a^2 - 1`` with ``a`` = the direct-arrow count,
    hard-gated to characteristic 0 (ALS Thm 6.5 hypothesis is ``k = C``).
  * ``toupie_block(A)`` -- the ``toupie`` algebra-only compute kind.

The ``a``-Kronecker ``Q_a`` (two vertices, ``a`` parallel arrows, no relations = the
length-1 toupie ``ToupieAlgebra([1]*a)``) has ``HH^* = [1, a^2-1, 0, ...]`` (hereditary,
char-independent, Euler ``chi = 2 - a^2``). The ``sl_a`` inclusion ``HH^1 ⊇ sl_a`` holds
over char 0 with ``a`` = the number of direct source->sink arrows; equality on the
``a``-Kronecker, a strict lower bound on a subdivided toupie.
"""
from __future__ import annotations

import re

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import NotFiniteDimensionalError, QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import _in_degree, _out_degree

_REFS = ("alsolotar_toupie", "assem_book")
_P_TOKEN = re.compile(r"\bp(\d+)\b")


def _free_dim(branches):
    """The relation-free toupie dimension = |Q_0| + #nonzero paths.

    |Q_0| = 2 + sum(l_i - 1) (source, sink, and l_i - 1 internal vertices per branch);
    #nonzero paths per branch of length l is l(l+1)/2 (all sub-paths of the length-l
    line), and the two length-0 idempotents at source/sink and the internal idempotents
    are already counted in |Q_0|. Verified: [1,1]->4, [1,1,1]->5, [2]->6, [1,2]->7,
    [2,2]->10.
    """
    verts = 2 + sum(l - 1 for l in branches)
    paths = sum(l * (l + 1) // 2 for l in branches)
    return verts + paths


def _build_quiver(branches):
    """(Quiver, branch_paths): source vertex ``1``, sink vertex ``2``; branch i a fresh
    directed path ``1 -> ... -> 2`` of length l (a length-1 branch is a single arrow
    ``1 -> 2``). Internal vertices are integers ``>= 3`` (so the quiver round-trips
    through the integer-vertex GUI). ``branch_paths[i]`` = the ordered arrow-name list."""
    source, sink = 1, 2
    verts = [source, sink]
    arrows = {}
    branch_paths = []
    nxt = 3
    for i, l in enumerate(branches):
        prev = source
        names = []
        for j in range(l):
            if j == l - 1:
                cur = sink
            else:
                cur = nxt
                nxt += 1
                verts.append(cur)
            name = f"a{i}_{j}"
            arrows[name] = (prev, cur)
            names.append(name)
            prev = cur
        branch_paths.append(names)
    return Quiver(verts, arrows), branch_paths


def _expand_relations(relations, branch_paths):
    """Expand the convenience branch-path tokens ``p0, p1, ...`` in each relation string
    to the arrow product of the corresponding branch (``p0 -> a0_0*a0_1*...``). Raw arrow
    products pass through unchanged. A token ``pK`` referencing a length-1 branch is a
    user error (a single source->sink arrow lies in rad, not rad^2, so a relation cannot
    involve it) and refuses loudly."""
    if not relations:
        return []
    out = []
    for rel in relations:
        def _sub(m):
            k = int(m.group(1))
            if k < 0 or k >= len(branch_paths):
                raise QuiverlabError(
                    f"ToupieAlgebra: relation references branch token p{k}, but there "
                    f"are only {len(branch_paths)} branches (p0..p{len(branch_paths)-1})",
                    hint="branch tokens are 0-indexed: p0 is the first branch")
            path = branch_paths[k]
            if len(path) < 2:
                raise QuiverlabError(
                    f"ToupieAlgebra: relation references branch token p{k}, whose branch "
                    "has length 1 (a single source->sink arrow); a relation must lie in "
                    "rad^2, but a length-1 branch path sits in rad, not rad^2",
                    hint="cross-branch relations are admissible only between branches of "
                         "length >= 2")
            return "*".join(path)
        out.append(_P_TOKEN.sub(_sub, rel))
    return out


def ToupieAlgebra(branches, relations=None, field=None, expected_dim=None):
    """Build the presented toupie ``kQ/I`` from a list of branch lengths.

    ``branches`` = ``[l_1, ..., l_a]`` (each ``l_i >= 1``). ``relations`` (optional) =
    grammar-token strings over the branch-path tokens ``p0..p_{a-1}`` (monomial
    truncations and/or NON-MONOMIAL cross-branch linear combinations) or over the raw
    arrow names ``a{i}_{j}``. ``expected_dim`` (optional) is checked loudly.

    Relation-free: certifies ``dim == 2 + sum(l_i - 1) + sum l_i(l_i+1)/2``. With
    relations: ``Quiver.algebra`` computes the dimension (raising if not finite), and
    ``expected_dim`` is checked when given. ``a``-Kronecker = ``ToupieAlgebra([1]*a)``
    (``dim = 2 + a``); a single branch of length ``n`` = ``kA_{n+1}``.
    """
    if field is None:
        field = QQ
    branches = list(branches)
    if len(branches) < 1:
        raise QuiverlabError(
            "ToupieAlgebra: need at least one branch",
            hint="ToupieAlgebra([1, 1]) is the 2-Kronecker; ToupieAlgebra([2]) is kA3")
    for i, l in enumerate(branches):
        if not isinstance(l, int) or isinstance(l, bool) or l < 1:
            raise QuiverlabError(
                f"ToupieAlgebra: branch {i} has length {l!r}; branch lengths must be "
                "integers >= 1",
                hint="a length-1 branch is a single source->sink arrow")

    Q, branch_paths = _build_quiver(branches)
    rel_strings = _expand_relations(relations, branch_paths)
    try:
        A = Q.algebra(relations=rel_strings, field=field)
    except NotFiniteDimensionalError as exc:
        raise QuiverlabError(
            f"ToupieAlgebra: the presentation kQ/I is not finite-dimensional ({exc})",
            hint="a genuine toupie is finite-dimensional; check the relations") from exc

    if not rel_strings:
        want = _free_dim(branches)
        if A.dim != want:
            raise QuiverlabError(
                f"ToupieAlgebra dimension certificate failed: relation-free dim should "
                f"be {want}, got {A.dim}",
                hint="please report this presentation -- the closed form "
                     "2 + sum(l_i - 1) + sum l_i(l_i+1)/2 was violated")
    if expected_dim is not None and A.dim != expected_dim:
        raise QuiverlabError(
            f"ToupieAlgebra dimension certificate failed: expected_dim={expected_dim}, "
            f"got dim {A.dim}",
            hint="the relations do not cut the algebra to the expected dimension")

    A._family_citations = _REFS
    return A


def _require_quiver(A, what):
    if A.quiver is None:
        raise QuiverlabError(
            f"{what}: needs a quiver-presented algebra",
            hint="the toupie recognizer walks the quiver shape; structure-constant "
                 "algebras carry no quiver")


def _endpoints(Q):
    """(source, sink, in_deg, out_deg) for a toupie shape, or None if the degree pattern
    fails (unique in-deg-0 source, unique out-deg-0 sink, every other vertex in/out 1)."""
    ind, outd = _in_degree(Q), _out_degree(Q)
    sources = [v for v in Q.vertices if ind[v] == 0]
    sinks = [v for v in Q.vertices if outd[v] == 0]
    if len(sources) != 1 or len(sinks) != 1:
        return None
    s, t = sources[0], sinks[0]
    if s == t:
        return None
    for v in Q.vertices:
        if v == s or v == t:
            continue
        if ind[v] != 1 or outd[v] != 1:
            return None
    return s, t, ind, outd


def is_toupie(A) -> bool:
    """The toupie graph-shape recognizer: the quiver is CONNECTED and ACYCLIC, has a
    UNIQUE source (in-deg 0), a UNIQUE sink (out-deg 0), and EVERY other vertex has
    in-deg 1 and out-deg 1. CONNECTED + ACYCLIC are load-bearing (the path U
    oriented-cycle quiver passes the degree checks but is neither). Refuses loudly on a
    presentation-less algebra."""
    _require_quiver(A, "is_toupie")
    Q = A.quiver
    ep = _endpoints(Q)
    if ep is None:
        return False
    if not Q.is_connected():
        return False
    if not Q.is_acyclic():
        return False
    return True


def toupie_branch_count(A) -> int:
    """The number of branches ``a = out_deg(source) == in_deg(sink)`` (asserts the
    equality). NOT the ``sl_a`` ``a`` -- see :func:`toupie_direct_arrow_count`."""
    if not is_toupie(A):
        raise QuiverlabError(
            "toupie_branch_count: A is not a toupie",
            hint="check is_toupie(A) first")
    s, t, ind, outd = _endpoints(A.quiver)
    a_out, a_in = outd[s], ind[t]
    if a_out != a_in:
        raise QuiverlabError(
            f"toupie_branch_count: out_deg(source)={a_out} != in_deg(sink)={a_in}",
            hint="a toupie is a vertex-disjoint bundle of branches, so the two must "
                 "agree -- please report this quiver")
    return a_out


def toupie_direct_arrow_count(A) -> int:
    """The ``sl_a`` ``a`` (ALS Thm 6.5): the number of arrows DIRECTLY source->sink
    (length-1 branches). ``a_direct <= a_branch``, equal iff every branch has length 1
    (the ``a``-Kronecker)."""
    if not is_toupie(A):
        raise QuiverlabError(
            "toupie_direct_arrow_count: A is not a toupie",
            hint="check is_toupie(A) first")
    s, t, _ind, _outd = _endpoints(A.quiver)
    return sum(1 for (src, tgt) in A.quiver.arrows.values()
               if src == s and tgt == t)


def toupie_sl_a_lower_bound(A) -> int:
    """``dim sl_a = a^2 - 1`` with ``a`` = the number of direct source->sink arrows
    (:func:`toupie_direct_arrow_count`). ALS Thm 6.5 gives ``HH^1(A) ⊇ sl_a`` over a
    field of characteristic 0 (``k = C``), so ``dim HH^1 >= a^2 - 1``. HARD-GATED to
    char 0 (the theorem's hypothesis); refuses loudly off it."""
    if A.domain.characteristic != 0:
        raise QuiverlabError(
            "toupie_sl_a_lower_bound: the sl_a inclusion (ALS Thm 6.5) requires "
            f"characteristic 0, got char {A.domain.characteristic}",
            hint="the theorem's hypothesis is k = C; over GF(p) the inclusion is not "
                 "claimed -- re-run over QQ")
    a = toupie_direct_arrow_count(A)
    return a * a - 1


def toupie_block(A, hh_top=4):
    """The ``toupie`` algebra-only compute kind (Plan 59 Task 5).

    A presentation-less algebra -> ``{"error": ..., "references": [...]}`` (never a 500).
    A non-toupie -> ``{"is_toupie": False, ...nulls..., "note": "not a toupie"}``. A
    toupie -> the branch/direct-arrow counts, ``HH^0..HH^{hh_top}`` (``engine="cs",
    auto_cs=True``), and the ``sl_a`` line (``dim sl_a = a^2 - 1``, ``a`` = direct-arrow
    count; ``char0`` flag + a note when the inclusion is not claimed).
    """
    refs = list(_REFS)
    try:
        it = is_toupie(A)
    except QuiverlabError as exc:
        return {"error": str(exc), "references": refs}
    if not it:
        return {"is_toupie": False, "branch_count": None, "direct_arrow_count": None,
                "hh": None, "hh_top": hh_top, "sl_a": None,
                "note": "not a toupie", "references": refs}

    bc = toupie_branch_count(A)
    dac = toupie_direct_arrow_count(A)
    hh = list(A.hochschild_cohomology(hh_top, engine="cs", auto_cs=True,
                                      verbose=False).dims)
    char0 = (A.domain.characteristic == 0)
    sla_dim = max(dac * dac - 1, 0)
    sl_a = {"a": dac, "dim": sla_dim, "char0": char0}
    note = None
    if not char0:
        note = ("sl_a ⊆ HH^1 (ALS Thm 6.5) is claimed only in characteristic 0; over "
                f"char {A.domain.characteristic} the abstract dim sl_a = {sla_dim} is "
                "reported without the inclusion")
    elif dac < 2:
        note = (f"a_direct = {dac} < 2: sl_a is trivial (dim {sla_dim}); the inclusion "
                "carries no content")
    return {"is_toupie": True, "branch_count": bc, "direct_arrow_count": dac,
            "hh": hh, "hh_top": hh_top, "sl_a": sl_a, "note": note, "references": refs}
