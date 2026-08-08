"""Skew-gentle module classification via special-string re-gluing (Plan 68; He-Zhou-Zhu
2004.11136, Garcia-Lavoue 2601.01744 Table 1, Crawley-Boevey / Geiss-de la Pena).

The indecomposable modules of a skew-gentle algebra ``A`` are the indecomposables of the
(isomorphic) split algebra ``A_hat = kQ_hat/I_hat`` (``split.py``).  This module offers
two views:

* ``classify`` runs the P46 string/band census on the **associated gentle algebra**
  ``A^g`` (nilpotent loops), keeps the ADMISSIBLE walks (those that do not use a special
  loop as a letter), and TYPES each ``(r, s) in {u, p}^2`` by whether its endpoints sit
  at a special vertex (Garcia-Lavoue sec 2.1).  A ``p``-endpoint re-glues into its two
  ``+/-`` forms -- which in the split model are literally the two split vertices
  ``i+``/``i-`` -- so a walk has ``2^{#special endpoints}`` split incarnations
  (characteristic-free, sidestepping the classical char != 2 ``k[T]/(T^2-1)`` split).

  (An admissible walk NEVER has a special vertex as an internal node: at a special
  vertex the loop occupies one in- and one out-slot, so there is at most one other
  in-arrow ``a`` and one other out-arrow ``b``, and gentleness forces ``a*b in I`` --
  a length-2 path that a valid walk cannot traverse.  Special vertices thus appear only
  as walk ENDPOINTS, and the ``+/-`` choice is made exactly there.)

* ``skew_gentle_indecomposables`` is the AUTHORITATIVE enumeration: the indecomposable
  modules of the split algebra itself, via the P41 AR quiver on rep-finite instances
  (count == AR vertex count -- the sufficiency oracle for the split relations).  This is
  a SUPERSET of the admissible-string re-gluings: modules on which the split idempotent
  mixes eigenvalues along a path (e.g. the projective ``P_1`` of the headline example)
  are genuine indecomposables that no single admissible string produces -- they are the
  loop-traversal re-gluings, and the AR route captures them all.

``skew_gentle_module`` materialises the string module of a chosen form ON the split
algebra (via the P46 ``_materialise`` self-certificate, which does not require the split
to be a string algebra), so it works on the non-monomial mesh split too.  Char caveat:
``is_indecomposable`` / ``is_isomorphic`` are rigorous over QQ / char > dim M, so the
batteries run over QQ.  Float-free / exact."""
from __future__ import annotations

import itertools
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.skewgentle.split import SkewGentleAlgebra, _copy_name
from quiverlab.skewgentle.triple import SkewGentleTriple, associated_gentle
from quiverlab.strings.modules import _materialise
from quiverlab.strings.walks import (_is_trivial, enumerate_strings, find_bands,
                                     letter_source, letter_target)


_SPLIT_CACHE = {}


def _split_algebra(triple, field):
    """The split algebra for ``(triple, field)``, memoised so repeated materialisations
    (e.g. the two forms of a special string) live over the SAME algebra object -- else
    ``is_isomorphic`` refuses them as "over different algebras" (mirrors the
    ``string_signs`` id-cache)."""
    key = (id(triple), repr(field))
    hit = _SPLIT_CACHE.get(key)
    if hit is not None and hit[0] is triple:
        return hit[1]
    A = SkewGentleAlgebra(triple, field=field)
    _SPLIT_CACHE[key] = (triple, A)
    return A


@dataclass(frozen=True)
class SkewGentleString:
    """A classified admissible walk of the associated gentle algebra.

    ``walk`` is the ``A^g`` walk; ``type`` is ``(r, s) in {u, p}^2`` (``p`` = endpoint
    at a special vertex); ``is_band`` flags a cyclic walk; ``forms`` lists the split
    incarnations (each a mapping ``{endpoint node index: '+'|'-'}`` over the special
    endpoints -- ``2^{#special endpoints}`` of them, one empty form for a ``(u,u)``
    string)."""
    walk: tuple
    type: tuple
    is_band: bool
    forms: tuple


def _uses_special_loop(walk, loop_names) -> bool:
    return any(nm in loop_names for (nm, _d) in walk)


def _node_vertices(Q_sp, walk):
    """The original vertices ``z_0..z_n`` visited by ``walk`` (a trivial walk yields the
    single vertex)."""
    if _is_trivial(walk):
        return [walk[0][1]]
    verts = [letter_source(Q_sp, walk[0])]
    for ell in walk:
        verts.append(letter_target(Q_sp, ell))
    return verts


def _endpoint_positions(walk):
    """The node indices that are walk endpoints: ``(0,)`` for a trivial walk, else the
    two extremes ``(0, len(walk))``."""
    return (0,) if _is_trivial(walk) else (0, len(walk))


def _classify_walk(Q_sp, walk, special, is_band):
    """Build a ``SkewGentleString``: type ``(r, s)`` from the endpoints, forms from the
    ``+/-`` choices at special endpoints."""
    node_verts = _node_vertices(Q_sp, walk)
    positions = _endpoint_positions(walk)
    # r from the first node, s from the last node (they coincide for a trivial walk).
    r = "p" if node_verts[positions[0]] in special else "u"
    s = "p" if node_verts[positions[-1]] in special else "u"
    special_positions = [p for p in positions if node_verts[p] in special]
    if special_positions:
        forms = tuple(dict(zip(special_positions, signs))
                      for signs in itertools.product(("+", "-"),
                                                     repeat=len(special_positions)))
    else:
        forms = ({},)                               # one form, no +/- choice
    return SkewGentleString(walk=tuple(walk), type=(r, s), is_band=is_band, forms=forms)


def classify(triple, max_length=8, budget=4096):
    """The admissible-walk census of the associated gentle algebra, typed and re-glued.

    Returns a list of ``SkewGentleString`` (strings first, then bands).  Completeness is
    inherited from ``enumerate_strings`` (complete iff rep-finite; see ``certificate``).
    """
    Ag = associated_gentle(triple)
    Q_sp = Ag.quiver
    loops = set(triple.loop_names.values())
    census = enumerate_strings(Ag, max_length=max_length, budget=budget)
    bands = find_bands(Ag, max_length=max_length)
    out = []
    for walk in census.walks:
        if _uses_special_loop(walk, loops):
            continue                                # not an admissible string
        out.append(_classify_walk(Q_sp, walk, triple.special, is_band=False))
    for band in bands:
        if _uses_special_loop(band, loops):
            continue
        out.append(_classify_walk(Q_sp, band, triple.special, is_band=True))
    return out


def _lift_walk(triple, walk, form):
    """Lift an ``A^g`` walk to a walk on the split quiver, choosing the ``+/-`` copy at
    each special endpoint per the ``form`` assignment (``{node index: sign}``)."""
    Q, Sp = triple.quiver, triple.special
    Q_sp = associated_gentle(triple).quiver
    node_verts = _node_vertices(Q_sp, walk)

    def node_copy(i):
        v = node_verts[i]
        if v not in Sp:
            return str(v)
        if i not in form:
            raise QuiverlabError(
                f"skew_gentle_module: special vertex {v!r} at node {i} is not an "
                "endpoint of this walk -- admissible walks touch special vertices only "
                "at their ends",
                hint="the walk crosses a special vertex internally; this should not "
                     "arise for an admissible (loop-free) string")
        return f"{v}{form[i]}"

    if _is_trivial(walk):
        return ((None, node_copy(0)),)

    split_walk = []
    for i, (a, d) in enumerate(walk):
        s0, t0 = Q.source(a), Q.target(a)
        if d > 0:                                   # direct: node i = source(a)
            csrc, ctgt = node_copy(i), node_copy(i + 1)
        else:                                       # inverse: node i = target(a)
            ctgt, csrc = node_copy(i), node_copy(i + 1)
        nm = _copy_name(a, s0, t0, csrc, ctgt, Sp)
        split_walk.append((nm, d))
    return tuple(split_walk)


def _materialise_split_walk(A_split, split_walk, name):
    """Materialise + self-certify the string module of ``split_walk`` on the split
    algebra (mirrors ``strings.string_module``'s core, but via ``_materialise`` so it
    works on the non-monomial mesh split too)."""
    Q, dom = A_split.quiver, A_split.domain
    if _is_trivial(split_walk):
        verts = [split_walk[0][1]]
    else:
        verts = [letter_source(Q, split_walk[0])]
        for ell in split_walk:
            verts.append(letter_target(Q, ell))
    n = len(verts)
    action = {a: lm.zeros(n, n, dom) for a in Q.arrows}
    if not _is_trivial(split_walk):
        for i, (nm, d) in enumerate(split_walk):
            if nm is None:
                continue
            if d > 0:
                action[nm][i + 1][i] = dom.one()
            else:
                action[nm][i][i + 1] = dom.one()
    return _materialise(A_split, verts, action, name)


def _form_name(walk, form):
    parts = [nm if d > 0 else f"{nm}^-1" for nm, d in walk]
    base = "M(e_%s)" % walk[0][1] if _is_trivial(walk) else "M(" + " ".join(parts) + ")"
    if form:
        tag = ",".join(f"{i}{s}" for i, s in sorted(form.items()))
        return f"{base}[{tag}]"
    return base


def skew_gentle_module(triple, sgstring, form=0, field=None):
    """Materialise the split-algebra module for the chosen ``form`` of a classified
    (non-band) string.  ``form`` is an integer index into ``sgstring.forms``.  Uses the
    split copies ``i+``/``i-`` as the two forms of a special string; self-certifies via
    ``check_module`` (inside ``_materialise``)."""
    if isinstance(triple, tuple):                   # (quiver, relations, special) tuple
        triple = SkewGentleTriple.make(*triple)
    if sgstring.is_band:
        raise QuiverlabError(
            "skew_gentle_module: band materialisation needs an eigenvalue + the special-"
            "band re-gluing (a backlog item); use band_module on the split directly",
            hint="pass a non-band SkewGentleString, or enumerate via "
                 "skew_gentle_indecomposables (AR route)")
    if not 0 <= form < len(sgstring.forms):
        raise QuiverlabError(
            f"skew_gentle_module: form index {form} out of range "
            f"(0..{len(sgstring.forms) - 1})",
            hint="a (u,u) string has one form; a special string has two per p-end")
    A = _split_algebra(triple, field)
    split_walk = _lift_walk(triple, sgstring.walk, sgstring.forms[form])
    return _materialise_split_walk(A, split_walk, _form_name(sgstring.walk,
                                                            sgstring.forms[form]))


def skew_gentle_indecomposables(triple, max_length=8, budget=4096, field=None):
    """All indecomposable modules of the skew-gentle algebra, materialised on the split
    algebra.

    AUTHORITATIVE route: the indecomposables of the split algebra via the P41 AR quiver
    (rep-finite instances) -- count == AR vertex count, each ``is_indecomposable``.  On
    a rep-infinite instance the AR quiver is not complete; we then fall back to the
    admissible-string form materialisations (a sound partial sample of the string-type
    indecomposables), flagged by the returned modules being a strict subset."""
    if isinstance(triple, tuple):
        triple = SkewGentleTriple.make(*triple)
    A = _split_algebra(triple, field)
    from quiverlab.modules.ar import knit_ar_quiver
    ar = knit_ar_quiver(A)
    if ar.is_complete:
        return [rec["module"] for rec in ar.vertices]
    # rep-infinite: no finite complete list; return the sound string-type sample.
    mods = []
    for s in classify(triple, max_length=max_length, budget=budget):
        if s.is_band:
            continue
        for f in range(len(s.forms)):
            mods.append(skew_gentle_module(triple, s, form=f, field=field))
    return mods
