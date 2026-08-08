"""The no-code ``skew_gentle`` compute block (Plan 68).

Mirrors ``strings_block`` / ``tau_tilting_block``: a single dict summarising the
skew-gentle world of a triple -- the recognizer verdict, the split-quiver shape, the
HZZ Lemma 1.5 dim law, the special-string classification counts, support tau-tilting
via the engine, and the brick-finite <=> rep-finite certificate.  Both runners call
this shared builder, then stamp ``block["citations"] = _citation_pairs(references)``.

Accepts a ``SkewGentleTriple``, a ``(quiver, relations, special)`` tuple, or a split
``Algebra`` carrying the ``_skew_gentle_triple`` construction marker (the dispatch path).
A bare presented algebra WITHOUT the marker reports ``is_skew_gentle: False`` + a note
(the ``BrauerGraphAlgebra`` honest-scope precedent -- recognizing an arbitrary algebra
as skew-gentle up to iso is the iso-problem, not attempted).  Everything is
JSON-serialisable (lists, ints, strings -- no tuples/sets), so it canonicalises through
the Plan-25 key and freezes as a golden.  Float-free."""
from __future__ import annotations

from quiverlab.skewgentle.certificate import (brick_finite_certificate,
                                              support_tau_tilting)
from quiverlab.skewgentle.modules import classify, skew_gentle_indecomposables
from quiverlab.skewgentle.split import SkewGentleAlgebra
from quiverlab.skewgentle.triple import (SkewGentleTriple, associated_gentle)

_REFERENCES = ["he_zhou_zhu", "chen_skew_gentle", "amiot_skew_gentle",
               "garcia_lavoue", "assem_book"]


def _triple_and_algebra(source, field):
    """Resolve ``source`` to a ``(triple, split_algebra)`` pair, or ``(None, None)`` for
    a bare algebra without the construction marker."""
    if isinstance(source, SkewGentleTriple):
        return source, SkewGentleAlgebra(source, field=field)
    if isinstance(source, tuple):
        t = SkewGentleTriple.make(*source)
        return t, SkewGentleAlgebra(t, field=field)
    triple = getattr(source, "_skew_gentle_triple", None)
    if triple is not None:
        return triple, source                       # already a split algebra
    return None, None                               # bare algebra: honest refusal below


def _triple_summary(triple):
    Q = triple.quiver
    return {"vertices": list(Q.vertices),
            "arrows": {name: [s, t] for name, (s, t) in Q.arrows.items()},
            "relations": list(triple.relations),
            "special": sorted(triple.special, key=repr)}


def skew_gentle_block(source, budget=512, field=None):
    """The ``skew_gentle`` block dict (see module docstring for the shape)."""
    triple, A = _triple_and_algebra(source, field)
    if triple is None:
        return {"kind": "skew_gentle", "is_skew_gentle": False,
                "note": "not a skew-gentle algebra: no construction marker "
                        "(_skew_gentle_triple). Recognizing an arbitrary presented "
                        "algebra as skew-gentle up to isomorphism is not attempted; "
                        "build via SkewGentleAlgebra(Q, I, Sp).",
                "references": _REFERENCES}

    Ag = associated_gentle(triple, field=field)
    strings = classify(triple)
    num_special = sum(1 for s in strings if "p" in s.type and not s.is_band)
    has_bands = any(s.is_band for s in strings)

    cert = brick_finite_certificate(triple, budget=budget, field=field)
    rep_finite = cert["rep_finite"]
    tau_finite = cert["tau_tilting_finite"]

    # Counts only when tau-tilting-finite (bounded + safe); rep-infinite reports None.
    num_indec, tau = None, {"num_pairs": None, "complete": False,
                            "status": cert["status"]}
    if tau_finite:
        tt = support_tau_tilting(triple, budget=budget, field=field)
        tau = {"num_pairs": tt["num_pairs"], "complete": tt["complete"],
               "status": tt["status"]}
        try:
            num_indec = len(skew_gentle_indecomposables(triple, field=field))
        except Exception:                            # honest omission, never a crash
            num_indec = None

    return {
        "kind": "skew_gentle",
        "is_skew_gentle": True,
        "triple": _triple_summary(triple),
        "split": {"num_vertices": len(list(A.quiver.vertices)),
                  "num_arrows": len(A.quiver.arrows), "dim": A.dim},
        "dim_law": {"split_dim": A.dim, "assoc_gentle_dim": Ag.dim,
                    "ok": A.dim == Ag.dim},          # HZZ Lemma 1.5 (NECESSARY check)
        "rank": len(list(A.quiver.vertices)),        # |Q_0| + |Sp|
        "classification": {"num_indecomposables": num_indec,
                           "num_special": num_special, "has_bands": has_bands,
                           "status": "complete" if tau_finite else "budget"},
        "tau_tilting": tau,
        "rep_type": {"rep_finite": rep_finite, "scope": cert["scope"],
                     "status": cert["status"]},
        "references": _REFERENCES,
    }
