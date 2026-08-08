"""The brick-finite <=> representation-finite certificate + support tau-tilting for
skew-gentle algebras (Plan 68).

Composed theorem: ``brick-finite <=> tau-tilting-finite`` (Demonet-Iyama-Jasso, any
f.d. algebra) composed with ``rep-finite <=> brick-finite`` (Garcia-Lavoue Thm 3.1,
skew-gentle, char != 2), so for a skew-gentle ``A`` over char != 2::

    rep-finite  <=>  tau-tilting-finite  <=>  exchange_graph(split).is_complete

Two routes:

* **route 1 (PRIMARY, rep-finite upgrade)** -- ``exchange_graph(split).is_complete`` on
  the char-FREE split model over QQ (the exchange-graph BFS leans on
  ``is_isomorphic``/``decompose``, rigorous only over char 0 / char > dim M, so it is
  always computed over QQ; the split construction is characteristic-free, so this is
  legitimate).  ``status == 'complete'`` => tau-tilting-finite; ``status == 'budget'``
  => tau-tilting-INfinite (rep-infinite over char != 2); any other status
  (e.g. an engine ``'error'``) leaves the verdict ``None`` (honest, never guessed).

* **route 2 (ONE-SIDED cross-check, rep-infinite downgrade only)** --
  ``find_bands(associated_gentle)``.  A band FOUND certifies rep-INFINITE (a real band
  is a real band -- infinitely many band modules), so it may only DOWNGRADE the verdict
  to rep-infinite; but ``find_bands(A^g)`` may MISS SPECIAL bands (bands touching Sp,
  which re-glue through the +/- split copies), so "no-bands" does NOT by itself certify
  rep-finite -- route 1 is authoritative for the upgrade (W4).  Because a found band is
  a fast, sound rep-infinite certificate, it short-circuits the (possibly expensive)
  route-1 BFS.

char 2 (M3): the rep-finite reading is WITHHELD (``rep_finite is None``) -- the
rep-finite equivalence is char != 2 (Garcia-Lavoue).  The tau-tilting-finiteness verdict
is still computed on the char-free split model over QQ (never over GF(2), never
``decompose`` over char 2) and reported, scope-flagged.  Char caveat: the batteries run
over QQ.  Float-free / exact."""
from __future__ import annotations

from quiverlab.skewgentle.split import SkewGentleAlgebra
from quiverlab.skewgentle.triple import (SkewGentleTriple, associated_gentle)
from quiverlab.fields import QQ


def _as_triple(triple):
    if isinstance(triple, tuple):
        return SkewGentleTriple.make(*triple)
    return triple


def _char_of(field):
    if field is None:
        return 0                                   # default domain is char-0 (QQ/CC)
    c = getattr(field, "characteristic", 0)
    return c() if callable(c) else c


def _analyze(triple, budget):
    """The char-free two-route analysis on the QQ split model.  Returns
    ``{band_route, tau_tilting_finite, num_bricks, status}``."""
    from quiverlab.strings.walks import find_bands

    Ag = associated_gentle(triple)
    loops = set(triple.loop_names.values())
    bands = [b for b in find_bands(Ag)
             if not any(nm in loops for (nm, _d) in b)]
    if bands:
        # route 2 sound downgrade: a band => rep-infinite => tau-tilting-infinite (the
        # char-free model has the same band); skip the expensive route-1 BFS.
        return {"band_route": "bands", "tau_tilting_finite": False,
                "num_bricks": None, "status": "band"}

    from quiverlab.tautilting.mutation import exchange_graph
    A = SkewGentleAlgebra(triple, field=QQ)         # char-free model (M3)
    eg = exchange_graph(A, budget_pairs=budget)
    if eg.status == "complete" and eg.is_complete:
        from quiverlab.tautilting.torsion import bricks
        return {"band_route": "no-bands", "tau_tilting_finite": True,
                "num_bricks": len(bricks(A, budget=budget)), "status": "complete"}
    if eg.status == "budget":
        return {"band_route": "no-bands", "tau_tilting_finite": False,
                "num_bricks": None, "status": "budget"}
    # engine error / genuinely undecided -- honest None, never a guess.
    return {"band_route": "no-bands", "tau_tilting_finite": None,
            "num_bricks": None, "status": eg.status}


def brick_finite_certificate(triple, budget=512, field=None):
    """The full two-route payload::

        {"tau_tilting_finite": bool|None, "num_bricks": int|None,
         "rep_finite": bool|None, "char": int,
         "scope": "char!=2"|"char==2 (narrowed)",
         "status": "complete"|"budget"|"band"|<engine status>,
         "band_route": "bands"|"no-bands"}
    """
    triple = _as_triple(triple)
    triple.validate()
    a = _analyze(triple, budget)
    char = _char_of(field)
    if char == 2:
        rep_finite, scope = None, "char==2 (narrowed)"       # M3: withhold the upgrade
    else:
        rep_finite, scope = a["tau_tilting_finite"], "char!=2"
    return {"tau_tilting_finite": a["tau_tilting_finite"],
            "num_bricks": a["num_bricks"],
            "rep_finite": rep_finite,
            "char": char, "scope": scope, "status": a["status"],
            "band_route": a["band_route"]}


def is_representation_finite(triple, budget=512, field=None):
    """``True`` / ``False`` (char != 2) or ``None`` (char 2 -- withheld, or a genuinely
    undecided route).  The record's certificate: DIJ ``brick-finite <=>
    tau-tilting-finite`` composed with Garcia-Lavoue Thm 3.1 ``rep-finite <=>
    brick-finite`` (char != 2)."""
    return brick_finite_certificate(triple, budget=budget, field=field)["rep_finite"]


def support_tau_tilting(triple, budget=512, field=None):
    """Support tau-tilting on the split algebra via P45 ``tau_tilting_block`` (the
    exchange graph, g-matrices, brick-labelled Hasse edges, counts; honest complete-iff
    contract).  Stamps ``n = |Q_0| + |Sp|`` (the split vertex count = the geometric
    rank, HZZ sec 6 ``|R| = |Q_0| + |Sp|``).

    Computed on the char-FREE split model over QQ, exactly as ``_analyze`` /
    ``brick_finite_certificate`` do: the exchange-graph BFS leans on ``is_isomorphic`` /
    ``decompose`` (rigorous only char 0 / char > dim M), and the split model is
    characteristic-free, so the pair count is a presentation invariant.  Running it over
    the caller's ``field`` (e.g. GF(2)) would compute over an unsound field -- M3; the
    ``field`` argument is accepted for API symmetry only and does not change the count."""
    from quiverlab.tautilting.block import tau_tilting_block

    triple = _as_triple(triple)
    A = SkewGentleAlgebra(triple, field=QQ)          # char-free model (M3), never GF(2)
    blk = tau_tilting_block(A, budget=budget)
    blk["n"] = len(list(A.quiver.vertices))         # = |Q_0| + |Sp|
    return blk
