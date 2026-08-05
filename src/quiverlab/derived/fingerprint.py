"""Derived fingerprint (Plan 43): a tuple of derived-invariants / necessary
conditions and an honest comparison. NEVER a derived-equivalence decider -- the
metaplan honest-scope rule: 'distinguished / not distinguished by these
invariants', never '(in)equivalent'.

The Coxeter polynomial and ``|det C|`` are derived invariants; HH/HC/center are
derived invariants (Rickard); the Cartan Smith factors are the two-sided
``GL_n(Z)``-equivalence class of the Cartan matrix -- an integer-equivalence
invariant, a NECESSARY condition, coarser than ``Z``-congruence (do not claim
congruence). Equal fingerprints do NOT imply derived equivalence (the cospectral
trees are the standing counterexample). A field that RAISES on this input
(singular Cartan -> coxeter/smith; presentation-less -> cartan) is captured as
``{"error": <msg>}`` per field, never a crash."""
from __future__ import annotations

import sympy as sp
from sympy.matrices.normalforms import invariant_factors

from quiverlab.errors import QuiverlabError


def _field(fn):
    # DepthLimitError is a QuiverlabError subclass, so a bar/tt blow-up (e.g. the
    # generic cyclic-homology mixed complex over a non-GF(p) field, which has no CS
    # route) is captured as an honest per-field error -- never a crash. Every field is
    # wrapped (the plan's block docstring promises this; the deep ones are not exempt).
    try:
        return fn()
    except QuiverlabError as exc:
        return {"error": str(exc)}


def derived_fingerprint(A, top=4):
    fp = {}
    fp["coxeter_polynomial"] = _field(lambda: str(A.coxeter_polynomial().as_expr()))
    fp["cartan_det"] = _field(lambda: int(sp.Matrix(A.cartan_matrix()).det()))
    fp["cartan_smith"] = _field(
        lambda: [int(x) for x in invariant_factors(sp.Matrix(A.cartan_matrix()))])
    fp["hh_cohomology_dims"] = _field(lambda: list(A.hochschild_cohomology(top).dims))
    fp["hh_homology_dims"] = _field(lambda: list(A.hochschild_homology(top).dims))
    fp["cyclic_dims"] = _field(lambda: list(A.cyclic_homology(top)))
    fp["center_dim"] = _field(lambda: A.center()[0])
    fp["gl_dim"] = _field(lambda: repr(A.global_dimension()))
    return fp


def compare_fingerprints(fa, fb):
    distinguished = []
    incomparable = []
    for key in fa:
        va, vb = fa.get(key), fb.get(key)
        if isinstance(va, dict) or isinstance(vb, dict):   # a field errored one side
            incomparable.append(key)                        # surfaced, never silent
            continue
        if va != vb:
            distinguished.append(key)
    verdict = ("distinguished" if distinguished
               else "not distinguished by these invariants")
    return {"distinguished_by": distinguished, "verdict": verdict,
            "incomparable_fields": incomparable}
