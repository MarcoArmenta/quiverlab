"""quiverlab GUI runner (Plan 10): executes landing-page GUI requests.

Runs IDENTICALLY under CPython (pytest, tests/gui/) and Pyodide (the docs-site
Web Worker) — this exact file is shipped into the built site and imported in
the browser. Import policy: the public quiverlab surface + three sanctioned
trace helpers (render_html / references_for / resolve_references); pinned by
tests/gui/test_interface_freshness.py. quiverlab.engine.* is forbidden.

All public functions take and return JSON STRINGS (postMessage-friendly)."""
import json
import os
import traceback

os.environ.setdefault("MPLBACKEND", "Agg")   # never let matplotlib probe for a display

import quiverlab

SCHEMA_VERSION = 1
MAX_DEGREE = 10

_state = {"algebra": None, "request": None, "events": None, "results": None}


class RequestError(Exception):
    """Invalid GUI request (schema violation) — reported like a library error."""


def _fail(exc):
    if isinstance(exc, (quiverlab.QuiverlabError, RequestError)):
        return {"ok": False,
                "error": {"type": type(exc).__name__, "message": str(exc)}}
    # Unexpected: generic message to the page, full traceback for the console.
    return {"ok": False,
            "error": {"type": "InternalError",
                      "message": "unexpected engine error — details in the browser console"},
            "detail": traceback.format_exc()}


def _field_from_spec(spec):
    kind = spec.get("kind") if isinstance(spec, dict) else None
    if kind == "CC":
        return quiverlab.CC
    if kind == "GF":
        p, n = spec.get("p"), spec.get("n", 1)
        if not (isinstance(p, int) and isinstance(n, int) and p >= 2 and n >= 1):
            raise RequestError("field GF needs integers p >= 2 and n >= 1")
        return quiverlab.GF(p ** n)   # FieldError (p not prime, ...) surfaces verbatim
    raise RequestError("unknown field kind %r (expected 'CC' or 'GF')" % (kind,))


def run_build(request_json):
    """Parse + validate a schema-1 request, build the algebra, reset all state."""
    _state.update(algebra=None, request=None, events=[], results=[])
    quiverlab.verbose = False   # the GUI renders its own report; never write trace files
    try:
        req = json.loads(request_json)
        if req.get("schema") != SCHEMA_VERSION:
            raise RequestError("unsupported schema %r (this GUI speaks schema 1)"
                               % (req.get("schema"),))
        alg = req.get("algebra") or {}
        kind = alg.get("kind")
        if kind == "family":
            raise RequestError("algebra kind 'family' is the server tier (Plan 09); "
                               "this GUI submits kind 'quiver' only")
        if kind != "quiver":
            raise RequestError("unknown algebra kind %r (expected 'quiver')" % (kind,))
        vertices = alg.get("vertices")
        if not (isinstance(vertices, list) and vertices
                and all(isinstance(v, int) for v in vertices)):
            raise RequestError("algebra.vertices must be a non-empty list of integers")
        arrows = alg.get("arrows")
        if not (isinstance(arrows, dict) and all(
                isinstance(st, list) and len(st) == 2
                and all(isinstance(x, int) for x in st)
                for st in arrows.values())):
            raise RequestError("algebra.arrows must map names to [source, target] pairs")
        relations = alg.get("relations", [])
        if not (isinstance(relations, list)
                and all(isinstance(r, str) for r in relations)):
            raise RequestError("algebra.relations must be a list of strings")
        field = _field_from_spec(alg.get("field"))
        Q = quiverlab.Quiver(vertices=vertices,
                             arrows={k: (s, t) for k, (s, t) in arrows.items()})
        A = Q.algebra(relations=relations, field=field)
        _state.update(algebra=A, request=req)
        out = {"ok": True, "dim": A.dim, "n_vertices": len(vertices),
               "n_arrows": len(arrows), "algebra": repr(A).splitlines()[0]}
    except Exception as exc:
        out = _fail(exc)
    return json.dumps(out)
