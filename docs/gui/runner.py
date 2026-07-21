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


def _parse_compute(spec):
    name, _, rng = spec.partition(":")
    top = None
    if rng:
        lo, _, hi = rng.partition("..")
        if lo != "0" or not hi.isdigit():
            raise RequestError("bad compute range %r (expected 'name:0..N')" % (spec,))
        top = int(hi)
        if top > MAX_DEGREE:
            raise RequestError("degree cap is %d (got %d)" % (MAX_DEGREE, top))
    return name, top


def _citation_pairs(keys):
    from quiverlab.trace.provenance import resolve_references
    return [list(p) for p in resolve_references(tuple(keys))]


def _latex_matrix(rows):
    body = r" \\ ".join(" & ".join(str(x) for x in row) for row in rows)
    return r"\begin{pmatrix} %s \end{pmatrix}" % body


def compute_one(spec):
    """Run ONE Plan-09 compute string against the built algebra."""
    A = _state["algebra"]
    try:
        if A is None:
            raise RequestError("no algebra built (run_build first)")
        name, top = _parse_compute(spec)
        if name in ("hh_cohomology", "hh_homology"):
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            method = (A.hochschild_cohomology if name == "hh_cohomology"
                      else A.hochschild_homology)
            table = method(top, verbose=False, trace=_state["events"])
            block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                     "engine": table.engine,
                     "citations": _citation_pairs(table.references)}
        elif name == "cartan":
            m = A.cartan_matrix()
            block = {"matrix": m, "latex": _latex_matrix(m),
                     "citations": _citation_pairs(A.citations())}
        elif name == "coxeter_polynomial":
            import sympy
            p = A.coxeter_polynomial()
            block = {"latex": sympy.latex(p.as_expr()), "text": str(p.as_expr()),
                     "citations": _citation_pairs(A.citations())}
        elif name == "global_dimension":
            g = A.global_dimension()
            block = {"text": str(g), "exact": g.exact, "value": g.value,
                     "citations": _citation_pairs(A.citations())}
        elif name == "center":
            dim_z, basis = A.center()
            # Basis entries are exact ints/rationals (sympy MPQ over CC) — not
            # JSON-serializable; ship them as exact strings.
            block = {"dim": dim_z,
                     "basis": [[str(x) for x in row] for row in basis],
                     "citations": _citation_pairs(A.citations())}
        else:
            raise RequestError("unknown invariant %r" % (name,))
        _state["results"].append(dict(block, invariant=spec))
        out = {"ok": True, "invariant": spec, "block": block}
    except Exception as exc:
        out = _fail(exc)
        out["invariant"] = spec
    return json.dumps(out)


def trace_html():
    """The worked-steps report as an HTML string ('' when nothing was traced)."""
    events = _state["events"] or []
    if not events:
        return ""
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_html import render_html
    title = "Worked steps — %s" % repr(_state["algebra"]).splitlines()[0]
    return render_html(list(events), title=title,
                       references=resolve_references(references_for(events)))


def tikz():
    return "" if _state["algebra"] is None else _state["algebra"].tikz()


def python_snippet():
    """Copy-paste reproduction of the GUI computation (the GUI-to-library bridge)."""
    req = _state["request"]
    if req is None:
        return ""
    alg = req["algebra"]
    f = alg["field"]
    if f["kind"] == "CC":
        field_name, field_expr = "CC", "CC"
    else:
        q = f["p"] ** f.get("n", 1)
        field_name, field_expr = "GF", "GF(%d)" % q
    arrows = ", ".join('"%s": (%d, %d)' % (k, s, t)
                       for k, (s, t) in alg["arrows"].items())
    lines = [
        "from quiverlab import Quiver, %s" % field_name,
        "",
        "Q = Quiver(vertices=%r, arrows={%s})" % (alg["vertices"], arrows),
        "A = Q.algebra(relations=%r, field=%s)" % (list(alg.get("relations", [])),
                                                   field_expr),
        "print(A.dim)",
    ]
    calls = {"hh_cohomology": "A.hochschild_cohomology(%d)",
             "hh_homology": "A.hochschild_homology(%d)",
             "cartan": "A.cartan_matrix()", "coxeter_polynomial": "A.coxeter_polynomial()",
             "global_dimension": "A.global_dimension()", "center": "A.center()"}
    for spec in req.get("compute", []):
        name, top = _parse_compute(spec)
        call = calls[name] % top if "%d" in calls[name] else calls[name]
        lines.append("print(%s)" % call)
    return "\n".join(lines) + "\n"


def result_bundle():
    return json.dumps({"schema": SCHEMA_VERSION, "request": _state["request"],
                       "quiverlab_version": quiverlab.__version__,
                       "results": _state["results"] or []}, indent=1)
