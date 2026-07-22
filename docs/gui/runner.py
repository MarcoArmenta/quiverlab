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
        # The downloadable bundle must not silently omit failures.
        # (No bundle exists until run_build resets results to a list.)
        if _state["results"] is not None:
            _state["results"].append({"invariant": spec, "error": out["error"]})
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


# --- Plan 11: wait-time estimation (pure arithmetic; spec 2026-07-22) --------

# Fitted on THIS machine (native, numba BLOCKED to match Pyodide's pure path),
# 2026-07-22, scripts/fit_eta_model.py. Worst-off factor on heavy (>0.3 s)
# grid points: bar 1.46x, fast 3.76x — inside one bucket width; the in-flight
# rescale absorbs the rest. Units are calibrated native-seconds; the browser
# factor comes from calibrate(). Do not hand-tune: rerun the fit script.
ETA_MODEL = {
    "bar":  {"alpha": 1.4622e-07, "p": 1.3},
    "fast": {"alpha": 5.3447e-07, "p": 1.1},
    "scalars": {"cartan": 0.01, "coxeter_polynomial": 0.2,
                "center": 0.05, "global_dimension": 0.5},
}
_MAX_CELLS = 4_000_000        # the library's bar guard (frozen contract)
_BUCKETS = (                  # (upper bound in seconds, id, label)
    (15.0, "seconds", "estimated: a few seconds"),
    (75.0, "minute", "estimated: under a minute"),
    (360.0, "minutes", "estimated: a few minutes"),
    (None, "long", "estimated: could be long — Cancel anytime"),
)


def _bar_guard_cells(m, n):
    """rows*cols of the bar coboundary d^n guard: (m(m-1)^{n+1}) * (m(m-1)^n)."""
    return (m * (m - 1) ** (n + 1)) * (m * (m - 1) ** n)


def _cap_degree(m, top):
    """First degree in 0..top whose guard exceeds the library's max_cells."""
    if m <= 2:
        return None               # (m-1) <= 1: sizes stay tiny forever
    for n in range(top + 1):
        if _bar_guard_cells(m, n) > _MAX_CELLS:
            return n
    return None


def _hh_units(m, top, route):
    mdl = ETA_MODEL[route]
    return mdl["alpha"] * sum(_bar_guard_cells(m, n) ** mdl["p"]
                              for n in range(top + 1))


def _units_for(dim, field_spec, compute):
    """(total units, per-invariant breakdown, cap info or None)."""
    route = ("fast" if field_spec.get("kind") == "GF"
             and field_spec.get("n", 1) == 1 else "bar")
    total, breakdown, cap = 0.0, [], None
    for spec in compute:
        name, top = _parse_compute(spec)
        if name in ("hh_cohomology", "hh_homology"):
            k = _cap_degree(dim, top)
            if k is not None and (cap is None or k < cap["degree"]):
                cap = {"degree": k, "invariant": spec}
            u = _hh_units(dim, min(top, (k - 1) if k is not None else top), route)
        else:
            u = ETA_MODEL["scalars"].get(name, 0.1)
        total += u
        breakdown.append({"invariant": spec, "units": u})
    return total, breakdown, cap


def bucket_for_seconds(seconds):
    for bound, bid, label in _BUCKETS:
        if bound is None or seconds < bound:
            return json.dumps({"bucket": bid, "label": label})
    raise AssertionError("unreachable")


def estimate(factor):
    """Estimate the CURRENT request against the just-built algebra."""
    try:
        A, req = _state["algebra"], _state["request"]
        if A is None or req is None:
            raise RequestError("no algebra built (run_build first)")
        units, breakdown, cap = _units_for(
            A.dim, req["algebra"]["field"], req.get("compute", []))
        seconds = units * float(factor)
        if cap is not None:
            bucket, label = "cap", ("will hit the engine's cell cap near "
                                    "degree %d" % cap["degree"])
        else:
            b = json.loads(bucket_for_seconds(seconds))
            bucket, label = b["bucket"], b["label"]
        out = {"ok": True, "dim": A.dim, "units": units, "seconds": seconds,
               "bucket": bucket, "label": label,
               "cap_degree": cap["degree"] if cap else None,
               "breakdown": breakdown}
    except Exception as exc:
        out = _fail(exc)
    return json.dumps(out)


_CAL_FIELD = {"kind": "GF", "p": 2, "n": 1}
_CAL_COMPUTE = ["hh_cohomology:0..6", "cartan", "center"]


def calibrate():
    """Time a fixed workload; factor = seconds per model unit on THIS machine.
    Builds locally (never via _state) so a visitor's probe state survives."""
    import time
    Q = quiverlab.Quiver(vertices=[1], arrows={"x": (1, 1)})
    t0 = time.monotonic()
    A = Q.algebra(relations=["x*x*x"], field=quiverlab.GF(2))
    A.hochschild_cohomology(6, verbose=False)
    A.cartan_matrix()
    A.center()
    seconds = time.monotonic() - t0
    units, _, _ = _units_for(3, _CAL_FIELD, _CAL_COMPUTE)
    return json.dumps({"seconds": seconds, "units": units,
                       "factor": seconds / units})
