"""quiverlab GUI runner (Plan 10): executes landing-page GUI requests.

Runs IDENTICALLY under CPython (pytest, tests/gui/) and Pyodide (the docs-site
Web Worker) — this exact file is shipped into the built site and imported in
the browser. Import policy: the public quiverlab surface + the sanctioned trace
helpers (render_html / render_json / references_for / resolve_references);
pinned by tests/gui/test_interface_freshness.py. quiverlab.engine.* is forbidden.

All public functions take and return JSON STRINGS (postMessage-friendly)."""
import json
import os
import traceback
from fractions import Fraction

os.environ.setdefault("MPLBACKEND", "Agg")   # never let matplotlib probe for a display

import quiverlab

SCHEMA_VERSION = 1
# The GUI tags requests carrying module / Ext / Tor blocks as schema 2 (the
# webapp validator's rule); this runner's dispatch handles both identically.
ACCEPTED_SCHEMAS = (1, 2)
MAX_DEGREE = 10
# Depth to which projective dimension is probed before reporting "infinite"
# (matches the library's injective_dimension(bound=32) default).
_PD_BOUND = 32

# Module compute kinds (Plan 26; Plan 30 adds `tor`/`decompose`). `ext`/`tor` also
# need a second module (`ext_target`/`tor_target`, the latter a LEFT A-module).
_MODULE_KINDS = frozenset({
    "dimension_vector", "rad_top_soc", "ext", "tor", "tau", "tau_minus",
    "projective_resolution", "injective_resolution",
    "projective_dimension", "injective_dimension", "decompose",
})

_state = {"algebra": None, "request": None, "events": None, "results": None,
          "module": None, "ext_target": None, "tor_target": None}


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
    _state.update(algebra=None, request=None, events=[], results=[],
                  module=None, ext_target=None, tor_target=None)
    quiverlab.verbose = False   # the GUI renders its own report; never write trace files
    try:
        req = json.loads(request_json)
        if req.get("schema") not in ACCEPTED_SCHEMAS:
            raise RequestError("unsupported schema %r (this GUI speaks schema 1/2)"
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
        # Module blocks (Plan 26) ride alongside the algebra; the module itself is
        # built lazily in compute_one, so a relation-violating matrix surfaces as a
        # per-computation error (rendered on the page), never a build crash.
        _state.update(algebra=A, request=req, module=req.get("module"),
                      ext_target=req.get("ext_target"),
                      tor_target=req.get("tor_target"))
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


# --- module block (Plan 26): build a module from the request + dispatch ------
# The webapp server tier (webapp/server/runner.py) carries the SAME logic; this
# is the Pyodide/client copy (the two runners already duplicate the algebra
# dispatch -- they cannot import each other).

def _parse_entry(x):
    """Parse an exact matrix entry to int/Fraction. Entries are DATA (never
    evaluated); a float is refused (exactness is the point)."""
    if isinstance(x, bool):
        raise RequestError("module entry %r is not a number" % (x,))
    if isinstance(x, int):
        return x
    if isinstance(x, float):
        raise RequestError("module entry %r is a float; entries must be exact "
                           "(integers or strings like '1/2')" % (x,))
    if isinstance(x, str):
        s = x.strip()
        try:
            return int(s)
        except ValueError:
            pass
        try:
            return Fraction(s)
        except (ValueError, ZeroDivisionError):
            raise RequestError("module entry %r is not an exact integer or "
                               "fraction" % (x,))
    raise RequestError("module entry %r is not a number" % (x,))


def _full_matrices(A, mspec):
    """Expand per-arrow BLOCK matrices (dim[target] x dim[source] in the
    representation quiver) into the full vertex-ordered arrow actions A.module
    consumes. side selects the representation quiver (right=A, left=A^op)."""
    rep = A if mspec.get("side", "right") == "right" else A.opposite()
    verts = list(rep.quiver.vertices)
    by_str = {str(v): v for v in verts}
    dimvec = {v: 0 for v in verts}
    for key, n in (mspec.get("dims") or {}).items():
        if key not in by_str:
            raise RequestError("module: no vertex %r in the algebra" % (key,))
        dimvec[by_str[key]] = int(n)
    start, off = {}, 0
    for v in verts:
        start[v] = off
        off += dimvec[v]
    n = off
    arrow_names = list(rep.quiver.arrows)
    for a in (mspec.get("maps") or {}):
        if a not in arrow_names:
            raise RequestError("module: no arrow %r in the algebra" % (a,))
    action = {}
    for a in arrow_names:
        s, t = rep.quiver.source(a), rep.quiver.target(a)
        rows, cols = dimvec[t], dimvec[s]
        full = [[0] * n for _ in range(n)]
        block = (mspec.get("maps") or {}).get(a)
        if block is not None:
            if len(block) != rows or any(len(r) != cols for r in block):
                raise RequestError("module map %r must be %dx%d (target x source "
                                   "dims)" % (a, rows, cols))
            for i in range(rows):
                for j in range(cols):
                    full[start[t] + i][start[s] + j] = _parse_entry(block[i][j])
        action[a] = full
    return dimvec, action


def _build_module(A, mspec, name):
    if not isinstance(mspec, dict):
        raise RequestError("this computation needs a module block")
    b = mspec.get("builtin")
    if b is not None:
        v, kind, side = b.get("vertex"), b.get("kind"), mspec.get("side", "right")
        builder = {"simple": A.simple, "projective": A.projective,
                   "injective": A.injective}.get(kind)
        if builder is None:
            raise RequestError("unknown builtin module kind %r" % (kind,))
        for vv in A.quiver.vertices:
            if vv == v or str(vv) == str(v):
                return builder(vv, side=side)
        raise RequestError("module: no vertex %r in the algebra" % (v,))
    dimvec, action = _full_matrices(A, mspec)
    return A.module(dimvec, action, side=mspec.get("side", "right"), name=name)


def _tor_target_spec():
    """The Tor second module N (a LEFT A-module). An omitted side defaults to
    ``"left"`` -- mirrors the schema/spec-core default -- so the GUI need not force
    the toggle. A builtin's own nested side is honored as-is."""
    ts = _state.get("tor_target")
    if isinstance(ts, dict):
        b = ts.get("builtin")
        has_side = "side" in ts or (isinstance(b, dict) and "side" in b)
        if not has_side:
            ts = dict(ts, side="left")
    return ts


def _dv(dimvec):
    return {str(v): int(n) for v, n in sorted(dimvec.items(), key=lambda kv: str(kv[0]))}


def _dv_latex(dimvec):
    d = _dv(dimvec)
    return "(" + ",\\, ".join(str(d[k]) for k in d) + ")" if d else "()"


def _homdim_latex(op, value):
    """Display latex for a homological-dimension block (``pd``/``id``), mirroring
    ``quiverlab.hpc.spec._homdim_latex``. An UNRESOLVED probe is not a proof of
    infinity -- the resolution merely did not terminate by ``_PD_BOUND`` -- so it
    states the certified lower bound, exactly as ``global_dimension`` does."""
    if value is None:
        return r"\operatorname{%s} M > %d" % (op, _PD_BOUND)
    return r"\operatorname{%s} M = %d" % (op, value)


_HOMDIM_UNRESOLVED = ("certified lower bound; the resolution did not terminate "
                      "within the probed depth %d" % _PD_BOUND)


def _mod_view(m):
    return {"dimvec": _dv(m.dimension_vector()), "dim": m.dim}


def _mod_repr(m):
    """A module as the no-code INPUT schema ``{"dims": {v: n}, "maps": {arrow: [[..]]}}``
    (Plan 34, Marco): the per-vertex dimension VECTOR + the exact per-arrow action
    matrices, redundant total dim dropped. Routed through the SAME library serializer as
    the server tier (``quiverlab.modules.qpa_module.module_blocks``) so the two runners
    cannot drift."""
    from quiverlab.modules.qpa_module import module_blocks
    return module_blocks(m)


_MOD_REFS = {
    "dimension_vector": ["assem_book"], "rad_top_soc": ["assem_book"],
    "tau": ["assem_book"], "tau_minus": ["assem_book"], "ext": ["module_ext"],
    "tor": ["minimal_resolution", "module_ext"],
    "projective_resolution": ["minimal_resolution"],
    "projective_dimension": ["minimal_resolution"],
    "injective_resolution": ["minimal_resolution", "assem_book"],
    "injective_dimension": ["minimal_resolution", "assem_book"],
    "decompose": ["assem_book"],
}


# The additive-tau note carried by an AR-translate block whose input decomposes.
# The block ships a stable KEY; each renderer localizes it (app.js via i18n, gui.js
# via its own English map) -- so the block stays language-neutral data.
_TAU_ADDITIVE_KEY = "mod.tau_additive"


def _summands_latex(vertices, letter):
    """A resolution term's summand multiset as LaTeX, e.g. ``P_{1}^{2} \\oplus P_{3}``
    (``letter`` = ``P`` projectives / ``I`` injectives). ``0`` for the zero term."""
    if not vertices:
        return "0"
    counts = {}
    for v in vertices:
        counts[v] = counts.get(v, 0) + 1

    def key(v):
        return (0, v) if isinstance(v, int) and not isinstance(v, bool) else (1, str(v))
    parts = []
    for v in sorted(counts, key=key):
        base = "%s_{%s}" % (letter, v)
        parts.append(base if counts[v] == 1 else "%s^{%d}" % (base, counts[v]))
    return " \\oplus ".join(parts)


# mirrors quiverlab.hpc.spec._differential_blocks (same cap, same shape)
_MAX_DIFF_CELLS = 250_000


def _differential_blocks(res, n_terms):
    """The resolution's maps as exact matrices (rows: target basis, columns:
    source basis). Projective: entry 0 = eps: Q_0 -> M, entry n = d_n: Q_n ->
    Q_{n-1}. Injective: entry 0 = iota: M -> E^0, entry n = d^n: E^{n-1} -> E^n."""
    out = []
    dmats = getattr(res, "dmats", None) or []
    for n in range(min(n_terms, len(dmats))):
        D = res.differential(n)
        nrows = len(D)
        ncols = len(D[0]) if nrows else 0
        if nrows * ncols > _MAX_DIFF_CELLS:
            out.append({"rows": nrows, "cols": ncols, "elided": True})
        else:
            out.append({"rows": nrows, "cols": ncols,
                        "matrix": [[str(x) for x in row] for row in D]})
    return out


def _term_basis_blocks(res, kind, M):
    """The ordered k-basis of each resolution term, as the concatenated path bases of
    its projective / injective summands (Plan 35 UNIT 2). Shape-identical to the server
    tier (``quiverlab.hpc.spec._term_basis_blocks``): ``term_basis[n]`` lists one label
    per basis vector of term n, so ``len(term_basis[n])`` = the term dimension (the
    differential's column count for a projective resolution, row count for an injective
    one). ``None`` (field omitted) when the algebra carries no path basis or the total
    is implausibly large -- renderers tolerate the absence."""
    try:
        from quiverlab.modules.builders import projective
        if kind == "projective_resolution":
            base, sym = M.algebra, "P"
        else:
            from quiverlab.modules.opposite import opposite_algebra
            base, sym = opposite_algebra(M.algebra), "I"
        cache = {}

        def paths(v):
            if v not in cache:
                cache[v] = list(projective(base, v)._pv_basis_labels)
            return cache[v]

        out, total = [], 0
        for n in range(len(res.terms)):
            labels = []
            for v in res.term(n):
                labels.extend("%s_%s: %s" % (sym, v, p) for p in paths(v))
            total += len(labels)
            if total > _MAX_DIFF_CELLS:               # implausible; omit (bloat guard)
                return None
            out.append(labels)
        return out
    except Exception:
        return None


def _summand_view(mod, mult):
    """One indecomposable summand: dimension vector, multiplicity, and either a
    STANDARD name (S_v / P_v / I_v) or its full per-arrow matrices (Marco
    2026-07-29). Same library serializer as the server tier -- no drift."""
    from quiverlab.modules.qpa_module import summand_blocks
    return summand_blocks(mod, mult)


def _decompose_engine():
    """(is_indecomposable, decompose) from the Plan-30 Part-A engine, or ``None``
    when not yet importable (so callers degrade to the pre-Plan-30 block shape)."""
    try:
        from quiverlab.modules.decompose import decompose, is_indecomposable
    except ImportError:
        return None
    return is_indecomposable, decompose


def _input_certificate(M):
    """Certify the INPUT of an AR translate (Marco #1): indecomposable, or its
    Krull-Schmidt decomposition + the additivity note. Best-effort + honest: ``{}``
    (no claim) when the decompose engine is unavailable OR cannot certify within
    budget (it is LOUD, e.g. char <= dim M) -- tau itself stays valid, so the block
    omits the certificate rather than assert what it could not prove."""
    eng = _decompose_engine()
    if eng is None:
        return {}
    is_indec, decompose = eng
    try:
        if is_indec(M):
            return {"indecomposable": True}
        return {"indecomposable": False,
                "decomposition": [_summand_view(s, m) for (s, m) in decompose(M)],
                "note_key": _TAU_ADDITIVE_KEY}
    except quiverlab.QuiverlabError:
        return {}


def _ar_translate(mod, kind, name):
    """One AR translate as a self-contained display payload: the symbol, its
    dimension-vector latex, the FULL representation ({dims, maps}) and the input's
    indecomposability certificate. Mirrors ``quiverlab.hpc.spec._ar_translate`` so
    the two runners ship the same shape."""
    t = mod.tau() if kind == "tau" else mod.tau_minus()
    sym = (r"\tau %s" if kind == "tau" else r"\tau^{-} %s") % name
    latex = ((sym + " = 0") if t.dim == 0
             else (r"\underline{\dim}\, " + sym + " = "
                   + _dv_latex(t.dimension_vector())))
    out = {"name": name, "side": t.side, "is_zero": t.dim == 0,
           "latex": latex, **_mod_view(t)}
    if t.dim > 0:
        out["repr"] = _mod_repr(t)
    out.update(_input_certificate(mod))
    return out


def _target_translates(A, kind):
    """The AR translates of the SECOND module(s) N the request names (the Ext/Tor
    argument), so a tau block covers every module in play, with full matrices
    (Marco, 2026-07-29). A loud refusal on one target becomes an honest error entry
    -- tau M is already computed and stays valid."""
    out = []
    for role, spec in (("ext_target", _state.get("ext_target")),
                       ("tor_target", _tor_target_spec())):
        if not isinstance(spec, dict):
            continue
        try:
            entry = _ar_translate(_build_module(A, spec, "N"), kind, "N")
        except (quiverlab.QuiverlabError, RequestError) as exc:
            entry = {"name": "N", "error": str(exc)}
        entry["role"] = role
        out.append(entry)
    return out


def _module_block(name, top):
    """Dispatch one module compute kind against the built module(s). Blocks carry
    a `latex` display and `citations`, like the algebra invariants."""
    A = _state["algebra"]
    M = _build_module(A, _state.get("module"), "M")
    keys = _MOD_REFS[name]
    cites = _citation_pairs(keys)
    if name == "dimension_vector":
        return {"kind": name, "side": M.side, "citations": cites, **_mod_view(M),
                "latex": r"\underline{\dim}\, M = " + _dv_latex(M.dimension_vector())}
    if name == "rad_top_soc":
        return {"kind": name, "side": M.side, "citations": cites,
                "radical": _mod_repr(M.radical()), "top": _mod_repr(M.top()),
                "socle": _mod_repr(M.socle())}
    if name in ("tau", "tau_minus"):
        # The translate ships as a full representation ({dims, maps}) -- mirrors the
        # hpc spec core, so both dispatches carry the AR translate's per-arrow
        # matrices -- together with the input certificate (Marco #1).
        entry = _ar_translate(M, name, "M")
        entry.pop("name", None)
        block = {"kind": name, "citations": cites, **entry}
        # ... and the same for the SECOND module N when the request names one
        # (Marco, 2026-07-29): a tau block covers every module in the request.
        targets = _target_translates(A, name)
        if targets:
            block["targets"] = targets
        return block
    if name == "decompose":
        eng = _decompose_engine()
        if eng is None:
            raise RequestError("module decomposition needs the Krull-Schmidt engine "
                               "(Plan 30 Part A); not available in this build")
        _, decompose = eng
        summands = [_summand_view(s, m) for (s, m) in decompose(M)]
        return {"kind": name, "side": M.side, "summands": summands,
                "iso_classes": len(summands), "citations": cites}
    if name == "ext":
        if top is None:
            raise RequestError("ext needs a range, e.g. 'ext:0..4'")
        from quiverlab.modules.ext import ext_dims
        N = _build_module(A, _state.get("ext_target"), "N")
        # Plan 35 wave 3a: explicit Ext cocycle representatives + self-cert data,
        # additive keys shared byte-for-byte with the hpc spec runner.
        # Plan 35 wave 3c: interpret=True captures the Yoneda exact sequence of each class.
        raw, reps = ext_dims(A, M, N, top, with_reps=True, interpret=True)
        block = {"kind": name, "top": top, "dims": [int(d) for d in raw],
                 "target": _mod_view(N), "citations": cites}
        block.update(reps)
        return block
    if name == "tor":
        if top is None:
            raise RequestError("tor needs a range, e.g. 'tor:0..4'")
        try:
            from quiverlab.modules.tor import tor_dims
        except ImportError:
            raise RequestError("Tor requires the Tor engine (Plan 29); not available "
                               "in this build")
        N = _build_module(A, _tor_target_spec(), "N")
        raw, reps = tor_dims(A, M, N, top, with_reps=True)
        block = {"kind": name, "top": top, "dims": [int(d) for d in raw],
                 "target": _mod_view(N), "citations": cites}
        block.update(reps)
        return block
    if name in ("projective_resolution", "injective_resolution"):
        if top is None:
            raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
        res = (M.projective_resolution(top) if name == "projective_resolution"
               else M.injective_resolution(top))
        letter = "P" if name == "projective_resolution" else "I"
        terms = [_dv(dv) for dv in res.dimension_vectors()]
        block = {"kind": name, "top": top, "terms": terms,
                 "betti": [res.betti(i) for i in range(len(terms))],
                 "summands": [_summands_latex(res.term(i), letter)
                              for i in range(len(terms))],
                 "differentials": _differential_blocks(res, len(terms)),
                 "citations": cites}
        term_basis = _term_basis_blocks(res, name, M)
        if term_basis is not None:
            block["term_basis"] = term_basis
        if name == "projective_resolution":
            block["pd"] = res.pd()
        else:
            block["injective_dimension"] = res.injective_dimension()
        return block
    if name == "projective_dimension":
        pd = M.projective_resolution(_PD_BOUND).pd()
        return {"kind": name, "value": pd, "finite": pd is not None, "citations": cites,
                "bound": _PD_BOUND, "latex": _homdim_latex("pd", pd),
                **({} if pd is not None else {"note": _HOMDIM_UNRESOLVED})}
    if name == "injective_dimension":
        idim = M.injective_dimension(bound=_PD_BOUND)
        return {"kind": name, "value": idim, "finite": idim is not None, "citations": cites,
                "bound": _PD_BOUND, "latex": _homdim_latex("id", idim),
                **({} if idim is not None else {"note": _HOMDIM_UNRESOLVED})}
    raise RequestError("unknown module invariant %r" % (name,))


# Module kinds with a Plan-30 Part-C worked-steps hook (quiverlab.trace.modules):
# when the report is requested (artifacts.pdf), a module compute also emits the
# exhaustive step events, so the HTML/.tex bundle covers modules like HH does.
_MODULE_TRACE_KINDS = frozenset({
    "projective_resolution", "injective_resolution", "ext", "tau", "tau_minus",
})


def _wants_trace():
    req = _state.get("request") or {}
    return bool((req.get("artifacts") or {}).get("pdf"))


def _emit_module_trace(name, top):
    """Append the worked-step events for a traceable module compute into
    ``_state['events']`` (only when the report is requested). Best-effort: the
    compute already succeeded, so a trace re-run failure silently skips the report
    enhancement rather than break the (already-rendered) result."""
    if name not in _MODULE_TRACE_KINDS or not _wants_trace():
        return
    A = _state["algebra"]
    try:
        from quiverlab.trace import modules as tm
        M = _build_module(A, _state.get("module"), "M")
        if name == "projective_resolution":
            events, _ = tm.trace_projective_resolution(M, top)
        elif name == "injective_resolution":
            events, _ = tm.trace_injective_resolution(M, top)
        elif name == "ext":
            N = _build_module(A, _state.get("ext_target"), "N")
            events, _ = tm.trace_ext(A, M, N, top)
        else:                                       # tau / tau_minus
            events, _ = tm.trace_tau(M, kind=name)
        if _state["events"] is not None:
            _state["events"].extend(events)
    except Exception:                               # best-effort report adornment
        pass


def compute_one(spec):
    """Run ONE Plan-09 compute string against the built algebra."""
    A = _state["algebra"]
    try:
        if A is None:
            raise RequestError("no algebra built (run_build first)")
        name, top = _parse_compute(spec)
        if name in _MODULE_KINDS:
            block = _module_block(name, top)
            _emit_module_trace(name, top)
        elif name in ("hh_cohomology", "hh_homology"):
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            method = (A.hochschild_cohomology if name == "hh_cohomology"
                      else A.hochschild_homology)
            table = method(top, verbose=False, trace=_state["events"])
            block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                     "engine": table.engine,
                     "citations": _citation_pairs(table.references)}
            # Plan 35 wave 3d: capture the explicit HH^n / HH_n representatives alongside
            # the dims (basis_classes / chain_basis / differentials / inner_dims per
            # degree) from the SAME dims path -- key-for-key identical to the server twin
            # (quiverlab.hpc.spec._dispatch), so the cross-runner contract holds. None
            # (dims-only) when no representative route applies.
            from quiverlab.hochschild.hh_reps import hh_reps_blocks
            try:                               # reps are ADDITIVE + best-effort: a
                reps = hh_reps_blocks(A, name, top, list(table.dims), table.engine)
            except Exception:                  # capture must NEVER break the dims block
                reps = None
            if reps:
                block.update(reps)
        elif name == "cyclic_homology":
            # Plan-35 follow-up: cyclic homology HC_0..HC_n (Connes (b, B) mixed
            # complex). Range kind; the block is key-for-key identical to the server
            # twin (quiverlab.hpc.spec._dispatch) -- kind/top/dims/engine/references/
            # citations -- so the cross-runner contract holds.
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            # Plan 35 wave 3b: capture the explicit HC representatives alongside the
            # dims (basis_classes / chain_basis / differentials / column_structure) from
            # the SAME (b, B) total complex -- key-for-key identical to the server twin
            # (quiverlab.hpc.spec._dispatch), so the cross-runner contract holds.
            table, reps = A.cyclic_homology(top, with_reps=True)
            keys = ["cyclic"]
            block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                     "engine": table.engine, "references": keys,
                     "citations": _citation_pairs(keys)}
            block.update(reps)
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
        elif name == "dimension":
            # Parity with the server/HPC runner (quiverlab.hpc.spec._dispatch),
            # which serves `dimension` = A.dim; same `value` semantics, GUI block
            # shape (value + citations, like global_dimension).
            block = {"value": A.dim, "citations": _citation_pairs(A.citations())}
        elif name in ("cup", "cap", "bracket", "connes_b"):
            # HH product surface (Plan 35): cup / cap / bracket / connes_b. Each
            # library method returns a frozen result whose .blocks() IS the block
            # dict (kind/top/engine + tables|matrices + references); we only add the
            # resolved citation pairs, exactly as the server twin does
            # (quiverlab.hpc.spec._dispatch). The block keeps `references` -- the
            # cross-runner contract asserts key-for-key equality with the server.
            if top is None:
                raise RequestError("%s needs a range, e.g. '%s:0..4'" % (name, name))
            method = {"cup": A.cup_products, "cap": A.cap_products,
                      "bracket": A.gerstenhaber_brackets,
                      "connes_b": A.connes_differentials}[name]
            block = method(top).blocks()
            block["citations"] = _citation_pairs(block["references"])
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


def _named_modules():
    """The modules the request named, as ``(label, Module)`` for the report's "The
    modules" section (mirrors ``quiverlab.hpc.spec._named``). Best-effort: a module
    that no longer builds is skipped -- the report describes, it does not verify."""
    A = _state.get("algebra")
    if A is None:
        return []
    specs = [("M", _state.get("module")), ("N", _state.get("ext_target")),
             ("N", _tor_target_spec())]
    named, seen_n = [], 0
    for label, mspec in specs:
        if not isinstance(mspec, dict):
            continue
        if label == "N":
            seen_n += 1
        try:
            named.append((label, _build_module(A, mspec, label)))
        except Exception:
            continue
    if seen_n > 1:                        # both an Ext and a Tor target: tell them apart
        roles = iter(("N (Ext target)", "N (Tor target)"))
        named = [(next(roles) if lbl == "N" else lbl, m) for lbl, m in named]
    return named


def trace_html():
    """The worked-steps report as an HTML string.

    Carries the session's COMPUTED RESULT BLOCKS as well as the worked steps, so the
    saved report is everything the page showed (Marco 2026-07-29) -- and is therefore
    produced even when nothing recorded worked steps (Cartan + centre, say).

    '' when the request did not ask for a report (``artifacts.pdf``), or when there
    is nothing at all to report."""
    events = _state["events"] or []
    results = _state["results"] or []
    if not events and not (results and _wants_trace()):
        return ""
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_html import render_html
    title = "Worked steps — %s" % repr(_state["algebra"]).splitlines()[0]
    return render_html(list(events), title=title, algebra=_state["algebra"],
                       references=resolve_references(references_for(events)),
                       results=results, modules=_named_modules())


def trace_json():
    """The worked-steps report as the JSON machine record ('' when nothing was traced).

    Plan 34: the third mandated artifact (PDF/HTML + this) -- the complete event
    stream, deterministic and schema-versioned. A GUI download button hooks this
    accessor beside trace_tex()/trace_html(). Uses the SAME (events, title,
    references) inputs as the other two renderers, so the byte-for-byte machine
    record matches what the writer persists for the same computation."""
    events = _state["events"] or []
    if not events:
        return ""
    from quiverlab.trace.provenance import references_for, resolve_references
    from quiverlab.trace.render_json import render_json
    title = "Worked steps — %s" % repr(_state["algebra"]).splitlines()[0]
    return render_json(list(events), title=title, algebra=_state["algebra"],
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
    if req.get("module"):
        lines += _module_snippet_lines(req["module"], "M")
    if req.get("ext_target"):
        lines += _module_snippet_lines(req["ext_target"], "N")
    if req.get("tor_target"):
        lines += _module_snippet_lines(_tor_target_spec(), "N")
    calls = {"hh_cohomology": "A.hochschild_cohomology(%d)",
             "hh_homology": "A.hochschild_homology(%d)",
             "cyclic_homology": "A.cyclic_homology(%d)",
             "cartan": "A.cartan_matrix()", "coxeter_polynomial": "A.coxeter_polynomial()",
             "global_dimension": "A.global_dimension()", "center": "A.center()",
             # `dimension` is a scalar invariant compute_one serves (A.dim) -- it MUST
             # have a snippet entry or python_snippet() KeyErrors when it is requested.
             "dimension": "A.dim",
             # HH product surface (Plan 35): same four calls as the server snippet
             # map (quiverlab.hpc.spec._snippet); each needs a range (%d = top).
             "cup": "A.cup_products(%d)", "cap": "A.cap_products(%d)",
             "bracket": "A.gerstenhaber_brackets(%d)",
             "connes_b": "A.connes_differentials(%d)",
             "dimension_vector": "M.dimension_vector()",
             "rad_top_soc": "(M.radical(), M.top(), M.socle())",
             "tau": "M.tau()", "tau_minus": "M.tau_minus()",
             "ext": "[A.ext(M, N, i) for i in range(%d + 1)]",
             "tor": "tor_dims(A, M, N, %d)  # from quiverlab.modules.tor",
             "decompose": "M.decompose()",
             "projective_resolution": "M.projective_resolution(%d).dimension_vectors()",
             "injective_resolution": "M.injective_resolution(%d).dimension_vectors()",
             "projective_dimension": "M.projective_resolution(%d).pd()" % _PD_BOUND,
             "injective_dimension": "M.injective_dimension()"}
    for spec in req.get("compute", []):
        name, top = _parse_compute(spec)
        tmpl = calls[name]
        call = tmpl % top if ("%d" in tmpl and top is not None) else tmpl
        lines.append("print(%s)" % call)
    return "\n".join(lines) + "\n"


def _fmt_scalar(x):
    if isinstance(x, Fraction):
        return (str(x.numerator) if x.denominator == 1
                else "Fraction(%d, %d)" % (x.numerator, x.denominator))
    return str(x)


def _pymat(mat):
    return ("[" + ", ".join("[" + ", ".join(_fmt_scalar(x) for x in row) + "]"
                            for row in mat) + "]")


def _module_snippet_lines(mspec, varname):
    """Runnable lines rebuilding `varname`: a builtin call, or the FULL exact
    arrow-action matrices A.module consumes."""
    b = mspec.get("builtin")
    if b is not None:
        return ['%s = A.%s(%r, side="%s")'
                % (varname, b.get("kind"), b.get("vertex"), mspec.get("side", "right"))]
    dimvec, action = _full_matrices(_state["algebra"], mspec)
    lines = []
    if any(isinstance(x, Fraction) for m in action.values() for row in m for x in row):
        lines.append("from fractions import Fraction")
    dv_lit = "{" + ", ".join("%r: %d" % (v, n) for v, n in dimvec.items()) + "}"
    maps_lit = "{" + ", ".join('"%s": %s' % (a, _pymat(m)) for a, m in action.items()) + "}"
    lines.append('%s = A.module(%s, %s, side="%s", name="%s")'
                 % (varname, dv_lit, maps_lit, mspec.get("side", "right"), varname))
    return lines


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
                "center": 0.05, "global_dimension": 0.5,
                # module kinds (Plan 26): cheap dim-vector reads up to
                # resolution/dimension probes that build syzygies to depth.
                "dimension_vector": 0.02, "rad_top_soc": 0.05,
                "tau": 0.1, "tau_minus": 0.1, "ext": 0.2, "tor": 0.2,
                "decompose": 0.3,
                "projective_resolution": 0.2, "injective_resolution": 0.2,
                "projective_dimension": 0.3, "injective_dimension": 0.3},
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
            if k is not None:
                # The homology engine's own DepthLimitError names b_{k+1} —
                # boundary indexing is shifted by one vs the coboundary guard —
                # so the SHOWN degree gets +1 for hh_homology (units keep raw k).
                shown = k + 1 if name == "hh_homology" else k
                if cap is None or shown < cap["degree"]:
                    cap = {"degree": shown, "invariant": spec}
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
