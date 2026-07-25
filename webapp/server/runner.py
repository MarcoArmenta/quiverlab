"""Executes a validated request against the public quiverlab API and writes
artifacts. No user code is executed: for a ``kind: "family"`` request the family
id and params are validated against the introspected catalog, then dispatched to
the matching top-level builder; for a ``kind: "quiver"`` request the algebra is
built from vertices/arrows/opaque relation strings, exactly as the Plan-10 GUI
runner does (``docs/gui/runner.py::run_build``). Relation strings are parsed by
the library's exact relation grammar (loud ``RelationError``).

Execution semantics are adopted from ``docs/gui/runner.py`` (the Plan-10
amendment names it authoritative): result-block shapes mirror its ``compute_one``
(HH tables read ``.kind``/``.dims``/``.engine``/``.references``; matrices via
``_latex_matrix``; ``coxeter_polynomial`` via ``sympy.latex``; ``center`` returns
``(dim, basis)`` with ``str()``-shipped entries; ``global_dimension`` reads
``.exact``/``.value``/``str()``), worked steps via ``quiverlab.trace.render_html``.

Never runs ``verbose=True``: the global ``quiverlab.verbose`` defaults to True and
would write worked-steps PDFs into a cwd-relative ``quiverlab_traces/`` -- the
server pins it to False (before building the algebra) for the whole run and passes
``verbose=False`` to every HH call that accepts it (amendment pt 6). Worked-steps
PDFs are instead produced deterministically via ``quiverlab.trace.writer.write_trace``
with an explicit ``out_dir`` (the job's artifact dir), so nothing is scanned from or
written to the cwd. QuiverlabError subclasses surface verbatim (type + message);
nothing else leaks a traceback into the payload."""
from __future__ import annotations

import json
import logging
from fractions import Fraction
from pathlib import Path
from typing import Callable

import quiverlab as ql
from quiverlab import errors as qerr

from webapp.server.catalog import validate_family, CatalogError, _iter_families
from webapp.server.references import resolve_references
from webapp.server.schema import (
    ComputeRequest, MODULE_KINDS, MODULE_RANGE_KINDS, parse_compute_item,
)

# Depth to which projective dimension is probed before reporting "infinite"
# (matches the library's `injective_dimension(bound=32)` default).
_PD_BOUND = 32

_log = logging.getLogger("quiverlab_web.runner")

# Honest labels for ``meta["pdf"]`` (finding 1). The HTML-fallback string is the
# longest candidate, so using it as the pre-cap provisional value keeps the
# result-size check conservative -- the value actually written after write_trace
# can only be the same length or shorter.
_PDF_OK = "trace.pdf"
_PDF_HTML_FALLBACK = ("PDF toolchain (pdflatex/tectonic) not found -- "
                      "worked steps in trace_steps.html")
_PDF_NO_HH = "no traced computation requested (PDF covers HH worked steps)"


class RunError(Exception):
    """A refusal the caller (worker/sync endpoint) turns into an honest error
    result. ``error_type`` is the library exception class name (verbatim) or one
    of the runner's own tags (``CatalogError``, ``ResultTooLarge``, ...)."""

    def __init__(self, error_type: str, message: str):
        super().__init__(f"{error_type}: {message}")
        self.error_type = error_type
        self.message = message


# --------------------------------------------------------------------------- #
# Field + builder resolution
# --------------------------------------------------------------------------- #

def _field(spec):
    """Build the exact field from the field spec (mirrors the GUI's
    ``_field_from_spec``). ``GF(p**n)`` covers both GF(p) and GF(p^n); a
    non-prime p surfaces as the library's ``FieldError``."""
    f = spec.field
    if f.kind == "CC":
        return ql.CC
    p, n = f.p, (f.n or 1)
    if not isinstance(p, int):
        raise RunError("FieldError", "field GF needs an integer p >= 2")
    return ql.GF(p ** n)


def _family_map() -> dict:
    """``{name: builder}`` for every buildable v1 family, reusing the catalog's
    own ``_iter_families`` so the app and the catalog agree on names by
    construction."""
    return {name: fn for name, fn in _iter_families()}


def build_algebra(spec):
    """Build a library ``Algebra`` from a validated ``AlgebraSpec`` (family or
    quiver). ``CatalogError`` / ``QuiverlabError`` propagate to the caller."""
    if spec.kind == "quiver":
        Q = ql.Quiver(vertices=list(spec.vertices),
                      arrows={k: tuple(v) for k, v in spec.arrows.items()})
        return Q.algebra(relations=list(spec.relations), field=_field(spec))
    # family
    validate_family(spec.family, spec.params)          # raises CatalogError
    builder = _family_map().get(spec.family)
    if builder is None:                                # validate_family passed
        raise RunError("CatalogError", f"no builder for family {spec.family!r}")
    return builder(field=_field(spec), **spec.params)


# --------------------------------------------------------------------------- #
# Top-level entry point
# --------------------------------------------------------------------------- #

def run_spec(req: ComputeRequest, artifact_dir,
             progress_cb: Callable[[dict], None] | None = None,
             result_max_bytes: int | None = None) -> dict:
    """Build the algebra, run each requested computation, and (once the result
    JSON is under the byte cap) write ``result.json`` plus the requested
    artifacts (the worked-steps ``trace.pdf`` -- or ``trace_steps.html`` when no
    LaTeX toolchain is present -- when ``pdf`` is set, ``tikz.tex`` when ``tikz``
    is set). Returns the result dict."""
    artifact_dir = Path(artifact_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    tikz_src = None
    hh_trace = None                    # (table, kind, top) of the LAST HH item
    A = None
    events: list = []
    try:
        items = [parse_compute_item(s) for s in req.compute]
        # Finding 3: two compute items of the same kind collide in the
        # kind-keyed ``results`` dict -- refuse loudly before computing anything.
        seen: set = set()
        for item in items:
            if item.kind in seen:
                raise RunError("DuplicateComputeItem",
                               f"{item.kind} requested twice; one range per "
                               "invariant kind")
            seen.add(item.kind)

        # Finding 4a: pin the global flag off BEFORE build_algebra (it defaults
        # True and would litter a cwd-relative quiverlab_traces/ with PDFs);
        # restore it in ``finally``.
        prev_verbose = getattr(ql, "verbose", None)
        if hasattr(ql, "verbose"):
            ql.verbose = False
        try:
            A = build_algebra(req.algebra)
            # Build the module(s) once if any module compute kind is requested. A
            # matrix that violates the relations raises the library's loud
            # QuiverlabError right here -- surfaced verbatim as a clean 4xx by the
            # caller (a QuiverlabError subclass name is a safe error type), never a
            # 500 and never a silent wrong answer.
            M = (_build_module(A, req.module, "M")
                 if any(it.kind in MODULE_KINDS for it in items) else None)
            N = (_build_module(A, req.ext_target, "N")
                 if any(it.kind == "ext" for it in items) else None)
            results: dict = {}
            for i, item in enumerate(items):
                if progress_cb:
                    progress_cb({"step": i, "of": len(items), "kind": item.kind})
                if item.kind in MODULE_KINDS:
                    results[item.kind] = _dispatch_module(A, item, M, N)
                    continue
                block, hh = _dispatch(A, item, events)
                results[item.kind] = block
                if hh is not None:
                    hh_trace = hh      # keep the last HH (table, kind, top)
        finally:
            if hasattr(ql, "verbose") and prev_verbose is not None:
                ql.verbose = prev_verbose

        if req.artifacts.tikz and hasattr(A, "tikz"):
            tikz_src = A.tikz()

        used_keys: list = []
        for payload in results.values():
            for k in payload.get("references", []):
                if k not in used_keys:
                    used_keys.append(k)

        meta: dict = {}
        if req.artifacts.pdf:
            # Provisional value for the cap check only: the no-HH case is final;
            # the HH case gets the (longest) fallback string here and its true
            # value after write_trace, below. result.json is written last, so
            # the cap-checked payload is >= what lands on disk.
            meta["pdf"] = _PDF_NO_HH if hh_trace is None else _PDF_HTML_FALLBACK

        result = {
            "quiverlab_version": getattr(ql, "__version__", "unknown"),
            "algebra": req.algebra.model_dump(),
            "results": results,
            "references": resolve_references(used_keys),
            "reproduce": _snippet(req, A),
            "meta": meta,
        }
    except CatalogError as exc:
        raise RunError("CatalogError", str(exc))
    except qerr.QuiverlabError as exc:
        raise RunError(type(exc).__name__, str(exc))
    except RunError:
        raise
    except Exception as exc:  # unexpected: full traceback server-side only
        _log.exception("run_spec failed for algebra=%s",
                       getattr(req.algebra, "family", req.algebra.kind))
        raise RunError(type(exc).__name__, str(exc))

    payload = json.dumps(result, indent=2, default=str)
    if result_max_bytes is not None and len(payload.encode("utf-8")) > result_max_bytes:
        raise RunError("ResultTooLarge",
                       f"result exceeds the {result_max_bytes}-byte cap; narrow "
                       "the degree range or run locally: pip install quiverlab")
    # --- write phase: nothing hits disk until the cap check has passed ---
    if req.artifacts.pdf and hh_trace is not None:
        table, kind, top = hh_trace
        meta["pdf"] = _write_worked_steps(events, table, A, kind, top,
                                          used_keys, artifact_dir)
        payload = json.dumps(result, indent=2, default=str)   # meta now honest
    (artifact_dir / "result.json").write_text(payload)
    if tikz_src is not None:
        (artifact_dir / "tikz.tex").write_text(tikz_src)
    return result


def _write_worked_steps(events, table, A, kind, top, used_keys, artifact_dir) -> str:
    """Render the LAST HH computation's worked steps into the artifact dir via
    the library's ``write_trace`` (PDF when pdflatex/tectonic is on PATH, else a
    self-contained no-JS HTML). Normalize the produced file to ``trace.pdf`` /
    ``trace_steps.html`` and return the honest ``meta["pdf"]`` label. ``out_dir``
    is explicit, so write_trace never reads or writes the cwd."""
    from quiverlab.trace.writer import write_trace
    produced = Path(write_trace(list(events), table, algebra=A, kind=kind, top=top,
                                references=_trace_references(used_keys, events),
                                out_dir=str(artifact_dir)))
    if produced.suffix == ".pdf":
        target, label = artifact_dir / "trace.pdf", _PDF_OK
    else:
        target, label = artifact_dir / "trace_steps.html", _PDF_HTML_FALLBACK
    if produced.resolve() != target.resolve():
        produced.replace(target)       # same dir -> atomic rename
    return label


def _trace_references(used_keys, events):
    """``(bibtex_key, formatted)`` pairs for the worked-steps bibliography,
    resolved through the library's provenance registry -- the same source a local
    ``verbose=True`` run uses. Prefers the run's aggregated keys; on a registry
    gap (KeyError) falls back to the trace-implied engine keys, then to ``()``."""
    from quiverlab.trace.provenance import references_for
    from quiverlab.trace.provenance import resolve_references as _rr
    try:
        return _rr(tuple(used_keys))
    except KeyError:
        try:
            return _rr(references_for(events))
        except KeyError:
            return ()


# --------------------------------------------------------------------------- #
# Per-invariant dispatch (block shapes mirror docs/gui/runner.py::compute_one)
# --------------------------------------------------------------------------- #

def _dispatch(A, item, events) -> tuple:
    """Returns ``(block, hh_capture)`` where ``hh_capture`` is
    ``(table, table.kind, top)`` for an HH item (so ``run_spec`` can render its
    worked-steps PDF from the LAST one) or ``None`` for any other invariant."""
    kind = item.kind
    if kind in ("hh_cohomology", "hh_homology"):
        top = item.hi
        if top is None:
            raise RunError("SchemaError",
                           f"{kind} needs a degree range, e.g. '{kind}:0..4'")
        method = (A.hochschild_cohomology if kind == "hh_cohomology"
                  else A.hochschild_homology)
        table = method(top, verbose=False, trace=events)
        keys = list(table.references)
        block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                 "engine": table.engine, "references": keys,
                 "citations": _citation_pairs(keys)}
        return block, (table, table.kind, top)
    if kind == "cartan":
        m = _rows(A.cartan_matrix())
        keys = list(A.citations())
        return {"matrix": m, "latex": _latex_matrix(m),
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "coxeter_polynomial":
        import sympy
        p = A.coxeter_polynomial()
        keys = list(A.citations())
        return {"latex": sympy.latex(p.as_expr()), "text": str(p.as_expr()),
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "global_dimension":
        g = A.global_dimension()
        keys = list(A.citations())
        return {"text": str(g), "exact": bool(g.exact), "value": g.value,
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "center":
        dim_z, basis = A.center()
        keys = list(A.citations())
        # Basis entries are exact ints/rationals (sympy MPQ over CC) -- not
        # JSON-serializable in general; ship them as exact strings.
        return {"dim": dim_z, "basis": [[str(x) for x in row] for row in basis],
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "dimension":
        keys = list(A.citations())
        return {"value": A.dim, "references": keys,
                "citations": _citation_pairs(keys)}, None
    raise RunError("SchemaError", f"unsupported computation {kind!r}")


def _citation_pairs(keys) -> list:
    """``[[bibtex_key, formatted], ...]`` via the library's provenance resolver
    (as ``compute_one`` does). Degrades to ``[]`` if a key is not in the registry
    -- a citation-registry gap must not abort a valid computation."""
    if not keys:
        return []
    try:
        from quiverlab.trace.provenance import resolve_references as _rr
        return [list(p) for p in _rr(tuple(keys))]
    except KeyError:                                    # unknown citation key
        # A citation-registry gap must not abort a valid computation; a real bug
        # (any non-KeyError) is left to surface loudly.
        _log.warning("citation resolution failed for keys=%s", keys, exc_info=True)
        return []


def _latex_matrix(rows) -> str:
    body = r" \\ ".join(" & ".join(str(x) for x in row) for row in rows)
    return r"\begin{pmatrix} %s \end{pmatrix}" % body


def _rows(mat) -> list:
    return [list(r) for r in (mat.tolist() if hasattr(mat, "tolist") else mat)]


# --------------------------------------------------------------------------- #
# Module block (Plan 26): build the module from the request, then dispatch the
# module-level invariants. All mathematics via the public quiverlab surface.
# --------------------------------------------------------------------------- #

def _match_vertex(algebra, label):
    """Resolve a JSON vertex label (a builtin's ``vertex``) to the actual vertex
    object of ``algebra`` (ints for quiver/family). Matched by string form so both
    ``1`` and ``"1"`` name vertex 1."""
    for v in algebra.quiver.vertices:
        if v == label or str(v) == str(label):
            return v
    raise RunError("SchemaError", f"module: no vertex {label!r} in the algebra")


def _parse_entry(x):
    """Parse a matrix entry to an EXACT python scalar (int or Fraction). Entries
    are DATA, never evaluated; a float is refused (exactness is the point). The
    library's field coercion turns these into GF(p)/CC elements."""
    if isinstance(x, bool):
        raise RunError("SchemaError", f"module entry {x!r} is not a number")
    if isinstance(x, int):
        return x
    if isinstance(x, float):
        raise RunError("ExactnessError",
                       f"module entry {x!r} is a float; entries must be exact "
                       "(integers or strings like '1/2')")
    if isinstance(x, str):
        s = x.strip()
        try:
            return int(s)
        except ValueError:
            pass
        try:
            return Fraction(s)
        except (ValueError, ZeroDivisionError):
            raise RunError("SchemaError",
                           f"module entry {x!r} is not an exact integer or fraction")
    raise RunError("SchemaError", f"module entry {x!r} is not a number")


def _full_matrices(algebra, mspec):
    """Turn the request's per-arrow BLOCK matrices (each sized dim[target] x
    dim[source] in the representation quiver) into the FULL n x n arrow actions
    the library's ``A.module`` consumes, in the vertex-ordered basis. Returns
    ``(dimension_vector, arrow_action)`` keyed by the actual vertex/arrow objects.
    ``side`` selects the representation quiver: right = A, left = A^op (a left
    A-module is a right A^op-module), so an arrow's source/target -- hence the
    block orientation -- come from that quiver, correct for both sides."""
    rep = algebra if mspec.side == "right" else algebra.opposite()
    verts = list(rep.quiver.vertices)
    by_str = {str(v): v for v in verts}
    dimvec = {v: 0 for v in verts}
    for key, n in (mspec.dims or {}).items():
        if key not in by_str:
            raise RunError("SchemaError", f"module: no vertex {key!r} in the algebra")
        dimvec[by_str[key]] = int(n)
    start, off = {}, 0
    for v in verts:
        start[v] = off
        off += dimvec[v]
    n = off
    arrow_names = list(rep.quiver.arrows)
    for a in (mspec.maps or {}):
        if a not in arrow_names:
            raise RunError("SchemaError", f"module: no arrow {a!r} in the algebra")
    action = {}
    for a in arrow_names:
        s, t = rep.quiver.source(a), rep.quiver.target(a)
        rows, cols = dimvec[t], dimvec[s]
        full = [[0] * n for _ in range(n)]
        block = (mspec.maps or {}).get(a)
        if block is not None:
            if len(block) != rows or any(len(r) != cols for r in block):
                got = f"{len(block)}x{len(block[0]) if block else 0}"
                raise RunError("SchemaError",
                               f"module map {a!r} must be {rows}x{cols} (target x "
                               f"source dims); got {got}")
            for i in range(rows):
                for j in range(cols):
                    full[start[t] + i][start[s] + j] = _parse_entry(block[i][j])
        action[a] = full
    return dimvec, action


def _build_module(algebra, mspec, name):
    if mspec is None:                          # defensive: the schema guarantees it
        raise RunError("SchemaError", f"module {name} block is required")
    if mspec.builtin is not None:
        b = mspec.builtin
        v = _match_vertex(algebra, b.vertex)
        builder = {"simple": algebra.simple, "projective": algebra.projective,
                   "injective": algebra.injective}[b.kind]
        return builder(v, side=mspec.side)
    dimvec, action = _full_matrices(algebra, mspec)
    return algebra.module(dimvec, action, side=mspec.side, name=name)


def _dv(dimvec) -> dict:
    """A JSON-safe dimension vector: string vertex keys (sorted), int values."""
    return {str(v): int(n) for v, n in sorted(dimvec.items(), key=lambda kv: str(kv[0]))}


def _mod_view(m) -> dict:
    return {"dimvec": _dv(m.dimension_vector()), "dim": m.dim}


# Citation step-ids (registry keys): GSZ2001 for the minimal engine / module Ext
# ("minimal_resolution" / "module_ext"), ASS2006 for the representation-theory
# constructions ("assem_book").
_MOD_REFS = {
    "dimension_vector": ["assem_book"],
    "rad_top_soc": ["assem_book"],
    "tau": ["assem_book"],
    "tau_minus": ["assem_book"],
    "ext": ["module_ext"],
    "projective_resolution": ["minimal_resolution"],
    "projective_dimension": ["minimal_resolution"],
    "injective_resolution": ["minimal_resolution", "assem_book"],
    "injective_dimension": ["minimal_resolution", "assem_book"],
}


def _with_refs(block: dict, kind: str) -> dict:
    keys = _MOD_REFS[kind]
    block["references"] = list(keys)
    block["citations"] = _citation_pairs(keys)
    return block


def _dispatch_module(A, item, M, N) -> dict:
    """Per-module-invariant dispatch. Every block carries ``references`` (registry
    step-ids) + resolved ``citations`` exactly like the algebra invariants."""
    kind = item.kind
    if kind == "dimension_vector":
        return _with_refs({"kind": "dimension_vector", "side": M.side,
                           **_mod_view(M)}, kind)
    if kind == "rad_top_soc":
        return _with_refs({"kind": "rad_top_soc", "side": M.side,
                           "radical": _mod_view(M.radical()),
                           "top": _mod_view(M.top()),
                           "socle": _mod_view(M.socle())}, kind)
    if kind in ("tau", "tau_minus"):
        t = M.tau() if kind == "tau" else M.tau_minus()
        return _with_refs({"kind": kind, "side": t.side, "is_zero": t.dim == 0,
                           **_mod_view(t)}, kind)
    if kind == "ext":
        top = item.hi
        if top is None:
            raise RunError("SchemaError", "ext needs a degree range, e.g. 'ext:0..4'")
        from quiverlab.modules.ext import ext_dims
        dims = ext_dims(A, M, N, top)          # loud if M, N are not comparable
        return _with_refs({"kind": "ext", "top": top,
                           "dims": [int(d) for d in dims],
                           "target": _mod_view(N)}, kind)
    if kind in ("projective_resolution", "injective_resolution"):
        top = item.hi
        if top is None:
            raise RunError("SchemaError", f"{kind} needs a degree range, e.g. "
                           f"'{kind}:0..4'")
        res = (M.projective_resolution(top) if kind == "projective_resolution"
               else M.injective_resolution(top))
        terms = [_dv(dv) for dv in res.dimension_vectors()]
        block = {"kind": kind, "top": top, "terms": terms,
                 "betti": [res.betti(i) for i in range(len(terms))]}
        if kind == "projective_resolution":
            block["pd"] = res.pd()
        else:
            block["injective_dimension"] = res.injective_dimension()
        return _with_refs(block, kind)
    if kind == "projective_dimension":
        pd = M.projective_resolution(_PD_BOUND).pd()
        return _with_refs({"kind": "projective_dimension", "value": pd,
                           "finite": pd is not None, "bound": _PD_BOUND}, kind)
    if kind == "injective_dimension":
        idim = M.injective_dimension(bound=_PD_BOUND)
        return _with_refs({"kind": "injective_dimension", "value": idim,
                           "finite": idim is not None, "bound": _PD_BOUND}, kind)
    raise RunError("SchemaError", f"unsupported module computation {kind!r}")


# --------------------------------------------------------------------------- #
# Copy-paste reproduction snippet
# --------------------------------------------------------------------------- #

def _fmt_scalar(x) -> str:
    if isinstance(x, Fraction):
        return (str(x.numerator) if x.denominator == 1
                else f"Fraction({x.numerator}, {x.denominator})")
    return str(x)


def _pymat(mat) -> str:
    return ("[" + ", ".join("[" + ", ".join(_fmt_scalar(x) for x in row) + "]"
                            for row in mat) + "]")


def _module_construction(mspec, varname, A) -> list:
    """Runnable lines building ``varname`` locally. Builtins reproduce the
    ``A.simple/projective/injective`` call; explicit modules reproduce the FULL
    arrow-action matrices ``A.module`` consumes (exact int/Fraction entries)."""
    if mspec.builtin is not None:
        b = mspec.builtin
        return [f'{varname} = A.{b.kind}({b.vertex!r}, side="{mspec.side}")']
    dimvec, action = _full_matrices(A, mspec)
    lines = []
    if any(isinstance(x, Fraction) for m in action.values() for row in m for x in row):
        lines.append("from fractions import Fraction")
    dv_lit = "{" + ", ".join(f"{v!r}: {n}" for v, n in dimvec.items()) + "}"
    maps_lit = "{" + ", ".join(f'"{a}": {_pymat(m)}' for a, m in action.items()) + "}"
    lines.append(f'{varname} = A.module({dv_lit}, {maps_lit}, '
                 f'side="{mspec.side}", name="{varname}")')
    return lines


def _snippet(req: ComputeRequest, A) -> str:
    f = req.algebra.field
    field_name = "CC" if f.kind == "CC" else "GF"
    field_expr = "CC" if f.kind == "CC" else f"GF({f.p ** (f.n or 1)})"
    if req.algebra.kind == "quiver":
        arrows = ", ".join('"%s": (%d, %d)' % (k, s, t)
                           for k, (s, t) in req.algebra.arrows.items())
        lines = [f"from quiverlab import Quiver, {field_name}", "",
                 f"Q = Quiver(vertices={list(req.algebra.vertices)!r}, "
                 f"arrows={{{arrows}}})",
                 f"A = Q.algebra(relations={list(req.algebra.relations)!r}, "
                 f"field={field_expr})"]
    else:
        params = "".join(f", {k}={v!r}" for k, v in req.algebra.params.items())
        lines = ["import quiverlab as ql",
                 f"A = ql.{req.algebra.family}(field=ql.{field_expr}{params})"]
    if req.module is not None:
        lines += _module_construction(req.module, "M", A)
    if req.ext_target is not None:
        lines += _module_construction(req.ext_target, "N", A)
    _snip = {"hh_cohomology": lambda it: f"A.hochschild_cohomology({it.hi})",
             "hh_homology": lambda it: f"A.hochschild_homology({it.hi})",
             "coxeter_polynomial": lambda it: "A.coxeter_polynomial()",
             "cartan": lambda it: "A.cartan_matrix()",
             "global_dimension": lambda it: "A.global_dimension()",
             "center": lambda it: "A.center()",
             "dimension": lambda it: "A.dim",
             "dimension_vector": lambda it: "M.dimension_vector()",
             "rad_top_soc": lambda it: "M.radical(), M.top(), M.socle()",
             "tau": lambda it: "M.tau()",
             "tau_minus": lambda it: "M.tau_minus()",
             "ext": lambda it: f"[A.ext(M, N, i) for i in range({it.hi} + 1)]",
             "projective_resolution":
                 lambda it: f"M.projective_resolution({it.hi}).dimension_vectors()",
             "injective_resolution":
                 lambda it: f"M.injective_resolution({it.hi}).dimension_vectors()",
             "projective_dimension":
                 lambda it: f"M.projective_resolution({_PD_BOUND}).pd()",
             "injective_dimension": lambda it: "M.injective_dimension()"}
    for s in req.compute:
        item = parse_compute_item(s)
        if item.kind in _snip:
            lines.append(_snip[item.kind](item))
    return "\n".join(lines)
