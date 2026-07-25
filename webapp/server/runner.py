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
from pathlib import Path
from typing import Callable

import quiverlab as ql
from quiverlab import errors as qerr

from webapp.server.catalog import validate_family, CatalogError, _iter_families
from webapp.server.references import resolve_references
from webapp.server.schema import ComputeRequest, parse_compute_item

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
            results: dict = {}
            for i, item in enumerate(items):
                if progress_cb:
                    progress_cb({"step": i, "of": len(items), "kind": item.kind})
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
            "reproduce": _snippet(req),
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
# Copy-paste reproduction snippet
# --------------------------------------------------------------------------- #

def _snippet(req: ComputeRequest) -> str:
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
    _snip = {"hh_cohomology": lambda it: f"A.hochschild_cohomology({it.hi})",
             "hh_homology": lambda it: f"A.hochschild_homology({it.hi})",
             "coxeter_polynomial": lambda it: "A.coxeter_polynomial()",
             "cartan": lambda it: "A.cartan_matrix()",
             "global_dimension": lambda it: "A.global_dimension()",
             "center": lambda it: "A.center()",
             "dimension": lambda it: "A.dim"}
    for s in req.compute:
        item = parse_compute_item(s)
        if item.kind in _snip:
            lines.append(_snip[item.kind](item))
    return "\n".join(lines)
