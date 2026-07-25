"""Plan 28 -- the HPC spec core: stdlib request parse/validate + the compute
dispatch, promoted out of ``webapp/server/runner.py`` so the wheel CLI and the
webapp share ONE implementation. No pydantic (the base wheel stays lean); no
``webapp`` import (the import-boundary test pins that). Everything mathematical
goes through the PUBLIC ``import quiverlab`` surface.

``webapp/server/runner.py`` now DELEGATES here: it validates at the HTTP boundary
with pydantic, then calls :func:`run` with ``req.model_dump(by_alias=True)`` and
returns the byte-identical result dict (pinned by
``tests/webapp/test_runner_delegation.py``). The request shape mirrors
``webapp/server/schema.py`` (schema 1/2; family|quiver algebra; ``compute`` items;
``artifacts``; ``module``/``ext_target`` blocks with exact-entry matrices) plus an
``hpc:`` block (``checkpoint_dir``, ``time_limit_s``, ``max_mem_bytes``,
``allow_large``, ``prime``, ``max_cells``, ``engine``) used only by the CLI.

The dispatch is ported near-verbatim from the runner (block shapes mirror
``docs/gui/runner.py::compute_one``); the only refusal class is
:class:`ComputeError` (``error_type`` + ``message``), which the webapp wrapper
re-raises as its ``RunError`` and the CLI maps to exit codes."""
from __future__ import annotations

import inspect
import json
import logging
import re
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

import quiverlab as ql
from quiverlab import errors as qerr

# The result envelope's schema version (Plan 28). ``run`` stamps it only when
# asked (the CLI does; the webapp does NOT, so its result dict is byte-unchanged).
RESULT_SCHEMA = 1

# Depth to which projective dimension is probed before reporting "infinite".
_PD_BOUND = 32
_DEEPEN_PRIME = 32003

_log = logging.getLogger("quiverlab.hpc")

SIDES = ("right", "left")

MODULE_KINDS = frozenset({
    "dimension_vector", "rad_top_soc", "ext", "tau", "tau_minus",
    "projective_resolution", "injective_resolution",
    "projective_dimension", "injective_dimension",
})
MODULE_RANGE_KINDS = frozenset({"ext", "projective_resolution", "injective_resolution"})

# Honest labels for ``meta["pdf"]`` (mirrors the runner verbatim).
_PDF_OK = "trace.pdf"
_PDF_HTML_FALLBACK = ("PDF toolchain (pdflatex/tectonic) not found -- "
                      "worked steps in trace_steps.html")
_PDF_NO_HH = "no traced computation requested (PDF covers HH worked steps)"

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


# --------------------------------------------------------------------------- #
# Errors
# --------------------------------------------------------------------------- #

class SpecError(qerr.QuiverlabError):
    """A malformed HPC/compute request (schema violation). A QuiverlabError so it
    is a client-safe error type (surfaced verbatim, never a 500)."""


class ComputeError(Exception):
    """A refusal the caller turns into an honest error. ``error_type`` is the
    library exception class name (verbatim) or one of the runner's own tags
    (``SchemaError``, ``CatalogError``, ``ResultTooLarge``, ...). The webapp maps
    this to its ``RunError``; the CLI maps ``error_type`` to an exit code."""

    def __init__(self, error_type: str, message: str):
        super().__init__(f"{error_type}: {message}")
        self.error_type = error_type
        self.message = message


class CheckpointStop(Exception):
    """The deepen path stopped cleanly at a checkpoint (time/memory budget) before
    completing the requested degree. The CLI turns this into exit 75 (resumable);
    the checkpoint on disk lets a rerun resume."""

    def __init__(self, degree: int, reason: str, reached: int):
        super().__init__(f"checkpoint stop at HH_{reached} (reason={reason}); "
                         f"resume to reach HH_{degree}")
        self.degree = degree
        self.reason = reason
        self.reached = reached


# --------------------------------------------------------------------------- #
# Request dataclasses (duck-type the pydantic models the runner read)
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class FieldSpec:
    kind: str
    p: int | None = None
    n: int = 1


@dataclass(frozen=True)
class FamilyAlgebra:
    kind: str
    family: str
    params: dict
    field: FieldSpec


@dataclass(frozen=True)
class QuiverAlgebra:
    kind: str
    vertices: list
    arrows: dict          # name -> (source, target) tuple
    relations: list
    field: FieldSpec


@dataclass(frozen=True)
class Artifacts:
    pdf: bool = False
    tikz: bool = False


@dataclass(frozen=True)
class BuiltinModule:
    kind: str
    vertex: Any


@dataclass(frozen=True)
class ModuleSpec:
    dims: dict | None = None
    maps: dict | None = None
    builtin: BuiltinModule | None = None
    side: str = "right"


@dataclass(frozen=True)
class HpcConfig:
    checkpoint_dir: str | None = None
    time_limit_s: int | None = None
    max_mem_bytes: int | None = None
    allow_large: bool = False
    prime: int | None = None      # None = use the algebra's own prime (deepen)
    max_cells: int | None = None
    engine: str | None = None


@dataclass(frozen=True)
class ComputeItem:
    kind: str
    lo: int | None
    hi: int | None


@dataclass(frozen=True)
class ComputeRequest:
    schema_version: int
    algebra: Any                 # FamilyAlgebra | QuiverAlgebra
    compute: list
    artifacts: Artifacts
    module: ModuleSpec | None
    ext_target: ModuleSpec | None
    hpc: HpcConfig | None
    raw_algebra: dict            # verbatim echo for the result envelope


# --------------------------------------------------------------------------- #
# Compute-item grammar (ported from webapp/server/schema.py)
# --------------------------------------------------------------------------- #

_RANGE = re.compile(r"^(?P<kind>[a-z_]+)(?::(?P<lo>\d+)\.\.(?P<hi>\d+))?$")


def parse_compute_item(s: str) -> ComputeItem:
    if not isinstance(s, str):
        raise SpecError(f"compute item {s!r} must be a string")
    m = _RANGE.match(s)
    if not m:
        raise SpecError(f"unparseable compute item {s!r}")
    lo = int(m["lo"]) if m["lo"] is not None else None
    hi = int(m["hi"]) if m["hi"] is not None else None
    if lo is not None and hi is not None and hi < lo:
        raise SpecError(f"empty degree range in {s!r}")
    return ComputeItem(kind=m["kind"], lo=lo, hi=hi)


# --------------------------------------------------------------------------- #
# Bibliography adapter (ported from webapp/server/references.py; library-sourced)
# --------------------------------------------------------------------------- #

def _bib_get(entry, name, default=None):
    if isinstance(entry, dict):
        return entry.get(name, default)
    return getattr(entry, name, default)


def _entry_view(entry) -> dict:
    key = _bib_get(entry, "key") or str(entry)
    formatted = _bib_get(entry, "formatted") or _bib_get(entry, "text") or key
    doi = _bib_get(entry, "doi")
    arxiv = _bib_get(entry, "arxiv") or _bib_get(entry, "arxiv_id")
    return {
        "key": key,
        "bibtex_key": _bib_get(entry, "bibtex_key"),
        "formatted": formatted,
        "doi": doi,
        "arxiv": arxiv,
        "doi_url": (f"https://doi.org/{doi}" if doi else None),
        "arxiv_url": (f"https://arxiv.org/abs/{arxiv}" if arxiv else None),
        "topic": _bib_get(entry, "topic") or "other",
        "annotation": _bib_get(entry, "annotation") or _bib_get(entry, "note") or "",
    }


@lru_cache(maxsize=1)
def _bib_index() -> dict:
    return {v["key"]: v for v in (_entry_view(e) for e in ql.bibliography())}


def resolve_references(keys) -> list:
    """Map used BibTeX keys to library entries, preserving order (byte-identical
    to the webapp's ``references.resolve_references``)."""
    idx = _bib_index()
    out = []
    for k in keys:
        out.append(idx.get(k, {"key": k, "bibtex_key": None, "formatted": k,
                               "doi": None, "arxiv": None, "doi_url": None,
                               "arxiv_url": None, "topic": "other",
                               "annotation": ""}))
    return out


# --------------------------------------------------------------------------- #
# Family catalog (ported from webapp/server/catalog.py; library-introspected)
# --------------------------------------------------------------------------- #

def _iter_families():
    for info in ql.families():
        name = info.name
        if name == "zoo":
            continue
        builder = getattr(ql, name, None)
        if builder is None:
            raise ComputeError("CatalogError",
                               f"families() lists {name!r} but quiverlab has no such "
                               "export; the family catalog has drifted")
        yield name, builder


@lru_cache(maxsize=1)
def _family_map() -> dict:
    return {name: fn for name, fn in _iter_families()}


def _param_names(builder) -> set:
    try:
        sig = inspect.signature(builder)
    except (TypeError, ValueError):
        return set()
    return {p for p in sig.parameters if p != "field"}


def validate_family(name: str, params: dict) -> None:
    fam = _family_map()
    if name not in fam:
        raise ComputeError("CatalogError", f"unknown family {name!r}")
    known = _param_names(fam[name])
    unknown = set(params) - known
    if unknown:
        raise ComputeError("CatalogError",
                           f"family {name!r} got unknown params {sorted(unknown)}")


# --------------------------------------------------------------------------- #
# Validation: dict -> ComputeRequest (stdlib only; mirrors schema.py semantics)
# --------------------------------------------------------------------------- #

def _valid_entry(x) -> bool:
    """A matrix entry is exact DATA: a JSON integer or an exact string literal;
    never a float, never a bool (mirrors schema.py::_valid_entry)."""
    if isinstance(x, bool):
        return False
    if isinstance(x, int):
        return True
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return False
        try:
            int(s)
            return True
        except ValueError:
            pass
        try:
            Fraction(s)
            return True
        except (ValueError, ZeroDivisionError):
            return False
    return False


def _parse_field(f) -> FieldSpec:
    if not isinstance(f, dict):
        raise SpecError("algebra.field must be an object with a 'kind'")
    kind = f.get("kind")
    if kind not in ("CC", "GF"):
        raise SpecError(f"unknown field kind {kind!r} (expected 'CC' or 'GF')")
    p = f.get("p")
    n = f.get("n", 1)
    if p is not None:
        if not isinstance(p, int) or isinstance(p, bool):
            raise SpecError("field GF needs an integer p")
        if p < 2:
            raise SpecError("GF(p): p must be a prime >= 2")
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise SpecError("field GF needs an integer n >= 1")
    return FieldSpec(kind=kind, p=p, n=n)


def _parse_algebra(a):
    if not isinstance(a, dict):
        raise SpecError("algebra must be an object")
    kind = a.get("kind")
    fld = _parse_field(a.get("field"))
    if kind == "family":
        family = a.get("family")
        if not isinstance(family, str):
            raise SpecError("algebra.family must be a string")
        params = a.get("params") or {}
        if not isinstance(params, dict):
            raise SpecError("algebra.params must be an object")
        return FamilyAlgebra("family", family, dict(params), fld)
    if kind == "quiver":
        vertices = a.get("vertices")
        if not (isinstance(vertices, list) and vertices
                and all(isinstance(v, int) and not isinstance(v, bool) for v in vertices)):
            raise SpecError("algebra.vertices must be a non-empty list of integers")
        arrows_in = a.get("arrows") or {}
        if not isinstance(arrows_in, dict):
            raise SpecError("algebra.arrows must map names to [source, target] pairs")
        arrows = {}
        for name, st in arrows_in.items():
            st = tuple(st) if isinstance(st, (list, tuple)) else st
            if (not isinstance(st, tuple) or len(st) != 2
                    or not all(isinstance(x, int) and not isinstance(x, bool) for x in st)):
                raise SpecError("algebra.arrows must map names to [source, target] pairs")
            arrows[name] = st
        relations = a.get("relations", [])
        if not (isinstance(relations, list) and all(isinstance(r, str) for r in relations)):
            raise SpecError("algebra.relations must be a list of strings")
        return QuiverAlgebra("quiver", list(vertices), arrows, list(relations), fld)
    raise SpecError(f"unknown algebra kind {kind!r} (expected 'family' or 'quiver')")


def _lift_builtin_side(data: dict) -> dict:
    """Lift a ``side`` nested inside ``builtin`` up to the top level (mirrors
    schema.py::_lift_builtin_side), so the canonical form carries ``side`` once."""
    if isinstance(data.get("builtin"), dict):
        b = dict(data["builtin"])
        if "side" in b:
            inner = b.pop("side")
            outer = data.get("side")
            if outer is not None and outer != inner:
                raise SpecError(f"conflicting side: builtin.side={inner!r} vs side={outer!r}")
            data = {**data, "builtin": b, "side": inner}
    return data


def _parse_module(data, what: str) -> ModuleSpec:
    if not isinstance(data, dict):
        raise SpecError(f"{what} must be an object")
    data = _lift_builtin_side(data)
    side = data.get("side", "right")
    if side not in SIDES:
        raise SpecError(f"{what}.side must be 'right' or 'left'")
    builtin = data.get("builtin")
    dims = data.get("dims")
    maps = data.get("maps")
    if builtin is not None:
        if dims is not None or maps is not None:
            raise SpecError(f"{what}: give either a builtin pick-list OR dims+maps, not both")
        if not isinstance(builtin, dict):
            raise SpecError(f"{what}.builtin must be an object")
        bkind = builtin.get("kind")
        if bkind not in ("simple", "projective", "injective"):
            raise SpecError(f"{what}.builtin.kind must be simple/projective/injective")
        return ModuleSpec(builtin=BuiltinModule(bkind, builtin.get("vertex")), side=side)
    if dims is None:
        raise SpecError(f"{what}: needs 'dims' (a dimension vector) or a 'builtin' pick-list")
    if not isinstance(dims, dict):
        raise SpecError(f"{what}.dims must be an object mapping vertex -> dimension")
    for v, n in dims.items():
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            raise SpecError(f"{what} dims[{v!r}] must be a non-negative integer")
    if maps is not None:
        if not isinstance(maps, dict):
            raise SpecError(f"{what}.maps must be an object mapping arrow -> matrix")
        for arrow, mat in maps.items():
            if not isinstance(mat, list):
                raise SpecError(f"{what} maps[{arrow!r}] must be a matrix (list of rows)")
            width = None
            for row in mat:
                if not isinstance(row, list):
                    raise SpecError(f"{what} maps[{arrow!r}] must be a matrix (list of rows)")
                if width is None:
                    width = len(row)
                elif len(row) != width:
                    raise SpecError(f"{what} maps[{arrow!r}] is not rectangular")
                for x in row:
                    if not _valid_entry(x):
                        raise SpecError(
                            f"{what} maps[{arrow!r}] has a non-exact entry {x!r}; entries "
                            "must be integers or exact strings like '1/2' (never floats)")
    return ModuleSpec(dims=dict(dims), maps=(dict(maps) if maps is not None else None),
                      side=side)


def _parse_hpc(data) -> HpcConfig:
    if not isinstance(data, dict):
        raise SpecError("hpc block must be an object")

    def _opt_int(key):
        v = data.get(key)
        if v is None:
            return None
        if not isinstance(v, int) or isinstance(v, bool):
            raise SpecError(f"hpc.{key} must be an integer")
        return v

    ckpt = data.get("checkpoint_dir")
    if ckpt is not None and not isinstance(ckpt, str):
        raise SpecError("hpc.checkpoint_dir must be a string path")
    prime = data.get("prime")
    if prime is not None and (not isinstance(prime, int) or isinstance(prime, bool)
                              or prime < 2):
        raise SpecError("hpc.prime must be a prime integer >= 2")
    engine = data.get("engine")
    if engine is not None and engine not in ("auto", "bar", "fast", "cs"):
        raise SpecError("hpc.engine must be one of auto/bar/fast/cs")
    return HpcConfig(
        checkpoint_dir=ckpt,
        time_limit_s=_opt_int("time_limit_s"),
        max_mem_bytes=_opt_int("max_mem_bytes"),
        allow_large=bool(data.get("allow_large", False)),
        prime=prime,
        max_cells=_opt_int("max_cells"),
        engine=engine,
    )


def load_config_text(text: str) -> dict:
    """Parse a config string as JSON (stdlib, always available) or, failing that,
    YAML (``pyyaml``, imported LAZILY so the base wheel works for JSON configs)."""
    text = text.strip()
    if not text:
        raise SpecError("empty config")
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml
        except ImportError:
            raise SpecError("config is not JSON and PyYAML is not installed; "
                            "install quiverlab[hpc] for YAML config support")
        try:
            data = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            raise SpecError(f"config is not valid YAML/JSON: {exc}")
    if not isinstance(data, dict):
        raise SpecError("config must be a mapping (object) at the top level")
    return data


def load_config(path) -> dict:
    """Read a config file (``-`` = stdin) and parse it (JSON or YAML)."""
    import sys
    if str(path) == "-":
        return load_config_text(sys.stdin.read())
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError as exc:
        raise SpecError(f"cannot read config {path}: {exc}")
    return load_config_text(text)


def parse_request(data) -> ComputeRequest:
    """Validate a request dict into a :class:`ComputeRequest`. Loud
    :class:`SpecError` on any violation. Accepts the webapp's pydantic
    ``model_dump(by_alias=True)`` form AND a hand-written YAML/JSON config."""
    if isinstance(data, ComputeRequest):
        return data
    if not isinstance(data, dict):
        raise SpecError("request must be a mapping (object)")
    schema_version = data.get("schema", 1)
    if schema_version not in (1, 2):
        raise SpecError(f"unsupported schema version {schema_version}; this tool speaks v1/v2")
    if "algebra" not in data:
        raise SpecError("request needs an 'algebra' block")
    algebra = _parse_algebra(data["algebra"])
    compute = data.get("compute")
    if not (isinstance(compute, list) and compute):
        raise SpecError("compute list must not be empty")
    items = [parse_compute_item(s) for s in compute]
    artifacts_in = data.get("artifacts") or {}
    if not isinstance(artifacts_in, dict):
        raise SpecError("artifacts must be an object")
    artifacts = Artifacts(bool(artifacts_in.get("pdf", False)),
                          bool(artifacts_in.get("tikz", False)))
    module = _parse_module(data["module"], "module") if data.get("module") is not None else None
    ext_target = (_parse_module(data["ext_target"], "ext_target")
                  if data.get("ext_target") is not None else None)
    hpc = _parse_hpc(data["hpc"]) if data.get("hpc") is not None else None

    if (module is not None or ext_target is not None) and schema_version != 2:
        raise SpecError("a 'module'/'ext_target' block requires schema 2")
    kinds = {it.kind for it in items}
    if kinds & MODULE_KINDS and module is None:
        need = sorted(kinds & MODULE_KINDS)
        raise SpecError(f"module compute kind(s) {need} require a 'module' block")
    if "ext" in kinds and ext_target is None:
        raise SpecError("Ext needs a second module 'ext_target' (the N in Ext^n(M, N))")

    return ComputeRequest(schema_version=schema_version, algebra=algebra,
                          compute=list(compute), artifacts=artifacts,
                          module=module, ext_target=ext_target, hpc=hpc,
                          raw_algebra=data["algebra"])


# --------------------------------------------------------------------------- #
# Algebra + field construction
# --------------------------------------------------------------------------- #

def _field(spec: FieldSpec):
    if spec.kind == "CC":
        return ql.CC
    p, n = spec.p, (spec.n or 1)
    if not isinstance(p, int):
        raise ComputeError("FieldError", "field GF needs an integer p >= 2")
    return ql.GF(p ** n)


def build_algebra(spec):
    """Build a library ``Algebra`` from a validated algebra spec (dataclass or
    dict). ``ComputeError`` (CatalogError/FieldError) or a raw ``QuiverlabError``
    (bad relation, non-prime field) propagate to the caller."""
    if isinstance(spec, dict):
        spec = _parse_algebra(spec)
    if isinstance(spec, QuiverAlgebra):
        Q = ql.Quiver(vertices=list(spec.vertices),
                      arrows={k: tuple(v) for k, v in spec.arrows.items()})
        return Q.algebra(relations=list(spec.relations), field=_field(spec.field))
    validate_family(spec.family, spec.params)
    builder = _family_map().get(spec.family)
    if builder is None:
        raise ComputeError("CatalogError", f"no builder for family {spec.family!r}")
    return builder(field=_field(spec.field), **spec.params)


# --------------------------------------------------------------------------- #
# Top-level entry point (ported from runner.run_spec, byte-stable)
# --------------------------------------------------------------------------- #

def run(req, artifact_dir, progress_cb: Callable[[dict], None] | None = None,
        result_max_bytes: int | None = None, *, result_schema: int | None = None,
        write_result: bool = True) -> dict:
    """Build the algebra, run each requested computation, and (once the result
    JSON is under the byte cap) write ``result.json`` plus the requested artifacts
    into ``artifact_dir``. Returns the result dict.

    ``result_schema`` (CLI only) stamps the envelope with ``result_schema``; the
    webapp passes None so the returned dict is byte-identical to the pre-Plan-28
    runner. ``write_result=False`` (CLI) writes only the sidecar artifacts and
    leaves the authoritative ``result.json`` write to the caller."""
    req = parse_request(req)
    artifact_dir = Path(artifact_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    tikz_src = None
    hh_trace = None
    A = None
    events: list = []
    hh_kwargs = _hh_kwargs(req.hpc)
    try:
        items = [parse_compute_item(s) for s in req.compute]
        seen: set = set()
        for item in items:
            if item.kind in seen:
                raise ComputeError("DuplicateComputeItem",
                                   f"{item.kind} requested twice; one range per "
                                   "invariant kind")
            seen.add(item.kind)

        prev_verbose = getattr(ql, "verbose", None)
        if hasattr(ql, "verbose"):
            ql.verbose = False
        try:
            A = build_algebra(req.algebra)
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
                if _deepen_applies(req.hpc, item, req.algebra):
                    results[item.kind] = _dispatch_deepen(A, item, req.hpc, progress_cb)
                    continue
                block, hh = _dispatch(A, item, events, hh_kwargs)
                results[item.kind] = block
                if hh is not None:
                    hh_trace = hh
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
            meta["pdf"] = _PDF_NO_HH if hh_trace is None else _PDF_HTML_FALLBACK

        result = {
            "quiverlab_version": getattr(ql, "__version__", "unknown"),
            "algebra": req.raw_algebra,
            "results": results,
            "references": resolve_references(used_keys),
            "reproduce": _snippet(req, A),
            "meta": meta,
        }
        if result_schema is not None:
            result["result_schema"] = result_schema
    except CheckpointStop:
        raise
    except ComputeError:
        raise
    except qerr.QuiverlabError as exc:
        raise ComputeError(type(exc).__name__, str(exc))
    except Exception as exc:  # unexpected: full traceback in the log only
        _log.exception("hpc.run failed for algebra=%s",
                       getattr(req.algebra, "family", req.algebra.kind))
        raise ComputeError(type(exc).__name__, str(exc))

    payload = json.dumps(result, indent=2, default=str)
    if result_max_bytes is not None and len(payload.encode("utf-8")) > result_max_bytes:
        raise ComputeError("ResultTooLarge",
                           f"result exceeds the {result_max_bytes}-byte cap; narrow "
                           "the degree range or run locally: pip install quiverlab")
    if req.artifacts.pdf and hh_trace is not None:
        table, kind, top = hh_trace
        meta["pdf"] = _write_worked_steps(events, table, A, kind, top,
                                          used_keys, artifact_dir)
        payload = json.dumps(result, indent=2, default=str)
    if write_result:
        (artifact_dir / "result.json").write_text(payload)
    if tikz_src is not None:
        (artifact_dir / "tikz.tex").write_text(tikz_src)
    return result


def _hh_kwargs(hpc: HpcConfig | None) -> dict:
    """Optional ``max_cells``/``engine`` threaded into ``hochschild_*`` from the
    hpc block. Empty when no hpc block (the webapp), so the HH call -- and hence
    the whole result -- is byte-identical to the pre-Plan-28 runner."""
    out: dict = {}
    if hpc is not None:
        if hpc.max_cells is not None:
            out["max_cells"] = hpc.max_cells
        if hpc.engine is not None:
            out["engine"] = hpc.engine
    return out


def _write_worked_steps(events, table, A, kind, top, used_keys, artifact_dir) -> str:
    from quiverlab.trace.writer import write_trace
    produced = Path(write_trace(list(events), table, algebra=A, kind=kind, top=top,
                                references=_trace_references(used_keys, events),
                                out_dir=str(artifact_dir)))
    if produced.suffix == ".pdf":
        target, label = artifact_dir / "trace.pdf", _PDF_OK
    else:
        target, label = artifact_dir / "trace_steps.html", _PDF_HTML_FALLBACK
    if produced.resolve() != target.resolve():
        produced.replace(target)
    return label


def _trace_references(used_keys, events):
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
# Deepen routing (CLI-only; sanctioned via deepen_runner)
# --------------------------------------------------------------------------- #

def _deepen_applies(hpc: HpcConfig | None, item: ComputeItem, algebra) -> bool:
    """The checkpointed minimal-A^e ``deepen`` path serves a big
    ``hh_homology`` over GF(p) when a checkpoint dir is configured; everything
    else takes the ordinary in-memory path."""
    return (hpc is not None and hpc.checkpoint_dir is not None
            and item.kind == "hh_homology" and item.hi is not None
            and algebra.field.kind == "GF")


def _dispatch_deepen(A, item: ComputeItem, hpc: HpcConfig, progress_cb) -> dict:
    from quiverlab.hpc.deepen_runner import run_deepen
    log = None
    if progress_cb is not None:
        log = lambda msg: progress_cb({"deepen": msg})   # noqa: E731
    dims, complete, reason = run_deepen(A, item.hi, hpc, log=log)
    if not complete:
        raise CheckpointStop(item.hi, reason, len(dims) - 1)
    keys = list(A.citations())
    return {"kind": "HH_", "top": item.hi, "dims": [int(d) for d in dims],
            "engine": "hanlab minimal A^e resolution (deepen/checkpointed)",
            "references": keys, "citations": _citation_pairs(keys)}


# --------------------------------------------------------------------------- #
# Per-invariant dispatch (block shapes mirror docs/gui/runner.py::compute_one)
# --------------------------------------------------------------------------- #

def _dispatch(A, item, events, hh_kwargs) -> tuple:
    kind = item.kind
    if kind in ("hh_cohomology", "hh_homology"):
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError",
                               f"{kind} needs a degree range, e.g. '{kind}:0..4'")
        method = (A.hochschild_cohomology if kind == "hh_cohomology"
                  else A.hochschild_homology)
        table = method(top, verbose=False, trace=events, **hh_kwargs)
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
        return {"dim": dim_z, "basis": [[str(x) for x in row] for row in basis],
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "dimension":
        keys = list(A.citations())
        return {"value": A.dim, "references": keys,
                "citations": _citation_pairs(keys)}, None
    raise ComputeError("SchemaError", f"unsupported computation {kind!r}")


def _citation_pairs(keys) -> list:
    if not keys:
        return []
    try:
        from quiverlab.trace.provenance import resolve_references as _rr
        return [list(p) for p in _rr(tuple(keys))]
    except KeyError:
        _log.warning("citation resolution failed for keys=%s", keys, exc_info=True)
        return []


def _latex_matrix(rows) -> str:
    body = r" \\ ".join(" & ".join(str(x) for x in row) for row in rows)
    return r"\begin{pmatrix} %s \end{pmatrix}" % body


def _rows(mat) -> list:
    return [list(r) for r in (mat.tolist() if hasattr(mat, "tolist") else mat)]


# --------------------------------------------------------------------------- #
# Module block (Plan 26)
# --------------------------------------------------------------------------- #

def _match_vertex(algebra, label):
    for v in algebra.quiver.vertices:
        if v == label or str(v) == str(label):
            return v
    raise ComputeError("SchemaError", f"module: no vertex {label!r} in the algebra")


def _parse_entry(x):
    if isinstance(x, bool):
        raise ComputeError("SchemaError", f"module entry {x!r} is not a number")
    if isinstance(x, int):
        return x
    if isinstance(x, float):
        raise ComputeError("ExactnessError",
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
            raise ComputeError("SchemaError",
                               f"module entry {x!r} is not an exact integer or fraction")
    raise ComputeError("SchemaError", f"module entry {x!r} is not a number")


def _full_matrices(algebra, mspec):
    rep = algebra if mspec.side == "right" else algebra.opposite()
    verts = list(rep.quiver.vertices)
    by_str = {str(v): v for v in verts}
    dimvec = {v: 0 for v in verts}
    for key, n in (mspec.dims or {}).items():
        if key not in by_str:
            raise ComputeError("SchemaError", f"module: no vertex {key!r} in the algebra")
        dimvec[by_str[key]] = int(n)
    start, off = {}, 0
    for v in verts:
        start[v] = off
        off += dimvec[v]
    n = off
    arrow_names = list(rep.quiver.arrows)
    for a in (mspec.maps or {}):
        if a not in arrow_names:
            raise ComputeError("SchemaError", f"module: no arrow {a!r} in the algebra")
    action = {}
    for a in arrow_names:
        s, t = rep.quiver.source(a), rep.quiver.target(a)
        rows, cols = dimvec[t], dimvec[s]
        full = [[0] * n for _ in range(n)]
        block = (mspec.maps or {}).get(a)
        if block is not None:
            if len(block) != rows or any(len(r) != cols for r in block):
                got = f"{len(block)}x{len(block[0]) if block else 0}"
                raise ComputeError("SchemaError",
                                   f"module map {a!r} must be {rows}x{cols} (target x "
                                   f"source dims); got {got}")
            for i in range(rows):
                for j in range(cols):
                    full[start[t] + i][start[s] + j] = _parse_entry(block[i][j])
        action[a] = full
    return dimvec, action


def _build_module(algebra, mspec, name):
    if mspec is None:
        raise ComputeError("SchemaError", f"module {name} block is required")
    if mspec.builtin is not None:
        b = mspec.builtin
        v = _match_vertex(algebra, b.vertex)
        builder = {"simple": algebra.simple, "projective": algebra.projective,
                   "injective": algebra.injective}[b.kind]
        return builder(v, side=mspec.side)
    dimvec, action = _full_matrices(algebra, mspec)
    return algebra.module(dimvec, action, side=mspec.side, name=name)


def _dv(dimvec) -> dict:
    return {str(v): int(n) for v, n in sorted(dimvec.items(), key=lambda kv: str(kv[0]))}


def _mod_view(m) -> dict:
    return {"dimvec": _dv(m.dimension_vector()), "dim": m.dim}


def _with_refs(block: dict, kind: str) -> dict:
    keys = _MOD_REFS[kind]
    block["references"] = list(keys)
    block["citations"] = _citation_pairs(keys)
    return block


def _dispatch_module(A, item, M, N) -> dict:
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
            raise ComputeError("SchemaError", "ext needs a degree range, e.g. 'ext:0..4'")
        from quiverlab.modules.ext import ext_dims
        dims = ext_dims(A, M, N, top)
        return _with_refs({"kind": "ext", "top": top,
                           "dims": [int(d) for d in dims],
                           "target": _mod_view(N)}, kind)
    if kind in ("projective_resolution", "injective_resolution"):
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError", f"{kind} needs a degree range, e.g. "
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
    raise ComputeError("SchemaError", f"unsupported module computation {kind!r}")


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
