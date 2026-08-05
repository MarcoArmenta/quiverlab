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
import sys
import time
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

# Absolute ceiling on a single module's total dimension (sum of its dimension
# vector), enforced at PARSE time -- before ``_full_matrices`` builds the
# per-arrow dense n x n action matrices (n = total module dim). This is the last
# line of defence for the WEBAPP entry paths (instant/queued/big, which all reach
# a module through :func:`parse_request` -> ``_parse_module`` and NEVER send an
# ``hpc`` block): without it a request like ``tor_target = {dims: {0: 10000}}``
# would allocate ~10^8 Python-int cells per arrow (multi-GB) regardless of tier.
# 2048 caps each arrow matrix at ~4.2M cells (~34 MB), and is ~5x above the largest
# module any test/legitimate request constructs (a dim-400 module, and only in a
# classification test that never builds it) -- so nothing legitimate is refused
# while a pathological module fails loudly instead of OOMing the box.
#
# The offline/cluster ``quiverlab-hpc run`` CLI is DIFFERENT: it is a local
# large-compute tool with its own memory budget. When its ``hpc`` block opts into
# large compute (``allow_large`` -- set by ``quiverlab-hpc run --allow-large``),
# this parse cap is BYPASSED (``parse_request`` passes ``max_total_dim=None``): the
# user has explicitly accepted the cost, and the real memory bound on that path is
# the CLI's own ``max_mem_bytes`` / RLIMIT deepen guard, not this webapp-shaped
# ceiling. The webapp never sets ``allow_large``, so it stays capped at 2048.
_MAX_MODULE_DIM = 2048

_log = logging.getLogger("quiverlab.hpc")

SIDES = ("right", "left")

MODULE_KINDS = frozenset({
    "dimension_vector", "rad_top_soc", "ext", "tor", "tau", "tau_minus",
    "projective_resolution", "injective_resolution",
    "projective_dimension", "injective_dimension", "decompose",
})
MODULE_RANGE_KINDS = frozenset({"ext", "tor", "projective_resolution",
                                "injective_resolution"})
# Module kinds with a Plan-30 Part-C worked-steps hook (quiverlab.trace.modules):
# an ``artifacts.pdf``-requesting module computation auto-emits the exhaustive
# worked-steps bundle, like HH. (decompose/tor have no hook yet -- follow-up.)
_MODULE_TRACE_KINDS = frozenset({
    "projective_resolution", "injective_resolution", "ext", "tau", "tau_minus",
})
# HH product surface (Plan 35): cup / cap / bracket / connes_b. Each backs a
# worked-steps chapter (quiverlab.trace.products) when a products request asks for
# ``artifacts.pdf`` and no HH/module trace already claimed the bundle.
PRODUCT_KINDS = frozenset({"cup", "cap", "bracket", "connes_b"})

# Honest labels for ``meta["pdf"]`` (the request flag is still named ``pdf``; the
# worked-steps report is now HTML + JSON, PDF/TeX output having been removed).
_WORKED_STEPS_OK = "worked steps in trace_steps.html"
_WORKED_STEPS_NO_HH = "no traced computation requested (the worked-steps report covers HH)"
_WORKED_STEPS_MODULE_FAIL = "worked-steps bundle could not be generated for this module computation"
_WORKED_STEPS_PRODUCT_FAIL = "worked-steps bundle could not be generated for this product computation"

_MOD_REFS = {
    "dimension_vector": ["assem_book"],
    "rad_top_soc": ["assem_book"],
    "tau": ["assem_book"],
    "tau_minus": ["assem_book"],
    "ext": ["module_ext"],
    "tor": ["minimal_resolution", "module_ext"],
    "projective_resolution": ["minimal_resolution"],
    "projective_dimension": ["minimal_resolution"],
    "injective_resolution": ["minimal_resolution", "assem_book"],
    "injective_dimension": ["minimal_resolution", "assem_book"],
    "decompose": ["assem_book"],
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
    tor_target: ModuleSpec | None    # Plan 30: the N in Tor^A_n(M, N), a left module
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
    # Degreewise dispatch always computes 0..hi (the engines return the full
    # 0..hi vector); a non-zero lower bound would be silently dropped, so reject it
    # loudly to match the GUI, which forbids lo != 0 (docs/gui/runner.py).
    if lo is not None and lo != 0:
        raise SpecError(
            f"compute range must start at 0 (got {s!r}); results are computed "
            f"0..N degreewise -- use '{m['kind']}:0..{hi}'")
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


# Families whose builders wrap OTHER algebras (no field= kwarg; parameters are
# base algebras given as PathAlgebra type strings) -- see build_algebra.
_WRAPPER_FAMILIES = frozenset({"TensorProduct", "TrivialExtension"})


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


def _parse_module(data, what: str, max_total_dim: int | None = _MAX_MODULE_DIM) -> ModuleSpec:
    """Validate a module block. ``max_total_dim`` caps the total module dimension
    at parse time (before any matrix is allocated); pass ``None`` to disable the cap
    -- the large-compute opt-in path (``hpc.allow_large``) does exactly that."""
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
    # Absolute allocation guard: the total module dimension drives the size of the
    # per-arrow dense action matrices ``_full_matrices`` builds (n x n each), so cap
    # it here -- before any matrix is allocated -- on the webapp entry paths. The
    # large-compute CLI opt-in (``hpc.allow_large``) passes ``max_total_dim=None`` to
    # bypass this ceiling (its own ``max_mem_bytes``/RLIMIT is the real bound).
    total = sum(int(n) for n in dims.values())
    if max_total_dim is not None and total > max_total_dim:
        raise SpecError(
            f"{what}: total module dimension {total} exceeds the {max_total_dim} "
            "cap; such a module would allocate a matrix too large to compute here -- "
            "narrow the dimension vector or run locally: pip install quiverlab")
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


def _parse_tor_target(data, max_total_dim: int | None = _MAX_MODULE_DIM) -> ModuleSpec | None:
    """Parse the ``tor_target`` (the N in Tor^A_n(M, N)) -- a LEFT A-module. An
    omitted side defaults to ``"left"`` (so it canonicalizes with an explicit
    ``"left"``); an explicit ``"right"`` is rejected loudly. Mirrors schema.py's
    ``_default_tor_side_left`` + the ``_module_rules`` side check. ``max_total_dim``
    is forwarded to :func:`_parse_module` (``None`` under ``hpc.allow_large``)."""
    if data is None:
        return None
    if not isinstance(data, dict):
        raise SpecError("tor_target must be an object")
    tt = dict(data)
    b = tt.get("builtin")
    has_side = "side" in tt or (isinstance(b, dict) and "side" in b)
    if not has_side:
        tt["side"] = "left"
    spec = _parse_module(tt, "tor_target", max_total_dim)
    if spec.side != "left":
        raise SpecError("Tor's second module 'tor_target' must be a LEFT A-module "
                        "(side='left'); Tor^A_n(M, N) pairs a right M with a left N")
    return spec


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
    # Parse the hpc block FIRST so the module-dimension cap can respect the
    # large-compute opt-in: ``allow_large`` (the offline/cluster CLI's --allow-large)
    # bypasses the 2048 parse ceiling. The webapp never sends an hpc block, so its
    # modules stay capped exactly as before.
    hpc = _parse_hpc(data["hpc"]) if data.get("hpc") is not None else None
    module_cap = None if (hpc is not None and hpc.allow_large) else _MAX_MODULE_DIM
    module = (_parse_module(data["module"], "module", module_cap)
              if data.get("module") is not None else None)
    ext_target = (_parse_module(data["ext_target"], "ext_target", module_cap)
                  if data.get("ext_target") is not None else None)
    tor_target = _parse_tor_target(data.get("tor_target"), module_cap)

    if (module is not None or ext_target is not None or tor_target is not None) \
            and schema_version != 2:
        raise SpecError("a 'module'/'ext_target'/'tor_target' block requires schema 2")
    kinds = {it.kind for it in items}
    if kinds & MODULE_KINDS and module is None:
        need = sorted(kinds & MODULE_KINDS)
        raise SpecError(f"module compute kind(s) {need} require a 'module' block")
    if "ext" in kinds and ext_target is None:
        raise SpecError("Ext needs a second module 'ext_target' (the N in Ext^n(M, N))")
    if "tor" in kinds and tor_target is None:
        raise SpecError("Tor needs a second module 'tor_target' (the N in "
                        "Tor^A_n(M, N), a LEFT A-module)")

    return ComputeRequest(schema_version=schema_version, algebra=algebra,
                          compute=list(compute), artifacts=artifacts,
                          module=module, ext_target=ext_target,
                          tor_target=tor_target, hpc=hpc,
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
    if spec.family in _WRAPPER_FAMILIES:
        # These builders take ALGEBRAS, not scalars, and no field= kwarg (the
        # result inherits the operands' field) -- the generic call below raised
        # a raw TypeError, so the catalog listed families nobody could build
        # (found live: TrivialExtension picked in the webapp form). A parameter
        # is a Dynkin/type STRING naming the base path algebra, built over the
        # request's field.
        f = _field(spec.field)
        args = {}
        for k, v in spec.params.items():
            if not isinstance(v, str) or not v.strip():
                raise ComputeError(
                    "CatalogError",
                    f"family {spec.family!r} parameter {k!r} names its base path "
                    f"algebra as a type string (e.g. \"A3\"); got {v!r}")
            args[k] = ql.PathAlgebra(v.strip(), field=f)
        return builder(**args)
    return builder(field=_field(spec.field), **spec.params)


# --------------------------------------------------------------------------- #
# Top-level entry point (ported from runner.run_spec, byte-stable)
# --------------------------------------------------------------------------- #

def run(req, artifact_dir, progress_cb: Callable[[dict], None] | None = None,
        result_max_bytes: int | None = None, *, result_schema: int | None = None,
        write_result: bool = True, capture_reps: bool = True) -> dict:
    """Build the algebra, run each requested computation, and (once the result
    JSON is under the byte cap) write ``result.json`` plus the requested artifacts
    into ``artifact_dir``. Returns the result dict.

    ``result_schema`` (CLI only) stamps the envelope with ``result_schema``; the
    webapp passes None so the returned dict is byte-identical to the pre-Plan-28
    runner. ``write_result=False`` (CLI) writes only the sidecar artifacts and
    leaves the authoritative ``result.json`` write to the caller.

    ``capture_reps=False`` skips the Plan-35 explicit-HH-representatives capture on the
    plain ``hh_cohomology`` / ``hh_homology`` dims blocks. The instant tier passes it:
    that capture runs the GF(p) bar route through ``engine.tt_calculus``, whose cold
    numba JIT (paid fresh in every spawned instant child) is tens of seconds -- far over
    the instant wall net -- while the reps only ever feed the report / GUI, which the
    instant tier discards. The dims block is byte-identical to before this wave when the
    flag is off, so a request that would have been instant stays instant."""
    t0_ns = time.monotonic_ns()
    req = parse_request(req)
    artifact_dir = Path(artifact_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    tikz_src = None
    hh_trace = None
    module_trace = None
    product_trace = None
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
            # ext/tor CONSUME their target, so those items build it unconditionally
            # (a missing block must still raise loudly there). A tau/tau^- request
            # additionally DISPLAYS the second module's translate (Marco,
            # 2026-07-29), which is opportunistic: built only when the request
            # actually supplies the target.
            def _target(spec, consumed_by):
                if any(it.kind == consumed_by for it in items):
                    return _build_module(A, spec, "N")
                if spec is not None and any(it.kind in ("tau", "tau_minus")
                                            for it in items):
                    return _build_module(A, spec, "N")
                return None

            N = _target(req.ext_target, "ext")
            T = _target(req.tor_target, "tor")
            results: dict = {}
            per_kind: dict = {}
            for i, item in enumerate(items):
                if progress_cb:
                    progress_cb({"step": i, "of": len(items), "kind": item.kind})
                t_item = time.monotonic_ns()
                try:
                    if item.kind in MODULE_KINDS:
                        results[item.kind] = _dispatch_module(A, item, M, N, T)
                        per_kind[item.kind] = _item_resources(t_item)
                        if (req.artifacts.pdf and module_trace is None
                                and item.kind in _MODULE_TRACE_KINDS):
                            # the first traceable module kind backs the worked-steps
                            # bundle (trace_steps.html); HH still takes precedence.
                            module_trace = (item.kind, item.hi, M, N)
                        continue
                    if _deepen_applies(req.hpc, item, req.algebra):
                        results[item.kind] = _dispatch_deepen(A, item, req.hpc, progress_cb)
                        per_kind[item.kind] = _item_resources(t_item)
                        continue
                    block, hh = _dispatch(A, item, events, hh_kwargs, capture_reps)
                    results[item.kind] = block
                    per_kind[item.kind] = _item_resources(t_item)
                    if hh is not None:
                        hh_trace = hh
                    elif (req.artifacts.pdf and product_trace is None
                          and item.kind in PRODUCT_KINDS):
                        # products are their own tables (no HH run): the first product
                        # kind backs the worked-steps bundle when nothing else claimed
                        # it. HH and module traces still take precedence at write time.
                        product_trace = (item.kind, item.hi)
                except qerr.DepthLimitError as exc:
                    # A single over-cap computation (a high-degree bar product/HH whose
                    # dense (co)chain matrices would exceed max_cells) must not sink the
                    # WHOLE request. Record an honest per-item error block and keep every
                    # other result + the worked-steps bundle -- parity with the Pyodide
                    # runner (docs/gui/runner.py), which already degrades per-computation.
                    # The report renders it "not computed -- DepthLimitError: ...", dims of
                    # the other invariants stand. Scoped to DepthLimitError ONLY (the
                    # cell-bound signal); schema/exactness/admissibility errors still
                    # propagate and fail the request, unchanged.
                    _log.info("compute item %s hit the cell bound; recording an honest "
                              "error block and continuing: %s", item.kind, exc)
                    results[item.kind] = {"error": {"type": type(exc).__name__,
                                                    "message": str(exc)},
                                          "references": []}
                    per_kind[item.kind] = _item_resources(t_item)
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
            meta["pdf"] = _WORKED_STEPS_NO_HH if hh_trace is None else _WORKED_STEPS_OK

        from quiverlab.trace.json_guide import build_json_guide
        result = {
            "quiverlab_version": getattr(ql, "__version__", "unknown"),
            "algebra": req.raw_algebra,
            "results": results,
            # Marco 2026-07-31 (ADDENDUM 2): a per-computation index of how to recover
            # every computed object from THIS result.json (generated from the actual
            # keys present, self-validated so it can never point at an absent object).
            "json_guide": build_json_guide(results),
            "references": resolve_references(used_keys),
            "reproduce": _snippet(req, A),
            "meta": meta,
        }
        if result_schema is not None:
            result["result_schema"] = result_schema
            # CLI-envelope-only echo of the computed-on modules as FULL
            # representations ({dims, maps} via the shared module_blocks
            # serializer), so the report can show the per-arrow action matrices.
            # Guarded by result_schema: the webapp passes None and its result
            # dict stays byte-identical (frozen goldens).
            for key, mod in (("module", M), ("ext_target", N), ("tor_target", T)):
                if mod is not None:
                    result[key] = {"side": mod.side, **_mod_repr(mod)}
            result["resources"] = _resources_used(t0_ns)
            result["resources"]["per_kind"] = per_kind
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
                                          used_keys, artifact_dir, meta,
                                          results=results, modules=_named(M, N, T))
        payload = json.dumps(result, indent=2, default=str)
    elif req.artifacts.pdf and module_trace is not None:
        m_kind, m_top, m_M, m_N = module_trace
        meta["pdf"] = _write_module_worked_steps(A, m_kind, m_top, m_M, m_N,
                                                 used_keys, artifact_dir, meta,
                                                 results=results,
                                                 modules=_named(M, N, T))
        payload = json.dumps(result, indent=2, default=str)
    elif req.artifacts.pdf and product_trace is not None:
        p_kind, p_top = product_trace
        meta["pdf"] = _write_product_worked_steps(A, p_kind, p_top, used_keys,
                                                  artifact_dir, meta,
                                                  results=results,
                                                  modules=_named(M, N, T))
        payload = json.dumps(result, indent=2, default=str)
    elif req.artifacts.pdf:
        # No traceable computation was requested (say: Cartan + centre only) -- the
        # report is still written, carrying the example and every computed result.
        # Without this the session's answers had nowhere to be saved (Marco
        # 2026-07-29).
        meta["pdf"] = _write_results_only(A, used_keys, artifact_dir, results,
                                          modules=_named(M, N, T))
        payload = json.dumps(result, indent=2, default=str)
    if write_result:
        # utf-8 explicitly: the default is the LOCALE codec (cp1252 on
        # Windows), which mangles any non-ASCII entry or citation.
        (artifact_dir / "result.json").write_text(payload, encoding="utf-8")
    if tikz_src is not None:
        (artifact_dir / "tikz.tex").write_text(tikz_src, encoding="utf-8")
    return result


def _item_resources(t0_ns: int) -> dict:
    """Per-computation footprint, exact ints: wall-clock ms for this item and
    the PROCESS peak RSS by the end of it (ru_maxrss is a high-water mark, so
    this is 'peak so far', not this item's own allocation)."""
    out = {"wall_ms": (time.monotonic_ns() - t0_ns) // 1_000_000}
    try:
        import resource as _res
        rss = _res.getrusage(_res.RUSAGE_SELF).ru_maxrss
        if sys.platform == "darwin":
            rss //= 1024
        out["peak_rss_mib"] = rss // 1024
    except Exception:
        pass
    return out


def _resources_used(t0_ns: int) -> dict:
    """Resources the run actually used / had available, as exact ints (the src/
    no-floats gate holds here): wall-clock ms, peak RSS MiB, detected usable
    cores and RAM. CLI-envelope-only (written beside ``result_schema``)."""
    out = {"wall_ms": (time.monotonic_ns() - t0_ns) // 1_000_000}
    try:
        import resource as _res
        rss = _res.getrusage(_res.RUSAGE_SELF).ru_maxrss   # KiB (Linux), B (macOS)
        if sys.platform == "darwin":
            rss //= 1024
        out["peak_rss_mib"] = rss // 1024
    except Exception:                                      # e.g. no resource module
        pass
    try:
        from quiverlab.hpc.resources import detect_resources
        det = detect_resources()
        out["cores_detected"] = det["cores"]
        if det.get("mem_bytes"):
            out["ram_mib_detected"] = det["mem_bytes"] // (1024 * 1024)
    except Exception:
        pass
    return out


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


def _promote_trace_artifacts(produced, artifact_dir) -> str:
    """Promote the writer's ``<stem>.html`` (+ sidecar ``.json``) to the fixed
    artifact names ``trace_steps.html`` and ``trace.json`` so ``pages.py`` can serve
    them (Plan 34: the .json machine record beside the print-ready HTML report).
    Returns the honest ``meta['pdf']`` label."""
    target = artifact_dir / "trace_steps.html"
    # The always-present JSON sidecar shares the produced stem (writer.py); promote it
    # to its stable name so every worked-steps run yields trace.json too.
    src = produced.with_suffix(".json")
    dst = artifact_dir / "trace.json"
    if src.exists() and src.resolve() != dst.resolve():
        src.replace(dst)
    if produced.resolve() != target.resolve():
        produced.replace(target)
    return _WORKED_STEPS_OK


def _named(M, N, T) -> list:
    """The modules the request named, as ``(label, Module)`` for the report's "The
    modules" section. ``N`` is the Ext target, ``T`` the Tor target; when both are
    present they are labelled apart so the reader knows which is which."""
    out = [("M", M)] if M is not None else []
    if N is not None and T is not None:
        out += [("N (Ext target)", N), ("N (Tor target)", T)]
    elif N is not None:
        out.append(("N", N))
    elif T is not None:
        out.append(("N", T))
    return out


def _write_worked_steps(events, table, A, kind, top, used_keys, artifact_dir,
                        meta=None, results=None, modules=()) -> str:
    from quiverlab.trace.writer import write_trace
    produced = Path(write_trace(list(events), table, algebra=A, kind=kind, top=top,
                                references=_trace_references(used_keys, events),
                                out_dir=str(artifact_dir), results=results,
                                modules=modules))
    return _promote_trace_artifacts(produced, artifact_dir)


def _write_results_only(A, used_keys, artifact_dir, results, modules=()) -> str:
    """The worked-steps bundle for a request with NO traceable computation: the
    example + every computed result block (Marco 2026-07-29 -- the report must save
    what the session computed, whatever was asked for). Best-effort, exactly like the
    module bundle: a rendering failure never loses the already-computed JSON."""
    from quiverlab.trace.writer import write_trace
    try:
        produced = Path(write_trace([], None, algebra=A, kind="results", top=0,
                                    references=_trace_references(used_keys, []),
                                    out_dir=str(artifact_dir), results=results,
                                    modules=modules))
        return _promote_trace_artifacts(produced, artifact_dir)
    except Exception:
        _log.exception("results-only worked-steps bundle failed")
        return _WORKED_STEPS_NO_HH


def _write_module_worked_steps(A, kind, top, M, N, used_keys, artifact_dir,
                               meta=None, results=None, modules=()) -> str:
    """Render the worked-steps bundle (HTML + JSON) for a MODULE computation via
    the Plan-30 Part-C trace hooks (``quiverlab.trace.modules``), mirroring the HH
    bundle. Best-effort: a trace failure never loses the already-computed JSON
    result -- it degrades to an honest note."""
    from quiverlab.trace import modules as _tm
    from quiverlab.trace.writer import write_trace
    try:
        if kind == "projective_resolution":
            events, _ = _tm.trace_projective_resolution(M, top)
        elif kind == "injective_resolution":
            events, _ = _tm.trace_injective_resolution(M, top)
        elif kind == "ext":
            events, _ = _tm.trace_ext(A, M, N, top)
        else:                                  # tau / tau_minus
            events, _ = _tm.trace_tau(M, kind=kind)
        produced = Path(write_trace(list(events), None, algebra=A, kind=kind,
                                    top=(top if top is not None else 0),
                                    references=_trace_references(used_keys, events),
                                    out_dir=str(artifact_dir), results=results,
                                    modules=modules))
        return _promote_trace_artifacts(produced, artifact_dir)
    except Exception:
        _log.exception("module worked-steps bundle failed for kind=%s", kind)
        return _WORKED_STEPS_MODULE_FAIL


def _product_object(A, kind, top):
    """Recompute the Task-1 product result object for the worked-steps chapter --
    the sibling of ``_write_module_worked_steps`` re-running the trace from its
    inputs (``A``, ``kind``, ``top``); the report re-derives, it does not cache."""
    method = {"cup": A.cup_products, "cap": A.cap_products,
              "bracket": A.gerstenhaber_brackets,
              "connes_b": A.connes_differentials}[kind]
    return method(top)


def _write_product_worked_steps(A, kind, top, used_keys, artifact_dir,
                                meta=None, results=None, modules=()) -> str:
    """Render the worked-steps bundle (HTML + JSON) for an HH-PRODUCT computation via
    the Plan-35 trace chapter (``quiverlab.trace.products``), mirroring the HH and
    module bundles. Best-effort: a trace failure never loses the already-computed
    JSON result -- it degrades to an honest note."""
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.writer import write_trace
    try:
        obj = _product_object(A, kind, top)
        events = products_chapter(A, kind, obj)
        produced = Path(write_trace(list(events), None, algebra=A, kind=kind,
                                    top=(top if top is not None else 0),
                                    references=_trace_references(used_keys, events),
                                    out_dir=str(artifact_dir), results=results,
                                    modules=modules))
        return _promote_trace_artifacts(produced, artifact_dir)
    except Exception:
        _log.exception("product worked-steps bundle failed for kind=%s", kind)
        return _WORKED_STEPS_PRODUCT_FAIL


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

def _dispatch(A, item, events, hh_kwargs, capture_reps=True) -> tuple:
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
        # Credit the resolution that actually computed the table: the engine
        # string names it, the references must too (Marco 2026-07-28: "if we
        # use Chouhy-Solotar we need to cite it").
        eng = str(table.engine).lower()
        for marker, ckey in (("chouhy", "chouhy_solotar"), ("bardzell", "bardzell")):
            if marker in eng and ckey not in keys:
                keys.append(ckey)
        block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                 "engine": table.engine, "references": keys,
                 "citations": _citation_pairs(keys)}
        # Plan 35 wave 3d: capture the explicit HH^n / HH_n representatives alongside the
        # dims (basis_classes / chain_basis / differentials / inner_dims per degree),
        # from the SAME dims path (GF(p) bar or Chouhy-Solotar). Additive block fields;
        # None (dims-only) when no representative route applies. Byte-identical Pyodide
        # twin (docs/gui/runner.py). The reader can read off HH^0's centre, HH^1's
        # derivations, HH^2's deformation cochain, HH_0's commutator residues.
        if capture_reps:                       # skipped by the instant tier (report-only
            from quiverlab.hochschild.hh_reps import hh_reps_blocks
            try:                               # data + a cold-JIT cost over its wall net)
                reps = hh_reps_blocks(A, kind, top, list(table.dims), table.engine)
            except Exception:                  # reps are ADDITIVE + best-effort: capture
                _log.warning("hh reps capture failed for %s; shipping dims only",
                             kind, exc_info=True)   # must NEVER break the dims block
                reps = None
            if reps:
                block.update(reps)
        return block, (table, table.kind, top)
    # Cyclic homology HC_0..HC_top (Connes (b, B) mixed complex). A range kind that
    # mirrors the hh_homology block (HHTable-based), but with NO worked-steps chapter
    # for now -- it is a dims table, so no hh_trace is returned (the None below).
    if kind == "cyclic_homology":
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError",
                               f"{kind} needs a degree range, e.g. '{kind}:0..4'")
        # Plan 35 wave 3b: capture the explicit HC representatives alongside the dims
        # (basis_classes / chain_basis / differentials / column_structure), from the
        # SAME (b, B) total complex -- additive block fields, byte-identical Pyodide twin.
        table, reps = A.cyclic_homology(top, with_reps=True)
        keys = ["cyclic"]
        block = {"kind": table.kind, "top": top, "dims": list(table.dims),
                 "engine": table.engine, "references": keys,
                 "citations": _citation_pairs(keys)}
        block.update(reps)
        return block, None
    # Per-invariant citation keys. NEVER A.citations() here: that set
    # ACCUMULATES across the run, so every block after (or beside) an HH
    # computation echoed the bar-resolution key -- the Cartan matrix was
    # citing Hochschild 1945 (Marco's report-example.pdf, 2026-07-28).
    if kind == "cartan":
        m = _rows(A.cartan_matrix())
        keys = ["assem_book"]
        return {"matrix": m, "latex": _latex_matrix(m),
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "coxeter_polynomial":
        import sympy
        p = A.coxeter_polynomial()
        keys = ["lenzing_delapena_spectral", "assem_book"]
        return {"latex": sympy.latex(p.as_expr()), "text": str(p.as_expr()),
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "global_dimension":
        g = A.global_dimension()
        keys = ["assem_book"]
        return {"text": str(g), "exact": bool(g.exact), "value": g.value,
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "center":
        dim_z, basis = A.center()
        keys = ["bar"]                     # Z(A) = HH^0(A) -- Hochschild's paper
        return {"dim": dim_z, "basis": [[str(x) for x in row] for row in basis],
                "references": keys, "citations": _citation_pairs(keys)}, None
    if kind == "dimension":
        keys = ["assem_book"]
        return {"value": A.dim, "references": keys,
                "citations": _citation_pairs(keys)}, None
    # Yoneda / Ext-algebra + Koszulity (Plan 38): a scalar kind on the algebra
    # block; the optional range gives the degree through which E(A) is computed
    # (default 6), read like the other range kinds. Both runners share the block
    # builder (modules.ext_algebra.ext_algebra_block), so they are byte-identical.
    if kind == "ext_algebra":
        top = item.hi if item.hi is not None else 6
        from quiverlab.modules.ext_algebra import ext_algebra_block
        block = ext_algebra_block(A, top)
        block["citations"] = _citation_pairs(block["references"])
        return block, None
    # Recognizer batch + type detection (Plan 38): a pure scalar kind on the
    # algebra block; per-flag honest errors, never a silent False.
    if kind == "recognizers":
        from quiverlab.invariants.recognizers import recognizers_block
        block = recognizers_block(A)
        block["citations"] = _citation_pairs(block["references"])
        return block, None
    # HH product surface (Plan 35): cup / cap / bracket / connes_b. Each library
    # method returns a frozen result object whose .blocks() IS the block dict
    # (kind/top/engine/basis/tables/window or hh_dims/matrices/ranks + references);
    # the only addition here is the resolved citation pairs, exactly as the other
    # kinds do above. No hh_trace (products are their own tables, not an HH run).
    if kind in ("cup", "cap", "bracket", "connes_b"):
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError",
                               f"{kind} needs a degree range, e.g. '{kind}:0..4'")
        method = {"cup": A.cup_products, "cap": A.cap_products,
                  "bracket": A.gerstenhaber_brackets,
                  "connes_b": A.connes_differentials}[kind]
        block = method(top).blocks()
        keys = list(block["references"])
        block["citations"] = _citation_pairs(keys)
        return block, None
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


def _dv_latex(dimvec) -> str:
    """Dimension vector as display latex, byte-identical to the Pyodide
    runner's composition (both feed the same block renderers)."""
    d = _dv(dimvec)
    return "(" + ",\\, ".join(str(d[k]) for k in d) + ")" if d else "()"


def _homdim_latex(op: str, value) -> str:
    """Display latex for a homological dimension block (``pd``/``id``).

    A finite value states the equality. An UNRESOLVED probe is not a proof of
    infinity -- the resolution merely did not terminate by ``_PD_BOUND`` -- so it
    states the certified lower bound ``> N``, the same honesty ``global_dimension``
    already uses ("certified lower bound; not resolved within depth 32"). Without
    this key the draw-page renderer typeset a literal "undefined" (Marco,
    2026-07-29)."""
    if value is None:
        return r"\operatorname{%s} M > %d" % (op, _PD_BOUND)
    return r"\operatorname{%s} M = %d" % (op, value)


_HOMDIM_UNRESOLVED = ("certified lower bound; the resolution did not terminate "
                      "within the probed depth %d" % _PD_BOUND)


def _mod_view(m) -> dict:
    return {"dimvec": _dv(m.dimension_vector()), "dim": m.dim}


def _mod_repr(m) -> dict:
    """A module as the no-code INPUT schema ``{"dims": {v: n}, "maps": {arrow: [[..]]}}``
    (Plan 34, Marco): the per-vertex dimension VECTOR and the exact per-arrow action
    matrices. The redundant total dim is intentionally dropped (the vector carries it),
    and the shape mirrors the module INPUT so a computed rad/top/soc feeds straight back
    into the panel. Both tiers route ``rad_top_soc`` through the SAME library serializer
    (``quiverlab.modules.qpa_module.module_blocks``), so the webapp and the Pyodide GUI
    cannot drift."""
    from quiverlab.modules.qpa_module import module_blocks
    return module_blocks(m)


def _with_refs(block: dict, kind: str) -> dict:
    keys = _MOD_REFS[kind]
    block["references"] = list(keys)
    block["citations"] = _citation_pairs(keys)
    return block


def _sort_key(v):
    """Deterministic vertex order for displayed summands (ints numerically, then
    everything else by string)."""
    return (0, v) if isinstance(v, int) and not isinstance(v, bool) else (1, str(v))


def _summands_latex(vertices, letter: str) -> str:
    """A resolution term's summand multiset as LaTeX, e.g. ``P_{1}^{2} \\oplus P_{3}``
    (``letter`` = ``P`` projectives / ``I`` injectives). ``0`` for the zero term."""
    if not vertices:
        return "0"
    counts: dict = {}
    for v in vertices:
        counts[v] = counts.get(v, 0) + 1
    parts = []
    for v in sorted(counts, key=_sort_key):
        base = f"{letter}_{{{v}}}"
        parts.append(base if counts[v] == 1 else f"{base}^{{{counts[v]}}}")
    return r" \oplus ".join(parts)


# A single differential past this many cells ships elided (shape only): the
# matrices are display data for the report, and one enormous syzygy map must not
# balloon result.json (the Plan-34 recorder backstop, same figure).
_MAX_DIFF_CELLS = 250_000


def _differential_blocks(res, n_terms) -> list:
    """The resolution's maps as exact matrices (rows: target basis, columns:
    source basis, vertex-ordered). Projective: entry 0 is the augmentation
    eps: Q_0 -> M, entry n is d_n: Q_n -> Q_{n-1}. Injective: entry 0 is
    iota: M -> E^0, entry n is d^n: E^{n-1} -> E^n."""
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
    its projective / injective summands (Plan 35 UNIT 2). ``term_basis[n]`` lists one
    label per basis vector of term n, in the SAME order the differentials use, so
    ``len(term_basis[n])`` equals the term's k-dimension (= the differential's column
    count for a projective resolution, row count for an injective one). A label is
    ``"P_v: <path>"`` (projective) / ``"I_v: <path>"`` (injective -- a dual basis over
    A^op). Returns ``None`` (caller omits the field) when the algebra carries no path
    basis (a structure-constants algebra) or the total basis is implausibly large; the
    renderers tolerate the field's absence. Shared shape with the Pyodide twin in
    ``docs/gui/runner.py`` so the tiers cannot drift."""
    try:
        from quiverlab.modules.builders import projective
        if kind == "projective_resolution":
            base, sym = M.algebra, "P"
        else:
            from quiverlab.modules.opposite import opposite_algebra
            base, sym = opposite_algebra(M.algebra), "I"
        cache: dict = {}

        def paths(v):
            if v not in cache:
                cache[v] = list(projective(base, v)._pv_basis_labels)
            return cache[v]

        out, total = [], 0
        for n in range(len(res.terms)):
            labels = []
            for v in res.term(n):
                labels.extend(f"{sym}_{v}: {p}" for p in paths(v))
            total += len(labels)
            if total > _MAX_DIFF_CELLS:               # implausible; omit (bloat guard)
                return None
            out.append(labels)
        return out
    except Exception:
        return None


def _summand_view(mod, mult) -> dict:
    """One indecomposable summand: its dimension vector, multiplicity, and either a
    STANDARD name (S_v / P_v / I_v) or its full per-arrow matrices -- a dimension
    vector alone does not say which module a non-standard summand is (Marco
    2026-07-29). Routed through the SAME library serializer as the Pyodide runner
    (``quiverlab.modules.qpa_module.summand_blocks``) so the tiers cannot drift."""
    from quiverlab.modules.qpa_module import summand_blocks
    return summand_blocks(mod, mult)


def _decompose(M):
    """(is_indecomposable, decompose) from the Plan-30 Part-A engine, or ``None``
    when it is not yet importable (so callers degrade to the pre-Plan-30 shape)."""
    try:
        from quiverlab.modules.decompose import decompose, is_indecomposable
    except ImportError:
        return None
    return is_indecomposable, decompose


def _input_certificate(M) -> dict:
    """Certify the INPUT of an AR translate (Marco #1): ``indecomposable: true`` or
    its Krull-Schmidt ``decomposition`` + the additivity note (tau is additive, so a
    decomposable input's tau is still exact -- computed summand-wise).

    Best-effort + honest: returns ``{}`` (no claim) when the decomposition engine is
    unavailable (Part A not merged) OR cannot certify within budget -- the engine is
    LOUD (e.g. char <= dim M), and tau itself stays valid, so the block simply omits
    the certificate rather than assert something it could not prove. The block is
    then byte-identical to the pre-Plan-30 shape (an explicit ``decompose`` request
    still surfaces the loud refusal for the user who asks for it)."""
    eng = _decompose(M)
    if eng is None:
        return {}
    is_indec, decompose = eng
    try:
        if is_indec(M):
            return {"indecomposable": True}
        return {"indecomposable": False,
                "decomposition": [_summand_view(s, m) for (s, m) in decompose(M)],
                "note_key": "mod.tau_additive"}
    except qerr.QuiverlabError:
        return {}


def _ar_translate(mod, kind: str, name: str) -> dict:
    """One AR translate as a self-contained display payload: the symbol, its
    dimension-vector latex, the FULL representation ({dims, maps}) and the input's
    indecomposability certificate. Shared by the ``M`` block and by the second
    module's entries below, so both carry the same fields."""
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


def _target_translates(kind: str, targets) -> list:
    """The AR translates of the SECOND module(s) N -- the Ext/Tor argument -- so a
    tau block covers every module the request names, with full matrices (Marco,
    2026-07-29). ``targets`` is a list of ``(role, module)``; an entry whose
    translate refuses LOUDLY is reported as an honest error entry rather than
    failing the whole block (tau M is already computed and valid)."""
    out = []
    for role, mod in targets:
        if mod is None:
            continue
        try:
            entry = _ar_translate(mod, kind, "N")
        except qerr.QuiverlabError as exc:
            entry = {"name": "N", "error": str(exc)}
        entry["role"] = role
        out.append(entry)
    return out


def _dispatch_module(A, item, M, N, T=None) -> dict:
    kind = item.kind
    if kind == "dimension_vector":
        # "latex" mirrors the Pyodide runner byte-for-byte: the draw page's
        # renderer typesets block.latex, and a missing key rendered a literal
        # "undefined" (Marco's report-example.pdf, 2026-07-28).
        return _with_refs({"kind": "dimension_vector", "side": M.side,
                           **_mod_view(M),
                           "latex": r"\underline{\dim}\, M = "
                                    + _dv_latex(M.dimension_vector())}, kind)
    if kind == "rad_top_soc":
        # "series" is the Loewy (radical) series top-to-bottom (Plan 37): each entry a
        # str-keyed composition-factor multiplicity dict (M.loewy_layers()), so the
        # report and both GUIs render the stacked Loewy diagram. Mirrors the Pyodide
        # runner byte-for-byte.
        return _with_refs({"kind": "rad_top_soc", "side": M.side,
                           "radical": _mod_repr(M.radical()),
                           "top": _mod_repr(M.top()),
                           "socle": _mod_repr(M.socle()),
                           "series": [dict(layer) for layer in M.loewy_layers()]}, kind)
    if kind in ("tau", "tau_minus"):
        # tau of a projective (dually tau^- of an injective) IS zero -- the shared
        # helper says so explicitly; renderers typeset block.latex (mirrors the
        # Pyodide runner). The AR translate IS a module, so it ships as a full
        # representation ({dims, maps}, like rad/top/soc) and reports/GUIs can show
        # the per-arrow action matrices (Marco, 2026-07-28).
        entry = _ar_translate(M, kind, "M")
        entry.pop("name", None)
        block = {"kind": kind, **entry}
        # ... and the same for the SECOND module N when the request names one
        # (Marco, 2026-07-29). Omitted entirely when there is no target, so a
        # single-module request's block is byte-identical.
        targets = _target_translates(kind, [("ext_target", N), ("tor_target", T)])
        if targets:
            block["targets"] = targets
        return _with_refs(block, kind)
    if kind == "decompose":
        eng = _decompose(M)
        if eng is None:
            raise SpecError("module decomposition requires the Krull-Schmidt engine "
                            "(Plan 30 Part A); not available in this build")
        _, decompose = eng
        summands = [_summand_view(s, m) for (s, m) in decompose(M)]
        return _with_refs({"kind": "decompose", "side": M.side,
                           "summands": summands, "iso_classes": len(summands)}, kind)
    if kind == "ext":
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError", "ext needs a degree range, e.g. 'ext:0..4'")
        from quiverlab.modules.ext import ext_dims
        # Plan 35 wave 3a: capture the explicit Ext cocycle representatives alongside
        # the dims (basis_classes / chain_basis / differentials per degree), from the
        # SAME Hom complex. Additive keys; renderers dispatch by kind.
        # Plan 35 wave 3c: interpret=True also captures the Yoneda exact sequence
        # 0 -> N -> Q -> ... -> M -> 0 realizing each class (an `interpretation` key).
        dims, reps = ext_dims(A, M, N, top, with_reps=True, interpret=True)
        block = {"kind": "ext", "top": top, "dims": [int(d) for d in dims],
                 "target": _mod_view(N),
                 # Marco 2026-08-03: say WHICH module the engine resolved and by
                 # WHICH resolution (renderers state it before the numbers).
                 "resolved": {"module": "M", "side": M.side,
                              "resolution": "minimal projective resolution"}}
        block.update(reps)
        return _with_refs(block, kind)
    if kind == "tor":
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError", "tor needs a degree range, e.g. 'tor:0..4'")
        try:
            from quiverlab.modules.tor import tor_dims
        except ImportError:
            raise SpecError("Tor requires the Tor engine (Plan 29); not available in "
                            "this build")
        # Plan 35 wave 3a: capture the explicit Tor cycle representatives alongside
        # the dims (Tor_0 = M (x)_A N as the cokernel), from the SAME tensor complex.
        dims, reps = tor_dims(A, M, T, top, with_reps=True)
        block = {"kind": "tor", "top": top, "dims": [int(d) for d in dims],
                 "target": _mod_view(T),
                 # Tor resolves the RIGHT module M and tensors with the left N
                 # (modules/tor.py, resolve="first").
                 "resolved": {"module": "M", "side": M.side,
                              "resolution": "minimal projective resolution"}}
        block.update(reps)
        return _with_refs(block, kind)
    if kind in ("projective_resolution", "injective_resolution"):
        top = item.hi
        if top is None:
            raise ComputeError("SchemaError", f"{kind} needs a degree range, e.g. "
                               f"'{kind}:0..4'")
        res = (M.projective_resolution(top) if kind == "projective_resolution"
               else M.injective_resolution(top))
        letter = "P" if kind == "projective_resolution" else "I"
        terms = [_dv(dv) for dv in res.dimension_vectors()]
        block = {"kind": kind, "top": top, "terms": terms,
                 "betti": [res.betti(i) for i in range(len(terms))],
                 "summands": [_summands_latex(res.term(i), letter)
                              for i in range(len(terms))],
                 "differentials": _differential_blocks(res, len(terms))}
        term_basis = _term_basis_blocks(res, kind, M)
        if term_basis is not None:
            block["term_basis"] = term_basis
        if kind == "projective_resolution":
            block["pd"] = res.pd()
        else:
            block["injective_dimension"] = res.injective_dimension()
        return _with_refs(block, kind)
    if kind == "projective_dimension":
        pd = M.projective_resolution(_PD_BOUND).pd()
        return _with_refs({"kind": "projective_dimension", "value": pd,
                           "finite": pd is not None, "bound": _PD_BOUND,
                           "latex": _homdim_latex("pd", pd),
                           **({} if pd is not None
                              else {"note": _HOMDIM_UNRESOLVED})}, kind)
    if kind == "injective_dimension":
        idim = M.injective_dimension(bound=_PD_BOUND)
        return _with_refs({"kind": "injective_dimension", "value": idim,
                           "finite": idim is not None, "bound": _PD_BOUND,
                           "latex": _homdim_latex("id", idim),
                           **({} if idim is not None
                              else {"note": _HOMDIM_UNRESOLVED})}, kind)
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
    if req.tor_target is not None:
        lines += _module_construction(req.tor_target, "N", A)
    _snip = {"hh_cohomology": lambda it: f"A.hochschild_cohomology({it.hi})",
             "hh_homology": lambda it: f"A.hochschild_homology({it.hi})",
             "cyclic_homology": lambda it: f"A.cyclic_homology({it.hi})",
             "coxeter_polynomial": lambda it: "A.coxeter_polynomial()",
             "cartan": lambda it: "A.cartan_matrix()",
             "global_dimension": lambda it: "A.global_dimension()",
             "center": lambda it: "A.center()",
             "dimension": lambda it: "A.dim",
             "ext_algebra":
                 lambda it: f"A.ext_algebra({it.hi if it.hi is not None else 6})",
             "recognizers": lambda it: ("[A.is_semisimple(), A.is_hereditary(), "
                                        "A.is_gentle(), A.dynkin_type(), "
                                        "A.form_type()]"),
             "cup": lambda it: f"A.cup_products({it.hi})",
             "cap": lambda it: f"A.cap_products({it.hi})",
             "bracket": lambda it: f"A.gerstenhaber_brackets({it.hi})",
             "connes_b": lambda it: f"A.connes_differentials({it.hi})",
             "dimension_vector": lambda it: "M.dimension_vector()",
             "rad_top_soc": lambda it: "M.radical(), M.top(), M.socle()",
             "tau": lambda it: "M.tau()",
             "tau_minus": lambda it: "M.tau_minus()",
             "ext": lambda it: f"[A.ext(M, N, i) for i in range({it.hi} + 1)]",
             "tor": lambda it: ("from quiverlab.modules.tor import tor_dims\n"
                                f"tor_dims(A, M, N, {it.hi})"),
             "decompose": lambda it: ("from quiverlab.modules.decompose import "
                                      "decompose\ndecompose(M)"),
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
