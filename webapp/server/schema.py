"""Versioned request schema. Everything is data — a family id with typed params,
or a quiver (vertices / arrows / opaque relation strings) — plus a field spec.
There is no code path that evaluates user strings; relation strings are parsed
later, loudly, by the library's exact relation grammar.

Schema v2 (Plan 26) adds an optional ``module`` block so a representation theorist
can specify a module with zero code: either explicitly (a dimension vector + one
exact-entry matrix per arrow) or via a zero-typing pick-list (a built-in
``simple`` / ``projective`` / ``injective`` at a vertex). Matrix entries stay pure
DATA (ints or exact strings like ``"1/2"``) — never evaluated, never floats; the
exact parse into the chosen field happens later, loudly, in the runner."""
from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction
from typing import Annotated, Any, Literal, Union

from pydantic import BaseModel, Field, field_validator, model_validator


class SchemaError(ValueError):
    pass


class FieldSpec(BaseModel):
    kind: Literal["CC", "GF"]
    p: int | None = None
    n: int = 1

    @field_validator("p")
    @classmethod
    def _p_positive(cls, v: int | None) -> int | None:
        if v is not None and v < 2:
            raise SchemaError("GF(p): p must be a prime >= 2")
        return v


def _reject_float_params(value: Any, path: str) -> None:
    """Recursively refuse any float/complex inside family params. The library is
    exact-only, so a float param can never build -- it would surface as a genericised
    500 in the runner. Vetting it at the boundary yields a clean typed 4xx instead
    (consistent with how matrix entries reject floats). Ints/strings/bools pass."""
    if isinstance(value, (float, complex)):
        raise SchemaError(
            f"family param {path} is a float ({value!r}); parameters must be exact "
            "(integers or strings, never floats)")
    if isinstance(value, dict):
        for k, v in value.items():
            _reject_float_params(v, f"{path}[{k!r}]")
    elif isinstance(value, (list, tuple)):
        for i, v in enumerate(value):
            _reject_float_params(v, f"{path}[{i}]")


class FamilyAlgebraSpec(BaseModel):
    kind: Literal["family"]
    family: str
    params: dict[str, Any] = Field(default_factory=dict)
    field: FieldSpec

    @field_validator("params")
    @classmethod
    def _params_exact(cls, v: dict[str, Any]) -> dict[str, Any]:
        for key, val in v.items():
            _reject_float_params(val, repr(key))
        return v


class QuiverAlgebraSpec(BaseModel):
    # The Plan-10 GUI (docs/gui/runner.py::run_build) emits this shape:
    # vertices, arrows {name: [source, target]}, opaque relation strings.
    kind: Literal["quiver"]
    vertices: list[int]
    arrows: dict[str, tuple[int, int]]
    relations: list[str] = Field(default_factory=list)
    field: FieldSpec

    @field_validator("vertices")
    @classmethod
    def _nonempty_vertices(cls, v: list[int]) -> list[int]:
        if not v:
            raise SchemaError("algebra.vertices must be a non-empty list of integers")
        return v


# Discriminated union: pydantic routes on `kind`, so a validation error names
# the kind rather than complaining about every member of the union.
AlgebraSpec = Annotated[
    Union[FamilyAlgebraSpec, QuiverAlgebraSpec],
    Field(discriminator="kind"),
]


class Artifacts(BaseModel):
    pdf: bool = False
    tikz: bool = False


# --------------------------------------------------------------------------- #
# Schema v2: the no-code module block (Plan 26)
# --------------------------------------------------------------------------- #

SIDES = ("right", "left")

# The module compute kinds served on top of a `module` block. `ext`/`tor` also need
# a second module (`ext_target`/`tor_target`); the rest act on the single `module`.
# Kept here so both the schema (routing rules) and the runner agree on the set by
# construction. `tor`/`decompose` are Plan 30 additions (Tor needs the Plan-29
# engine; decompose is the Krull-Schmidt splitter).
MODULE_KINDS = frozenset({
    "dimension_vector", "rad_top_soc", "ext", "tor", "tau", "tau_minus",
    "projective_resolution", "injective_resolution",
    "projective_dimension", "injective_dimension", "decompose",
    "tilting_check",
})
# Module kinds that consume a degree range (`kind:0..n`); the rest are scalars.
MODULE_RANGE_KINDS = frozenset({"ext", "tor", "projective_resolution",
                                "injective_resolution"})


def _valid_entry(x: Any) -> bool:
    """A matrix entry is exact DATA: a JSON integer or an exact string literal
    (e.g. ``"3"``, ``"-2"``, ``"1/2"``). Never a float (exactness is the whole
    point) and never a bool (``True``/``False`` are not field entries). The
    string is NOT parsed here -- the exact parse into the chosen field happens in
    the runner, where the domain (GF(p) vs CC) is known -- but its lexical shape
    is vetted so a malformed entry fails loudly at request time, not with a 500."""
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
            Fraction(s)                       # accepts "1/2"; also exact decimals
            return True
        except (ValueError, ZeroDivisionError):
            return False
    return False


class BuiltinModule(BaseModel):
    """A zero-typing pick-list module: the simple / projective / injective at a
    vertex (the library's ``A.simple/projective/injective(v, side=)`` builders)."""
    kind: Literal["simple", "projective", "injective"]
    vertex: Any                               # a vertex label (int for quiver/family)


class ModuleSpec(BaseModel):
    """A module, either specified explicitly (``dims`` + one exact-entry matrix per
    arrow in ``maps``) or via a ``builtin`` pick-list. ``side`` defaults to
    ``"right"`` and is always emitted, so an omitted side and an explicit
    ``"right"`` canonicalize to the SAME cache key (Plan 25)."""
    dims: dict[str, int] | None = None
    maps: dict[str, list[list[Any]]] | None = None
    builtin: BuiltinModule | None = None
    side: Literal["right", "left"] = "right"

    @model_validator(mode="before")
    @classmethod
    def _lift_builtin_side(cls, data):
        """Accept the task's ``{"builtin": {"kind", "vertex", "side"}}`` shape by
        lifting a ``side`` nested inside ``builtin`` up to the top level, so the
        canonical form carries ``side`` in exactly one place (a builtin and an
        explicit module both dump ``side`` as a sibling field)."""
        if isinstance(data, dict) and isinstance(data.get("builtin"), dict):
            b = dict(data["builtin"])
            if "side" in b:
                inner = b.pop("side")
                outer = data.get("side")
                if outer is not None and outer != inner:
                    raise SchemaError(
                        f"conflicting side: builtin.side={inner!r} vs side={outer!r}")
                data = {**data, "builtin": b, "side": inner}
        return data

    @model_validator(mode="after")
    def _one_form(self):
        if self.builtin is not None:
            if self.dims is not None or self.maps is not None:
                raise SchemaError("module: give either a builtin pick-list OR "
                                  "dims+maps, not both")
            return self
        if self.dims is None:
            raise SchemaError("module: needs 'dims' (a dimension vector) or a "
                              "'builtin' pick-list")
        for v, n in self.dims.items():
            if not isinstance(n, int) or isinstance(n, bool) or n < 0:
                raise SchemaError(f"module dims[{v!r}] must be a non-negative integer")
        for arrow, mat in (self.maps or {}).items():
            width = None
            for row in mat:
                if not isinstance(row, list):
                    raise SchemaError(f"module maps[{arrow!r}] must be a matrix "
                                      "(list of rows)")
                if width is None:
                    width = len(row)
                elif len(row) != width:
                    raise SchemaError(f"module maps[{arrow!r}] is not rectangular")
                for x in row:
                    if not _valid_entry(x):
                        raise SchemaError(
                            f"module maps[{arrow!r}] has a non-exact entry {x!r}; "
                            "entries must be integers or exact strings like '1/2' "
                            "(never floats)")
        return self


class ComputeRequest(BaseModel):
    schema_version: int = Field(1, alias="schema")
    algebra: AlgebraSpec
    compute: list[str]
    artifacts: Artifacts = Field(default_factory=Artifacts)
    module: ModuleSpec | None = None          # v2 (Plan 26)
    ext_target: ModuleSpec | None = None      # v2: the N in Ext^n(M, N), a RIGHT A-module
    tor_target: ModuleSpec | None = None      # v2 (Plan 30): the N in Tor^A_n(M, N)

    @model_validator(mode="before")
    @classmethod
    def _default_tor_side_left(cls, data):
        """A ``tor_target`` (the N in Tor^A_n(M, N)) is a LEFT A-module. Default an
        omitted side to ``"left"`` BEFORE the ModuleSpec validates (whose own default
        is ``"right"``), so an omitted side and an explicit ``"left"`` canonicalize to
        the SAME cache key (mirrors the Plan-25/26 side-default discipline). An
        explicit ``"right"`` is left in place and rejected loudly after validation."""
        if isinstance(data, dict) and isinstance(data.get("tor_target"), dict):
            tt = dict(data["tor_target"])
            b = tt.get("builtin")
            has_side = "side" in tt or (isinstance(b, dict) and "side" in b)
            if not has_side:
                tt["side"] = "left"
                data = {**data, "tor_target": tt}
        return data

    @field_validator("schema_version")
    @classmethod
    def _schema_known(cls, v: int) -> int:
        if v not in (1, 2):
            raise SchemaError(f"unsupported schema version {v}; this server speaks v1/v2")
        return v

    @field_validator("compute")
    @classmethod
    def _nonempty(cls, v: list[str]) -> list[str]:
        if not v:
            raise SchemaError("compute list must not be empty")
        for item in v:
            parse_compute_item(item)          # validates each entry
        return v

    @model_validator(mode="after")
    def _module_rules(self):
        """The module block is a v2 feature; and any module compute kind needs a
        ``module`` (and, for ``ext``/``tor``, an ``ext_target``/``tor_target``)."""
        if (self.module is not None or self.ext_target is not None
                or self.tor_target is not None) and self.schema_version != 2:
            raise SchemaError("a 'module'/'ext_target'/'tor_target' block requires "
                              "schema 2")
        kinds = {parse_compute_item(s).kind for s in self.compute}
        if kinds & MODULE_KINDS and self.module is None:
            need = sorted(kinds & MODULE_KINDS)
            raise SchemaError(f"module compute kind(s) {need} require a 'module' block")
        if "ext" in kinds and self.ext_target is None:
            raise SchemaError("Ext needs a second module 'ext_target' (the N in "
                              "Ext^n(M, N))")
        if "tor" in kinds and self.tor_target is None:
            raise SchemaError("Tor needs a second module 'tor_target' (the N in "
                              "Tor^A_n(M, N), a LEFT A-module)")
        if self.ext_target is not None and self.ext_target.side != "right":
            raise SchemaError("Ext's second module 'ext_target' must be a RIGHT "
                              "A-module (side='right'); Ext^n(M, N) pairs a right M "
                              "with a right N")
        if self.tor_target is not None and self.tor_target.side != "left":
            raise SchemaError("Tor's second module 'tor_target' must be a LEFT "
                              "A-module (side='left'); Tor^A_n(M, N) pairs a right M "
                              "with a left N")
        return self

    def model_dump(self, *args, **kwargs):
        """Drop ``module``/``ext_target``/``tor_target`` when absent so a non-module
        request canonicalizes byte-identically to the pre-Plan-26 shape (the cache
        key of every existing family/quiver request -- and every Plan-26 ext request
        -- is unchanged; only genuine Tor requests carry the extra block)."""
        d = super().model_dump(*args, **kwargs)
        for k in ("module", "ext_target", "tor_target"):
            if d.get(k) is None:
                d.pop(k, None)
        return d


@dataclass(frozen=True)
class ComputeItem:
    kind: str
    lo: int | None
    hi: int | None


_RANGE = re.compile(r"^(?P<kind>[a-z_]+)(?::(?P<lo>\d+)\.\.(?P<hi>\d+))?$")


def parse_compute_item(s: str) -> ComputeItem:
    m = _RANGE.match(s)
    if not m:
        raise SchemaError(f"unparseable compute item {s!r}")
    lo = int(m["lo"]) if m["lo"] is not None else None
    hi = int(m["hi"]) if m["hi"] is not None else None
    if lo is not None and hi is not None and hi < lo:
        raise SchemaError(f"empty degree range in {s!r}")
    # Degreewise dispatch always computes 0..hi; a non-zero lower bound would be
    # silently dropped by the runner, so reject it here to match the GUI (which
    # forbids lo != 0, docs/gui/runner.py) -- server and GUI agree.
    if lo is not None and lo != 0:
        raise SchemaError(
            f"compute range must start at 0 (got {s!r}); results are computed "
            f"0..N degreewise -- use '{m['kind']}:0..{hi}'")
    return ComputeItem(kind=m["kind"], lo=lo, hi=hi)
