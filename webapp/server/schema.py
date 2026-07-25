"""Versioned request schema. Everything is data — a family id with typed params,
or a quiver (vertices / arrows / opaque relation strings) — plus a field spec.
There is no code path that evaluates user strings; relation strings are parsed
later, loudly, by the library's exact relation grammar."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Annotated, Any, Literal, Union

from pydantic import BaseModel, Field, field_validator


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


class FamilyAlgebraSpec(BaseModel):
    kind: Literal["family"]
    family: str
    params: dict[str, Any] = Field(default_factory=dict)
    field: FieldSpec


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


class ComputeRequest(BaseModel):
    schema_version: int = Field(1, alias="schema")
    algebra: AlgebraSpec
    compute: list[str]
    artifacts: Artifacts = Field(default_factory=Artifacts)

    @field_validator("schema_version")
    @classmethod
    def _schema_one(cls, v: int) -> int:
        if v != 1:
            raise SchemaError(f"unsupported schema version {v}; this server speaks v1")
        return v

    @field_validator("compute")
    @classmethod
    def _nonempty(cls, v: list[str]) -> list[str]:
        if not v:
            raise SchemaError("compute list must not be empty")
        for item in v:
            parse_compute_item(item)          # validates each entry
        return v


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
    return ComputeItem(kind=m["kind"], lo=lo, hi=hi)
