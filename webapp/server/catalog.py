"""Catalog built by introspecting ``quiverlab.families()``. Generated, not
hand-maintained, so it cannot diverge from the library's real offering.

The real ``families()`` (confirmed by probing the installed library at Plans
10-19 state, per the 2026-07-24 amendment) returns a ``FamilyListing`` of
``FamilyInfo`` records with fields ``name``, ``signature``, ``route``,
``citations``, ``summary`` -- NOT ``(name, builder)`` pairs. Each family name is
a top-level ``quiverlab`` export, so the builder is ``getattr(quiverlab,
info.name)`` and its parameters come from ``inspect.signature(builder)``. The
``zoo`` entry is excluded: it yields a batch of algebras, not a single buildable
family. This module imports only the PUBLIC ``quiverlab`` surface."""
from __future__ import annotations

import inspect
from functools import lru_cache

import quiverlab as ql


class CatalogError(ValueError):
    pass


def _iter_families():
    """Yield ``(name, builder)`` for each buildable v1 family.

    ``ql.families()`` yields ``FamilyInfo`` records; the builder is the
    top-level ``quiverlab`` export named ``info.name``. ``zoo`` is skipped (it
    is a batch iterator, not a single family). A missing export is drift -- we
    raise loudly rather than silently drop it."""
    for info in ql.families():
        name = info.name
        if name == "zoo":
            continue
        builder = getattr(ql, name, None)
        if builder is None:
            raise CatalogError(
                f"families() lists {name!r} but quiverlab has no such export; "
                "the family catalog has drifted from the library surface"
            )
        yield name, builder


def _params_of(builder) -> list[dict]:
    out = []
    try:
        sig = inspect.signature(builder)
    except (TypeError, ValueError):
        return out
    for pname, p in sig.parameters.items():
        if pname == "field":
            continue
        default = None if p.default is inspect.Parameter.empty else p.default
        # bool is an int subclass -- test it FIRST (amendment pt 2).
        if isinstance(default, bool):
            kind = "bool"
        elif isinstance(default, int):
            kind = "int"
        else:
            kind = "str"
        out.append({"name": pname, "kind": kind, "default": default})
    return out


@lru_cache(maxsize=1)
def build_catalog() -> dict:
    families = [{"name": name, "params": _params_of(builder),
                 "fields": ["CC", "GF"]}
                for name, builder in _iter_families()]
    return {
        "families": sorted(families, key=lambda f: f["name"]),
        "fields": {"CC": {"needs": []}, "GF": {"needs": ["p"], "optional": ["n"]}},
    }


@lru_cache(maxsize=1)
def family_names() -> frozenset[str]:
    return frozenset(f["name"] for f in build_catalog()["families"])


def validate_family(name: str, params: dict) -> None:
    if name not in family_names():
        raise CatalogError(f"unknown family {name!r}; see GET /api/catalog")
    known = {p["name"] for f in build_catalog()["families"]
             if f["name"] == name for p in f["params"]}
    unknown = set(params) - known
    if unknown:
        raise CatalogError(f"family {name!r} got unknown params {sorted(unknown)}")
