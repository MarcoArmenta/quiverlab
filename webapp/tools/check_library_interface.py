"""Fail-loud verification that the released quiverlab surface matches what
webapp/ assumes. Run at the start of Plan 09; STOP if it reports drift.

The reference algebra is built the SAME way the runner builds it -- by looking
up a builder in families() and calling it by keyword -- so this gate fails
exactly where the runner would, not on a parallel code path.

Written against the REAL installed surface (see the 2026-07-24 plan amendment):
families() yields FamilyInfo records (the builder is a top-level export named
info.name); the reference family is QuantumCI (truncated_polynomial is not in
families()); the algebra exposes A.dim as an attribute; HH tables carry
.dims/.kind/.engine/.references but NOT .latex(); citation pairs come from
quiverlab.trace.provenance.resolve_references; the bibliography is a
Bibliography object with .groups/.keys/.bibtex()."""
from __future__ import annotations

import importlib.util

# Version floor this plan was written against (the installed 0.1.0.dev0 library).
# Bump when the plan is re-synced to a newer library; a lower installed version
# STOPS the plan.
MIN_QUIVERLAB_VERSION = "0.1.0"


def _version_tuple(v: str) -> tuple:
    parts = []
    for chunk in v.split(".")[:3]:
        num = "".join(ch for ch in chunk if ch.isdigit())
        parts.append(int(num) if num else 0)
    while len(parts) < 3:
        parts.append(0)
    return tuple(parts)


def _families_map(ql) -> dict:
    """Mirror of the catalog's family resolution: families() yields FamilyInfo
    records whose builder is the top-level export named info.name."""
    out = {}
    for info in ql.families():
        out[info.name] = getattr(ql, info.name, None)
    return out


def check() -> list[str]:
    problems: list[str] = []
    try:
        import quiverlab as ql
    except Exception as exc:  # pragma: no cover - import failure is the report
        return [f"cannot import quiverlab: {exc!r}"]

    for name in ("families", "zoo", "GF", "CC", "__version__", "bibliography",
                 "verbose"):
        if not hasattr(ql, name):
            problems.append(f"quiverlab.{name} missing")

    # Bibliography subsystem (spec 3.9): the /literature page renders .groups,
    # the download button serves .bibtex().
    if hasattr(ql, "bibliography"):
        try:
            bib = ql.bibliography()
            if not getattr(bib, "groups", None):
                problems.append("quiverlab.bibliography().groups is empty")
            if not hasattr(bib, "keys"):
                problems.append("quiverlab.bibliography() lacks .keys")
            if not hasattr(bib, "bibtex"):
                problems.append("quiverlab.bibliography() lacks .bibtex()")
            elif not bib.bibtex():
                problems.append("quiverlab.bibliography().bibtex() returned nothing")
        except Exception as exc:
            problems.append(f"quiverlab.bibliography() raised {exc!r}")

    # Version floor: STOP if the installed library predates this plan.
    ver = getattr(ql, "__version__", "0")
    if _version_tuple(ver) < _version_tuple(MIN_QUIVERLAB_VERSION):
        problems.append(
            f"quiverlab {ver} < required floor {MIN_QUIVERLAB_VERSION}; the plan "
            "was written against a newer library -- re-sync before executing")

    # The [fast] extra must be present so prod runs the fast GF(p) engine path.
    if importlib.util.find_spec("numba") is None:
        problems.append("numba not importable -- install the [fast] extra")

    # families() must enumerate the FamilyInfo records the catalog resolves by name.
    fam_map: dict = {}
    try:
        fam_map = _families_map(ql)
        if not fam_map:
            problems.append("quiverlab.families() returned nothing")
    except Exception as exc:
        problems.append(f"quiverlab.families() raised {exc!r}")

    # Build the reference algebra exactly as the runner will: by keyword.
    # QuantumCI(q=1) stays nondegenerate in char 2 (q=2 would degenerate the relation).
    A = None
    builder = fam_map.get("QuantumCI")
    if builder is None:
        problems.append("families() lacks 'QuantumCI' (runner reference family)")
    else:
        try:
            A = builder(q=1, field=ql.GF(2))
        except Exception as exc:
            problems.append(f"QuantumCI(q=1, field=GF(2)) raised {exc!r}")

    if A is not None:
        for attr in ("dim", "hochschild_cohomology", "hochschild_homology",
                     "coxeter_polynomial", "cartan_matrix", "global_dimension",
                     "center", "citations"):
            if not hasattr(A, attr):
                problems.append(f"algebra.{attr} missing")
        # verbose=False is REQUIRED: the global verbose flag defaults to True and
        # HH calls would otherwise emit worked-steps trace PDFs into quiverlab_traces/.
        for hh in ("hochschild_cohomology", "hochschild_homology"):
            if not hasattr(A, hh):
                continue
            try:
                table = getattr(A, hh)(3, verbose=False)
                for field in ("dims", "kind", "engine", "references"):
                    if not hasattr(table, field):
                        problems.append(f"{hh} table lacks .{field}")
                if hasattr(table, "dims") and not isinstance(table.dims, (list, tuple)):
                    problems.append(f"{hh} table .dims is not a sequence")
                if hasattr(table, "references") and not isinstance(
                        table.references, (list, tuple)):
                    problems.append(f"{hh} table .references is not a sequence")
            except Exception as exc:
                problems.append(f"{hh}(3, verbose=False) raised {exc!r}")

    # Citation-pair helper the worker uses to turn keys into (key, formatted) pairs.
    try:
        from quiverlab.trace.provenance import resolve_references
        pairs = resolve_references(("bar",))
        if not pairs:
            problems.append("resolve_references(('bar',)) returned nothing")
        elif not all(isinstance(p, (list, tuple)) and len(p) == 2 for p in pairs):
            problems.append("resolve_references(('bar',)) is not a sequence of pairs")
    except Exception as exc:
        problems.append(
            f"cannot use quiverlab.trace.provenance.resolve_references: {exc!r}")

    # Error hierarchy the app surfaces verbatim.
    try:
        from quiverlab import errors as qerr
        for err in ("QuiverlabError", "ExactnessError", "RelationError",
                    "DepthLimitError", "FieldError"):
            if not hasattr(qerr, err):
                problems.append(f"quiverlab.errors.{err} missing")
    except Exception as exc:
        problems.append(f"cannot import quiverlab.errors: {exc!r}")

    return problems
