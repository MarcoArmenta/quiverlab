"""Plan 08 freshness gate: refuse to build release infrastructure against a
half-built library. Prints a report; exits nonzero (STOP) on any drift.

Run:  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 \
      /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python scripts/release_freshness.py
"""
from __future__ import annotations

import importlib
import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent


def check() -> list[str]:
    """Return a list of PREREQUISITE drift messages (Plans 03-07 library surfaces +
    docs sources); empty list == fresh. The deprecated-license line is NOT included
    here -- it is an informational Task-3 TODO returned by license_todo() -- so a
    not-yet-SPDX license never forces a STOP."""
    drift: list[str] = []

    # --- prerequisite library surfaces (Plans 03-07) -----------------------
    import quiverlab
    from quiverlab import Quiver, CC

    # Plan 03: general (non-monomial) kQ/I lowers, groebner.system present.
    try:
        importlib.import_module("quiverlab.groebner.system")
        A = Quiver([1, 2, 3, 4],
                   {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                   ).algebra(relations=["a*b - c*d"], field=CC)
        assert A.dim == 9
    except Exception as e:  # noqa: BLE001
        drift.append(f"Plan 03 (general kQ/I / groebner.system) incomplete: {e!r}")

    # Plan 04: Chouhy-Solotar resolution package.
    if importlib.util.find_spec("quiverlab.resolutions_cs") is None:
        drift.append("Plan 04 (quiverlab.resolutions_cs) not present")

    # Plan 05: modules + Ext + gl.dim on the public Algebra.
    for name in ("ext", "simple", "projective", "global_dimension"):
        if not hasattr(quiverlab.Algebra, name):
            drift.append(f"Plan 05 surface Algebra.{name} missing")

    # Plan 06: families catalog + bibliography() + the PACKAGED citations registry.
    for name in ("zoo", "families", "bibliography"):
        if not hasattr(quiverlab, name):
            drift.append(f"Plan 06 surface quiverlab.{name}() missing")
    try:
        from quiverlab import citations
        if not citations.references_bib_path().is_file():
            drift.append("Plan 06 packaged citations.references_bib_path() does not "
                         "resolve to a file (JOSS paper.bib + docs + web /literature depend on it)")
    except Exception as e:  # noqa: BLE001
        drift.append(f"Plan 06 quiverlab.citations registry missing: {e!r}")

    # Marker policy: heavy test dirs must be in the deep bucket. Any NEW top-level
    # tests/ dir not covered by tests/conftest.py buckets is triaged loudly by
    # tests/release/test_markers.py::test_buckets_partition_the_suite; this note is
    # the standing reminder to add Plan 05/06/07 heavy dirs to _DEEP_DIRS.

    # Plan 07: viz + trace on the public Algebra.
    for name in ("draw", "tikz"):
        if not hasattr(quiverlab.Algebra, name):
            drift.append(f"Plan 07 surface Algebra.{name} missing")

    # --- docs sources this plan will publish -------------------------------
    tut = sorted((ROOT / "docs" / "tutorials").glob("*.ipynb"))
    if len(tut) < 3:
        drift.append(f"expected >=3 tutorial notebooks, found {len(tut)}")
    internals = sorted((ROOT / "docs" / "internals").glob("[0-9][0-9]-*.md"))
    if len(internals) < 7:
        drift.append(f"expected >=7 internals chapters, found {len(internals)}")

    # The deprecated license={text=...} table is NOT prerequisite drift -- it is a
    # Task-3 fix-it TODO, reported separately by license_todo() as an INFORMATIONAL
    # note and never a STOP. Keeping it out of `drift` is what lets this gate exit 0
    # on current `main` (where the license is not SPDX-fixed until Task 3 runs), so
    # Step 2's "exit=0" and Task 13's "release_freshness.py exits 0" both hold.
    return drift


def license_todo() -> list[str]:
    """Informational only (never a STOP): the deprecated PEP-639 license={text=...}
    table that Task 3 replaces with the SPDX `license = "MIT"` string. main() prints
    these notes but they do NOT affect the exit code."""
    notes: list[str] = []
    pp = (ROOT / "pyproject.toml").read_text()
    if 'license = { text = "MIT" }' in pp or 'license = {text = "MIT"}' in pp:
        notes.append("pyproject still uses the deprecated license={text=...} table "
                     "(Task 3 replaces it with the PEP 639 SPDX string)")
    return notes


def main() -> int:
    drift = check()
    for note in license_todo():
        print(f"PLAN 08 FRESHNESS GATE: NOTE (informational, not a STOP) -- {note}")
    if drift:
        print("PLAN 08 FRESHNESS GATE: STOP -- prerequisites drifted:\n")
        for d in drift:
            print(f"  - {d}")
        print("\nComplete the named prerequisite plan(s) before running Plan 08.")
        return 1
    print("PLAN 08 FRESHNESS GATE: OK -- library surfaces, docs sources, and "
          "pyproject baseline are as expected.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
