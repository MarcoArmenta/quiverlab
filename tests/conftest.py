"""Global test fixtures + CI-bucket auto-assignment.

Plan 07: the worked-steps trace subsystem is ON by default (quiverlab.verbose =
True, spec D9); force it OFF for the whole suite so that unrelated computations do
not write ./quiverlab_traces/ files. Trace tests that need it opt back in explicitly
(set quiverlab.verbose = True or pass verbose=True / trace=[...]).

Plan 08 Task 2: auto-assign CI buckets by test location so the suite can be split
across the GitHub Actions matrix. Explicit bucket markers on a test win.
Buckets (exactly one per test; disjoint + exhaustive, enforced by the partition test):
  qpa  -- needs the [qpa] extra (passagemath-gap + QPA); CI QPA job only.
  deep -- heavy engine / resolution / module / families / batch suites; one Linux leg.
  fast -- everything else; runs on every OS x Python matrix cell.
Orthogonal sub-tag:
  slow -- an individually long test; IMPLIES deep (never runs in the fast matrix).
"""
import pytest

import quiverlab


# --- Plan 07 (PRESERVE VERBATIM) -------------------------------------------
def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "verbose_default: observe the shipped quiverlab.verbose default "
        "(opt out of the quiet-suite fixture)",
    )


@pytest.fixture(autouse=True)
def _quiet_traces(request):
    # Always capture + restore, so even a `verbose_default`-marked test that mutates
    # quiverlab.verbose cannot leak state into the next test (review hardening).
    prev = getattr(quiverlab, "verbose", True)
    if request.node.get_closest_marker("verbose_default") is None:
        # Force the trace subsystem OFF so unrelated tests don't write ./quiverlab_traces/.
        quiverlab.verbose = False
    # else: leave the shipped default (quiverlab.verbose = True) observable.
    try:
        yield
    finally:
        quiverlab.verbose = prev


# --- Plan 08 Task 2 (ADD) --------------------------------------------------
# Top-level dirs (relative to tests/) whose tests are heavy -> the deep leg.
# (resolutions_cs is the Chouhy-Solotar suite; there is NO tests/chouhy_solotar dir.)
_DEEP_DIRS = ("engine", "resolutions_cs", "modules", "families", "batch")
# Individually heavy files that may live outside the deep dirs.
_DEEP_FILES = ("test_complete.py", "test_deepen.py", "test_properties.py",
               "test_acceptance.py", "test_cs_", "test_bardzell", "test_minimal")
_BUCKETS = {"fast", "deep", "qpa"}


def _bucket(nodeid: str) -> str:
    # nodeid looks like "tests/engine/test_foo.py::test_bar"
    parts = nodeid.replace("\\", "/").split("/")
    top = parts[1] if len(parts) > 1 else ""
    fname = parts[-1].split("::")[0]
    if top == "qpa":
        return "qpa"
    if top in _DEEP_DIRS or any(fname.startswith(f) or fname == f for f in _DEEP_FILES):
        return "deep"
    return "fast"


def pytest_collection_modifyitems(config, items):
    for item in items:
        names = {m.name for m in item.iter_markers()}
        # `slow` always rides the deep leg (unless it is a qpa test).
        if "slow" in names and not (names & {"deep", "qpa"}):
            item.add_marker(pytest.mark.deep)
            names.add("deep")
        if names & _BUCKETS:          # explicit bucket (fast/deep/qpa) wins -> no auto-mark
            continue
        item.add_marker(getattr(pytest.mark, _bucket(item.nodeid)))
