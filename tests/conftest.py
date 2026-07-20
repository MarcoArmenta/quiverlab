"""Global test fixtures. The worked-steps trace subsystem is ON by default
(quiverlab.verbose = True, spec D9); force it OFF for the whole suite so that
unrelated computations do not write ./quiverlab_traces/ files. Trace tests that
need it opt back in explicitly (set quiverlab.verbose = True or pass
verbose=True / trace=[...])."""
import pytest

import quiverlab


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
