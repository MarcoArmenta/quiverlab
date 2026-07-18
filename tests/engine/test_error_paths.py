"""Item 7 (TODO P3): the negative / error branches that the rest of the pipeline relies on.

The scan drivers depend on `check_associative` to reject malformed structure constants,
and `ChouhySolotarResolution` depends on `_classify` to refuse algebras whose CS closed-form
differentials were never derived. Both guards are only ever hit by the happy path elsewhere;
here we exercise the rejection branches directly.
"""
import numpy as np

from quiverlab.engine.hh_engine import check_associative, Algebra
from quiverlab.engine.scan3 import quantum_ci

# not ported: test_cs_classify_rejects_unsupported_reduction_system,
# test_cs_classify_accepts_supported_family -- both construct ChouhySolotarResolution from
# resolutions_cs, which is EXCLUDED from Plan 02. The check_associative guard is kept.


def test_check_associative_flags_nonassociative_with_triple():
    # unit e_0; x*x = y, x*y = x, y*x = y, y*y = 0:
    #   (x*x)*y = y*y = 0   but   x*(x*y) = x*x = y   -> non-associative at (1,1,1).
    m = 3
    T = np.zeros((m, m, m), dtype=np.int64)
    for j in range(m):
        T[0, j, j] = 1
        T[j, 0, j] = 1
    T[1, 1, 2] = 1   # x*x = y
    T[1, 2, 1] = 1   # x*y = x
    T[2, 1, 2] = 1   # y*x = y
    A = Algebra(m, T, np.array([1, 0, 0], dtype=np.int64), name="nonassoc")
    ok, bad = check_associative(A)
    assert ok is False
    assert bad == (1, 1, 1)


def test_check_associative_accepts_known_good():
    ok, bad = check_associative(quantum_ci(2))
    assert ok is True and bad is None
