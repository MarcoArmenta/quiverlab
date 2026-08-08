"""Plan 61 x Plan 56 interplay: on ada algebras over CC, the Theorem-B HH^1 verdict AGREES
with P56's is_simply_connected verdict wherever P56 is definite; where P56 returns None
(Adian-Rabin undecidable), Theorem B RESOLVES it. Auto-flips to a real assert when P56
lands."""
import pytest

pytest.importorskip("quiverlab.invariants.coverings")   # P56 alignment point
from quiverlab.invariants.coverings import is_simply_connected  # noqa: E402
from quiverlab import Quiver, RadicalSquareZero  # noqa: E402
from quiverlab.fields import CC  # noqa: E402
from quiverlab.modules.recognizers_ladder import recognizer_ladder  # noqa: E402

pytestmark = pytest.mark.oracle_crossengine


def _ada_cc_members():
    # two ada algebras over CC (both simply connected here): the tree fixture + the comm-square.
    Q5 = Quiver([1, 2, 3, 4, 5],
                {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    comm = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                  ).algebra(relations=["a*b - c*d"], field=CC)
    return [RadicalSquareZero(Q5, field=CC), comm]


@pytest.mark.parametrize("A", _ada_cc_members())
def test_theorem_B_agrees_with_p56_on_ada_cc_members(A):
    # ada/CC members: Theorem-B HH^1 verdict AGREES with P56 wherever P56 is definite;
    # where P56 returns None (Adian-Rabin), Theorem B is the complete resolution (headline).
    L = recognizer_ladder(A)
    assert L.verdict("ada") is True and L.simple_connectedness.applicable is True
    p56 = is_simply_connected(A).verdict
    if p56 is not None:                                     # P56 definite => must match
        assert p56 == L.simple_connectedness.verdict
    else:                                                   # P56 undecided => P61 resolves it
        assert L.simple_connectedness.verdict is not None
