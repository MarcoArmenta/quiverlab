"""Live QPA crosscheck: IsGentleAlgebra / IsSpecialBiserialAlgebra."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.invariants import recognizers as rec
from quiverlab.qpa import scripts, session

pytestmark = [pytest.mark.qpa,
              pytest.mark.skipif(session.should_skip_qpa(),
                                 reason="[qpa] backend not installed")]


CASES = [
    ("gentle_a3", lambda: Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))),
    ("square", lambda: Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4),
                                             "c": (1, 3), "d": (3, 4)}).algebra(
        relations=["a*b - c*d"], field=GF(5))),
    ("three_star", lambda: Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3),
                                                 "c": (1, 4)}).algebra(field=GF(5))),
]


@pytest.mark.parametrize("name,build", CASES, ids=[c[0] for c in CASES])
def test_gentle_and_sb_match_qpa(name, build):
    A = build()
    decl = scripts.quiver_and_algebra_script(A)      # binds GAP variable `A`
    ours_sb = rec.is_special_biserial(A)
    ours_g = rec.is_gentle(A)
    qpa_sb = bool(session.run(decl + "\nIsSpecialBiserialAlgebra(A);"))
    qpa_g = bool(session.run(decl + "\nIsGentleAlgebra(A);"))
    assert (ours_sb, ours_g) == (qpa_sb, qpa_g)
