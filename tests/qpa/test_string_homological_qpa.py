"""QPA scope for R34 (Plan 59). QPA recognizes special-biserial/gentle but has NO
homological string test (no extension/middle-term verb). Crosschecks: a DIRECT-session
IsSpecialBiserialAlgebra parity call (crosscheck.py has no recognizer verb) + a standing
guard that FAILS if such a surface ever appears."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_special_biserial
from quiverlab.qpa import scripts, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_homological_string_surface():
    lg = session.libgap_handle()
    for name in ("MiddleTermsOfExtensions", "HomologicalStringTest"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), \
            f"QPA now ships {name} -- add a real crosscheck (honest scope changed)"


@pytest.mark.parametrize("A, expected", [
    # gentle kA3/(ab): special-biserial + string; both should say True.
    (Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ),
     True),
    # kD4 (subspace): a valency-3 vertex -> NOT special-biserial (the R34 non-string side).
    (Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}).algebra(
        relations=[], field=QQ), False),
])
def test_special_biserial_parity_direct_session(A, expected):
    session.require_gap()
    lg = session.libgap_handle()
    if not bool(lg.eval('IsBoundGlobal("IsSpecialBiserialAlgebra")')):
        pytest.skip("this QPA version has no IsSpecialBiserialAlgebra")
    base = scripts.quiver_and_algebra_script(A)            # binds `A` in the GAP session
    qpa = bool(session.run(base + "\nIsSpecialBiserialAlgebra(A);"))
    assert qpa == is_special_biserial(A) == expected
