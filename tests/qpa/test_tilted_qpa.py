"""QPA cross-check for the tilted-algebra recognizer (Plan 60 / R17).

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.

QPA 1.37 has NO tilted-algebra recognizer (no ``IsTiltedAlgebra`` / slice search) -- this test
DOCUMENTS that as a fail-if-appears trip-wire (Plan 60's honest-scope note goes stale the day QPA
ships one). What QPA CAN corroborate is the Ringel-slice HALF of the certificate: the slice module
``S = (+) Sigma`` returned on a ``tilted`` verdict is a TILTING module, which QPA's COMPUTATIONAL
``TiltingModule(S, 1)`` confirms (returns a ``[true, coresolution]`` list, not ``false``) -- the
P44 tilting crosscheck precedent (``tests/qpa/test_tilting_qpa.py``)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero
from quiverlab.fields import QQ
from quiverlab.modules.qpa_module import graded_form
from quiverlab.qpa import session
from quiverlab.qpa.scripts import module_decl, quiver_and_algebra_script

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _kA3_radsq():
    return RadicalSquareZero(Quiver([1, 2, 3], {"a1": (2, 1), "a2": (3, 2)}), field=QQ)


def _qpa_tilting(A, T, n, var):
    dv, arr = graded_form(T)
    scr = quiver_and_algebra_script(A) + "\n" + module_decl(A, dv, arr, var)
    scr += "\nt := TiltingModule(%s, %d);" % (var, n)
    return bool(session.run(scr))                     # false -> not tilting; list -> tilting


def test_qpa_has_no_tilted_recognizer():
    lg = session.libgap_handle()
    for name in ("IsTiltedAlgebra", "IsTilted", "TiltedAlgebra", "SliceInModuleCategory"):
        assert not bool(lg.eval('IsBoundGlobal("%s")' % name)), (
            "QPA now exposes %s -- wire a real crosscheck (this assert is the trip-wire "
            "that Plan 60's QPA scope note is stale)" % name)


def test_slice_module_is_tilting_via_qpa():
    A = _kA3_radsq()
    rep = A.tilted_check()
    assert rep.verdict == "tilted"
    # the slice module S = (+) Sigma is a tilting A-module (Ringel Thm 1.9(2), tilting half);
    # QPA's computational TiltingModule(S, 1) confirms it (non-false).
    assert _qpa_tilting(A, rep.slice_module, 1, "SS") is True
