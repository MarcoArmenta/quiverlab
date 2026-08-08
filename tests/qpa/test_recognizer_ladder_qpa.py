"""Plan 61 / R18 QPA cross-oracle (probe-first, honest fallback).

QPA 1.37 has NO quasi-tilted / shod / weakly-shod / laura / ada recognizer -- a live
NamesGVars() + IsBoundGlobal() sweep proves it and FAILS if a future QPA ever grows one
(the Plan-35 fail-if-appears pattern, so this honest-scope note is revisited). What QPA
DOES expose corroborates the POINTWISE ingredients behind the ladder: GlobalDimensionOfAlgebra
(the gl.dim theorem gates shod => gl.dim <= 3, quasi-tilted => gl.dim <= 2) and
InjDimensionOfModule (the pd/id data behind the parts -- e.g. S_3's id = 2 is exactly why
the rad^2=0 A5 is not shod). Recognizer verdicts run over QQ (M1 char-0-decisive
identification); the gl.dim/id crosschecks run over GF(5) on the SAME combinatorial algebra."""
import pytest

from quiverlab import GF, Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.recognizers_ladder import recognizer_ladder
from quiverlab.qpa import session
from quiverlab.qpa.crosscheck import (crosscheck_global_dimension,
                                      crosscheck_inj_dimension)

pytestmark = [pytest.mark.qpa,
              pytest.mark.skipif(session.should_skip_qpa(),
                                 reason="[qpa] backend not installed")]

# The names QPA might plausibly use for a recognizer of any ladder rung.
_RECOGNIZER_NAMES = ("IsQuasiTiltedAlgebra", "IsShodAlgebra", "IsWeaklyShodAlgebra",
                     "IsLauraAlgebra", "IsAdaAlgebra")


def _radsq_a5(field):
    return RadicalSquareZero(
        Quiver([1, 2, 3, 4, 5], {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)}),
        field=field)


def test_qpa_exposes_no_recognizer_ladder_surface():
    lg = session.libgap_handle()
    # (1) Scan the global name table FIRST -- IsBoundGlobal("Foo") REGISTERS "Foo" into
    #     NamesGVars(), so a scan after the queries would echo them back (Plan-35 note).
    #     "ada" is too generic a substring to scan; the distinctive rung names suffice.
    gvar_names = [str(n).lower() for n in lg.eval("NamesGVars()")]
    ladder_like = sorted({n for n in gvar_names
                          if "shod" in n or "laura" in n or "quasitilt" in n})
    # (2) None of the named recognizer entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _RECOGNIZER_NAMES}
    present = [name for name, ok in bound.items() if ok]
    if not present and not ladder_like:
        pytest.skip(
            "QPA 1.37 exposes no quasi-tilted/shod/weakly-shod/laura/ada recognizer "
            "(no Is{QuasiTilted,Shod,WeaklyShod,Laura,Ada}Algebra; NamesGVars has no "
            "shod/laura/quasitilt name). Covering oracle: the Plan-61 literature + "
            "self-cert batteries in tests/modules/; the gl.dim/id ingredients ARE "
            "QPA-checked below.")
    pytest.fail(
        "QPA now exposes a recognizer-ladder surface -- wire a real crosscheck against "
        f"the Plan-61 verdicts. Found: present={present}, ladder_like={ladder_like}.")


@pytest.mark.parametrize("name,build,gld", [
    ("kA3", lambda f: linear_path_algebra(3, field=f), 1),
    ("radsq_a5", _radsq_a5, 4),
])
def test_gldim_gate_corroborated_by_qpa(name, build, gld):
    # The gl.dim theorem gates are QPA-checkable bounds: quiverlab's gl.dim (used inside the
    # ladder over QQ) equals QPA's GlobalDimensionOfAlgebra on the same algebra over GF(5).
    crosscheck_global_dimension(build(GF(5))).assert_agree()
    L = recognizer_ladder(build(QQ))
    assert L.is_complete and L.gldim == gld
    if L.verdict("shod"):
        assert L.gldim <= 3            # Coelho-Lanzilotta
    if L.verdict("quasi_tilted"):
        assert L.gldim <= 2            # HRS QT1


def test_pointwise_id_behind_shod_refusal():
    # The rad^2=0 A5 is NOT shod because S_3 has both pd > 1 and id > 1. QPA corroborates the
    # POINTWISE id datum: InjDimensionOfModule(S_3) = 2 (> 1). This is the ingredient the
    # complement/parts computation reads; QPA has no shod verdict, but it agrees on id.
    A_gf = _radsq_a5(GF(5))
    crosscheck_inj_dimension(A_gf, A_gf.simple(3), bound=6).assert_agree()
    L = recognizer_ladder(_radsq_a5(QQ))
    assert L.verdict("shod") is False and L.rungs["shod"].witness["module"] == "S_3"
    assert L.rungs["shod"].witness["id_le_1"] is False
