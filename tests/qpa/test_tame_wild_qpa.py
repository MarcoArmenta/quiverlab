"""QPA (GAP) cross-check for the Tits-form FORM layer + honest scope on the
representation-type VERDICT layer (Plan 62 / R19).

Live discovery (QPA 1.37, verified 2026-08-08): QPA DOES ship a combinatorial
Tits-form + weak-positivity/nonnegativity surface --
``TitsUnitFormOfAlgebra(A)``, ``IsWeaklyPositiveUnitForm``,
``IsWeaklyNonnegativeUnitForm`` -- so the Plan-62 FORM layer gets a genuine
two-implementation cross-engine oracle (this file), NOT an honest skip. (The plan's
draft assumed no such surface; adjusting to reality per the plan's own "cross-engine
+ QPA agreement wherever QPA implements the feature" rule.)

QPA has NO representation-type (tame/wild) VERDICT verb (no IsTameAlgebra /
IsWildAlgebra / RepresentationType bound; verified via NamesGVars). So the Plan-62
VERDICT layer -- rep-finite / tame / wild gated on strong simple connectivity --
has no QPA oracle; its covering oracles are the literature/theory pins (Bongartz;
Brustle-de la Pena-Skowronski). ``test_qpa_has_no_representation_type_verb`` records
this honest scope and FAILS if QPA ever ships such a verb.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.qpa import crosscheck, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _star(arms, field):
    verts, ars, nxt = [0], {}, 1
    for ai, length in enumerate(arms):
        prev = 0
        for _ in range(length):
            verts.append(nxt)
            ars[f"e{ai}_{nxt}"] = (prev, nxt)
            prev = nxt
            nxt += 1
    return Quiver(verts, ars).algebra(field=field)


def _kron(m, field):
    return Quiver([1, 2], {f"a{k}": (1, 2) for k in range(m)}).algebra(field=field)


# The form is field-free; GF(p) is a field QPA supports exactly.
_F = GF(7)
_CASES = {
    "A5": linear_path_algebra(5, field=_F),          # Dynkin: wp True, wnn True
    "E8": _star([1, 2, 4], _F),                      # Dynkin: wp True, wnn True
    "~E8": _star([1, 2, 5], _F),                     # Euclidean: wp False, wnn True
    "T237": _star([1, 2, 6], _F),                    # wild (list-decided): wp/wnn False
    "K2": _kron(2, _F),                              # ~A1: wp False, wnn True
    "K3": _kron(3, _F),                              # wild: wp False, wnn False
    # a relation algebra (nonzero r_ij): kA3 with a*b=0
    "kA3rel": Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=_F),
}


@pytest.mark.parametrize("name", list(_CASES))
def test_weak_pos_nonneg_agree_with_qpa(name):
    # our weak positivity / weak nonnegativity of q_A == QPA's
    # IsWeaklyPositiveUnitForm / IsWeaklyNonnegativeUnitForm of TitsUnitFormOfAlgebra.
    crosscheck(_CASES[name], "tits_weak").assert_agree()


def test_T237_wnn_false_confirmed_by_qpa():
    # the T_{2,3,7} 10-variable case: our LIST decides weakly nonnegative = False; QPA's
    # INDEPENDENT IsWeaklyNonnegativeUnitForm confirms it (its sincere defect is outside
    # any box, so this is a real validation of the classified-list route).
    rep = crosscheck(_CASES["T237"], "tits_weak")
    assert rep.ours == [False, False] and rep.qpa == [False, False]


def test_qpa_has_no_representation_type_verb():
    # honest scope: QPA has the FORM surface but NO tame/wild VERDICT verb. FAIL if it
    # ever ships one (then wire a real verdict-layer crosscheck).
    lg = session.libgap_handle()
    names = ("IsTameAlgebra", "IsWildAlgebra", "RepresentationType",
             "IsOfFiniteRepresentationType", "RepresentationTypeOfAlgebra")
    bound = [n for n in names if bool(lg.eval(f'IsBoundGlobal("{n}")'))]
    assert not bound, (
        f"QPA now exposes a representation-type verdict verb {bound} -- wire a "
        "real crosscheck against the Plan-62 rep-finite/tame/wild verdicts.")
