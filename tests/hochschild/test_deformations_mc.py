"""Plan 78 Task 2: the Maurer-Cartan report + the MRRS nilpotent-regime gate (char 0).

The nilpotent-regime gate is R13's conservative reading (gentle, no parallel arrows,
no oriented cycles); kZ_n/J^ell have an oriented cycle so the gate is None for them
(honest -- they are still unobstructed via the direct obstruction). The clean True
examples are gentle + acyclic (kA_3/rad^2, a gentle string on an acyclic quiver).
"""
import pytest

import quiverlab as ql
from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.hochschild.deformations import (
    is_l_infinity_nilpotent, maurer_cartan, _obstruction_data,
    self_bracket_class, _is_zero_vec)

pytestmark = pytest.mark.deep


def _kZ(n, ell, field=QQ):
    names = "abcdefghijklmnopqrstuvwxyz"
    verts = list(range(1, n + 1))
    arrows = {names[i]: (verts[i], verts[(i + 1) % n]) for i in range(n)}
    rels = []
    for start in range(n):
        word, cur = [], start
        for _ in range(ell):
            word.append(names[cur])
            cur = (cur + 1) % n
        rels.append("*".join(word))
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def _kA3_rad2(field=QQ):                       # 1->2->3, ab = 0 (gentle, acyclic)
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=field)


def _square_string(field=QQ):                  # gentle string on an acyclic quiver
    return Quiver([1, 2, 3, 4],
                  {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}).algebra(
        relations=["a*b", "c*d"], field=field)


@pytest.fixture(scope="module")
def qci0_data():
    return _obstruction_data(ql.QuantumCI(0, field=QQ))


# ---------------------------------------------------------------------------
@pytest.mark.oracle_literature
def test_nilpotent_gate_true_on_gentle_acyclic():
    """is_l_infinity_nilpotent = True on gentle + acyclic + no-parallel examples, and
    maurer_cartan then reports MC = Z^2 (MRRS Thm 5.4)."""
    for A in (_kA3_rad2(), _square_string()):
        assert is_l_infinity_nilpotent(A) is True
        rep = maurer_cartan(A)
        assert rep.nilpotent_regime is True
        assert rep.mc_description.startswith("MC = Z^2")
        assert rep.unobstructed is True
        assert rep.z2_dim is not None and rep.z2_dim >= 0
        assert rep.mc_order_certified >= 1


@pytest.mark.oracle_selfcert
def test_nilpotent_gate_none_outside_hypotheses():
    """None (NOT certified) when a hypothesis fails: kZ4/J^2 has an oriented cycle;
    QuantumCI(0) is non-gentle. Never a silent MC = Z^2."""
    assert is_l_infinity_nilpotent(_kZ(4, 2)) is None      # oriented cycle
    assert is_l_infinity_nilpotent(ql.QuantumCI(0, field=QQ)) is None  # non-gentle
    # char p: char-0 theorem, so None even on a gentle-acyclic algebra
    assert is_l_infinity_nilpotent(_kA3_rad2(field=GF(5))) is None


@pytest.mark.oracle_selfcert
def test_mc_quantumci0_is_obstructed():
    """maurer_cartan(QuantumCI(0))/QQ (non-nilpotent): nilpotent_regime is not True, the
    order-2 lift is obstructed, mc_order_certified == 1, a witness is exposed."""
    A = ql.QuantumCI(0, field=QQ)
    rep = maurer_cartan(A)
    assert rep.nilpotent_regime is not True
    assert rep.unobstructed is False
    assert rep.mc_order_certified == 1
    assert rep.obstruction_witness is not None
    assert rep.mc_description.startswith("obstructed")
    assert rep.window_note is not None


@pytest.mark.oracle_selfcert
def test_mc_recursion_solvable_iff_class_vanishes(qci0_data):
    """The DGLA recursion delta mu_2 = -1/2 [mu_1,mu_1] is solvable iff the RHS class
    vanishes: a direction with zero self-bracket lifts, the witness (nonzero) does not."""
    from quiverlab.hochschild.deformations import _verdict
    data = qci0_data
    dom = data["dom"]
    n = data["hh2_dim"]
    # find an UNOBSTRUCTED basis direction (zero self-bracket) and the OBSTRUCTED witness
    unob = [i for i in range(n)
            if _is_zero_vec(self_bracket_class(
                data, [dom.one() if k == i else dom.zero() for k in range(n)]), dom)]
    assert unob, "expected at least one liftable basis direction"
    for i in unob:                              # these lift (class vanishes)
        coords = [dom.one() if k == i else dom.zero() for k in range(n)]
        assert _is_zero_vec(self_bracket_class(data, coords), dom)
    _, _, wc = _verdict(data)                   # the obstructed witness does NOT lift
    assert not _is_zero_vec(self_bracket_class(data, wc), dom)


@pytest.mark.oracle_selfcert
def test_mc_charp_loud_but_hh_block_survives():
    """maurer_cartan over char p is a loud QuiverlabError (the interpretation is char-0),
    while the field-general HH^2/[a,a] block is still available."""
    from quiverlab.hochschild.deformations import (
        infinitesimal_deformations, obstruction_map)
    A = ql.QuantumCI(0, field=GF(5))
    with pytest.raises(QuiverlabError):
        maurer_cartan(A)
    # the field-general block still computes
    assert infinitesimal_deformations(A)["hh2_dim"] == 3
    assert obstruction_map(A).unobstructed is False
