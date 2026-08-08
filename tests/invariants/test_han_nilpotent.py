"""Leg (i) of the Han bounded-extension recognizer -- tensor nilpotency of ``A/B``
over ``B`` (CLMS ``2101.02597`` Def. 2.3, Sec. 5). Plan 73 / R7, Task group I.

The two literature oracles: Ex. 5.3 (``A/B`` tensor-nilpotent, index 2, every
relative cycle ``J``-interrupted) and Ex. 5.5 (NOT nilpotent, witness ``d*c`` -- a
relative cycle with no ``J``-interrupter). Composition left-to-right.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import QQ
from quiverlab.families.extension import arrow_removal_subalgebra

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def build_ex53(field=QQ):
    return Quiver(
        [1, 2, 3, 4, 5],
        {"al": (5, 1), "be": (1, 5), "d": (1, 4), "c": (4, 3), "b": (3, 2), "a": (2, 1)},
    ).algebra(relations=["al*be", "d*c*b*a - be*al"], field=field)


def build_ex55(field=QQ):
    return Quiver(
        [1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3), "d": (3, 1)}
    ).algebra(relations=["a*b", "d*a", "b*d", "d*c*d"], field=field)


# --------------------------------------------------------------------------- #
# Task I1: relative cycles + J-interrupters (CLMS Def. 5.11 / 5.13)
# --------------------------------------------------------------------------- #
@lit
def test_ex53_relative_cycles_all_J_interrupted():
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    cycles = ext.relative_cycles(cap=8)
    assert cycles                                   # Ex 5.3 HAS relative cycles
    # every relative cycle carries a J-interrupter (the join collapses into B) ->
    # A/B is tensor-nilpotent (CLMS Thm 5.14/5.16)
    assert all(ext.has_J_interrupter(w) for w in cycles)


@lit
def test_ex55_dc_cycle_has_no_J_interrupter():
    ext = arrow_removal_subalgebra(build_ex55(), ("d",))
    cycles = ext.relative_cycles(cap=8)
    assert "d*c" in cycles                           # the cycle at vertex 3 through d
    # no J-interrupter -> the cycle survives every tensor power (CLMS Ex. 5.5)
    assert ext.has_J_interrupter("d*c") is False


@selfcert
def test_relative_cycles_are_vertex_cycles():
    ext = arrow_removal_subalgebra(build_ex55(), ("d",))
    Q = ext.A.quiver
    for w in ext.relative_cycles():
        toks = tuple(w.split("*"))
        assert Q.word_source(toks) == Q.word_target(toks)
