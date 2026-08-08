"""Legs (ii)+(iii) of the Han bounded-extension recognizer -- finite ``pd_{B^e}(A/B)``
(CLMS Ex. 6.1 gl.dim route + enveloping fallback) and one-sided ``B``-projectivity
(CLMS Thm 5.20). Plan 73 / R7, Task group II.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import QQ
from quiverlab.families.extension import arrow_removal_subalgebra
from quiverlab.invariants.han import finite_pd_Be

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
# Task II1: the gl.dim B < inf PRIMARY route (CLMS Ex. 6.1) -- free, no B^e
# --------------------------------------------------------------------------- #
@lit
def test_ex53_pd_Be_finite_via_gldim():
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    pd = finite_pd_Be(ext)
    assert pd["route"] == "gldim" and pd["status"] == "finite"
    assert pd["gldim_B"] == 2                        # gl.dim B = 2 (live-verified)


@lit
def test_ex55_pd_Be_finite_via_gldim():
    ext = arrow_removal_subalgebra(build_ex55(), ("d",))
    pd = finite_pd_Be(ext)
    assert pd["route"] == "gldim" and pd["status"] == "finite"
    assert pd["gldim_B"] == 2
