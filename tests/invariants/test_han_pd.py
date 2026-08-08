"""Legs (ii)+(iii) of the Han bounded-extension recognizer -- finite ``pd_{B^e}(A/B)``
(CLMS Ex. 6.1 gl.dim route + enveloping fallback) and one-sided ``B``-projectivity
(CLMS Thm 5.20). Plan 73 / R7, Task group II.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import QQ
from quiverlab.families.extension import arrow_removal_subalgebra
from quiverlab.invariants.han import finite_pd_Be, one_sided_projective

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


# --------------------------------------------------------------------------- #
# Task II3: one-sided B-projectivity (leg iii, CLMS Thm 5.20)
# --------------------------------------------------------------------------- #
@lit
def test_ex53_one_sided_projective():
    """Ex. 5.3: A/B is one-sided B-projective (pd_B(A/B) = 0). In quiverlab's
    left-to-right convention the projective side is RIGHT (CLMS right-to-left
    'left'-bounded => right-projective here)."""
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    osp = one_sided_projective(ext, side="auto")
    assert osp.projective is True
    assert osp.side == "right"                        # convention-translated (see docstring)
    assert osp.pd_right == 0


@selfcert
def test_ex53_explicit_right_projective():
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    osp = one_sided_projective(ext, side="right")
    assert osp.projective is True and osp.pd_right == 0


@lit
def test_ex55_not_one_sided_projective():
    """Ex. 5.5 (not bounded): A/B is NOT projective on either side."""
    ext = arrow_removal_subalgebra(build_ex55(), ("d",))
    osp = one_sided_projective(ext, side="auto")
    assert osp.projective is False and osp.side is None
    assert osp.pd_left != 0 and osp.pd_right != 0
