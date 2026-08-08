"""The Jacobi-Zariski tool -- relative Hochschild homology HH_*(A|B) via the CLMS
normalized relative bar complex (Plan 73 / R7, Task III1). Computed over GF(p) for
speed (the relative complex term C_1 blows up over QQ).
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import GF, QQ
from quiverlab.families.extension import arrow_removal_subalgebra
from quiverlab.hochschild.jacobi_zariski import (
    jacobi_zariski_sequence, relative_homology)

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


def build_ex53(field):
    return Quiver(
        [1, 2, 3, 4, 5],
        {"al": (5, 1), "be": (1, 5), "d": (1, 4), "c": (4, 3), "b": (3, 2), "a": (2, 1)},
    ).algebra(relations=["al*be", "d*c*b*a - be*al"], field=field)


def build_ex55(field):
    return Quiver(
        [1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3), "d": (3, 1)}
    ).algebra(relations=["a*b", "d*a", "b*d", "d*c*d"], field=field)


@lit
def test_ex53_relative_homology_vanishes_above_index():
    """Cor. 2.4: HH_*(A|B) = 0 for * >= n (n = 2) on the nilpotent Ex. 5.3."""
    ext = arrow_removal_subalgebra(build_ex53(GF(32003)), ("a",))
    rh = relative_homology(ext, 6).dims
    assert rh == [5, 0, 0, 0, 0, 0, 0]                # HH_0 = 5, HH_{>=1} = 0


@selfcert
def test_hh0_relative_is_coinvariants():
    """HH_0(A|B, A) = dim of the B-coinvariants A_B = A (x)_{B^e} B."""
    ext = arrow_removal_subalgebra(build_ex53(GF(32003)), ("a",))
    rh = relative_homology(ext, 3).dims
    assert rh[0] == 5


@selfcert
def test_not_nilpotent_relative_homology_refuses():
    """Ex. 5.5 is not tensor-nilpotent => the relative complex is infinite => loud."""
    ext = arrow_removal_subalgebra(build_ex55(GF(32003)), ("d",))
    with pytest.raises(QuiverlabError):
        relative_homology(ext, 6)


@selfcert
def test_jz_exact_twice_in_three():
    """The JZ self-cert: exactness at the HH_*(B) spot (injection bound) and the
    HH_*(A|B) spot (Cor. 2.4 vanishing) both hold on Ex. 5.3."""
    ext = arrow_removal_subalgebra(build_ex53(GF(32003)), ("a",))
    jz = jacobi_zariski_sequence(ext, 6)
    assert all(jz.exact_at_B) and all(jz.exact_at_rel)
    assert jz.nilpotency_index == 2
    assert jz.hh_rel == [5, 0, 0, 0, 0, 0, 0]


@xeng
def test_jz_dims_from_shipped_engines():
    """The JZ sequence assembles HH_*(A), HH_*(B) from the shipped engines; on the
    bounded Ex. 5.3 they AGREE (iso), the injection bound is an equality."""
    ext = arrow_removal_subalgebra(build_ex53(GF(32003)), ("a",))
    jz = jacobi_zariski_sequence(ext, 6)
    assert jz.hh_A == jz.hh_B == [5, 0, 0, 0, 0, 0, 0]
    assert all(g == 0 for g in jz.gap_at_A)           # equality => zero gap here
