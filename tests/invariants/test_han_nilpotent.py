"""Leg (i) of the Han bounded-extension recognizer -- tensor nilpotency of ``A/B``
over ``B`` (CLMS ``2101.02597`` Def. 2.3, Sec. 5). Plan 73 / R7, Task group I.

The two literature oracles: Ex. 5.3 (``A/B`` tensor-nilpotent, index 2, every
relative cycle ``J``-interrupted) and Ex. 5.5 (NOT nilpotent, witness ``d*c`` -- a
relative cycle with no ``J``-interrupter). Composition left-to-right.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import GF, QQ
from quiverlab.families import truncated_polynomial
from quiverlab.families.extension import arrow_removal_subalgebra
from quiverlab.invariants.han import is_tensor_nilpotent

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
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


# --------------------------------------------------------------------------- #
# Task I2: is_tensor_nilpotent -- capped semi-decision + certificate + witness
# --------------------------------------------------------------------------- #
@lit
def test_ex53_tensor_nilpotent_index_2():
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    tn = is_tensor_nilpotent(ext, cap=8)
    assert tn.status == "nilpotent" and tn.index == 2
    assert tn.route in ("J_interrupter", "no_relative_cycles", "direct_cap")


@lit
def test_ex55_not_tensor_nilpotent_witness():
    ext = arrow_removal_subalgebra(build_ex55(), ("d",))
    tn = is_tensor_nilpotent(ext, cap=8)
    assert tn.status == "not_nilpotent"
    assert tn.witness is not None                    # the surviving relative cycle


@selfcert
def test_certificate_forces_vanishing():
    """J-interrupter route => (A/B)^{ox index} recomputed directly == 0."""
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    tn = is_tensor_nilpotent(ext, cap=8)
    assert ext.tensor_power_dim(tn.index) == 0 and ext.tensor_power_dim(tn.index - 1) != 0


@selfcert
def test_cap_returns_exactly_undecided_not_a_lie():
    """The cap-only route sees 3 nonzero powers and CANNOT conclude non-nilpotency
    (powers may grow before dying) -- it must return EXACTLY 'undecided'."""
    ext = arrow_removal_subalgebra(build_ex55(), ("d",))
    tn = is_tensor_nilpotent(ext, cap=3, use_certificate=False)
    assert tn.status == "undecided"
    assert tn.route == "direct_cap"


@selfcert
def test_cap_route_finds_nilpotent_when_a_power_vanishes():
    """The naive cap route (no certificate) still concludes 'nilpotent' on Ex 5.3."""
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    tn = is_tensor_nilpotent(ext, cap=8, use_certificate=False)
    assert tn.status == "nilpotent" and tn.index == 2 and tn.route == "direct_cap"


# --------------------------------------------------------------------------- #
# Task III3: CLMS Ex. 5.4 -- rad is E-tensor-nilpotent iff Q is acyclic
# (E = kQ_0 vertex subalgebra: new_arrows = ALL arrows, so B = kQ_0, A/B = rad A)
# --------------------------------------------------------------------------- #
@lit
def test_ex54_acyclic_rad_is_E_nilpotent():
    """Linear A_3 is ACYCLIC: rad is E-tensor-nilpotent (CLMS Ex. 5.4)."""
    A3 = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)
    ext = arrow_removal_subalgebra(A3, ("a", "b"))    # B = kQ_0
    assert ext.dim_B == 3                             # E = k^3 (the vertices)
    tn = is_tensor_nilpotent(ext, cap=8)
    assert tn.status == "nilpotent"


@lit
def test_ex54_cyclic_rad_is_not_E_nilpotent():
    """A loop (oriented cycle) is CYCLIC: rad is NOT E-tensor-nilpotent -- the loop
    survives every tensor power (CLMS Ex. 5.4)."""
    loop = truncated_polynomial(2, field=QQ)          # k[x]/(x^2), one loop x
    ext = arrow_removal_subalgebra(loop, ("x",))      # B = kQ_0 = k
    assert ext.dim_B == 1
    tn = is_tensor_nilpotent(ext, cap=8)
    assert tn.status == "not_nilpotent"


@xeng
def test_relative_homology_equals_absolute_over_kQ0():
    """E = kQ_0 is separable, so HH_*(A|E) = absolute HH_*(A) (cross-checks the CLMS
    relative bar complex against the shipped HH engine on acyclic A_3)."""
    from quiverlab.hochschild.jacobi_zariski import relative_homology
    A3 = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=GF(32003))
    ext = arrow_removal_subalgebra(A3, ("a", "b"))
    rh = relative_homology(ext, 4).dims
    abs_hh = A3.hochschild_homology(4, verbose=False).dims
    # HH_0(kA_3) = k^3 (the three vertices; no oriented cycles). Three independent
    # routes agree: the fast engine, relative.py's E-relative route, and my complex.
    assert rh == abs_hh == [3, 0, 0, 0, 0]
