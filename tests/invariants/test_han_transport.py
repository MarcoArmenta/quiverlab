"""Han transport -- the injection/isomorphism LADDER (CLMS Thm 3.1/4.6, Plan 73 /
R7, Task III2). Ex. 5.3 bounded (iso), Ex. 5.5 not-bounded; the ladder LOGIC is
verified independently of any algebra.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import GF, QQ
from quiverlab.invariants.han import transport_verdict

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


@lit
def test_ex53_bounded_transport():
    ht = build_ex53(QQ).han_transport(("a",), hh_top=9)
    assert ht.transport == "bounded"                 # all three legs (a side) -> ISO
    assert ht.certificate.tensor_nilpotent.index == 2
    assert ht.certificate.pd_Be["route"] == "gldim"  # finite pd_{B^e} FREE via gl.dim B = 2
    assert ht.injection_from is not None             # certified LOWER BOUND (H3)
    assert ht.injection_bound_ok is True             # dim HH_m(B) <= dim HH_m(A); EQUALITY here
    assert "clms_bounded_extensions" in ht.references


@lit
def test_ex55_not_bounded():
    ht = build_ex55(QQ).han_transport(("d",), hh_top=9)
    assert ht.transport == "not_bounded"             # A/B not tensor-nilpotent (leg i False)


@xeng
def test_injection_bound_self_cert_gate():
    """The honest JZ consequence on the SHIPPED HH engines is the INJECTION BOUND
    dim HH_m(B) <= dim HH_m(A). Ex 5.3 is BOUNDED so it is an equality; the
    not-bounded Ex 5.5 sequences differ."""
    A53, B53 = build_ex53(GF(32003)), build_ex53(GF(32003)).han_transport(("a",)).certificate
    a53 = build_ex53(GF(32003)).hochschild_homology(9, verbose=False).dims
    from quiverlab.families.extension import arrow_removal_subalgebra
    b53 = arrow_removal_subalgebra(build_ex53(GF(32003)), ("a",)).B \
        .hochschild_homology(9, verbose=False).dims
    assert all(b <= a for a, b in zip(a53, b53))
    assert a53 == b53                                # EQUALITY only because 5.3 is bounded
    a55 = build_ex55(GF(32003)).hochschild_homology(9, verbose=False).dims
    b55 = arrow_removal_subalgebra(build_ex55(GF(32003)), ("d",)).B \
        .hochschild_homology(9, verbose=False).dims
    assert a55 != b55                                # not bounded -> differ


@selfcert
def test_injection_ladder_logic():
    """W3(b): the LADDER LOGIC itself (H1), independent of any concrete algebra."""
    assert transport_verdict(nilp=True, pd_finite=True, one_sided=True) == "bounded"
    assert transport_verdict(nilp=True, pd_finite=True, one_sided=False) == "pd_injection"
    assert transport_verdict(nilp=True, pd_finite=None, one_sided=None) == "nilpotent_injection"
    assert transport_verdict(nilp=False, pd_finite=None, one_sided=None) == "not_bounded"
    assert transport_verdict(nilp=None, pd_finite=None, one_sided=None) == "undecided"
    # there is NO iso-from-leg-(i): pd finite but one-sided undecided stays an injection
    assert transport_verdict(nilp=True, pd_finite=True, one_sided=None) == "pd_injection"


@lit
def test_ex53_bounded_certificate_legs():
    cert = build_ex53(QQ).bounded_extension(("a",))
    assert cert.bounded is True and cert.side == "right"
    assert cert.tensor_nilpotent.status == "nilpotent"
    assert cert.one_sided.projective is True
    assert cert.pd_Be["status"] == "finite"


@lit
def test_ex55_not_bounded_certificate():
    cert = build_ex55(QQ).bounded_extension(("d",))
    assert cert.bounded is False
    assert cert.tensor_nilpotent.status == "not_nilpotent"


# --------------------------------------------------------------------------- #
# Task III4: the gl.dim B = infinity fallback + the strict-injection case
# --------------------------------------------------------------------------- #
def _loop_ext(field=QQ):
    """A = k<x,a>/(x^2) on Q = {x:1->1 loop, a:1->2}, F = {a}. B = k[x]/(x^2) (+) k
    has gl.dim B = infinity (dual-numbers factor) and nonzero higher HH."""
    A = Quiver([1, 2], {"x": (1, 1), "a": (1, 2)}).algebra(relations=["x*x"], field=field)
    return A


@selfcert
def test_enveloping_fallback_exercised():
    """gl.dim B = infinity routes finite_pd_Be to the ENVELOPING fallback (the primary
    oracles never reach it): A/B tensor-nilpotent, pd_{B^e} finite via B^e, nonzero
    higher HH_*(B). The transport is 'bounded' (A/B here is one-sided projective)."""
    from quiverlab.families.extension import arrow_removal_subalgebra
    from quiverlab.invariants.han import finite_pd_Be
    A = _loop_ext(QQ)
    ext = arrow_removal_subalgebra(A, ("a",))
    assert ext.B.global_dimension().exact is False           # gl.dim B = infinity
    pd = finite_pd_Be(ext, pd_cap=12)
    assert pd["route"] == "enveloping" and pd["status"] == "finite"
    ht = A.han_transport(("a",), hh_top=8)
    assert ht.certificate.tensor_nilpotent.status == "nilpotent"
    assert ht.transport == "bounded"
    hhB = ext.B.hochschild_homology(8, engine="cs", verbose=False).dims
    assert hhB[4] > 0 and hhB[6] > 0                          # nonzero HIGHER HH_*(B)


@selfcert
def test_strict_injection_unexemplified_ladder_logic_verified():
    """H2/W3c honest scope: a STRICT-injection instance (tensor-nilpotent-but-UNBOUNDED,
    giving dim HH_m(B) < dim HH_m(A) under a theorem-guaranteed injection) is
    UNEXEMPLIFIED within the arrow-extension model on the small quivers surveyed at
    implementation: tensor-nilpotency forces every relative cycle to be J-interrupted,
    which (on one/two-vertex loop extensions of dual-numbers B) simultaneously makes
    A/B one-sided B-projective -> bounded (equality), not a proper injection. The
    'pd_injection'/'nilpotent_injection' rows and the STRICT injection are therefore
    LOGIC-verified (the ladder function), not instantiated by a concrete algebra."""
    assert transport_verdict(nilp=True, pd_finite=True, one_sided=False) == "pd_injection"
    assert transport_verdict(nilp=True, pd_finite=None, one_sided=False) == "nilpotent_injection"
    # the strict-injection case is theory-only; the ladder rows above encode it.
