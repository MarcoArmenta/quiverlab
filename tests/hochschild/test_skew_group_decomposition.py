"""Plan 74 Task 4: stefan_decomposition assembly + the DIRECT cross-check +
orbit / free-orbit / non-abelian / modular oracles."""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families.basic import truncated_polynomial
from quiverlab.families.skew_group import (
    GroupAction, QuiverAutomorphism, is_free_action, skew_group_algebra)
from quiverlab.fields import GF, QQ
from quiverlab.hochschild.skew_group import stefan_decomposition

pytestmark = [pytest.mark.oracle_crossengine, pytest.mark.oracle_literature,
              pytest.mark.oracle_selfcert]


def _dual_z2():
    A = truncated_polynomial(2, field=QQ)
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": -1})
    return A, GroupAction.cyclic(sig, 2)


def _two_cycle_nakayama(field=QQ):
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2), "b": (2, 1)})
    return Q.algebra(relations=["a*b", "b*a"], field=field)


# ------------------------------------------------------------------ R8 oracle
@pytest.mark.oracle_crossengine
def test_r8_oracle_dual_numbers_both_sides():
    A, ga = _dual_z2()
    coh = stefan_decomposition(A, ga, 3, side="coh", verify_direct=True)
    assert coh.dims == [1, 1, 1, 1]
    assert coh.direct_dims == [1, 1, 1, 1]
    assert coh.agrees is True
    assert coh.status == "complete"
    hom = stefan_decomposition(A, ga, 3, side="hom", verify_direct=True)
    assert hom.dims == [2, 1, 1, 1] == hom.direct_dims
    assert hom.agrees is True
    # per-summand split pinned (findings 1a/2b): NOT [1, 0] at n=1
    assert hom.summands[0]["inv"] == [1, 0, 0, 0]
    assert hom.summands[1]["inv"] == [1, 1, 1, 1]


# ------------------------------------------------------- free-orbit oracle 2d
@pytest.mark.oracle_crossengine
def test_free_swap_orbit_is_truncated_polynomial():
    A = _two_cycle_nakayama(QQ)
    swap = QuiverAutomorphism({1: 2, 2: 1}, {"a": "b", "b": "a"})
    ga = GroupAction.cyclic(swap, 2)
    assert is_free_action(A, ga) is True
    S = skew_group_algebra(A, ga)
    assert S.dim == 8
    pf = S.presented_form()
    assert sorted(pf.quiver.vertices) == [1]
    assert len(pf.quiver.arrows) == 1
    loop = next(iter(pf.quiver.arrows))
    assert sorted(str(r) for r in pf.relations) == [f"{loop}*{loop}"]
    assert pf.is_selfinjective() is True
    ref = truncated_polynomial(2, field=QQ)
    assert S.hochschild_cohomology(2).dims == ref.hochschild_cohomology(2).dims == [2, 1, 1]


# ------------------------------------------------------- Z/3 non-involution 2a
@pytest.mark.oracle_crossengine
def test_z3_on_cubic_gf7():
    A = truncated_polynomial(3, field=GF(7))
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": 2})   # 2^3 = 1 mod 7
    ga = GroupAction.cyclic(sig, 3)
    coh = stefan_decomposition(A, ga, 2, side="coh", verify_direct=True)
    assert coh.dims == [1, 1, 1] == coh.direct_dims
    assert coh.agrees is True
    hom = stefan_decomposition(A, ga, 2, side="hom", verify_direct=True)
    assert hom.dims == [3, 2, 2] == hom.direct_dims
    # HH_0 splits 1 + 1 + 1 across the three summands
    assert [s["inv"][0] for s in hom.summands] == [1, 1, 1]


# ------------------------------------------------- non-abelian S3 on k^3 (2c)
@pytest.mark.oracle_literature
def test_s3_on_k3_non_abelian():
    A = Quiver(vertices=[1, 2, 3], arrows={}).algebra(relations=[], field=GF(7))
    transp = QuiverAutomorphism({1: 2, 2: 1, 3: 3}, {})
    cyc = QuiverAutomorphism({1: 2, 2: 3, 3: 1}, {})
    ga = GroupAction.from_generators([transp, cyc], A)
    assert ga.order == 6
    S = skew_group_algebra(A, ga)
    assert S.dim == 18                              # dim law
    rep = stefan_decomposition(A, ga, 0, side="coh", verify_direct=True)
    assert rep.dims == [2]                          # HH^0 = 2 (Morita k x k)
    assert rep.direct_dims == [2]
    assert rep.agrees is True
    # three conjugacy classes with proper centralizers Z(g) = S3 / Z2 / Z3
    zsizes = sorted(len(ga.centralizer(c[0])) for c in ga.conjugacy_classes)
    assert zsizes == [2, 3, 6]


# ----------------------------------- identity-summand monomorphism (selfcert)
@pytest.mark.oracle_selfcert
def test_identity_summand_monomorphism():
    A, ga = _dual_z2()
    for side in ("coh", "hom"):
        rep = stefan_decomposition(A, ga, 3, side=side)
        ident = rep.summands[0]              # the g=e summand = HH^n(A)^G
        for n in range(len(rep.dims)):
            assert ident["inv"][n] <= rep.dims[n]


# ------------------------------------------------------ modular loud refusal
@pytest.mark.oracle_selfcert
def test_modular_refusal_char_divides_order():
    A = _two_cycle_nakayama(GF(2))
    swap = QuiverAutomorphism({1: 2, 2: 1}, {"a": "b", "b": "a"})
    ga = GroupAction.cyclic(swap, 2)             # char 2 | |G| = 2
    with pytest.raises(QuiverlabError, match="modular|char"):
        stefan_decomposition(A, ga, 2, side="coh")
    # the constructor + DIRECT HH stay available even in the modular case
    S = skew_group_algebra(A, ga)
    assert S.dim == 8
