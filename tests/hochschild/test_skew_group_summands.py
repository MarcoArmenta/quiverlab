"""Plan 74 Task 3: the twisted summand + Z(g)-action transport + Reynolds
invariants (the P52 bridge). The bridge-non-triviality pins (findings 1a/2b/2e):
a NONZERO twisted invariant at homology n=1, hand-derived two ways."""
import pytest

from quiverlab.families.basic import truncated_polynomial
from quiverlab.families.skew_group import GroupAction, QuiverAutomorphism
from quiverlab.fields import QQ
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.hochschild.skew_group import stefan_decomposition

pytestmark = [pytest.mark.oracle_crossengine, pytest.mark.oracle_selfcert]


def _flagship():
    A = truncated_polynomial(2, field=QQ)
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": -1})
    return A, GroupAction.cyclic(sig, 2)


@pytest.mark.oracle_crossengine
def test_twisted_summand_and_regular_reuse():
    A, _ = _flagship()
    dom = A.domain
    sig_mat = [[dom.one(), dom.zero()], [dom.zero(), dom.coerce(-1)]]
    Msig = Bimodule.twisted(A, phi=sig_mat)
    assert Msig.check() is True
    assert A.hochschild_cohomology(3, coefficients=Msig).dims == [1, 1, 1, 1]
    assert A.hochschild_homology(3, coefficients=Msig).dims == [1, 1, 1, 1]
    # P52 M=A oracle: regular bimodule reproduces the ordinary HH
    reg = Bimodule.regular(A)
    assert A.hochschild_cohomology(3, coefficients=reg).dims == \
        A.hochschild_cohomology(3).dims == [2, 1, 1, 1]


@pytest.mark.oracle_selfcert
def test_cohomology_degree0_hand_pins():
    A, ga = _flagship()
    rep = stefan_decomposition(A, ga, 3, side="coh")
    # summands ordered [identity, sigma]
    assert rep.summands[0]["g"] == "e"
    assert rep.summands[0]["inv"][0] == 1        # HH^0(A,A)^G = span{1}
    assert rep.summands[1]["inv"][0] == 0        # HH^0(A,{}_sigma A)^G = 0
    # full cohomology invariant profile
    assert rep.summands[0]["inv"] == [1, 1, 1, 1]
    assert rep.summands[1]["inv"] == [0, 0, 0, 0]
    assert rep.dims == [1, 1, 1, 1]


@pytest.mark.oracle_crossengine
def test_homology_per_summand_split_nonzero_twisted_at_n1():
    A, ga = _flagship()
    rep = stefan_decomposition(A, ga, 3, side="hom")
    # n=0: identity coinv 1 AND twisted coinv 1 (finding 2e)
    assert rep.summands[0]["inv"][0] == 1
    assert rep.summands[1]["inv"][0] == 1
    # n=1: identity coinv 0 (sigma = -1 on HH_1 = Omega^1 = k.dx), twisted coinv 1
    # (NONZERO, carries the whole degree) -- findings 1a/2b, NOT [1, 0].
    assert rep.summands[0]["inv"][1] == 0
    assert rep.summands[1]["inv"][1] == 1
    assert rep.summands[0]["inv"] == [1, 0, 0, 0]
    assert rep.summands[1]["inv"] == [1, 1, 1, 1]
    assert rep.dims == [2, 1, 1, 1]
