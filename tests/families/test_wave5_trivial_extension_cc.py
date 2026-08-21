"""Wave 5 (v1.0.1): ``TrivialExtension(A)`` over the DEFAULT field CC.

Finding 5c: Plan 31 delivers a presented ``T(A)`` "over QQ/GF(p)", but on CC
(the default field, a SympyExactDomain that computes exactly in QQ) it silently
dropped to the structure-constant fallback where CS refuses and bar blows up, so
``T.is_symmetric()`` / ``T.hochschild_cohomology`` hit a wall. The CC-over-QQ
path now routes through the rational presentation (the coeff emitter renders exact
rationals), re-tagged to CC, and is certified per instance by ``dim kQ_T/I_T ==
2*dim A`` (the Plan-31 gate). GF(p)/QQ stay numerically identical.
"""
import pytest

import quiverlab
from quiverlab import fields

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
def test_trivial_extension_cc_is_presented():
    T = quiverlab.TrivialExtension(quiverlab.linear_path_algebra(2))  # CC default
    assert T.quiver is not None                       # presented, not structure-constant
    assert T.dim == 2 * 3                             # dim T = 2 dim A (kA_2 has dim 3)


@selfcert
def test_trivial_extension_cc_is_symmetric():
    T = quiverlab.TrivialExtension(quiverlab.linear_path_algebra(2))  # CC default
    assert T.is_symmetric() is True                  # every T(A) is symmetric


@xeng
def test_trivial_extension_cc_hochschild_matches_qq():
    # HH^*(T(kA_2)) must agree between the CC presentation and the QQ one
    # (char-0 field-independent), and equal the Plan-31 pin [3,1,1,1,1].
    Tcc = quiverlab.TrivialExtension(quiverlab.linear_path_algebra(2))
    Tqq = quiverlab.TrivialExtension(quiverlab.linear_path_algebra(2, field=fields.QQ))
    coh_cc = Tcc.hochschild_cohomology(4).dims
    coh_qq = Tqq.hochschild_cohomology(4).dims
    assert coh_cc == coh_qq == [3, 1, 1, 1, 1]


@lit
def test_trivial_extension_qq_and_gfp_unchanged():
    # DO NOT BREAK: the literature-pinned QQ/GF(p) builds stay identical.
    Tqq = quiverlab.TrivialExtension(quiverlab.linear_path_algebra(2, field=fields.QQ))
    assert Tqq.dim == 6 and Tqq.quiver is not None and Tqq.is_symmetric() is True
    Tp = quiverlab.TrivialExtension(quiverlab.linear_path_algebra(2, field=quiverlab.GF(5)))
    assert Tp.dim == 6 and Tp.quiver is not None
