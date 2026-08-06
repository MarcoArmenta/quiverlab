"""Live M2 SpectralSequences crosscheck (commutative-only, m2-gated). Compares the
convention-robust E_inf totals (== total-complex homology) of a small commutative
Koszul double complex both systems build.

ARBITRATED / RECORDED (Plan-42, live probe under M2 1.26.06): M2's
``SpectralSequences`` package could NOT be scripted -- it rides the ``ChainComplex``
type removed in the Complexes-based core (``chainComplex`` errors as an undefined
Symbol), so its page objects do not reach ``--script``/sentinel stdout the way
``Complexes`` does. Per the plan's sanctioned fallback, the compared quantity is the
``E_inf`` totals == M2's ``Complexes`` homology of the SAME total complex (strong
convergence makes this convention-robust: both systems must agree on the
convergence target). The ``E_2`` cell grid is NOT compared (M2's page indexing was
not confirmable live). The example is the Koszul complex ``K(x) (x) K(y)`` over
``k[x,y]/(x^2, xy, y^2)`` -- genuinely non-regular, so the homology ``[1, 3, 2]`` is
nonzero in every degree and the SS does real work (it degenerates at ``E_2``, not
``E_1``)."""
import pytest

from quiverlab.m2 import crosscheck as cc
from quiverlab.m2 import session

pytestmark = pytest.mark.skipif(session.should_skip_m2(),
                                reason="Macaulay2 not installed")


def test_commutative_ss_einf_totals_match_m2():
    # a small filtered/double complex over ZZ/7; compare E_inf totals (robust).
    cc.crosscheck_commutative_ss(p=7, top=4).assert_agree()
