"""Omega/tau-periodicity certificates (is_isomorphic-certified). Cyclic Nakayama
kZ_n/J^2 has Omega-periodic simples of period n (from the Kupisch series);
hereditary => not Omega-periodic (finite pd)."""
import pytest

from quiverlab import NakayamaAlgebra, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.homdims import omega_periodicity, tau_periodicity

pytestmark = pytest.mark.oracle_literature


def _cyclic_rad2(n):
    # kZ_n / J^2: Kupisch [2]*n, self-injective; Omega(S_i) = S_{i-1}, period n.
    return NakayamaAlgebra(n=n, l=2, cyclic=True, field=QQ)


def test_cyclic_nakayama_omega_period_is_n():
    A = _cyclic_rad2(3)
    S = A.simple(list(A.quiver.vertices)[0])
    p = omega_periodicity(S)
    assert p == 3
    # certificate self-check: Omega^p S ~ S
    Om = S
    for _ in range(p):
        Om = Om.syzygy()
    from quiverlab.modules.hom import is_isomorphic
    assert is_isomorphic(Om, S)


def test_hereditary_not_omega_periodic():
    A = linear_path_algebra(3, field=QQ)               # hereditary: pd finite
    assert omega_periodicity(A.simple(1)) is None


def test_self_injective_tau_periodic():
    A = _cyclic_rad2(3)
    S = A.simple(list(A.quiver.vertices)[0])
    p = tau_periodicity(S)
    assert p is not None and 1 <= p <= 3
    T = S
    for _ in range(p):
        T = T.tau()
    from quiverlab.modules.hom import is_isomorphic
    assert is_isomorphic(T, S)


@pytest.mark.oracle_selfcert
def test_projective_not_omega_periodic():
    A = linear_path_algebra(3, field=QQ)
    assert omega_periodicity(A.projective(1)) is None  # Omega P = 0
