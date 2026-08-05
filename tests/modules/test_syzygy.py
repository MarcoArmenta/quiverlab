"""Public syzygy / cosyzygy. Self-certifying: the syzygy is the kernel of the
projective cover (the exact submodule the minimal resolution already builds)."""
import pytest

from quiverlab import GF, Quiver, truncated_polynomial
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.resolution import cosyzygy, syzygy

pytestmark = pytest.mark.oracle_selfcert


def _a3rel():
    # kA3 with a*b = 0 (Quiver 1->2->3, radical-square-zero on the 1->3 path).
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_syzygy_of_simple_is_radical_of_cover():
    A = _a3rel()
    S1 = A.simple(1)
    Om = syzygy(S1)                          # Omega S1 = rad P1 = S2 here
    assert Om.dimension_vector() == A.simple(2).dimension_vector()
    assert is_isomorphic(Om, A.simple(2))


def test_syzygy_of_projective_is_zero():
    A = _a3rel()
    assert syzygy(A.projective(2)).dim == 0


def test_truncated_polynomial_syzygy_shift():
    # k[x]/(x^3): Omega(k[x]/(x^i)) = k[x]/(x^{3-i}); Omega S = M_2 (dim 2).
    A = truncated_polynomial(3, field=GF(7))
    S = A.simple(1)
    Om = syzygy(S)
    assert Om.dim == 2
    assert syzygy(Om).dim == 1               # Omega^2 S = M_1 = S again (dim 1)


def test_cosyzygy_matches_dual_route():
    # k[x]/(x^3) is self-injective: cosyzygy(S) = Omega^{-1} S = M_2 (dim 2).
    A = truncated_polynomial(3, field=GF(7))
    S = A.simple(1)
    co = cosyzygy(S)
    assert co.dim == 2
    from quiverlab.modules.duality import dualize
    assert co.dimension_vector() == dualize(syzygy(dualize(S))).dimension_vector()


def test_module_method_delegates():
    A = _a3rel()
    assert A.simple(1).syzygy().dim == syzygy(A.simple(1)).dim
    assert A.simple(1).cosyzygy().dim == cosyzygy(A.simple(1)).dim
