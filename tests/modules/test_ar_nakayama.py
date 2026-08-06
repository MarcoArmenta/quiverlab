"""Nakayama functor. Self-cert: nu(P_v) ~ I_v, nu^-(I_v) ~ P_v, and the ker of
the induced injective map ~ tau M (ties nu to the trusted translate). Additive."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.modules.ar import nakayama_functor, nakayama_functor_minus
from quiverlab.modules.hom import is_isomorphic

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    from quiverlab.fields import QQ
    return linear_path_algebra(3, field=QQ)


def test_nu_of_projective_is_injective():
    A = _a3()
    for v in (1, 2, 3):
        assert is_isomorphic(nakayama_functor(A.projective(v)), A.injective(v))


def test_nu_minus_of_injective_is_projective():
    A = _a3()
    for v in (1, 2, 3):
        assert is_isomorphic(nakayama_functor_minus(A.injective(v)), A.projective(v))


def test_nu_ties_to_tau_on_nonprojectives():
    # SHARPENED (Step 3): the internal certificate is ker(nu P_1 -> nu P_0) ~ tau M;
    # `_return_ker=True` exposes that kernel so the test pins the identity exactly
    # against the module's own trusted tau() (the arbiter). S1 over kA3 is
    # non-projective (S1 != P1 = [1,2,3]), so tau S1 is nonzero.
    A = _a3()
    S1 = A.simple(1)
    assert not is_isomorphic(S1, A.projective(1))       # S1 non-projective
    _nuS1, ker_g = nakayama_functor(S1, _return_ker=True)
    assert is_isomorphic(ker_g, S1.tau())               # ker g == tau S1, exactly


def test_nu_additive():
    A = _a3()
    from quiverlab.modules.morphism import direct_sum
    D, _, _ = direct_sum(A.simple(1), A.simple(2))
    nD = nakayama_functor(D)
    n1, n2 = nakayama_functor(A.simple(1)), nakayama_functor(A.simple(2))
    assert nD.dimension_vector() == {v: n1.dimension_vector().get(v, 0)
                                     + n2.dimension_vector().get(v, 0)
                                     for v in A.quiver.vertices}
