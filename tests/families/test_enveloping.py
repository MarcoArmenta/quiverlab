"""``B^e = B (x) B^op`` as a first-class bound quiver algebra (Plan 73 / R7,
Task II2) -- the ``gl.dim B = infinity`` fallback substrate for leg (ii). The
primary oracles (Ex. 5.3/5.5) do NOT exercise this (both have ``gl.dim B = 2``).
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import GF, QQ
from quiverlab.families import truncated_polynomial
from quiverlab.families.extension import enveloping_algebra
from quiverlab.invariants.han import _bimodule_to_Be_module

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA(n, field=QQ):
    verts = list(range(1, n + 1))
    arrows = {f"a{i}": (i, i + 1) for i in range(1, n)}
    return Quiver(verts, arrows).algebra(relations=[], field=field)


@selfcert
def test_enveloping_dim_and_quiver():
    B = _kA(2)                                       # kA2, dim 3
    Be = enveloping_algebra(B)
    assert Be.quiver is not None
    assert Be.dim == B.dim ** 2 == 9                 # certified iso to B (x) B^op
    # product-quiver vertices = (Q_B)_0 x (Q_B)_0
    assert len(list(Be.quiver.vertices)) == len(list(B.quiver.vertices)) ** 2


@selfcert
def test_enveloping_dual_numbers():
    B = truncated_polynomial(2, field=QQ)            # k[x]/(x^2), dim 2
    Be = enveloping_algebra(B)
    assert Be.dim == 4 and Be.quiver is not None
    assert len(list(Be.quiver.vertices)) == 1


@xeng
def test_pd_Be_regular_matches_hochschild_vanishing():
    """pd_{B^e}(B) via the enveloping algebra vs the SHIPPED public HH engine:
    ``pd_{B^e}(B)`` is the length of the minimal bimodule resolution, so
    ``HH^n(B, B) = 0`` for ``n > pd_{B^e}(B)``. Hereditary ``kA_n`` -> ``pd = 1``
    (HH^{>=2} = 0); ``k[x]/(x^2)`` -> infinite (HH never vanishes)."""
    from quiverlab.hochschild.coefficients import Bimodule

    for n in (2, 3):
        B = _kA(n)
        Be = enveloping_algebra(B)
        N = _bimodule_to_Be_module(Bimodule.regular(B), B, Be)
        pd_env = N.projective_resolution(16).pd()
        assert pd_env == 1                            # hereditary -> Hochschild dim 1
        # cross-engine anchor: HH^n(B) vanishes above pd_{B^e}
        coh = B.hochschild_cohomology(4, engine="cs").dims
        assert coh[2:] == [0, 0, 0]                   # HH^{>=2} = 0, consistent with pd 1

    # dual numbers: infinite Hochschild dimension (unresolved + HH never dies)
    B = truncated_polynomial(2, field=QQ)
    Be = enveloping_algebra(B)
    N = _bimodule_to_Be_module(Bimodule.regular(B), B, Be)
    assert N.projective_resolution(12).pd() is None   # not resolved -> infinite
    coh = B.hochschild_cohomology(5, engine="cs").dims
    assert coh[4] > 0 and coh[5] > 0                  # HH^n keeps up -> infinite pd
