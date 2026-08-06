"""Irreducible maps = dim rad(M,N)/rad^2(M,N). Literature: kA_n has exactly the
mesh arrows (each interior irreducible-map multiplicity 1)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import irreducible_maps

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def _indecs(A):
    # kA3 indecomposables = the 6 interval modules; realise the small ones here.
    return [A.simple(v) for v in (1, 2, 3)] + [A.projective(v) for v in (1, 2)]


@selfcert
def test_no_irreducible_self_map_for_bricks():
    A = _a3()
    S2 = A.simple(2)
    within = _indecs(A)
    assert irreducible_maps(S2, S2, within) == 0     # rad End(S2) = 0 (a brick)


@lit
def test_ka3_mesh_arrow_multiplicities():
    A = _a3()
    within = _indecs(A)
    S2, S3, P1, P2 = A.simple(2), A.simple(3), A.projective(1), A.projective(2)
    # ZA3 mesh (see test_ar_almost_split): [3] --1--> [2,3] is an irreducible map, so
    # irreducible_maps(S3, P2) = 1; there is NO irreducible map [2] -> [1,2,3], so
    # irreducible_maps(S2, P1) = 0. (The arrow [1,2,3] -> [1] factors through the
    # missing interval [1,2]; over this partial `within` it is not yet resolved, which
    # is the honest "relative to within" contract -- exact only on a closed set.)
    assert irreducible_maps(S3, P2, within) == 1     # mesh edge [3] -> [2,3]
    assert irreducible_maps(S2, P1, within) == 0     # no mesh edge [2] -> [1,2,3]
    assert all(irreducible_maps(X, X, within) == 0 for X in within)  # no loops
