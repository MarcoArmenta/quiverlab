"""P55/P60 tie: the left/right support algebras A_lambda, A_rho are PRODUCTS OF TILTED algebras
(ACLV Thm A for ada). On the shared ACLV 2.2(b) rad^2=0 A5 fixture, A_lambda (verts {1,2,3}) and
A_rho (verts {3,4,5}) are each kA3/rad^2; EACH connected component is certified tilted by the
Plan-60 recognizer (this is the content P55 PIN'd for P60). QQ-scope."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

pytestmark = pytest.mark.oracle_literature


def _radsq_a5():
    return RadicalSquareZero(Quiver([1, 2, 3, 4, 5],
        {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)}), field=QQ)


def test_hereditary_support_components_tilted():
    atlas = left_right_parts(linear_path_algebra(3, field=QQ))
    for comp in atlas.left_support.components:
        assert comp["algebra"].is_tilted() is True       # hereditary factor => tilted


def test_aclv_22b_support_components_tilted():
    atlas = left_right_parts(_radsq_a5())
    for support in (atlas.left_support, atlas.right_support):
        assert len(support.components) == 1               # each connected
        for comp in support.components:
            assert comp["algebra"].is_tilted() is True    # kA3/rad^2 factor => tilted (type A3)
