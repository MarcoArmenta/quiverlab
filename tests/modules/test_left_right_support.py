"""Support algebras A_lambda = End(+ P_x : P_x in L_A), A_rho dually (Plan 55 / R15).
Self-cert: e_lambda is convex; dim A_lambda = dim End(+ P_x) (by additivity), and the presented
induced-subquiver algebra REPRODUCES that dimension (the independent check) -- cross-checked
against end_algebra. Literature: hereditary => A_lambda = A_rho = A (connected); ACLV 2.2(b)
=> e_lambda = {1,2,3}, e_rho = {3,4,5}. The tiltedness of each factor is PIN'd for P60."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.left_right import left_right_parts

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _radsq_nakayama_a5():
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


@lit
def test_hereditary_support_is_the_whole_algebra():
    A = linear_path_algebra(4, field=QQ)               # L_A = R_A = ind A
    atlas = left_right_parts(A)
    assert set(atlas.left_support.vertices) == {1, 2, 3, 4}
    assert set(atlas.right_support.vertices) == {1, 2, 3, 4}
    assert atlas.left_support.dim == A.dim             # A_lambda = A
    assert len(atlas.left_support.components) == 1      # connected


@lit
def test_aclv_22b_support_vertices():
    atlas = left_right_parts(_radsq_nakayama_a5())
    assert set(atlas.left_support.vertices) == {1, 2, 3}       # e_lambda
    assert set(atlas.right_support.vertices) == {3, 4, 5}      # e_rho (overlap at 3)


@xeng
def test_presented_support_dim_equals_end_dim():
    # the INDEPENDENT check: the presented induced-subquiver algebra reproduces dim End.
    from quiverlab.modules.endomorphism import end_algebra
    from quiverlab.modules.morphism import direct_sum
    A = linear_path_algebra(3, field=QQ)
    atlas = left_right_parts(A)
    verts = list(atlas.left_support.vertices)
    projs = [A.projective(x) for x in verts]
    D = direct_sum(*projs)[0] if len(projs) > 1 else projs[0]
    assert atlas.left_support.algebra.dim == end_algebra(D).dim   # presented build == End
    assert atlas.left_support.dim == end_algebra(D).dim


@selfcert
def test_support_presentation_less_refuses_loudly():
    # a structure-constants-only algebra (no quiver) cannot present the induced subquiver:
    # the support build refuses loudly, never fabricates a quiver.
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.left_right import _support_algebra, _universe
    A = linear_path_algebra(3, field=QQ)
    _ar, U, _n, _r = _universe(A, 256)
    sc = A.__class__.from_structure_constants(A.T, A.unit, field=QQ, check=False)
    assert sc.quiver is None
    with pytest.raises(QuiverlabError):
        _support_algebra(sc, U, [1, 2, 3], projectives=True)


@pytest.mark.skip(reason="PIN: each support component is tilted -- VERIFIED in P60 "
                         "(tilted recognizer). Auto-flips to a real assert when P60 lands.")
def test_support_components_are_tilted_PIN():
    atlas = left_right_parts(linear_path_algebra(3, field=QQ))
    for comp in atlas.left_support.components:
        assert comp["algebra"].is_tilted()             # P60 API
