"""Q(T)/W(T)/Jac(T) (Plan 48, deep). ARBITER: fan of the (n+3)-gon -> linear A_n quiver
EXACTLY (fixes the angle orientation). Literature: hexagon-internal-triangle Jac dim 6.
Cross-subsystem (P44+P46+P48): is_gentle(jacobian_of(T)) True across the surface zoo.
Loud: punctures / closed surfaces are out of v1 scope."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.invariants.recognizers import is_gentle
from quiverlab.surfaces.qp import jacobian_of, potential_of, quiver_of
from quiverlab.surfaces.triangulation import (annulus_triangulation, fan_triangulation,
                                             hexagon_with_internal_triangle,
                                             once_punctured_torus, Triangulation)

arb = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _hexagon_internal_triangle():
    return hexagon_with_internal_triangle()


def _linear_A(n):
    arrows = {f"a{i}": (i, i + 1) for i in range(1, n)}
    return Quiver(list(range(1, n + 1)), arrows)


@arb
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_disc_fan_is_linear_A_n(n):
    # THE ARBITER: fan of the (n+3)-gon -> linear A_n quiver 1->2->...->n.
    Q = quiver_of(fan_triangulation(n + 3))
    verts = list(Q.vertices)
    assert len(verts) == n
    # underlying arrows == a linear chain on the arc ids (up to the arc-id bijection):
    ends = sorted((Q.source(a), Q.target(a)) for a in Q.arrows)
    chain = sorted((v, v + 1) for v in verts[:-1]) if n >= 2 else []
    assert len(ends) == max(n - 1, 0)
    # the arc ids 1..n are consecutive and each consecutive pair carries exactly one arrow
    assert {tuple(sorted(e)) for e in ends} == {tuple(sorted(c)) for c in chain}


@arb
def test_disc_fan_is_linear_A_n_directed():
    # the arbiter pins the DIRECTION too: 1 -> 2 -> ... -> n (not the reverse).
    Q = quiver_of(fan_triangulation(7))                  # A_4
    ends = sorted((Q.source(a), Q.target(a)) for a in Q.arrows)
    assert ends == [(1, 2), (2, 3), (3, 4)]


@lit
def test_hexagon_internal_triangle_jacobian_dim_6():
    # hexagon (6 marked pts) triangulated with the central triangle {13,35,51}: one
    # internal triangle -> one 3-cycle -> Jac = triangle algebra dim 6 (P44's pin).
    T = _hexagon_internal_triangle()
    J = jacobian_of(T, field=QQ)
    assert J.dim == 6
    assert is_gentle(J)


@lit
def test_disc_fan_jacobian_is_kA_n():
    # fan => no internal triangle => W = 0 => Jac = kA_n (hereditary). D_5 (pentagon) -> kA_2.
    J = jacobian_of(fan_triangulation(5), field=QQ)      # pentagon
    assert J.dim == 3                                    # kA_2: 2 vertices + 1 arrow
    assert is_gentle(J)


@arb
@pytest.mark.parametrize("factory", [
    lambda: fan_triangulation(5), lambda: fan_triangulation(7),
    lambda: annulus_triangulation(2, 1), lambda: annulus_triangulation(2, 2),
    lambda: _hexagon_internal_triangle(),
])
def test_surface_jacobians_are_gentle(factory):
    # ABCP/LFS: unpunctured-with-boundary surface Jacobians are gentle. P44+P46+P48 tie.
    J = jacobian_of(factory(), field=QQ)
    assert is_gentle(J)


@lit
def test_annulus_C21_is_affine_A2_shape():
    # smallest non-degenerate annulus: C(2,1) -> 3 arcs, affine A~_2 (acyclic), hereditary,
    # gentle. Pin the exact arrow count (arbitrated by the disc-fixed convention +
    # finiteness).
    Q = quiver_of(annulus_triangulation(2, 1))
    assert len(list(Q.vertices)) == 3
    assert len(Q.arrows) == 3                            # affine A~_2 cycle, acyclically oriented
    assert Q.is_acyclic()                                # finiteness doubly-guards orientation
    J = jacobian_of(annulus_triangulation(2, 1), field=QQ)
    assert is_gentle(J)


@selfcert
def test_punctured_and_closed_refused():
    with pytest.raises(QuiverlabError):
        quiver_of(once_punctured_torus())                # p>0 out of v1 scope
