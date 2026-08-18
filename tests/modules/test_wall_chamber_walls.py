"""The wall D(B) of a brick as an EXACT inequality system (Plan 63 / R25). D(B) =
{theta : theta.dim B = 0, theta.dim N <= 0 for every submodule N <= B}. Literature (BST
2019): D(B) is a rational polyhedral cone of codim >= 1; a simple brick gives the full
hyperplane, a non-simple brick a proper face. THE PIN: over kA2, D(S1)/D(S2) are full lines
but D(P1) is a RAY (the record's 'the P1 wall is a ray')."""
import pytest
from fractions import Fraction

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.tautilting.stability import is_theta_semistable
from quiverlab.tautilting.wallchamber import Wall, wall_of_brick

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@lit
def test_simple_brick_wall_is_full_hyperplane():
    A = _kA2()
    w1 = wall_of_brick(A.simple(1))
    assert w1.brick_dimvec == (1, 0) and w1.equality == (1, 0)
    assert w1.is_full_hyperplane is True and w1.inequalities == ()
    assert w1.codim == 1 and w1.brick_name == "S1"
    # full line {theta_1 = 0}: two opposite rays along the theta_2 axis
    assert set(map(tuple, [[Fraction(r[0]), Fraction(r[1])] for r in w1.rays])) == {
        (Fraction(0), Fraction(1)), (Fraction(0), Fraction(-1))}


@lit
def test_P1_wall_is_a_ray():
    # THE record pin. P1 = (1,1), unique proper submodule S2 = (0,1); D(P1) = {theta_1 +
    # theta_2 = 0, theta_2 <= 0} = the single ray (1,-1). NOT the full line.
    A = _kA2()
    w = wall_of_brick(A.projective(1))
    assert w.brick_dimvec == (1, 1) and w.equality == (1, 1)
    assert w.is_full_hyperplane is False           # a proper face -- a RAY
    assert (0, 1) in w.inequalities                # the S2 submodule constraint theta_2 <= 0
    rays = [(Fraction(r[0]), Fraction(r[1])) for r in w.rays]
    assert len(rays) == 1                          # a ray, not a line (one extreme ray)
    d = rays[0]
    assert d[0] > 0 and d[1] < 0 and d[0] + d[1] == 0   # direction (1,-1) up to scale


@selfcert
def test_wall_points_are_theta_semistable():
    # every theta in the relative interior of D(B) makes B theta-semistable (King) -- the
    # cross-check tying the inequality system to the shipped stability predicate.
    A = _kA2()
    B = A.projective(1)
    w = wall_of_brick(B)
    # a point on the ray (1,-1): theta = (1,-1)
    assert is_theta_semistable(B, [Fraction(1), Fraction(-1)]) is True
    # a point on the OPPOSITE ray (-1,1) is NOT in D(P1) (theta_2 = 1 > 0): not semistable
    assert is_theta_semistable(B, [Fraction(-1), Fraction(1)]) is False
    # the equality holds on every ray; every inequality holds on every ray
    for r in w.rays:
        theta = [Fraction(r[0]), Fraction(r[1])]
        assert sum(theta[i] * w.equality[i] for i in range(2)) == 0
        for d in w.inequalities:
            assert sum(theta[i] * d[i] for i in range(2)) <= 0


@selfcert
def test_codim_and_equality_normal():
    A = _kA2()
    for M in (A.simple(1), A.simple(2), A.projective(1)):
        w = wall_of_brick(M)
        dv = M.dimension_vector()
        assert w.equality == tuple(dv[v] for v in (1, 2))
        assert w.codim == 1
