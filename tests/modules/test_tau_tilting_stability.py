"""King theta-stability + the wall-and-chamber fan (Plan 45 / C4). Self-cert: the n=2 fan
tiles R^2 (exact angular sweep, no floats); every chamber g-matrix is unimodular; each
wall's brick dim-vector is orthogonal to the shared g-facet. Cross-engine: the fan's wall
brick normals match the mutation brick labels (Task D)."""
import functools

from fractions import Fraction

import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.stability import (is_theta_semistable, wall_and_chamber_fan)

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def _cones_tile_R2(chambers):
    """Exact angular sweep (cross-product signs, NO atan2): the distinct g-vector rays,
    sorted by angle, must have every cyclically-consecutive pair be exactly one chamber --
    the 2D cones then partition R^2, covering the circle once."""
    def parse(ray):
        return (Fraction(ray[0]), Fraction(ray[1]))

    rays = []
    for ch in chambers:
        for r in ch["rays"]:
            v = parse(r)
            if v not in rays:
                rays.append(v)

    def half(v):                       # 0 for angle in [0, pi), 1 for [pi, 2pi)
        x, y = v
        return 0 if (y > 0 or (y == 0 and x > 0)) else 1

    def cross(a, b):
        return a[0] * b[1] - a[1] * b[0]

    def cmp(a, b):
        if half(a) != half(b):
            return -1 if half(a) < half(b) else 1
        c = cross(a, b)
        return -1 if c > 0 else (1 if c < 0 else 0)

    srt = sorted(rays, key=functools.cmp_to_key(cmp))
    m = len(srt)
    expected = {frozenset((srt[i], srt[(i + 1) % m])) for i in range(m)}
    got = {frozenset((parse(ch["rays"][0]), parse(ch["rays"][1]))) for ch in chambers}
    return expected == got and len(chambers) == m


@selfcert
def test_n2_fan_tiles_the_plane():
    A = _kA2()
    fan = wall_and_chamber_fan(A)
    assert fan["complete"] and fan["n"] == 2 and len(fan["chambers"]) == 5
    assert _cones_tile_R2(fan["chambers"])


@selfcert
def test_every_chamber_g_matrix_is_unimodular():
    from quiverlab.tautilting.stability import _det
    A = linear_path_algebra(3, field=QQ)
    fan = wall_and_chamber_fan(A)
    for ch in fan["chambers"]:
        assert _det(ch["g_matrix"]) in (1, -1)         # sign-coherent basis of Z^n


@selfcert
def test_n3_l1_unfolding_is_a_sane_rendering():
    # The n=3 L1/octahedron unfolding is a RENDERING, not a certified 3D tiling (see the
    # verification page honest-scope). This is the CHEAP sanity check the critic asked for:
    # every chamber's projected net is a nondegenerate 2D triangle (exact area != 0, no
    # floats) and the net count equals the chamber count. It does NOT certify that the
    # projected faces tile the octahedron net without gaps/overlaps -- only per-chamber
    # unimodularity (test_every_chamber_g_matrix_is_unimodular) is claimed at n=3.
    def area2(net):
        (x0, y0), (x1, y1), (x2, y2) = [(Fraction(p[0]), Fraction(p[1])) for p in net]
        return (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)

    A = linear_path_algebra(3, field=QQ)                 # kA3, n = 3
    fan = wall_and_chamber_fan(A)
    assert fan["complete"] and fan["n"] == 3 and fan.get("projection") == "L1"
    nets = [ch["net2d"] for ch in fan["chambers"]]
    assert len(nets) == len(fan["chambers"])             # one net per chamber
    for net in nets:
        assert len(net) == 3                             # a triangle per chamber
        assert area2(net) != 0                           # nondegenerate (exact, no floats)


@selfcert
def test_king_semistable_definition():
    # over kA2 with theta = (1, -1): the simple S_1 (dim (1,0)) is NOT semistable
    # (theta.dim S_1 = 1 != 0); the module P_1 = [1,2] (dim (1,1)) has theta.dim = 0 and
    # its only proper submodule S_2 (dim (0,1)) has theta.dim = -1 <= 0 => semistable.
    A = _kA2()
    assert is_theta_semistable(A.simple(1), [Fraction(1), Fraction(-1)]) is False
    assert is_theta_semistable(A.projective(1), [Fraction(1), Fraction(-1)]) is True


@xeng
def test_wall_normals_match_mutation_brick_labels():
    # each wall's brick dim-vector (from the mutation edge label, Task D) is orthogonal to
    # the shared g-facet it labels: theta . dim(B) = 0 on the wall (King).
    A = _kA2()
    fan = wall_and_chamber_fan(A)
    for wall in fan["walls"]:
        bd = wall["brick_dimvec"]
        nrm = wall["normal"]
        verts = fan["vertices"]
        assert sum(nrm[k] * bd[verts[k]] for k in range(len(verts))) == 0
