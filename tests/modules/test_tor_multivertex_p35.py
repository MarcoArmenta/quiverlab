"""Plan 35 wave 3a review carry-over -- a MULTI-VERTEX Tor self-cert that exercises the
nonzero ``_tor_boundary`` collapse.

On a 3-vertex cyclic Nakayama with rad^2 = 0 (kZ_3 / J^2, infinite global dimension) the
minimal resolution of a simple cycles through the three vertices, so ``d_n (x) 1`` has to
COLLAPSE a genuinely nonzero contribution ``(g_l . p) (x) y = g_l (x) (p . y)`` back into
a target summand (the branch that solves into ``e_w N``). This pins:
  * the captured Tor dims equal the engine ``tor_dims`` (the collapse is rank-correct);
  * at least one shipped differential carries a NONZERO entry (the collapse path fires);
  * every shipped class vector is annihilated by its shipped differential (self-cert).
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.complex_reps import tor_reps
from quiverlab.modules.tor import tor_dims


def _nakayama_rad2():
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (3, 1)})
    return Q.algebra(relations=["a*b", "b*c", "c*a"], field=GF(5))


def _dense(vector, ncols, dom):
    v = [dom.zero()] * ncols
    for idx, c in vector:
        v[idx] = dom.coerce(int(c)) if str(c).lstrip("-").isdigit() else dom.coerce(c)
    return v


@pytest.mark.oracle_crossengine
def test_multivertex_tor_boundary_collapse():
    A = _nakayama_rad2()
    M = A.simple(2)                                  # right module S_2
    N = A.injective(1, side="left")                  # left module with nonzero rad action
    dom = A.domain
    top = 4
    dims, pl = tor_reps(A, M, N, top)
    # (a) dims match the engine -- the collapse is rank-correct
    assert dims == tor_dims(A, M, N, top)
    # (b) at least one differential has a nonzero entry (the collapse path fired)
    nonzero_degrees = []
    for k, d in pl["differentials"].items():
        rows = d.get("rows")
        if rows and any(str(x) != "0" for row in rows for x in row):
            nonzero_degrees.append(int(k))
    assert nonzero_degrees, "expected a nonzero _tor_boundary on this Nakayama example"
    # (c) self-cert: every shipped class vector is annihilated by its differential
    certified = 0
    for deg, classes in pl["basis_classes"].items():
        d = pl["differentials"][deg]
        assert not d.get("elided")
        rows, ncols = d.get("rows") or [], d["shape"][1]
        for cl in classes:
            v = _dense(cl["vector"], ncols, dom)
            if rows:
                R = [[dom.coerce(int(x)) if str(x).lstrip("-").isdigit() else dom.coerce(x)
                      for x in row] for row in rows]
                assert all(dom.is_zero(p) for p in lm.matvec(R, v, dom))
            certified += 1
    assert certified == sum(dims)
