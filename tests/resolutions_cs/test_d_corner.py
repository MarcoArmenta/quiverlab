"""Plan 51 Task A.0: the per-corner resolution-differential matrix D_corner
(the new primitive underlying the homotopy-lifting solve).

Deep bucket (tests/resolutions_cs -> deep).  D_corner(res, m, o, t) is the matrix
of d_P: P_m -> P_{m-1} restricted to the (o,t)-corner (columns = the corner
k-basis of P_m, rows = the corner k-basis of P_{m-1}), assembled off res.d_terms
through pelt.apply_lower -- the single-complex sibling of
diagonal.TensorComplex.tensor_matrix, NOT res.matrix (the Hom-dual coboundary).

THE correctness gate for the primitive is d^2 = 0:
    D_corner(m-1, o, t) . D_corner(m, o, t) == 0  on every corner.
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.homotopy_lifting import D_corner, _corner_complex

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_selfcert]


def _res(rels, arrows, verts, md, p=5):
    f = GF(p)
    Q = Quiver(verts, arrows)
    A = Q.algebra(relations=rels, field=f)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, f), max_degree=md)


def _kx2(md=6):
    return _res(["x*x"], {"x": (1, 1)}, [1], md)


def _straddle(md=4):
    return _res(["x*x", "y*y", "x*y*x"], {"x": (1, 1), "y": (1, 1)}, [1], md)


def _qci(md=4):
    return _res(["x*x", "y*y", "y*x - 2*x*y"], {"x": (1, 1), "y": (1, 1)}, [1], md)


def _comm_square(md=3):
    # kQ/(a*b - c*d): 1->2->4, 1->3->4; gldim 2 (S(3) empty).
    return _res(["a*b - c*d"],
                {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}, [1, 2, 3, 4], md)


def _matmul(dom, A_, B_):
    if not A_ or not B_:
        return []
    inner, rows, cols = len(B_), len(A_), len(B_[0])
    out = [[dom.zero()] * cols for _ in range(rows)]
    for i in range(rows):
        for k in range(inner):
            if dom.is_zero(A_[i][k]):
                continue
            for j in range(cols):
                out[i][j] = dom.add(out[i][j], dom.mul(A_[i][k], B_[k][j]))
    return out


def _corners(res):
    """All (o, t) vertex pairs of the underlying quiver."""
    verts = list(res.rs.quiver.vertices)
    return [(o, t) for o in verts for t in verts]


# --------------------------------------------------------------------------- #
# A.0 -- shape + the d^2 = 0 self-certification                                #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk,name,top", [
    (_kx2, "kx2", 5), (_straddle, "straddle", 3),
    (_qci, "qci", 3), (_comm_square, "comm_square", 2)])
def test_d_corner_dd_zero(mk, name, top):
    res = mk()
    dom = res.dom
    saw = False
    for (o, t) in _corners(res):
        for m in range(2, top + 1):
            lo = D_corner(res, m - 1, o, t)
            hi = D_corner(res, m, o, t)
            prod = _matmul(dom, lo, hi)
            assert not any(not dom.is_zero(x) for row in prod for x in row), \
                f"d^2 != 0 on {name} corner ({o},{t}) at m={m}"
            if lo and hi and lo[0] and hi[0]:
                saw = True
    assert saw, f"no non-trivial D_corner product exercised on {name}"


def test_d_corner_shapes_and_basis_order_kx2():
    """columns index P_m's corner k-basis, rows index P_{m-1}'s; the basis ordering
    is the deterministic (S-sequence, ascending corner indices) order, and the
    empty branch (P_{-1} = 0) makes D_corner(0) have no rows."""
    res = _kx2()
    cc = _corner_complex(res)
    for m in range(1, 5):
        M = D_corner(res, m, 1, 1)
        cols = cc.p_basis(m, 1, 1)
        rows = cc.p_basis(m - 1, 1, 1)
        assert len(M) == len(rows) and (not M or len(M[0]) == len(cols)), \
            f"D_corner({m}) shape {len(M)}x{len(M[0]) if M else 0} != rows/cols"
    # base: P_{-1} = 0 -> no rows (the "no equations" branch of the psi-solve).
    assert D_corner(res, 0, 1, 1) == []
    # p_basis is the deterministic single-PELT enumeration (degree-0 chain word ()).
    assert cc.p_basis(0, 1, 1) == [(0, (), 0), (0, (), 1), (1, (), 0), (1, (), 1)]


def test_d_corner_is_not_the_hom_coboundary_kx2():
    """D_corner is the RESOLUTION differential d_P (single-complex), a different
    object from res.matrix(n,'coh') (the Hom-dual coboundary delta^n: C^n->C^{n+1}):
    different shape and direction."""
    res = _kx2()
    d1 = D_corner(res, 1, 1, 1)               # P_1 -> P_0: 4 x 4 (both corner bases)
    coh1 = res.matrix(1, "coh")               # C^1 -> C^2: dim_C(2) x dim_C(1)
    # they are genuinely different matrices (shape or content); assert not equal.
    assert d1 != coh1
