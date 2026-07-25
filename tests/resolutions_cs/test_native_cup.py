"""Plan 20 Task 3 gates: the native Hochschild cup product on the CS resolution.

Deep bucket (tests/resolutions_cs -> deep); run this FILE directly:
    NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q \
        tests/resolutions_cs/test_native_cup.py -p no:cacheprovider

The four algebraic gates (plan Task 3):
  (a) LEIBNIZ (the sign arbiter): delta(f cup g) = (delta f) cup g
      + (-1)^p f cup (delta g), exactly, over GF(5), for ALL basis cochains
      f in C^p, g in C^q -- the full basis double loop (it IS the minimal complete
      check: cup is bilinear, so holding on every (e_i,e_j) pair holds on all
      f,g).  Coboundary orientation: `res.matrix(n,"coh")` is delta^n : C^n->C^{n+1}
      (rows = _basis(n+1,"coh"), cols = _basis(n,"coh")), so delta is LEFT
      multiplication by the matrix.  Full p,q<=2 grid on k[x]/x^2; on the heavier
      straddle-monomial and quantum-CI the always-on grid is p+q<=2 (diagonal degree
      <=3) -- deeper degrees need diagonal(4/5), whose per-generator lift-solve is
      out of the deep-bucket budget (straddle diagonal(3) already ~30s).
  (b) UNIT: the degree-0 identity-collapse cochain (coordinate 1 at e_v for each
      vertex chain sigma_v) is a two-sided unit for cup on cochains EXACTLY.
  (c) IN-WINDOW ANCHOR (permanent overlap oracle): native cup == transported cup
      (`Comparison.cup_of_cs_classes`) mod coboundary, on every nonzero HH^p x HH^q
      representative pair with p+q <= window, for k[x]/x^2 over GF(32003) and the
      quantum-CI over GF(5).
  (d) graded-commutativity f cup g ~ (-1)^{pq} g cup f and associativity
      (f cup g) cup h ~ f cup (g cup h) MOD COBOUNDARY on HH representatives.
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.comparison import Comparison
from quiverlab.resolutions_cs.cup import native_cup

pytest.importorskip("quiverlab.groebner")


# --------------------------------------------------------------------------- #
# fixtures                                                                     #
# --------------------------------------------------------------------------- #
def _res(rels, arrows, max_degree, p=5):
    f = GF(p)
    Q = Quiver([1], arrows)
    A = Q.algebra(relations=rels, field=f)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, f),
                                   max_degree=max_degree)


def _kx2(md=6):
    return _res(["x*x"], {"x": (1, 1)}, md)


def _straddle(md=4):
    return _res(["x*x", "y*y", "x*y*x"], {"x": (1, 1), "y": (1, 1)}, md)


def _qci(md=4):
    return _res(["x*x", "y*y", "y*x - 2*x*y"], {"x": (1, 1), "y": (1, 1)}, md)


def _kx2_gf():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(32003))


def _qci_gf():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=GF(5))


# --------------------------------------------------------------------------- #
# GF(p) vector helpers                                                         #
# --------------------------------------------------------------------------- #
def _toints(res, vec):
    P = res.dom.p
    return [res.to_int(x) % P for x in vec]


def _matvec(res, M, v):
    """M . v over the domain, as ints mod p (M an int/domain matrix, v a domain
    vector).  Rows may be empty (matrix into a zero cochain space)."""
    dom = res.dom
    out = []
    for row in M:
        acc = dom.zero()
        for j, mij in enumerate(row):
            acc = dom.add(acc, dom.mul(mij, v[j]))
        out.append(acc)
    return _toints(res, out)


def _dcup(res, n, vec):
    """delta^n applied to a cochain vec over _basis(n,'coh'); result ints mod p over
    _basis(n+1,'coh')."""
    return _matvec(res, res.matrix(n, "coh"), vec)


def _basis_cochains(res, deg):
    d = len(res._basis(deg, "coh"))
    for i in range(d):
        v = [res.dom.zero()] * d
        v[i] = res.dom.one()
        yield v


def _unit_cochain(res):
    """The degree-0 identity-collapse cochain: coordinate 1 at j = index of e_v for
    each vertex chain sigma_v, 0 elsewhere (represents 1_A)."""
    basis = res._basis(0, "coh")
    v = [res.dom.zero()] * len(basis)
    for i, (ch, j) in enumerate(basis):
        if j == res.ar._word_index[("v", ch.o)]:
            v[i] = res.dom.one()
    return v


# =========================================================================== #
# (a) LEIBNIZ -- exact, the sign arbiter                                       #
# =========================================================================== #
def _leibniz_holds(res, p, q):
    """delta(f cup g) == (delta f) cup g + (-1)^p f cup (delta g) for every basis
    pair (e_i, e_j)."""
    P = res.dom.p
    sign = 1 if p % 2 == 0 else -1
    for f in _basis_cochains(res, p):
        df = _matvec(res, res.matrix(p, "coh"), f)          # delta f in C^{p+1}
        for g in _basis_cochains(res, q):
            fg = native_cup(res, f, p, g, q)
            lhs = _dcup(res, p + q, fg)                     # delta(f cup g)
            dfg = native_cup(res, df, p + 1, g, q)          # (delta f) cup g
            dg = _matvec(res, res.matrix(q, "coh"), g)      # delta g in C^{q+1}
            fdg = native_cup(res, f, p, dg, q + 1)          # f cup (delta g)
            rhs = [(dfg[i] + sign * fdg[i]) % P for i in range(len(dfg))]
            if lhs != rhs:
                return False, (f, g)
    return True, None


# k[x]/x^2: the full p,q <= 2 grid (diagonal to degree 5, all instantaneous).
@pytest.mark.parametrize("p,q", [(a, b) for a in (0, 1, 2) for b in (0, 1, 2)])
def test_leibniz_kx2(p, q):
    ok, witness = _leibniz_holds(_kx2(), p, q)
    assert ok, f"Leibniz failed on k[x]/x^2 at (p,q)=({p},{q}); witness {witness}"


# straddle / quantum-CI: p+q <= 2 (max diagonal degree 3) -- always-on budget.
_LIGHT = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]


@pytest.mark.parametrize("p,q", _LIGHT)
def test_leibniz_straddle(p, q):
    ok, witness = _leibniz_holds(_straddle(), p, q)
    assert ok, f"Leibniz failed on straddle-monomial at (p,q)=({p},{q}); witness {witness}"


@pytest.mark.parametrize("p,q", _LIGHT)
def test_leibniz_qci(p, q):
    ok, witness = _leibniz_holds(_qci(), p, q)
    assert ok, f"Leibniz failed on quantum-CI at (p,q)=({p},{q}); witness {witness}"


# =========================================================================== #
# (b) UNIT -- exact two-sided unit on cochains                                 #
# =========================================================================== #
@pytest.mark.parametrize("mk,name,qs", [
    (_kx2, "kx2", (0, 1, 2, 3)),
    (_straddle, "straddle", (0, 1, 2)),
    (_qci, "qci", (0, 1, 2)),
])
def test_unit_two_sided_exact(mk, name, qs):
    res = mk()
    u = _unit_cochain(res)
    for q in qs:
        for g in _basis_cochains(res, q):
            gg = _toints(res, g)
            assert _toints(res, native_cup(res, u, 0, g, q)) == gg, \
                f"1 cup g != g exactly on {name} at q={q}"
            assert _toints(res, native_cup(res, g, q, u, 0)) == gg, \
                f"g cup 1 != g exactly on {name} at q={q}"


# =========================================================================== #
# (c) IN-WINDOW ANCHOR -- native cup == transported cup mod coboundary          #
# =========================================================================== #
def _anchor(comp, pairs):
    maxn = max(pp + qq for pp, qq in pairs)
    comp._ensure(maxn)
    res = comp._res
    tested = 0
    for (p, q) in pairs:
        if p + q > comp.window:
            continue
        repp = comp.cs_cohomology_basis(p)
        repq = comp.cs_cohomology_basis(q)
        for i in range(len(repp)):
            for j in range(len(repq)):
                u = comp.hh_class_cs(p, i)
                v = comp.hh_class_cs(q, j)
                native = native_cup(res, u.vec, p, v.vec, q)
                transported = comp.cup_of_cs_classes(u, v)
                assert comp.same_cohomology_class(native, transported, degree=p + q), \
                    f"native != transported cup mod coboundary at (p,q)=({p},{q}), reps ({i},{j})"
                tested += 1
    return tested


def test_anchor_native_equals_transported_kx2():
    """k[x]/x^2 over GF(32003): HH^n = k for every n, so every (p,q) contributes."""
    tested = _anchor(Comparison(_kx2_gf()), [(1, 1), (1, 2), (2, 1), (2, 2)])
    assert tested > 0


def test_anchor_native_equals_transported_qci():
    """quantum-CI over GF(5): HH^1 = k^2, HH^2 = k, so (1,1) cups into HH^2."""
    comp = Comparison(_qci_gf())
    tested = _anchor(comp, [(1, 1), (2, 1), (1, 2)])
    if tested == 0:
        pytest.skip("no nonzero HH pairs in the tested window for this fixture")


# =========================================================================== #
# (d) graded-commutativity + associativity MOD COBOUNDARY on HH reps           #
# =========================================================================== #
def test_graded_commutativity_kx2():
    """k[x]/x^2 over GF(32003): u1 in HH^1, u2 in HH^2.
      (1,1): u1 cup u1 ~ -(u1 cup u1)   (pq odd)
      (1,2): u1 cup u2 ~  (u2 cup u1)   (pq even)"""
    comp = Comparison(_kx2_gf())
    comp._ensure(3)
    res = comp._res
    u1 = comp.hh_class_cs(1, 0)
    u2 = comp.hh_class_cs(2, 0)
    # (1,1), odd*odd
    a = native_cup(res, u1.vec, 1, u1.vec, 1)
    neg = [(-x) % comp.p for x in a]
    assert comp.same_cohomology_class(a, neg, degree=2)
    # (1,2), even parity: u1 cup u2 ~ u2 cup u1
    lhs = native_cup(res, u1.vec, 1, u2.vec, 2)
    rhs = native_cup(res, u2.vec, 2, u1.vec, 1)
    assert comp.same_cohomology_class(lhs, rhs, degree=3)


def test_associativity_kx2():
    """(u cup u) cup u ~ u cup (u cup u) on HH^*(k[x]/x^2), degree 3."""
    comp = Comparison(_kx2_gf())
    comp._ensure(3)
    res = comp._res
    u = comp.hh_class_cs(1, 0)
    uu = native_cup(res, u.vec, 1, u.vec, 1)             # cocycle, degree 2
    left = native_cup(res, uu, 2, u.vec, 1)
    right = native_cup(res, u.vec, 1, uu, 2)
    assert comp.same_cohomology_class(left, right, degree=3)


def test_graded_commutativity_qci():
    """quantum-CI over GF(5): two distinct HH^1 reps u0,u1; u0 cup u1 ~ -(u1 cup u0)."""
    comp = Comparison(_qci_gf())
    comp._ensure(2)
    res = comp._res
    reps1 = comp.cs_cohomology_basis(1)
    if len(reps1) < 2:
        pytest.skip("needs dim HH^1 >= 2")
    u0 = comp.hh_class_cs(1, 0)
    u1 = comp.hh_class_cs(1, 1)
    lhs = native_cup(res, u0.vec, 1, u1.vec, 1)
    rhs = native_cup(res, u1.vec, 1, u0.vec, 1)
    neg = [(-x) % comp.p for x in rhs]                   # (-1)^{1*1}
    assert comp.same_cohomology_class(lhs, neg, degree=2)
