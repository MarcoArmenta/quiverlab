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


def _kx2_gf5():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(5))


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


# =========================================================================== #
# (Task 4) PAST-WINDOW DELIVERY -- the headline: cup routes to the native CS   #
# diagonal once max(p,q) exceeds the (deliberately tiny) comparison window.    #
# =========================================================================== #
def test_native_cup_past_window_kx2():
    """Past the window, `Comparison.cup_of_cs_classes` no longer refuses -- it routes
    to the native CS diagonal (Plan 20 Task 4).  k[x]/x^2 over GF(5), char 5 != 2:
    the HH ring is k[alpha]/alpha^2 (x) k[beta], |alpha|=1 (odd), |beta|=2 (even), so

        alpha cup alpha  (deg 2, ODD square)  ~ 0   (a coboundary; forced by graded
                                                     commutativity in char != 2)
        beta  cup beta   (deg 4, EVEN square) !~ 0  (nonzero -- beta is the period-2
                                                     periodicity element)

    Both cups have max(p,q) > window == 0, so today's transport raises
    NotImplementedError; after wiring they compute natively.  Class dims are pinned
    via same_cohomology_class against zero (never bytes)."""
    comp = Comparison(_kx2_gf5(), max_cells=8)
    # _default_window: eff=max(m-1,2)=2, largest n with m*eff^(2n+1) <= 8:
    #   n=0 -> 2*2 = 4 <= 8 (ok); n=1 -> 2*8 = 16 > 8 -> window 0.
    assert comp.window == 0, f"expected tiny window 0, got {comp.window}"
    u1 = comp.hh_class_cs(1, 0)
    u2 = comp.hh_class_cs(2, 0)

    a = comp.cup_of_cs_classes(u1, u1)                     # deg 2 > window: native now
    zero2 = [0] * len(comp._res._basis(2, "coh"))
    assert comp.same_cohomology_class(a, zero2, degree=2), \
        "alpha cup alpha must be a coboundary (odd square, char != 2)"

    b = comp.cup_of_cs_classes(u2, u2)                     # deg 4 > window: native now
    zero4 = [0] * len(comp._res._basis(4, "coh"))
    assert not comp.same_cohomology_class(b, zero4, degree=4), \
        "beta cup beta must be a nonzero class (even periodicity square)"


def test_native_cup_bridge_to_longer_transport_kx2():
    """Bridge oracle: a cup past a TINY window computed NATIVELY equals -- mod
    coboundary -- the SAME cup computed by TRANSPORT on a wider-window instance
    (Plan 20 Task 4).  k[x]/x^2 over GF(5): u1 in HH^1, u2 in HH^2, and
    u1 cup u2 = alpha*beta in HH^3 is nonzero, so the match is not a vacuous
    coboundary coincidence.

    The two instances build SEPARATE CS resolutions; Plan-17 canonicalization makes
    their CS bases identical, so the coordinate vectors are directly comparable --
    asserted loudly (equal basis lengths) before the comparison."""
    small = Comparison(_kx2_gf5(), max_cells=8)            # window 0 (native past it)
    big = Comparison(_kx2_gf5())                           # default window 9 (transport)
    assert small.window == 0 and big.window >= 3, (small.window, big.window)

    u1s, u2s = small.hh_class_cs(1, 0), small.hh_class_cs(2, 0)
    u1b, u2b = big.hh_class_cs(1, 0), big.hh_class_cs(2, 0)

    native = small.cup_of_cs_classes(u1s, u2s)            # max(1,2)=2 > 0 -> native, deg 3
    transported = big.cup_of_cs_classes(u1b, u2b)         # max(1,2)=2 <= 9 -> transport

    small_basis = [(ch.word, j) for ch, j in small._res._basis(3, "coh")]
    big_basis = [(ch.word, j) for ch, j in big._res._basis(3, "coh")]
    assert small_basis == big_basis, \
        "Plan-17 canonicalization must give ELEMENT-WISE identical CS bases across " \
        "instances (chain word + corner index), so the coordinate vectors are directly " \
        "comparable"
    assert big.same_cohomology_class(native, transported, degree=3), \
        "native (tiny window) != longer transport (wide window) mod coboundary"
    zero3 = [0] * len(big._res._basis(3, "coh"))
    assert not big.same_cohomology_class(transported, zero3, degree=3), \
        "the bridged class (alpha*beta) must be nonzero -- not a vacuous match"


def test_cup_engine_selector_kx2():
    """`engine=` forces a route: 'native' computes at any degree; 'transport' keeps
    the window refusal; an invalid engine is a ValueError naming the three options."""
    comp = Comparison(_kx2_gf5(), max_cells=8)             # window 0
    u1 = comp.hh_class_cs(1, 0)
    # in-window is degree <= 0 for a cup, so (1,1) is past the window either way.
    forced_native = comp.cup_of_cs_classes(u1, u1, engine="native")
    assert isinstance(forced_native, list)
    with pytest.raises(NotImplementedError):
        comp.cup_of_cs_classes(u1, u1, engine="transport")   # transport still refuses
    with pytest.raises(ValueError):
        comp.cup_of_cs_classes(u1, u1, engine="bogus")


# =========================================================================== #
# (Task 5) DEEP PINS -- the past-window headline at real depth, the first      #
# genuinely multi-vertex diagonal, and the domain-generic (QQ) smoke.          #
#                                                                              #
# Controller-adjudicated small-window pattern (Task-4 style): the plan's "2-3  #
# degrees past each fixture's DEFAULT window" is infeasible (default windows    #
# are ~9-10; Delta(4) costs ~25-30s on QCI and is out of budget on the         #
# straddle).  Instead a TINY/ZERO comparison window makes "past-window" land at #
# LOW absolute degree, so the heavy Delta build is the only cost and stays      #
# inside the deep bucket.                                                       #
# =========================================================================== #
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace
from quiverlab.resolutions_cs.diagonal import TensorComplex, diagonal
from quiverlab.resolutions_cs.comparison import CSClass
from quiverlab.families.zoo import load_catalog, build_from_record


def _cw(ch):
    """Storage word of a chain (== d_terms targets / TensorComplex._chain_word)."""
    return ("__v__", ch.o) if ch.degree == 0 else ch.word


def _comm_square(field):
    """The Plan-18 multi-vertex zoo record kQ/(cd-ab): commutative square
    1->2->4, 1->3->4 (arrows a,b,c,d), gldim 2, dim 9."""
    rec = next(r for r in load_catalog() if r["name"] == "comm_square")
    return build_from_record(rec, field=field)


def test_deep_qci_past_window_leibniz_and_nonzero_class():
    """QCI over GF(5), PAST a zero window at real depth (needs Delta_4, ~28s --
    deep-bucket-acceptable).  window=0 forces every cup (max(p,q) >= 1 > 0) onto the
    native route, so "past-window" lands at low absolute degree.

    Two pins Task 3 budgeted out (they need Delta_4):
      * LEIBNIZ at (p,q) = (2,1) and (1,2) -- the sign arbiter at total degree 3,
        exact over GF(5) on the full basis grid (reuses the cached Delta_4).
      * one NONZERO product class: cup(HH^1[0], HH^1[1]) is a nonzero HH^2 class
        (the off-diagonal (1,1)-cup; the diagonal squares vanish by graded
        commutativity, char 5 != 2).

    Honest scope note: HH^*(QCI/GF5) = (2,2,1,0,2,4,...), so HH^3 = 0 -- there is NO
    nonzero degree-3 class; the nonzero product class lives at degree 2 while the
    Leibniz pins are the degree-3 (Delta_4) ones.  A degree-3 native cup is
    additionally shown to COMPUTE past the window (transport would raise) although
    its class is forced to zero."""
    comp = Comparison(_qci_gf(), window=0)
    comp._ensure(4)                        # CS resolution to max_degree 5
    res = comp._res
    diagonal(res, 4)                       # the ~28s build; cached for every cup below

    # (i) Leibniz at (2,1) and (1,2) -- reuse the cached Delta_4.
    ok, w = _leibniz_holds(res, 2, 1)
    assert ok, f"Leibniz failed on QCI at (2,1); witness {w}"
    ok, w = _leibniz_holds(res, 1, 2)
    assert ok, f"Leibniz failed on QCI at (1,2); witness {w}"

    # (ii) one nonzero product class, PAST window=0 (auto routes native).
    reps1 = comp.cs_cohomology_basis(1)
    assert len(reps1) >= 2, "QCI HH^1 should be 2-dimensional"
    u0, u1 = comp.hh_class_cs(1, 0), comp.hh_class_cs(1, 1)
    prod = comp.cup_of_cs_classes(u0, u1)              # max(1,1)=1 > 0 -> native, deg 2
    zero2 = [0] * len(res._basis(2, "coh"))
    assert not comp.same_cohomology_class(prod, zero2, degree=2), \
        "cup(HH^1[0], HH^1[1]) must be a nonzero HH^2 class"

    # (iii) a degree-3 native cup COMPUTES past the window (HH^3 = 0 -> zero class).
    u2 = comp.hh_class_cs(2, 0)
    c3 = comp.cup_of_cs_classes(u0, u2)                # deg 3 > 0 -> native (no refusal)
    zero3 = [0] * len(res._basis(3, "coh"))
    assert comp.same_cohomology_class(c3, zero3, degree=3)


def test_deep_comm_square_multivertex_diagonal():
    """comm_square (Plan-18 multi-vertex zoo, GF(5)): the FIRST genuine multi-vertex
    exercise of the lifted diagonal -- the empty-corner branch a single-vertex
    algebra never reaches (a Task-2 review carry-forward).

    A = kQ/(cd - ab), Q the commutative square 1->2->4, 1->3->4; gldim 2, so the CS
    resolution is finite: S(0..2) nonempty, S(3) = empty.  Two validations:

      * Delta EXISTS and CHAIN-MAPS at every built degree (n = 0,1,2): for every
        sigma in S(n), apply_tensor_d(n, Delta_n(sigma)) == zeta(sigma), and zeta is a
        d^{PxP}-cycle for n >= 2 -- the multi-vertex three-slot corner bookkeeping
        with empty corners, plus the Delta_0 base shape (one generator per vertex).
      * native cup == transported cup at degree 2, EXACTLY as cochains (engine
        "native" vs "transport" on the SAME in-window instance; window >= 2 so both
        routes run).  HH^2(comm_square) = 0 collapses all of C^2 onto coboundaries,
        so the class-level test would be VACUOUS here -- the EXACT cochain match plus
        a nonzero native cochain are the substantive multi-vertex cross-check.
    """
    A = _comm_square(GF(5))
    comp = Comparison(A)
    assert comp.window >= 2, f"expected a wide window, got {comp.window}"
    comp._ensure(3)                                   # resolution to max_degree 4
    res = comp._res
    assert len(res.ss.S(3)) == 0, "comm_square (gldim 2) should have S(3) empty"

    tc = TensorComplex(res)
    d0 = tc.diagonal(0)
    assert set(d0) == {("__v__", v) for v in res.rs.quiver.vertices}, \
        "Delta_0 must have exactly one generator per vertex (multi-vertex base shape)"
    for n in range(1, 3):                             # top nonvanishing degree is 2
        diag_n = tc.diagonal(n)
        prev = tc.diagonal(n - 1)
        for sigma in res.ss.S(n):
            zeta = tc._zeta(n, sigma, prev)
            assert tc.apply_tensor_d(n, diag_n[_cw(sigma)]) == zeta, \
                f"comm_square chain-map failed at n={n}, sigma={sigma.word}"
            if n >= 2:
                assert tc.apply_tensor_d(n - 1, zeta) == {}, \
                    f"zeta not a d^(PxP)-cycle at n={n}, sigma={sigma.word}"

    # native cup == transported cup (multi-vertex), exact cochains at degree 2.
    dom = comp.dom
    d1 = comp._dcs(1)
    Z1 = [[int(x) % comp.p for x in v]
          for v in nullspace([[dom.coerce(x) for x in row] for row in d1], dom)]
    assert len(Z1) >= 1, "comm_square should have degree-1 cocycles"
    saw_nonzero = False
    for i in range(len(Z1)):
        for j in range(len(Z1)):
            u, v = CSClass(1, Z1[i]), CSClass(1, Z1[j])
            nat = [x % comp.p for x in comp.cup_of_cs_classes(u, v, engine="native")]
            tr = [x % comp.p for x in comp.cup_of_cs_classes(u, v, engine="transport")]
            assert nat == tr, \
                f"native != transported cup on comm_square at (z{i},z{j}): {nat} vs {tr}"
            if any(nat):
                saw_nonzero = True
    assert saw_nonzero, \
        "every native cup vanished -- the multi-vertex diagonal produced no content"


def test_deep_qq_diagonal_smoke():
    """The diagonal is DOMAIN-GENERIC.  On k[x]/x^2 over QQ (exact rationals;
    Comparison is GF(p)-gated and NOT used here) build Delta to n=3 directly and
    assert the chain-map identity and the zeta-cycle property EXACTLY over QQ."""
    Q = Quiver([1], {"x": (1, 1)})
    A = Q.algebra(relations=["x*x"], field=QQ)
    rs = build_reduction_system(Q, ["x*x"], QQ)
    res = ChouhySolotarResolution(A, rs, max_degree=5)
    assert res.dom.name == "QQ", f"expected the QQ domain, got {res.dom.name}"
    tc = TensorComplex(res)
    for n in range(4):
        diag_n = tc.diagonal(n)
        prev = tc.diagonal(n - 1) if n >= 1 else None
        for sigma in res.ss.S(n):
            dn = diag_n[_cw(sigma)]
            if n >= 1:
                zeta = tc._zeta(n, sigma, prev)
                assert tc.apply_tensor_d(n, dn) == zeta, \
                    f"QQ chain-map identity failed at n={n}, sigma={sigma.word}"
                if n >= 2:
                    assert tc.apply_tensor_d(n - 1, zeta) == {}, \
                        f"QQ zeta not a d^(PxP)-cycle at n={n}, sigma={sigma.word}"
