"""Plan 21 gates: the native Hochschild cap product on the CS resolution.

The cap is the HOMOLOGY-side `b·w·a` collapse of the SAME lifted diagonal Δ Plan 20
built for the cup (`resolutions_cs/diagonal.py`); no new construction. A cochain
f ∈ C^p capped with a chain z ∈ C_n lands in C_{n-p} via the sign-free op-twisted
collapse `b_c · x · b_a · f(τ) · b_mid` — f eats the degree-p FIRST factor τ of the
(p,n-p)-split of Δ_n(σ), ρ (degree n-p) survives, and the value is read against the
`hom` corner of ρ. Notation matches `cap_of_cs_classes(f, z)` (cohomology first).

Deep bucket (tests/resolutions_cs → deep); run this FILE directly:
    NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q \
        tests/resolutions_cs/test_native_cap.py -p no:cacheprovider

Gates (per the repo's verification standard; sources cited, no invented theorems):
  (a) UNIT cap  z ∩ 1 = z: `native_cap(ε, z) = z` for ε the degree-0 augmentation
      cochain (value e_v at each σ_v).  Δ's (0,n)-component is the standard AW term
      `e⊗σ_v⊗e⊗σ⊗e`, so the unit cap is the identity EXACTLY on chains (not merely mod
      boundary).  Conventions per Tamarkin–Tsygan / engine.tt_calculus:29-33
      (`1_A ∩ z = z`); the CS diagonal realises it.
  (b) IN-WINDOW ANCHOR (permanent overlap oracle): native cap == the transported
      `Comparison.cap_of_cs_classes` (Plan 14) mod boundary, on every nonzero
      HH^p × HH_n class pair with max(p,n) <= window.  The transported cap is the
      in-window ground truth (nose-tested tt_calculus cap), exactly as the transported
      cup anchored Plan 20; this PINS the sign-free `b·w·a` collapse.
  (c) MODULE IDENTITY  (z ∩ f) ∩ g ~ z ∩ (f ∪ g): HH_• is a module over the cup ring
      HH^•.  In `cap_of_cs_classes(f,z)` notation this is
      `native_cap(g, native_cap(f, z)) ~ native_cap(native_cup(f,g), z)` mod boundary —
      the native cup carries f∪g PAST the window.  (tt_calculus:32-33 module law
      `f ∩ (g ∩ z) = (g ⌣ f) ∩ z`, transcribed to the CS side.)
  (d) CAP-LEIBNIZ (exact, the sign check independent of Comparison):
      `b(f ∩ z) = (-1)^{p+1}(δf ∩ z) + (-1)^p (f ∩ b z)` exactly over GF(5), on the full
      basis grid — the homology mirror of the cup's Leibniz gate and the arbiter that
      the collapse is sign-free (all Koszul sign already lives in Δ).  Signs match
      tt_calculus test_cap_leibniz.
  (e) DEGREE EDGES: n = p lands in C_0; p > n raises ValueError (bar convention,
      tt_calculus.cap_cochain — a cap into negative degree is refused, not silently 0).
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.comparison import Comparison, CSClass
from quiverlab.resolutions_cs.cap import native_cap
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
    dom = res.dom
    out = []
    for row in M:
        acc = dom.zero()
        for j, mij in enumerate(row):
            acc = dom.add(acc, dom.mul(mij, v[j]))
        out.append(acc)
    return out


def _basis_cochains(res, deg):
    d = len(res._basis(deg, "coh"))
    for i in range(d):
        v = [res.dom.zero()] * d
        v[i] = res.dom.one()
        yield v


def _basis_chains(res, deg):
    d = len(res._basis(deg, "hom"))
    for i in range(d):
        v = [res.dom.zero()] * d
        v[i] = res.dom.one()
        yield v


def _unit_cochain(res):
    """The degree-0 augmentation cochain ε: value e_v at each vertex chain σ_v
    (represents 1_A ∈ C^0)."""
    basis = res._basis(0, "coh")
    v = [res.dom.zero()] * len(basis)
    for i, (ch, j) in enumerate(basis):
        if j == res.ar._word_index[("v", ch.o)]:
            v[i] = res.dom.one()
    return v


# =========================================================================== #
# (a) UNIT cap -- exact identity on chains                                     #
# =========================================================================== #
@pytest.mark.parametrize("mk,name,ns", [
    (_kx2, "kx2", (0, 1, 2, 3)),
    (_straddle, "straddle", (0, 1, 2)),
    (_qci, "qci", (0, 1, 2, 3)),
])
def test_unit_cap_exact(mk, name, ns):
    res = mk()
    u = _unit_cochain(res)
    for n in ns:
        for z in _basis_chains(res, n):
            capped = native_cap(res, u, 0, z, n)
            assert _toints(res, capped) == _toints(res, z), \
                f"1 ∩ z != z exactly on {name} at n={n}"


# =========================================================================== #
# (d) CAP-LEIBNIZ -- exact, the sign arbiter                                   #
# =========================================================================== #
def _cap_leibniz_holds(res, p, n):
    """b(f ∩ z) == (-1)^{p+1}(δf ∩ z) + (-1)^p (f ∩ b z) on every basis (f, z)."""
    P = res.dom.p
    s1 = 1 if (p + 1) % 2 == 0 else -1
    s2 = 1 if p % 2 == 0 else -1
    b_out = res.matrix(n - p, "hom")            # b_{n-p}: C_{n-p} -> C_{n-p-1}
    dcoh = res.matrix(p, "coh")                 # δ^p: C^p -> C^{p+1}
    bhom = res.matrix(n, "hom")                 # b_n: C_n -> C_{n-1}
    for f in _basis_cochains(res, p):
        df = _matvec(res, dcoh, f)
        for z in _basis_chains(res, n):
            lhs = _toints(res, _matvec(res, b_out, native_cap(res, f, p, z, n)))
            t1 = _toints(res, native_cap(res, df, p + 1, z, n))
            t2 = _toints(res, native_cap(res, f, p, _matvec(res, bhom, z), n - 1))
            rhs = [(s1 * a + s2 * b) % P for a, b in zip(t1, t2)]
            if lhs != rhs:
                return False, (f, z)
    return True, None


@pytest.mark.parametrize("p,n", [(1, 2), (1, 3), (2, 3), (0, 2), (0, 3)])
def test_cap_leibniz_kx2(p, n):
    ok, w = _cap_leibniz_holds(_kx2(), p, n)
    assert ok, f"cap-Leibniz failed on k[x]/x^2 at (p,n)=({p},{n}); witness {w}"


@pytest.mark.parametrize("p,n", [(1, 2), (0, 2)])
def test_cap_leibniz_straddle(p, n):
    ok, w = _cap_leibniz_holds(_straddle(), p, n)
    assert ok, f"cap-Leibniz failed on straddle at (p,n)=({p},{n}); witness {w}"


@pytest.mark.parametrize("p,n", [(1, 2), (0, 2), (1, 3)])
def test_cap_leibniz_qci(p, n):
    ok, w = _cap_leibniz_holds(_qci(), p, n)
    assert ok, f"cap-Leibniz failed on quantum-CI at (p,n)=({p},{n}); witness {w}"


# =========================================================================== #
# (e) DEGREE EDGES                                                             #
# =========================================================================== #
def test_cap_degree_edge_n_equals_p():
    """n = p lands in C_0 (a genuine degree-0 chain), never refused."""
    res = _kx2()
    u = _unit_cochain(res)
    # 1 ∩ z for z ∈ C_2, p=0, n=2 -> C_2 (already covered); here the p=n case:
    # a degree-2 cochain capped with a degree-2 chain -> C_0.
    for f in _basis_cochains(res, 2):
        for z in _basis_chains(res, 2):
            out = native_cap(res, f, 2, z, 2)
            assert len(out) == len(res._basis(0, "hom"))


def test_cap_degree_edge_p_gt_n_raises():
    """p > n (cap into negative degree) raises ValueError (bar convention,
    tt_calculus.cap_cochain / test_cap_degree_guard)."""
    res = _kx2()
    f = [res.dom.zero()] * len(res._basis(2, "coh"))
    z = [res.dom.zero()] * len(res._basis(1, "hom"))
    with pytest.raises(ValueError):
        native_cap(res, f, 2, z, 1)


# =========================================================================== #
# (b) IN-WINDOW ANCHOR -- native cap == transported cap mod boundary           #
# =========================================================================== #
def _anchor(comp, pairs):
    """native_cap == cap_of_cs_classes(transport) mod boundary on nonzero class pairs.
    Returns (#tested, #nonzero-transported) so callers can demand a substantive run."""
    maxn = max(n for _, n in pairs)
    comp._ensure(maxn + 1)
    res = comp._res
    tested = nonzero = 0
    for (p, n) in pairs:
        if max(p, n) > comp.window:
            continue
        freps = comp.cs_cohomology_basis(p)
        zreps = comp.cs_homology_basis(n)
        for i in range(len(freps)):
            for j in range(len(zreps)):
                f = comp.hh_class_cs(p, i)
                z = comp.hh_class_cs_hom(n, j)
                native = native_cap(res, f.vec, p, z.vec, n)
                transported = comp.cap_of_cs_classes(f, z)
                assert comp.same_homology_class(native, transported, degree=n - p), \
                    f"native != transported cap mod boundary at (p,n)=({p},{n}), reps ({i},{j})"
                tested += 1
                if any(x % comp.p for x in transported):
                    nonzero += 1
    return tested, nonzero


def test_anchor_native_equals_transported_kx2():
    """k[x]/x^2 over GF(32003): HH_n = k for every n, so the caps are substantive."""
    tested, nonzero = _anchor(Comparison(_kx2_gf()), [(1, 1), (1, 2), (1, 3), (2, 3)])
    assert tested > 0 and nonzero > 0


def test_anchor_native_equals_transported_qci():
    """quantum-CI over GF(5), NON-COMMUTATIVE (yx = 2xy): the distinguisher that pins
    the op-twisted `b·w·a` collapse (the sign-free `a·w·b` cohomology order would
    disagree here)."""
    comp = Comparison(_qci_gf())
    tested, nonzero = _anchor(comp, [(1, 1), (1, 3), (2, 3)])
    assert tested > 0 and nonzero > 0, "expected substantive nonzero QCI caps"


# =========================================================================== #
# (c) MODULE IDENTITY -- (z ∩ f) ∩ g ~ z ∩ (f ∪ g), NATIVE cup + cap           #
# =========================================================================== #
def _module_identity(comp, triples):
    comp._ensure(max(n for _, _, n in triples) + 1)
    res = comp._res
    tested = 0
    for (p, q, n) in triples:
        freps = comp.cs_cohomology_basis(p)
        greps = comp.cs_cohomology_basis(q)
        zreps = comp.cs_homology_basis(n)
        for fr in freps:
            for gr in greps:
                for zr in zreps:
                    f, g, z = CSClass(p, fr), CSClass(q, gr), CSClass(n, zr)
                    inner = native_cap(res, f.vec, p, z.vec, n)             # deg n-p
                    lhs = native_cap(res, g.vec, q, inner, n - p)           # deg n-p-q
                    fg = native_cup(res, f.vec, p, g.vec, q)                # deg p+q
                    rhs = native_cap(res, fg, p + q, z.vec, n)              # deg n-p-q
                    assert comp.same_homology_class(lhs, rhs, degree=n - p - q), \
                        f"(z∩f)∩g !~ z∩(f∪g) at (p,q,n)=({p},{q},{n})"
                    tested += 1
    return tested


def test_module_identity_qci():
    """quantum-CI over GF(5): HH_• a module over HH^• via the NATIVE cup and cap."""
    tested = _module_identity(Comparison(_qci_gf()),
                              [(1, 1, 2), (1, 1, 3), (2, 1, 3)])
    assert tested > 0


def test_module_identity_kx2():
    tested = _module_identity(Comparison(_kx2_gf()),
                              [(1, 1, 2), (1, 1, 3), (2, 1, 3)])
    assert tested > 0


# =========================================================================== #
# multi-vertex: comm_square (the plan-named case) + cn_3_2 (nonzero at depth)  #
# =========================================================================== #
from quiverlab.families.zoo import load_catalog, build_from_record


def _record(name, field):
    rec = next(r for r in load_catalog() if r["name"] == name)
    return build_from_record(rec, field=field)


def test_multivertex_comm_square_unit_cap():
    """comm_square kQ/(cd-ab), 1->2->4, 1->3->4 (gldim 2): the quiver is acyclic, so
    C_n^hom = 0 for n >= 1 and the only homology is HH_0 = k^4.  The multi-vertex UNIT
    cap `ε ∩ z = z` (four vertices, four-corner collapse) is the substantive check —
    native == z EXACTLY on every C_0 basis chain, and == the transported cap on the
    HH^0 = center = k·1 class."""
    comp = Comparison(_record("comm_square", GF(5)))
    comp._ensure(1)
    res = comp._res
    assert len(res._basis(1, "hom")) == 0, "comm_square is acyclic: C_1^hom must be 0"
    u = _unit_cochain(res)
    d0 = len(res._basis(0, "hom"))
    assert d0 == 4, f"HH_0 = C_0 should be k^4 for comm_square, got {d0}"
    for z in _basis_chains(res, 0):
        assert _toints(res, native_cap(res, u, 0, z, 0)) == _toints(res, z), \
            "multi-vertex unit cap must be the identity on C_0"
    # the HH^0 generator is the center = the unit; native cap == transported cap.
    f0 = comp.hh_class_cs(0, 0)
    z0 = comp.hh_class_cs_hom(0, 0)
    native = native_cap(res, f0.vec, 0, z0.vec, 0)
    transported = comp.cap_of_cs_classes(f0, z0)
    assert comp.same_homology_class(native, transported, degree=0)


def test_multivertex_cn_3_2_cap_leibniz():
    """cn_3_2 (cyclic Nakayama, 3 vertices, HH_ = (3,0,1,1)): the genuinely
    non-vanishing multi-vertex cap.  The exact cap-Leibniz gate at (p,n) = (1,3) and
    (0,2) — where native caps are nonzero (empty-corner + cross-vertex bridges) —
    substantively exercises the multi-vertex `b·w·a` collapse."""
    res = _record_res("cn_3_2", GF(5), md=4)
    for (p, n) in [(1, 3), (0, 2)]:
        ok, w = _cap_leibniz_holds(res, p, n)
        assert ok, f"cn_3_2 cap-Leibniz failed at (p,n)=({p},{n}); witness {w}"
    # substantive: at least one native cap is a nonzero chain
    saw = False
    for f in _basis_cochains(res, 1):
        for z in _basis_chains(res, 3):
            if any(_toints(res, native_cap(res, f, 1, z, 3))):
                saw = True
    assert saw, "expected a nonzero multi-vertex cap on cn_3_2"


def _record_res(name, field, md):
    A = _record(name, field)
    from quiverlab.resolutions_cs.build import reduction_system_of
    return ChouhySolotarResolution(A, reduction_system_of(A), max_degree=md)
