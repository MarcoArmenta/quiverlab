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
# MULTI-PRIME + BYTE-REPRODUCIBILITY                                           #
# =========================================================================== #
@pytest.mark.parametrize("prime", [2, 3, 5, 32003])
def test_cap_multiprime_kx2(prime):
    """The cap is domain-generic and prime-independent: the exact unit cap and
    cap-Leibniz hold on k[x]/x^2 over several primes -- including char 2, where the
    Koszul signs (−1)^p collapse to +1 (a genuine edge for the sign bookkeeping)."""
    res = _res(["x*x"], {"x": (1, 1)}, 6, p=prime)
    u = _unit_cochain(res)
    for n in (1, 2, 3):
        for z in _basis_chains(res, n):
            assert _toints(res, native_cap(res, u, 0, z, n)) == _toints(res, z), \
                f"unit cap != z on kx2/GF({prime}) at n={n}"
    for (p, n) in [(1, 2), (1, 3), (2, 3)]:
        ok, w = _cap_leibniz_holds(res, p, n)
        assert ok, f"cap-Leibniz failed on kx2/GF({prime}) at ({p},{n}); witness {w}"


def test_cap_byte_reproducible_qci():
    """native_cap is byte-reproducible: two FRESH CS resolutions give identical cap
    vectors (follows from Δ's byte-reproducibility, Plan 17; pinned here for the cap).
    Uses quantum-CI at n=3, where the diagonal's correction solve has genuine nullspace
    freedom -- the case the Plan-17 canonicalization exists for -- and asserts at least
    one nonzero cap so the match is substantive."""
    ra, rb = _qci(), _qci()
    saw = False
    for fi in range(len(ra._basis(1, "coh"))):
        f = [ra.dom.zero()] * len(ra._basis(1, "coh"))
        f[fi] = ra.dom.one()
        for zi in range(len(ra._basis(3, "hom"))):
            z = [ra.dom.zero()] * len(ra._basis(3, "hom"))
            z[zi] = ra.dom.one()
            va = _toints(ra, native_cap(ra, f, 1, z, 3))
            vb = _toints(rb, native_cap(rb, f, 1, z, 3))
            assert va == vb, f"native_cap not byte-reproducible at (f{fi}, z{zi})"
            if any(va):
                saw = True
    assert saw, "byte-repro check was vacuous -- every cap vanished"


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


# =========================================================================== #
# PAST-WINDOW DELIVERY -- the headline: cap routes to the native CS diagonal   #
# once max(p,n) exceeds the (deliberately tiny) comparison window.             #
# =========================================================================== #
def test_native_cap_past_window_kx2():
    """Past the window `Comparison.cap_of_cs_classes` no longer refuses -- it routes to
    the native CS diagonal (`engine='auto'`).  k[x]/x^2 over GF(5), window shrunk to 0:
    f ∈ HH^1 capped with z ∈ HH_1 (max(1,1) = 1 > 0) computes natively and is a nonzero
    HH_0 class; today's transport route would raise NotImplementedError."""
    comp = Comparison(_kx2_gf5(), max_cells=8)
    assert comp.window == 0, f"expected tiny window 0, got {comp.window}"
    f = comp.hh_class_cs(1, 0)
    z = comp.hh_class_cs_hom(1, 0)
    capped = comp.cap_of_cs_classes(f, z)                  # deg (1,1) > window: native
    zero0 = [0] * len(comp._res._basis(0, "hom"))
    assert not comp.same_homology_class(capped, zero0, degree=0), \
        "f ∩ z (HH^1 ∩ HH_1) must be a nonzero HH_0 class"


def test_native_cap_bridge_to_longer_transport_kx2():
    """Bridge oracle: a cap past a TINY window computed NATIVELY equals -- mod
    boundary -- the SAME cap computed by TRANSPORT on a wider-window instance.
    k[x]/x^2 over GF(5): f ∈ HH^1 ∩ z ∈ HH_3 = a nonzero HH_2 class, so the match is not
    a vacuous boundary coincidence.  Plan-17 canonicalization makes the two instances'
    CS bases identical -- asserted ELEMENT-WISE (chain word + corner index), not merely
    equal length."""
    small = Comparison(_kx2_gf5(), max_cells=8)            # window 0 (native past it)
    big = Comparison(_kx2_gf5())                           # default window (transport)
    assert small.window == 0 and big.window >= 3, (small.window, big.window)

    fs, zs = small.hh_class_cs(1, 0), small.hh_class_cs_hom(3, 0)
    fb, zb = big.hh_class_cs(1, 0), big.hh_class_cs_hom(3, 0)

    native = small.cap_of_cs_classes(fs, zs)              # max(1,3)=3 > 0 -> native, deg 2
    transported = big.cap_of_cs_classes(fb, zb)           # max(1,3)=3 <= window -> transport

    small_basis = [(ch.word, j) for ch, j in small._res._basis(2, "hom")]
    big_basis = [(ch.word, j) for ch, j in big._res._basis(2, "hom")]
    assert small_basis == big_basis, \
        "Plan-17 canonicalization must give element-wise identical CS hom bases"
    assert big.same_homology_class(native, transported, degree=2), \
        "native (tiny window) != longer transport (wide window) mod boundary"
    zero2 = [0] * len(big._res._basis(2, "hom"))
    assert not big.same_homology_class(transported, zero2, degree=2), \
        "the bridged class must be nonzero -- not a vacuous match"


def test_cap_engine_selector_kx2():
    """`engine=` forces a route: 'native' computes at any degree; 'transport' keeps the
    window refusal; an invalid engine is a ValueError naming the three options."""
    comp = Comparison(_kx2_gf5(), max_cells=8)             # window 0
    f = comp.hh_class_cs(1, 0)
    z = comp.hh_class_cs_hom(1, 0)
    forced_native = comp.cap_of_cs_classes(f, z, engine="native")
    assert isinstance(forced_native, list)
    with pytest.raises(NotImplementedError):
        comp.cap_of_cs_classes(f, z, engine="transport")   # transport still refuses
    with pytest.raises(ValueError):
        comp.cap_of_cs_classes(f, z, engine="bogus")


def test_cap_auto_in_window_is_transport_byte_unchanged():
    """`engine='auto'` in-window is byte-for-byte the transported route (Plan-14
    behavior), and both equal `engine='transport'` there.  k[x]/x^2 over GF(32003),
    wide window; (1,1) is well inside it."""
    comp = Comparison(_kx2_gf())
    f = comp.hh_class_cs(1, 0)
    z = comp.hh_class_cs_hom(1, 0)
    auto = comp.cap_of_cs_classes(f, z, engine="auto")
    transport = comp.cap_of_cs_classes(f, z, engine="transport")
    assert auto == transport, "auto in-window must be byte-identical to transport"


# =========================================================================== #
# DEEP PINS -- the past-window cap headline at real depth (shared Delta_4      #
# fixture), and the domain-generic (QQ) smoke.                                 #
# =========================================================================== #
def test_deep_qci_cap_past_window(qci_gf5_diag4):
    """QCI over GF(5), PAST a zero window at real depth, reusing the SESSION-scoped
    Delta_4 (`qci_gf5_diag4`, shared with the cup deep pin -- Plan 21 review item).
    window=0 forces every cap onto the native route.  HH_•(QCI/GF5) = (3,2,3,4,3,2).

    Three pins:
      * exact CAP-LEIBNIZ at (p,n) = (1,3) and (2,3) -- the sign arbiter at real depth,
        past the window (reuses the cached Delta_3 inside Delta_4);
      * one NONZERO cap class: HH^1 ∩ HH_3 -> HH_2 (HH_2 = k^3), routed native;
      * the MODULE IDENTITY (z∩f)∩g ~ z∩(f∪g) with the NATIVE cup for f∪g, all past
        window=0 -- HH_• as a module over HH^• with no bar object anywhere."""
    comp = qci_gf5_diag4
    res = comp._res

    # (i) exact cap-Leibniz at (1,3), (2,3) -- past window=0.
    for (p, n) in [(1, 3), (2, 3)]:
        ok, w = _cap_leibniz_holds(res, p, n)
        assert ok, f"cap-Leibniz failed on QCI at (p,n)=({p},{n}); witness {w}"

    # (ii) one nonzero cap class HH^1 ∩ HH_3 -> HH_2 (past window=0 -> native).
    reps1 = comp.cs_cohomology_basis(1)
    reps3 = comp.cs_homology_basis(3)
    assert len(reps1) >= 1 and len(reps3) >= 1
    zero2 = [0] * len(res._basis(2, "hom"))
    saw = False
    for i in range(len(reps1)):
        for j in range(len(reps3)):
            capped = comp.cap_of_cs_classes(comp.hh_class_cs(1, i),
                                            comp.hh_class_cs_hom(3, j))
            if not comp.same_homology_class(capped, zero2, degree=2):
                saw = True
    assert saw, "expected a nonzero HH^1 ∩ HH_3 -> HH_2 cap class past the window"

    # (iii) module identity via the NATIVE cup, past window=0.
    f = comp.hh_class_cs(1, 0)
    g = comp.hh_class_cs(1, 1) if len(reps1) >= 2 else comp.hh_class_cs(1, 0)
    z = comp.hh_class_cs_hom(3, 0)
    inner = native_cap(res, f.vec, 1, z.vec, 3)                     # deg 2
    lhs = native_cap(res, g.vec, 1, inner, 2)                       # deg 1
    fg = native_cup(res, f.vec, 1, g.vec, 1)                        # deg 2 (native)
    rhs = native_cap(res, fg, 2, z.vec, 3)                          # deg 1
    assert comp.same_homology_class(lhs, rhs, degree=1), \
        "(z∩f)∩g !~ z∩(f∪g) past the window (native cup + cap)"


def _dom_eq(res, a, b):
    """Exact equality of two domain vectors (a - b == 0 componentwise)."""
    dom = res.dom
    return len(a) == len(b) and all(dom.is_zero(dom.sub(x, y)) for x, y in zip(a, b))


def _dom_axpy(res, s1, v1, s2, v2):
    """s1*v1 + s2*v2 in the domain, s1,s2 in {+1,-1}."""
    dom = res.dom
    one, none = dom.one(), dom.neg(dom.one())
    c1, c2 = (one if s1 == 1 else none), (one if s2 == 1 else none)
    return [dom.add(dom.mul(c1, a), dom.mul(c2, b)) for a, b in zip(v1, v2)]


def test_deep_qq_cap_smoke():
    """The native cap is DOMAIN-GENERIC.  On k[x]/x^2 over QQ (exact rationals;
    Comparison is GF(p)-gated and NOT used) assert the two exact chain-level identities
    directly over the QQ domain: the unit cap `1 ∩ z = z` and the cap-Leibniz
    `b(f∩z) = (-1)^{p+1}(δf∩z) + (-1)^p(f∩bz)`."""
    from quiverlab.fields import QQ
    Q = Quiver([1], {"x": (1, 1)})
    A = Q.algebra(relations=["x*x"], field=QQ)
    rs = build_reduction_system(Q, ["x*x"], QQ)
    res = ChouhySolotarResolution(A, rs, max_degree=5)
    assert res.dom.name == "QQ", f"expected the QQ domain, got {res.dom.name}"

    # unit cap exact over QQ
    u = _unit_cochain(res)
    for n in (1, 2, 3):
        for z in _basis_chains(res, n):
            assert _dom_eq(res, native_cap(res, u, 0, z, n), z), \
                f"QQ unit cap 1 ∩ z != z at n={n}"

    # cap-Leibniz exact over QQ (s1 = (-1)^{p+1}, s2 = (-1)^p)
    for (p, n) in [(1, 2), (1, 3), (2, 3)]:
        s1 = 1 if (p + 1) % 2 == 0 else -1
        s2 = 1 if p % 2 == 0 else -1
        b_out = res.matrix(n - p, "hom")
        dcoh = res.matrix(p, "coh")
        bhom = res.matrix(n, "hom")
        for f in _basis_cochains(res, p):
            df = _matvec(res, dcoh, f)
            for z in _basis_chains(res, n):
                lhs = _matvec(res, b_out, native_cap(res, f, p, z, n))
                t1 = native_cap(res, df, p + 1, z, n)
                t2 = native_cap(res, f, p, _matvec(res, bhom, z), n - 1)
                assert _dom_eq(res, lhs, _dom_axpy(res, s1, t1, s2, t2)), \
                    f"QQ cap-Leibniz failed at (p,n)=({p},{n})"
