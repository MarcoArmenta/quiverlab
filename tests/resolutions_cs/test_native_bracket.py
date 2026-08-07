"""Plan 51 Task B/C: native_bracket -- the CS Gerstenhaber bracket past the bar
window, any Domain -- and Comparison.bracket_of_cs_classes engine dispatch.

Deep bucket (tests/resolutions_cs -> deep).  The arbiter chain (Plan 51 DD4):
  * DESCENT (self-cert): [f,g] of cocycles is a cocycle (delta[f,g]=0 exactly).
  * ANTISYMMETRY (consistency check, tautological for the relative minus -- NOT a
    sign oracle by itself): [f,g] ~ -(-1)^{(p-1)(q-1)}[g,f] mod coboundary.
  * THE ANCHOR (crossengine, THE sign-fixing oracle): in-window, over GF(p),
    native_bracket == the transported Comparison.bracket_of_cs_classes mod
    coboundary, on every nonzero HH^p x HH^q pair.  The transported route is the
    classical Gerstenhaber bar bracket (an INDEPENDENT method), so agreement fixes
    the implemented sign; a FLIPPED sign breaks it (asserted -> non-vacuity).

PLAN CORRECTION (folded into research R1, dated 2026-08-07): the plan's B.2/MINOR-a
names QuantumCI (2,2) as the odd-exponent discriminator, but HH^3(QCI/GF5)=0
(dims [2,2,1,0,2,4]), so the (2,2) bracket lands in the ZERO space and is VACUOUS
for the sign.  The genuine non-vacuous odd-exponent anchor is **QCI (2,4) -> HH^5**
(dim 4, (p-1)(q-1)=1*3=3 ODD): the bracket class is nonzero, native == transported,
and the flipped sign disagrees (test_odd_exponent_sign_anchor_qci_2_4).
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.comparison import Comparison
from quiverlab.resolutions_cs.bracket import native_bracket, _eval_cochain_on_pelt
from quiverlab.resolutions_cs.homotopy_lifting import homotopy_lifting
from quiverlab.resolutions_cs.cup import _cochain_evaluator

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_selfcert]


# --------------------------------------------------------------------------- #
# fixtures                                                                     #
# --------------------------------------------------------------------------- #
def _res(rels, arrows, md, field=None):
    f = GF(5) if field is None else field
    Q = Quiver([1], arrows)
    A = Q.algebra(relations=rels, field=f)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, f), max_degree=md)


def _kx2(md=6, field=None):
    return _res(["x*x"], {"x": (1, 1)}, md, field)


def _qci(md=4):
    return _res(["x*x", "y*y", "y*x - 2*x*y"], {"x": (1, 1), "y": (1, 1)}, md)


def _kx2_gf():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(32003))


def _qci_gf():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=GF(5))


def _cocycles(res, n):
    dom = res.dom
    M = res.matrix(n, "coh")
    if not M:
        d = len(res._basis(n, "coh"))
        return [[dom.one() if i == j else dom.zero() for j in range(d)] for i in range(d)]
    return nullspace([[dom.coerce(x) for x in row] for row in M], dom)


def _dcup(res, n, vec):
    dom = res.dom
    M = res.matrix(n, "coh")
    return [sum((M[r][j] * vec[j] for j in range(len(vec))), dom.zero())
            for r in range(len(M))] if M else []


def _bracket_flipped(res, f, p, g, q):
    """The SAME assembly as native_bracket but with the bracket sign FLIPPED
    (multiply the second term's sign by -1) -- the sign-discrimination control."""
    ar, dom = res.ar, res.dom
    n = p + q - 1
    psi_f = homotopy_lifting(res, f, p)
    psi_g = homotopy_lifting(res, g, q)
    f_at = _cochain_evaluator(res, f, p)
    g_at = _cochain_evaluator(res, g, q)
    # native uses -(-1)^{(p-1)(q-1)}; flip: -(-1)^{(p-1)(q-1)+1}
    sign = dom.neg(dom.one()) if ((p - 1) * (q - 1)) % 2 == 0 else dom.one()
    ob = res._basis(n, "coh")
    pos = {(ch.word, j): i for i, (ch, j) in enumerate(ob)}
    out = [dom.zero()] * len(ob)
    for sigma in res.ss.S(n):
        vfg = _eval_cochain_on_pelt(res, f_at, psi_g.apply(n, sigma))
        vgf = _eval_cochain_on_pelt(res, g_at, psi_f.apply(n, sigma))
        for j in ar.corner(sigma.o, sigma.t, "coh"):
            v = dom.add(vfg[j], dom.neg(dom.mul(sign, vgf[j])))
            if not dom.is_zero(v):
                out[pos[(sigma.word, j)]] = v
    return out


# =========================================================================== #
# B.1 -- descent, antisymmetry, and the p,q>=1 guard                           #
# =========================================================================== #
@pytest.mark.parametrize("mk,name,pairs", [
    (_kx2, "kx2", [(1, 1), (1, 2), (2, 1), (2, 2)]),
    (_qci, "qci", [(1, 1), (1, 2), (2, 1)])])
def test_bracket_descent(mk, name, pairs):
    """[f,g] of cocycles is a cocycle: delta[f,g] = 0 exactly (over GF(5))."""
    res = mk()
    dom = res.dom
    for (p, q) in pairs:
        for f in _cocycles(res, p):
            for g in _cocycles(res, q):
                br = native_bracket(res, f, p, g, q)
                db = _dcup(res, p + q - 1, br)
                assert all(dom.is_zero(dom.coerce(x)) for x in db), \
                    f"descent failed on {name} at (p,q)=({p},{q})"


def test_bracket_antisymmetry_mod_coboundary():
    """[f,g] ~ -(-1)^{(p-1)(q-1)}[g,f] mod coboundary (a consistency check; the
    RELATIVE minus is forced by antisymmetry -- the flipped assembly is NOT
    antisymmetric, asserted in test_flip_breaks_antisymmetry)."""
    comp = Comparison(_kx2_gf())
    comp._ensure(4)
    res = comp._res
    for (p, q) in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        for f in _cocycles(res, p):
            for g in _cocycles(res, q):
                fg = native_bracket(res, f, p, g, q)
                gf = native_bracket(res, g, q, f, p)
                s = 1 if ((p - 1) * (q - 1)) % 2 == 0 else -1
                neg = [(-s * x) % comp.p for x in gf]
                assert comp.same_cohomology_class(fg, neg, degree=p + q - 1), \
                    f"antisymmetry failed at (p,q)=({p},{q})"


def test_bracket_degree_zero_out_of_scope():
    res = _kx2()
    d0 = len(res._basis(0, "coh"))
    d1 = len(res._basis(1, "coh"))
    f0 = [res.dom.one()] + [res.dom.zero()] * (d0 - 1)
    g1 = [res.dom.one()] + [res.dom.zero()] * (d1 - 1)
    with pytest.raises(ValueError):
        native_bracket(res, f0, 0, g1, 1)
    with pytest.raises(ValueError):
        native_bracket(res, g1, 1, f0, 0)


# =========================================================================== #
# B.2 -- the IN-WINDOW ANCHOR (crossengine, the sign-fixing oracle)            #
# =========================================================================== #
pytestmark_anchor = [pytest.mark.oracle_crossengine, pytest.mark.oracle_selfcert]


def _anchor(comp, pairs):
    comp._ensure(max(p + q for p, q in pairs) + 2)
    res = comp._res
    tested = 0
    for (p, q) in pairs:
        if max(p, q) > comp.window:
            continue
        repp = comp.cs_cohomology_basis(p)
        repq = comp.cs_cohomology_basis(q)
        for i in range(len(repp)):
            for j in range(len(repq)):
                u, v = comp.hh_class_cs(p, i), comp.hh_class_cs(q, j)
                native = native_bracket(res, u.vec, p, v.vec, q)
                transported = comp.bracket_of_cs_classes(u, v)
                assert comp.same_cohomology_class(native, transported, degree=p + q - 1), \
                    f"native != transported bracket at (p,q)=({p},{q}), reps ({i},{j})"
                tested += 1
    return tested


@pytest.mark.oracle_crossengine
def test_anchor_native_equals_transported_kx2():
    """k[x]/x^2 over GF(32003): HH^n = k for every n, so every (p,q) contributes."""
    tested = _anchor(Comparison(_kx2_gf()), [(1, 1), (1, 2), (2, 1), (2, 2)])
    assert tested > 0


@pytest.mark.oracle_crossengine
def test_anchor_native_equals_transported_qci():
    """quantum-CI over GF(5): HH = [2,2,1,0,...]; (1,1)->HH^1, (1,2)/(2,1)->HH^2."""
    comp = Comparison(_qci_gf())
    tested = _anchor(comp, [(1, 1), (1, 2), (2, 1)])
    assert tested > 0


@pytest.mark.oracle_crossengine
def test_flip_breaks_anchor_kx2():
    """NON-VACUITY of the anchor as a sign oracle: at k[x]/x^2 (1,2) (even exponent,
    a nonzero HH^2 class) the CORRECT-sign native bracket matches transport while the
    FLIPPED-sign assembly does NOT -- so the anchor genuinely constrains the sign."""
    comp = Comparison(_kx2_gf())
    comp._ensure(4)
    res = comp._res
    u, v = comp.hh_class_cs(1, 0), comp.hh_class_cs(2, 0)
    tr = comp.bracket_of_cs_classes(u, v)                 # (1,2) -> HH^2
    correct = native_bracket(res, u.vec, 1, v.vec, 2)
    flipped = _bracket_flipped(res, u.vec, 1, v.vec, 2)
    zero = [0] * len(res._basis(2, "coh"))
    assert not comp.same_cohomology_class(tr, zero, degree=2), "expected a nonzero class"
    assert comp.same_cohomology_class(correct, tr, degree=2)
    assert not comp.same_cohomology_class(flipped, tr, degree=2), \
        "FLIPPED sign must disagree -- else the anchor is sign-blind here"


# The genuine ODD-exponent non-vacuous discriminator (the QCI(2,2) substitute).
# Heavy: needs Delta_5 for QCI (~2 min).  Module-scoped so it builds once.
@pytest.fixture(scope="module")
def _qci_diag5():
    from quiverlab.resolutions_cs.diagonal import diagonal
    comp = Comparison(_qci_gf())
    comp._ensure(6)
    diagonal(comp._res, 5)
    return comp


@pytest.mark.oracle_crossengine
@pytest.mark.slow
def test_odd_exponent_sign_anchor_qci_2_4(_qci_diag5):
    """THE odd-exponent sign anchor (Plan 51 DD4, correcting the vacuous QCI(2,2)):
    QCI/GF5 (2,4) -> HH^5 (dim 4), (p-1)(q-1)=3 ODD.  The bracket class is NONZERO,
    native == transported, and the FLIPPED sign disagrees -- the genuine crossengine
    fix of the (-1)^{(p-1)(q-1)} factor at an odd exponent."""
    comp = _qci_diag5
    res = comp._res
    p, q, n = 2, 4, 5
    assert 4 <= comp.window, "transport needs window >= 4 for (2,4)"
    repp, repq = comp.cs_cohomology_basis(p), comp.cs_cohomology_basis(q)
    saw_nonzero = False
    for i in range(len(repp)):
        for j in range(len(repq)):
            u, v = comp.hh_class_cs(p, i), comp.hh_class_cs(q, j)
            tr = comp.bracket_of_cs_classes(u, v)
            native = native_bracket(res, u.vec, p, v.vec, q)
            flipped = _bracket_flipped(res, u.vec, p, v.vec, q)
            zero = [0] * len(res._basis(n, "coh"))
            assert comp.same_cohomology_class(native, tr, degree=n), \
                f"native != transported at (2,4) reps ({i},{j})"
            if not comp.same_cohomology_class(tr, zero, degree=n):
                saw_nonzero = True
                assert not comp.same_cohomology_class(flipped, tr, degree=n), \
                    "FLIPPED sign must disagree on a nonzero odd-exponent class"
    assert saw_nonzero, "expected a nonzero (2,4) bracket class (HH^5 dim 4)"


def test_native_bracket_over_qq_computes():
    """Domain-generic: native_bracket computes over QQ on k[x]/x^2 (past the bar
    world entirely -- Comparison is GF(p)-gated and unused)."""
    res = _kx2(md=6, field=QQ)
    assert res.dom.name == "QQ"
    for f in _cocycles(res, 1):
        for g in _cocycles(res, 2):
            br = native_bracket(res, f, 1, g, 2)
            assert len(br) == len(res._basis(2, "coh"))
