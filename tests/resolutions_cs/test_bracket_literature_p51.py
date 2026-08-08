"""Plan 51 Task F: literature oracles for the native Gerstenhaber bracket.

Deep bucket.
  F.1 -- k[x]/(x^n) off GF(p) (the "any exact field" headline): over QQ the classical
         HH*(k[x]/x^2) Gerstenhaber structure -- the record's k[x]/(x^n) ZERO self-
         brackets, plus the one NONZERO mixed bracket, class-level; char-2 native route.
  F.2 -- QuantumCI ties: (in-window) CS-native == Plan-35 transported over GF(5)/GF(3)
         [crossengine]; (past-window) the CS-native bracket COMPUTES and its table DIMS
         line up with the QCI HH dims [selfcert -- dims, NOT a literature bracket value].
  F.3 -- Oke arXiv:2103.12331 section 7 bracket tables: BLOCKED-until-transcribed xfail
         fence (never fabricate; open the PDF, transcribe verbatim, then flip).
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace, solve
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.comparison import Comparison
from quiverlab.resolutions_cs.bracket import native_bracket

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_literature]


def _res(rels, arrows, md, field):
    Q = Quiver([1], arrows)
    A = Q.algebra(relations=rels, field=field)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, field), max_degree=md)


def _cocycles(res, n):
    dom = res.dom
    M = res.matrix(n, "coh")
    if not M:
        d = len(res._basis(n, "coh"))
        return [[dom.one() if i == j else dom.zero() for j in range(d)] for i in range(d)]
    return nullspace([[dom.coerce(x) for x in row] for row in M], dom)


def _is_coboundary(res, vec, n):
    dom = res.dom
    for row in res.matrix(n, "coh"):
        acc = dom.zero()
        for j, m in enumerate(row):
            acc = dom.add(acc, dom.mul(m, dom.coerce(vec[j])))
        if not dom.is_zero(acc):
            return False                    # not a cocycle -> not a coboundary
    if n == 0:
        return all(dom.is_zero(dom.coerce(v)) for v in vec)
    Md = res.matrix(n - 1, "coh")
    if not (Md and Md[0]):
        return all(dom.is_zero(dom.coerce(v)) for v in vec)
    x = solve([[dom.coerce(c) for c in row] for row in Md],
              [dom.coerce(v) for v in vec], dom)
    return x is not None


def _nonzero_rep(res, n):
    for z in _cocycles(res, n):
        if not _is_coboundary(res, z, n):
            return z
    return None


# =========================================================================== #
# F.1 -- k[x]/x^2 off GF(p): the classical bracket zeros + one nonzero          #
# =========================================================================== #
def test_kx2_bracket_zeros_and_nonzero_over_qq():
    """k[x]/x^2 over QQ (char 0): [alpha,alpha] = [beta,beta] = 0 (the record's
    k[x]/(x^n) zero self-brackets) and [alpha,beta] != 0 (the nontrivial Lie bracket).
    Class-level (mod coboundary), never bytes."""
    res = _res(["x*x"], {"x": (1, 1)}, 6, QQ)
    assert res.dom.name == "QQ"
    a = _nonzero_rep(res, 1)                                 # HH^1 = k.alpha
    b = _nonzero_rep(res, 2)                                 # HH^2 = k.beta
    assert a is not None and b is not None
    assert _is_coboundary(res, native_bracket(res, a, 1, a, 1), 1), "[a,a] must vanish"
    assert _is_coboundary(res, native_bracket(res, b, 2, b, 2), 3), "[b,b] must vanish"
    ab = native_bracket(res, a, 1, b, 2)                     # HH^2
    assert not _is_coboundary(res, ab, 2), "[alpha,beta] must be a NONZERO HH^2 class"


def test_kx2_bracket_char2_native_matches_transport():
    """The char-2 native route: over GF(2), the CS-native bracket computes and equals
    the transported bracket mod coboundary on every in-window nonzero pair."""
    comp = Comparison(Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(2)))
    comp._ensure(4)
    res = comp._res
    tested = 0
    for (p, q) in [(1, 1), (1, 2), (2, 1)]:
        if max(p, q) > comp.window:
            continue
        for i in range(len(comp.cs_cohomology_basis(p))):
            for j in range(len(comp.cs_cohomology_basis(q))):
                u, v = comp.hh_class_cs(p, i), comp.hh_class_cs(q, j)
                native = native_bracket(res, u.vec, p, v.vec, q)
                assert comp.same_cohomology_class(
                    native, comp.bracket_of_cs_classes(u, v), degree=p + q - 1)
                tested += 1
    assert tested > 0


# =========================================================================== #
# F.2 -- QuantumCI ties                                                         #
# =========================================================================== #
def _qci_comp(prime):
    return Comparison(Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=GF(prime)))


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("prime", (5, 3))
def test_qci_bracket_in_window_native_equals_transported(prime):
    """IN-WINDOW: CS-native bracket == Plan-35 GF(p) transported bracket, mod
    coboundary, over GF(5)/GF(3) (a real crossengine comparator)."""
    comp = _qci_comp(prime)
    comp._ensure(4)
    res = comp._res
    tested = 0
    for (p, q) in [(1, 1), (1, 2), (2, 1)]:
        if max(p, q) > comp.window:
            continue
        for i in range(len(comp.cs_cohomology_basis(p))):
            for j in range(len(comp.cs_cohomology_basis(q))):
                u, v = comp.hh_class_cs(p, i), comp.hh_class_cs(q, j)
                native = native_bracket(res, u.vec, p, v.vec, q)
                assert comp.same_cohomology_class(
                    native, comp.bracket_of_cs_classes(u, v), degree=p + q - 1)
                tested += 1
    assert tested > 0


@pytest.mark.oracle_selfcert
def test_qci_bracket_table_dims_line_up_with_hh_dims():
    """PAST-WINDOW: the CS-native bracket table COMPUTES and its per-bidegree dims are
    (dim HH^p, dim HH^q, dim HH^{p+q-1}), matching the QCI HH dims [2,2,1,0,...]
    (BGMS/Bergh-Erdmann shape).  This certifies self-consistency + correct DIMS, NOT a
    literature bracket VALUE (there is no transport comparator or published table past
    the window; the Oke section-7 values are F.3, blocked)."""
    from quiverlab.resolutions_cs.homology import cs_cohomology_dims
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - 2*x*y"], field=GF(5))
    hh = list(cs_cohomology_dims(A, 3).dims)
    assert hh[:4] == [2, 2, 1, 0], f"QCI/GF5 HH dims drifted: {hh[:4]}"
    hb = A.gerstenhaber_brackets(2, engine="cs")
    assert hb.window is None and hb.basis == "cs/GF(5)"
    for (p, q), t in hb.tables.items():
        assert t.dims == (hh[p], hh[q], hh[p + q - 1]), \
            f"bracket table dims at {(p, q)} != HH dims"


# =========================================================================== #
# F.3 -- Oke section 7 tables: BLOCKED-until-transcribed                        #
# =========================================================================== #
@pytest.mark.xfail(reason="Oke arXiv:2103.12331 section 7 quiver/relations/bracket "
                          "values not yet transcribed from the PDF -- DO NOT fabricate; "
                          "open the PDF, transcribe verbatim with proposition/equation "
                          "numbers, then flip to a real assert (Plan 51 F.3 / section 2).",
                   strict=True)
def test_oke_section7_koszul_bracket_tables():
    # BLOCKED: the Koszul quiver-algebra family of Oke arXiv:2103.12331 section 7 (7.1
    # deg-2 cocycles, 7.2 deg-1 cocycles, 7.3 derivation operators) with its explicit
    # homotopy liftings and bracket values must be transcribed VERBATIM (quiver,
    # relations, each value + equation number) before pinning.  Until then this is an
    # honest strict-xfail fence (the Plan-31 auto-flip precedent), never a fabricated
    # value.  When transcribed: build the section-7 algebra, compute native_bracket on
    # the section-7.1/7.2 cocycles, and assert the section-7 values class-level.
    raise AssertionError("Oke section 7 not yet transcribed")
