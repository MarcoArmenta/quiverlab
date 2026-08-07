"""Plan 51 Task E: the native-bracket identity batteries -- graded Jacobi, the
native Poisson/cup-Leibniz law, and off-GF(p) antisymmetry -- class-level (mod
coboundary), over GF primes AND over QQ, on the CS resolution.

Deep bucket.  All checks are field-agnostic (a coboundary membership solve over
res.dom), so the SAME identity is exercised over GF(p) and over QQ with no
Comparison (past the bar world).  The bracket sign is FIXED by the Task-B anchor;
these are Lie-structure consistency + genuine sign constraints (Poisson).

DEVIATION NOTE (Plan E.1): delivered class-level (mod coboundary) over GF primes +
QQ rather than by composing the structure-constant TABLES -- class-level is the
mathematically meaningful form and the table constants are independently pinned by
the D.1 bar==CS cross-engine gate.  Non-vacuity is asserted (a nonzero intermediate
bracket/cup class must occur), documented per test.
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace, solve
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.bracket import native_bracket
from quiverlab.resolutions_cs.cup import native_cup

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_selfcert]

PRIMES = (32003, 5, 3)


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


def _is_cocycle(res, vec, n):
    dom = res.dom
    for row in res.matrix(n, "coh"):
        acc = dom.zero()
        for j, m in enumerate(row):
            acc = dom.add(acc, dom.mul(m, dom.coerce(vec[j])))
        if not dom.is_zero(acc):
            return False
    return True


def _is_coboundary(res, vec, n):
    """vec is a coboundary in C^n: a cocycle AND in colspan(delta^{n-1})."""
    dom = res.dom
    if not _is_cocycle(res, vec, n):
        return False
    if n == 0:
        return all(dom.is_zero(dom.coerce(v)) for v in vec)
    Md = res.matrix(n - 1, "coh")           # rows C^n, cols C^{n-1}; columns span im
    if not (Md and Md[0]):
        return all(dom.is_zero(dom.coerce(v)) for v in vec)
    x = solve([[dom.coerce(c) for c in row] for row in Md],
              [dom.coerce(v) for v in vec], dom)
    return x is not None


def _lin(res, *terms):
    """sum_k s_k * v_k with s_k in {+1,-1}, over res.dom."""
    dom = res.dom
    n = len(terms[0][1])
    out = [dom.zero()] * n
    for s, v in terms:
        ss = dom.one() if s > 0 else dom.neg(dom.one())
        for i, x in enumerate(v):
            out[i] = dom.add(out[i], dom.mul(ss, dom.coerce(x)))
    return out


def _sgn(e):
    return 1 if e % 2 == 0 else -1


# =========================================================================== #
# E.1 -- graded Jacobi                                                         #
# =========================================================================== #
def _jacobi_holds(res, triples):
    """(-1)^{(p-1)(r-1)}[[f,g],h] + cyclic ~ 0 mod coboundary, for each degree
    triple (p,q,r) on all basis-cocycle triples. Returns (ok, saw_nonzero_bracket)."""
    saw_nonzero = False
    for (p, q, r) in triples:
        Zp, Zq, Zr = _cocycles(res, p), _cocycles(res, q), _cocycles(res, r)
        out_deg = p + q + r - 2
        for f in Zp:
            for g in Zq:
                for h in Zr:
                    fg = native_bracket(res, f, p, g, q)
                    gh = native_bracket(res, g, q, h, r)
                    hf = native_bracket(res, h, r, f, p)
                    if not _is_coboundary(res, fg, p + q - 1):
                        saw_nonzero = True
                    t1 = native_bracket(res, fg, p + q - 1, h, r)      # [[f,g],h]
                    t2 = native_bracket(res, gh, q + r - 1, f, p)      # [[g,h],f]
                    t3 = native_bracket(res, hf, r + p - 1, g, q)      # [[h,f],g]
                    jac = _lin(res,
                               (_sgn((p - 1) * (r - 1)), t1),
                               (_sgn((q - 1) * (p - 1)), t2),
                               (_sgn((r - 1) * (q - 1)), t3))
                    if not _is_coboundary(res, jac, out_deg):
                        return False, saw_nonzero
    return True, saw_nonzero


@pytest.mark.parametrize("prime", PRIMES)
def test_graded_jacobi_kx2_gfp(prime):
    res = _res(["x*x"], {"x": (1, 1)}, 5, GF(prime))
    ok, nz = _jacobi_holds(res, [(1, 1, 1), (1, 1, 2), (2, 1, 1)])
    assert ok, f"graded Jacobi failed on k[x]/x^2 over GF({prime})"
    assert nz, "vacuous: every intermediate bracket was a coboundary"


def test_graded_jacobi_kx2_qq():
    res = _res(["x*x"], {"x": (1, 1)}, 5, QQ)
    ok, nz = _jacobi_holds(res, [(1, 1, 1), (1, 1, 2), (2, 1, 1)])
    assert ok, "graded Jacobi failed on k[x]/x^2 over QQ"
    assert nz


def test_graded_jacobi_qci_hh1_lie():
    """QCI/GF5 at (1,1,1): a multi-generator Jacobi on the 2-dimensional HH^1.  Here
    the outer-derivation Lie algebra is ABELIAN ([HH^1,HH^1] = 0 as classes), so
    Jacobi holds trivially -- the content-bearing (non-vacuous) Jacobi is the k[x]/x^2
    battery above (its (2,1,1) triple has a nonzero degree-2 bracket)."""
    res = _res(["x*x", "y*y", "y*x - 2*x*y"], {"x": (1, 1), "y": (1, 1)}, 4, GF(5))
    assert len(_cocycles(res, 1)) >= 2, "QCI HH^1 should be >= 2-dimensional"
    ok, _ = _jacobi_holds(res, [(1, 1, 1)])
    assert ok, "graded Jacobi failed on QCI/GF5 HH^1"


# =========================================================================== #
# E.2 -- native Poisson / cup-Leibniz                                          #
# =========================================================================== #
def _poisson_holds(res, triples):
    """[f, g cup h] = [f,g] cup h + (-1)^{(p-1)q} g cup [f,h] mod coboundary."""
    saw_nonzero = False
    for (p, q, r) in triples:
        Zp, Zq, Zr = _cocycles(res, p), _cocycles(res, q), _cocycles(res, r)
        out_deg = p + q + r - 1
        for f in Zp:
            for g in Zq:
                for h in Zr:
                    gh = native_cup(res, g, q, h, r)                    # C^{q+r}
                    lhs = native_bracket(res, f, p, gh, q + r)          # C^{p+q+r-1}
                    fg = native_bracket(res, f, p, g, q)
                    fh = native_bracket(res, f, p, h, r)
                    if not _is_coboundary(res, lhs, out_deg):
                        saw_nonzero = True
                    t1 = native_cup(res, fg, p + q - 1, h, r)
                    t2 = native_cup(res, g, q, fh, p + r - 1)
                    rhs = _lin(res, (1, t1), (_sgn((p - 1) * q), t2))
                    diff = _lin(res, (1, lhs), (-1, rhs))
                    if not _is_coboundary(res, diff, out_deg):
                        return False, saw_nonzero
    return True, saw_nonzero


def test_poisson_cup_leibniz_kx2_qq():
    res = _res(["x*x"], {"x": (1, 1)}, 6, QQ)
    ok, nz = _poisson_holds(res, [(1, 1, 1), (1, 1, 2), (2, 1, 1)])
    assert ok, "native Poisson/cup-Leibniz failed on k[x]/x^2 over QQ"
    assert nz, "vacuous: [f, g cup h] was always a coboundary"


def test_poisson_cup_leibniz_qci_gf5():
    """QCI/GF5 at (1,1,1) -- native cup + native bracket, off any transport window."""
    res = _res(["x*x", "y*y", "y*x - 2*x*y"], {"x": (1, 1), "y": (1, 1)}, 4, GF(5))
    ok, _ = _poisson_holds(res, [(1, 1, 1)])
    assert ok, "native Poisson/cup-Leibniz failed on QCI over GF(5)"


# =========================================================================== #
# E.3 -- antisymmetry, native / off-GF(p) row                                  #
# =========================================================================== #
@pytest.mark.parametrize("field,name", [(QQ, "QQ"), (GF(5), "GF5")])
def test_antisymmetry_native(field, name):
    """[f,g] ~ -(-1)^{(p-1)(q-1)} [g,f] mod coboundary, class-level over QQ and GF5
    (extends the table-level GF(p) test_bracket_antisymmetry to the native route)."""
    res = _res(["x*x"], {"x": (1, 1)}, 6, field)
    for (p, q) in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        for f in _cocycles(res, p):
            for g in _cocycles(res, q):
                fg = native_bracket(res, f, p, g, q)
                gf = native_bracket(res, g, q, f, p)
                s = _sgn((p - 1) * (q - 1))
                diff = _lin(res, (1, fg), (s, gf))       # fg + s*gf ~ 0  (since fg ~ -s gf)
                assert _is_coboundary(res, diff, p + q - 1), \
                    f"antisymmetry failed over {name} at (p,q)=({p},{q})"
