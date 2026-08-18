"""Plan 51 Task A.1/A.1a/A.2: the homotopy-lifting tower psi_eta on the CS
resolution (Negron-Witherspoon / Volkov), the past-window bracket's carrier.

Deep bucket (tests/resolutions_cs -> deep).  The DEFINING equation (self-cert):

    d_P(psi(sigma)) - (-1)^{n-1} psi(d sigma) == (eta (x) 1 - 1 (x) eta) Delta_k(sigma)

exactly, as PELT dicts, on every generator sigma in S_k.  Both sides are recomputed
INDEPENDENTLY here (d_P via pelt.apply_lower; the half-collapse from diagonal();
psi(d sigma) by the outer bimodule action) so the test does not merely echo the
solver.  Plus: psi-solve CONSISTENCY at every built degree (Volkov Lemma 2 -- never
the defensive assertion) and the CocycleError guard (a non-cocycle raises a DISTINCT
typed error, never Delta's NotImplementedError).
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.fields import QQ
from quiverlab.fields.linalg import nullspace
from quiverlab.groebner import build_reduction_system
from quiverlab.errors import CocycleError
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.diagonal import diagonal
from quiverlab.resolutions_cs.cup import _cochain_evaluator
from quiverlab.resolutions_cs.pelt import _accum, _resolve_chain, _vecs, apply_lower
from quiverlab.resolutions_cs.homotopy_lifting import homotopy_lifting

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_selfcert]


def _res(rels, arrows, md, field=None):
    f = GF(5) if field is None else field
    Q = Quiver([1], arrows)
    A = Q.algebra(relations=rels, field=f)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, f), max_degree=md)


def _kx2(md=6, field=None):
    return _res(["x*x"], {"x": (1, 1)}, md, field)


def _straddle(md=4):
    return _res(["x*x", "y*y", "x*y*x"], {"x": (1, 1), "y": (1, 1)}, md)


def _qci(md=4):
    return _res(["x*x", "y*y", "y*x - 2*x*y"], {"x": (1, 1), "y": (1, 1)}, md)


def _cocycles(res, n):
    dom = res.dom
    M = res.matrix(n, "coh")
    if not M:
        d = len(res._basis(n, "coh"))
        return [[dom.one() if i == j else dom.zero() for j in range(d)] for i in range(d)]
    return nullspace([[dom.coerce(x) for x in row] for row in M], dom)


def _cw(ch):
    return ("__v__", ch.o) if ch.degree == 0 else ch.word


# --------------------------------------------------------------------------- #
# INDEPENDENT recomputation of both sides of (star)                            #
# --------------------------------------------------------------------------- #
def _apply_dP(res, m, pelt):
    """d_m: P_m -> P_{m-1} on a single-PELT, via apply_lower (d_{N-1} on degree N-1)."""
    out = apply_lower(res, m + 1, pelt)
    return {k: v for k, v in out.items() if not res.dom.is_zero(v)}


def _half_collapse(res, eta_at, n, k, sigma):
    """(eta (x) 1 - 1 (x) eta) Delta_k(sigma), single-PELT over P_{k-n}. Independent
    re-implementation of the module's construction (with the Koszul sign on 1(x)eta)."""
    ar, dom = res.ar, res.dom
    tc = res._tensor_complex
    dpelt = diagonal(res, k).get(_cw(sigma), {})
    out = {}
    for (ai, tau_w, mi, rho_w, ci), coeff in dpelt.items():
        if dom.is_zero(coeff):
            continue
        tau, rho = tc._chain(tau_w), tc._chain(rho_w)
        b_a, b_mid, b_c = ar.A._basis_vec(ai), ar.A._basis_vec(mi), ar.A._basis_vec(ci)
        if tau.degree == n:
            ev = eta_at(tau)
            if not all(dom.is_zero(v) for v in ev):
                lvec = ar.mul(b_a, ar.mul(ev, b_mid))
                for aj, av in enumerate(lvec):
                    if not dom.is_zero(av):
                        _accum(out, (aj, rho.word, ci), dom.mul(coeff, av), dom)
        if rho.degree == n:
            ev = eta_at(rho)
            if not all(dom.is_zero(v) for v in ev):
                rvec = ar.mul(b_mid, ar.mul(ev, b_c))
                ks = dom.one() if (n * tau.degree) % 2 == 0 else dom.neg(dom.one())
                for cj, cv in enumerate(rvec):
                    if not dom.is_zero(cv):
                        _accum(out, (ai, tau.word, cj),
                               dom.neg(dom.mul(dom.mul(coeff, ks), cv)), dom)
    return out


def _psi_on_diff(res, psi, n, k, sigma):
    """psi(d_k sigma), single-PELT over P_{k-n}, by the outer bimodule action."""
    ar, dom = res.ar, res.dom
    tc = res._tensor_complex
    out = {}
    for (coeff, a_word, tw, c_word) in res.d_terms(k, sigma):
        if dom.is_zero(coeff):
            continue
        a_vec, c_vec = _vecs(res, sigma, a_word, c_word)
        target = _resolve_chain(res, tw)
        psi_tw = psi.apply(k - 1, target)
        for (ai, chw, ci), kappa in psi_tw.items():
            base = dom.mul(coeff, kappa)
            if dom.is_zero(base):
                continue
            left = ar.mul(a_vec, ar.A._basis_vec(ai))
            right = ar.mul(ar.A._basis_vec(ci), c_vec)
            for aj, av in enumerate(left):
                if dom.is_zero(av):
                    continue
                for cj, cv in enumerate(right):
                    if dom.is_zero(cv):
                        continue
                    _accum(out, (aj, chw, cj), dom.mul(base, dom.mul(av, cv)), dom)
    return out


def _star_holds(res, eta, n, kmax):
    """(star) exactly on every generator sigma in S_k for n <= k <= kmax."""
    dom = res.dom
    psi = homotopy_lifting(res, eta, n)
    eta_at = _cochain_evaluator(res, eta, n)
    sign = dom.one() if (n - 1) % 2 == 0 else dom.neg(dom.one())
    for k in range(n, kmax + 1):
        for sigma in res.ss.S(k):
            m = k - n + 1
            lhs = _apply_dP(res, m, psi.apply(k, sigma))
            rhs = _half_collapse(res, eta_at, n, k, sigma)
            for key, v in _psi_on_diff(res, psi, n, k, sigma).items():
                _accum(rhs, key, dom.mul(sign, v), dom)
            rhs = {kk: v for kk, v in rhs.items() if not dom.is_zero(v)}
            if lhs != rhs:
                return False, (k, sigma.word, lhs, rhs)
    return True, None


# --------------------------------------------------------------------------- #
# A.1 -- the defining equation holds exactly                                   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk,name,n,kmax", [
    (_kx2, "kx2", 1, 5), (_kx2, "kx2", 2, 5),
    (_straddle, "straddle", 1, 3), (_straddle, "straddle", 2, 3),
    (_qci, "qci", 1, 3), (_qci, "qci", 2, 3)])
def test_star_equation(mk, name, n, kmax):
    res = mk()
    for eta in _cocycles(res, n):
        ok, w = _star_holds(res, eta, n, kmax)
        assert ok, f"(star) failed on {name} at n={n}: {w}"


# --------------------------------------------------------------------------- #
# A.1a -- psi consistency at every built degree + the CocycleError guard       #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk,name,kmax", [
    (_straddle, "straddle", 3), (_qci, "qci", 3)])
def test_psi_consistent_every_degree(mk, name, kmax):
    """For every basis cocycle eta the psi-solve is CONSISTENT (never the defensive
    'inconsistent' AssertionError) at every degree Delta was built -- the concrete
    witness of the Lemma-2 guarantee (Plan 51 DD1a)."""
    res = mk()
    for n in (1, 2):
        for eta in _cocycles(res, n):
            psi = homotopy_lifting(res, eta, n)
            for k in range(kmax + 1):
                psi.tower(k)          # raises AssertionError if inconsistent -- must not


def test_cocycle_guard_is_distinct_typed_error():
    """A non-cocycle raises the DISTINCT CocycleError, never Delta's
    NotImplementedError (Plan 51 MINOR c)."""
    res = _kx2()
    dom = res.dom
    for n in (1, 2, 3):
        M = res.matrix(n, "coh")
        d = len(res._basis(n, "coh"))
        for i in range(d):
            e = [dom.one() if k == i else dom.zero() for k in range(d)]
            img = [sum((M[r][j] * e[j] for j in range(d)), dom.zero())
                   for r in range(len(M))] if M else []
            if any(not dom.is_zero(dom.coerce(x)) for x in img):
                with pytest.raises(CocycleError):
                    homotopy_lifting(res, e, n)
                return
    pytest.fail("no non-cocycle basis vector found to exercise the guard")


def test_homotopy_lifting_rejects_degree_zero():
    res = _kx2()
    with pytest.raises(ValueError):
        homotopy_lifting(res, [res.dom.one()], 0)


# --------------------------------------------------------------------------- #
# A.2 -- Domain-generic smoke over QQ                                          #
# --------------------------------------------------------------------------- #
def test_star_equation_over_qq():
    """The lifting is DOMAIN-GENERIC: (star) holds exactly over QQ on k[x]/x^2
    (Comparison is GF(p)-gated and not used -- drive the resolution directly)."""
    res = _kx2(md=6, field=QQ)
    assert res.dom.name == "QQ"
    for n in (1, 2):
        for eta in _cocycles(res, n):
            ok, w = _star_holds(res, eta, n, 5)
            assert ok, f"(star) failed over QQ at n={n}: {w}"
