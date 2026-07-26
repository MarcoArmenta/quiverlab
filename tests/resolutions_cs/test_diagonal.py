"""Plan 20 Task 1 gates: the double-PELT model of P_p ⊗_A P_q and its Koszul-signed
tensor differential d^{P⊗P} = d_p ⊗ 1 + (−1)^p 1 ⊗ d_q.

Deep bucket (tests/resolutions_cs → deep): run this FILE directly, e.g.
    NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q \
        tests/resolutions_cs/test_diagonal.py -p no:cacheprovider

Gates
-----
* (d^{P⊗P})² = 0 on the zoo — k[x]/x² (quadratic), the straddle-monomial
  k⟨x,y⟩/(xx,yy,xyx), and the quantum-CI k⟨x,y⟩/(xx,yy,yx−2xy) — over GF(5),
  degrees n = 1..4, every corner.  d^{P⊗P} is an A^e-module (bimodule) map, so
  d²=0 on the whole free module ⟺ d²=0 on each free A^e-generator (an element
  with vertex idempotents in the two outer slots).  The GENERATOR form of the
  gate (apply_tensor_d twice on every generator, over every corner) is therefore
  the exact matrix identity in its faithful, materialisable form — it is not a
  weaker check.  The literal dense `tensor_matrix(n-1)·tensor_matrix(n)==0` is
  ALSO asserted on the fixtures/degrees where the full (o,t)-corner k-space is
  small enough to materialise, to exercise `tensor_matrix` directly.
* Koszul-sign unit pin: on k[x]/x³ at p=q=1, one hand-computed entry confirms the
  −1 sits on the 1⊗d summand (and the paired d⊗1 entry is +1).
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.diagonal import TensorComplex

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_selfcert]


# --------------------------------------------------------------------------- #
# fixtures — single/double-loop quivers matching the mathematical intent       #
# (k[x]/x² is a genuine single loop, NOT k⟨x,y⟩/x²).                            #
# --------------------------------------------------------------------------- #
def _res(rels, arrows, max_degree=6, p=5):
    f = GF(p)
    Q = Quiver([1], arrows)
    A = Q.algebra(relations=rels, field=f)
    return ChouhySolotarResolution(A, build_reduction_system(Q, rels, f),
                                   max_degree=max_degree)


def _kx2():
    return _res(["x*x"], {"x": (1, 1)})


def _straddle():
    return _res(["x*x", "y*y", "x*y*x"], {"x": (1, 1), "y": (1, 1)})


def _qci():
    return _res(["x*x", "y*y", "y*x - 2*x*y"], {"x": (1, 1), "y": (1, 1)})


def _kx3():
    return _res(["x*x*x"], {"x": (1, 1)})


ZOO = [(_kx2, "kx2"), (_straddle, "straddle-monomial"), (_qci, "qci")]


# --------------------------------------------------------------------------- #
# helpers                                                                      #
# --------------------------------------------------------------------------- #
def _dd_on_generators(tc, n):
    """(d^{P⊗P})² = 0 checked on every free A^e-generator over every corner —
    the faithful, materialisable form of the matrix identity."""
    res, dom = tc.res, tc.res.dom
    verts = list(res.rs.quiver.vertices)
    for o in verts:
        for t in verts:
            for gen in tc.generators(n, o, t):
                img = tc.apply_tensor_d(n, {gen: dom.one()})
                dd = tc.apply_tensor_d(n - 1, img)
                assert not dd, f"d²≠0 at n={n}, generator {gen}: {dd}"


def _matrix_dd_zero(tc, n, o, t):
    """literal dense identity tensor_matrix(n-1,·)·tensor_matrix(n,·) == 0."""
    res, dom = tc.res, tc.res.dom
    prod = res._matmul(tc.tensor_matrix(n - 1, o, t), tc.tensor_matrix(n, o, t))
    assert all(dom.is_zero(x) for row in prod for x in row), \
        f"tensor_matrix d²≠0 at n={n}, corner ({o},{t})"


# --------------------------------------------------------------------------- #
# (d^{P⊗P})² = 0 — generator form, all three algebras, n=1..4                  #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk,name", ZOO, ids=[n for _, n in ZOO])
def test_tensor_dd_zero_generators(mk, name):
    tc = TensorComplex(mk())
    for n in range(1, 5):
        _dd_on_generators(tc, n)


# --------------------------------------------------------------------------- #
# (d^{P⊗P})² = 0 — literal dense matrix identity where the corner is small      #
# --------------------------------------------------------------------------- #
def test_tensor_matrix_dd_zero_kx2():
    tc = TensorComplex(_kx2())
    for o in tc.res.rs.quiver.vertices:
        for t in tc.res.rs.quiver.vertices:
            for n in range(1, 5):
                _matrix_dd_zero(tc, n, o, t)


def test_tensor_matrix_dd_zero_qci():
    tc = TensorComplex(_qci())
    for n in range(1, 5):
        _matrix_dd_zero(tc, n, 1, 1)


def test_tensor_matrix_dd_zero_straddle_lowdeg():
    # full-corner space blows up at n=4 (dim³·Σ|S_p||S_q| ≈ 10⁴); the generator
    # gate covers n=4 for the straddle. Here we exercise the dense matrix at n≤2.
    tc = TensorComplex(_straddle())
    for n in (1, 2):
        _matrix_dd_zero(tc, n, 1, 1)


# --------------------------------------------------------------------------- #
# Koszul-sign unit pin (k[x]/x³, p=q=1): the −1 sits on 1⊗d                     #
# --------------------------------------------------------------------------- #
def test_koszul_sign_pin_kx3():
    """Hand computation of d^{P⊗P}_2 on the generator g = 1 ⊗ x ⊗ 1 ⊗ x ⊗ 1 of
    P_1 ⊗_A P_1 over A = k[x]/x³.  With d_1(1⊗x⊗1) = x⊗e⊗1 − 1⊗e⊗x,

        d(g) = (d_1⊗1)(g)  +  (−1)^1 (1⊗d_1)(g)
             = [ x⊗e ⊗ e⊗x⊗e  −  e⊗e ⊗ x⊗x⊗e ]        # d⊗1, sign +1
               − [ e⊗x ⊗ x⊗e⊗e  −  e⊗x ⊗ e⊗e⊗x ]      # 1⊗d, Koszul −1

    so the FIRST 1⊗d term e⊗x⊗x⊗e⊗e carries −1 while the paired d⊗1 term
    x⊗e⊗e⊗x⊗e carries +1 — the Koszul sign lives on the second summand only.
    Basis indices (kx3): e_1 → 0, path x → 1."""
    tc = TensorComplex(_kx3())
    res, dom = tc.res, tc.res.dom
    ev = res.ar._word_index[("v", 1)]          # 0
    xi = res.ar._word_index[("p", ("x",))]     # 1
    gen = (ev, ("x",), ev, ("x",), ev)
    img = {k: res.to_int(v) for k, v in tc.apply_tensor_d(2, {gen: dom.one()}).items()}

    key_1d = (ev, ("x",), xi, ("__v__", 1), ev)   # 1⊗d first term
    key_d1 = (xi, ("__v__", 1), ev, ("x",), ev)   # d⊗1 first term
    assert img[key_1d] == (-1) % 5, "the −1 must sit on the 1⊗d summand"
    assert img[key_d1] == 1, "the paired d⊗1 entry is sign-free"

    assert img == {
        key_d1:                              1,      #  x ⊗ e ⊗ e ⊗ x ⊗ e
        (ev, ("__v__", 1), xi, ("x",), ev): (-1) % 5,  # −e ⊗ e ⊗ x ⊗ x ⊗ e
        key_1d:                             (-1) % 5,  # −e ⊗ x ⊗ x ⊗ e ⊗ e  (Koszul)
        (ev, ("x",), ev, ("__v__", 1), xi):  1,     #  e ⊗ x ⊗ e ⊗ e ⊗ x
    }


# --------------------------------------------------------------------------- #
# Task 2: the comparison-lifted diagonal Δ: P → P ⊗_A P                        #
#   Δ is the chain-map lift of id_A through the lifting equation               #
#     d^{P⊗P}_n · Δ_n(σ) = Δ_{n−1}(d_n σ)                                       #
#   solved per generator (fields.linalg.solve + reduce_mod_nullspace), the     #
#   structural clone of resolution.py::_d_general's correction solve.          #
# --------------------------------------------------------------------------- #
from quiverlab.resolutions_cs.diagonal import diagonal as diagonal_fn


def _cw(ch):
    """Storage word of a chain (== TensorComplex._chain_word == d_terms targets)."""
    return ("__v__", ch.o) if ch.degree == 0 else ch.word


def _assert_pq_structure(tc, n, sigma, dpelt):
    """(p,q)-component sanity: every key (a, τ, mid, ρ, c) of Δ_n(σ) has
    deg τ + deg ρ == n and endpoint-compatible corners
    a ∈ e_{o(σ)}Ae_{o(τ)},  mid ∈ e_{t(τ)}Ae_{o(ρ)},  c ∈ e_{t(ρ)}Ae_{t(σ)}."""
    ar = tc.ar
    for (ai, tau_w, mi, rho_w, ci) in dpelt:
        tau, rho = tc._chain(tau_w), tc._chain(rho_w)
        assert tau.degree + rho.degree == n, \
            f"key {(tau_w, rho_w)}: deg {tau.degree}+{rho.degree} != {n}"
        assert ai in ar.corner(sigma.o, tau.o, "coh"), f"a-slot corner: {(ai, tau_w)}"
        assert mi in ar.corner(tau.t, rho.o, "coh"), f"mid-slot corner: {(mi, tau_w, rho_w)}"
        assert ci in ar.corner(rho.t, sigma.t, "coh"), f"c-slot corner: {(ci, rho_w)}"


def _canon(diag, res):
    """int-normalised copy of a Δ dict for exact cross-run (byte) comparison."""
    return {w: {k: res.to_int(v) for k, v in dp.items()} for w, dp in diag.items()}


# --------------------------------------------------------------------------- #
# (a) chain-map identity  apply_tensor_d(n, Δ_n(σ)) == Δ_{n−1}(d_n σ)          #
#     + (d) (p,q)-component structural sanity, folded in per generator.        #
#     kx2/qci to n=4; straddle to n=3 ONLY — the n=4 straddle system is        #
#     ≈4752×10152 dense over GF(5) (pure-python solve infeasible), budgeted out.#
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk,name,maxn", [
    (_kx2, "kx2", 4),
    (_qci, "qci", 4),
    (_straddle, "straddle-monomial", 3),
])
def test_diagonal_chain_map(mk, name, maxn):
    """Δ is a chain map. Two validations beyond the by-construction solve:
      * the lift is NEVER inconsistent on these in-CS-scope fixtures (a wrong ζ
        gluing would generically make ζ a non-cycle, and the solve would raise
        the loud NotImplementedError);
      * ζ = Δ_{n−1}(d_n σ) is itself a d^{P⊗P}-cycle (apply_tensor_d(n−1, ζ) == {})
        for n ≥ 2 — the propagated chain-map property, checked independently of
        the lift-solve."""
    tc = TensorComplex(mk())
    res = tc.res
    for n in range(1, maxn + 1):
        diag_n = tc.diagonal(n)
        prev = tc.diagonal(n - 1)
        for sigma in res.ss.S(n):
            zeta = tc._zeta(n, sigma, prev)
            dn = diag_n[_cw(sigma)]
            assert tc.apply_tensor_d(n, dn) == zeta, \
                f"chain-map identity failed at {name} n={n}, σ={sigma.word}"
            if n >= 2:
                assert tc.apply_tensor_d(n - 1, zeta) == {}, \
                    f"ζ(σ) is not a d^(P⊗P)-cycle at {name} n={n}, σ={sigma.word}"
            _assert_pq_structure(tc, n, sigma, dn)


# --------------------------------------------------------------------------- #
# (b) byte-reproducibility: two fresh resolutions → identical Δ dicts          #
# --------------------------------------------------------------------------- #
def test_diagonal_byte_reproducible_qci():
    """Plan-17 law: Δ (via solve + reduce_mod_nullspace) is byte-reproducible."""
    ra, rb = _qci(), _qci()
    for n in range(4):
        da = diagonal_fn(ra, n)
        db = diagonal_fn(rb, n)
        assert _canon(da, ra) == _canon(db, rb), f"Δ_{n} not byte-reproducible"


# --------------------------------------------------------------------------- #
# (c) Δ_0 base shape: Δ_0(σ_v) = e_v ⊗ σ_v ⊗ e_v ⊗ σ_v ⊗ e_v                   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk,name", ZOO, ids=[n for _, n in ZOO])
def test_diagonal_base_shape(mk, name):
    tc = TensorComplex(mk())
    res, dom = tc.res, tc.res.dom
    d0 = tc.diagonal(0)
    verts = list(res.rs.quiver.vertices)
    assert set(d0) == {("__v__", v) for v in verts}
    for v in verts:
        ev = res.ar._word_index[("v", v)]
        w = ("__v__", v)
        assert d0[w] == {(ev, w, ev, w, ev): dom.one()}, f"Δ_0 base shape wrong at v={v}"


# --------------------------------------------------------------------------- #
# (d) (p,q)-component structure smoke — a dedicated, lightweight gate          #
# --------------------------------------------------------------------------- #
def test_diagonal_pq_structure_smoke():
    """Every Δ_n(σ) key splits as τ.degree + ρ.degree == n with endpoint-compatible
    corners (checked cheaply on kx2 to n=4 and qci to n=2)."""
    for mk, maxn in ((_kx2, 4), (_qci, 2)):
        tc = TensorComplex(mk())
        for n in range(maxn + 1):
            diag_n = tc.diagonal(n)
            for sigma in tc.res.ss.S(n):
                _assert_pq_structure(tc, n, sigma, diag_n[_cw(sigma)])


# --------------------------------------------------------------------------- #
# (e) the scope edge: an inconsistent lift-solve raises the loud               #
#     NotImplementedError -- never a silent fallback (Plan 21 review item).    #
# --------------------------------------------------------------------------- #
def test_diagonal_inconsistent_lift_raises_loudly(monkeypatch):
    """The diagonal's lift-solve `d^{P⊗P}·Δ_n(σ) = ζ(σ)` shares the exact scope edge of
    resolution.py::_d_general: an inconsistent solve (no lift) is refused with a loud
    NotImplementedError ("higher CS homotopy correction … spec §6 risk register"),
    NEVER a silent fallback.  No in-CS-scope algebra reaches it (the lift always closes
    on the certified fixtures), so the edge is exercised by a CONSTRUCTED refusal:
    force `fields.linalg.solve` (as imported in diagonal.py) to report inconsistency by
    returning None, and demand the loud raise at the first degree that actually solves.

    Δ_0 (base case, no solve) is unaffected; degree 1 hits the solve branch and must
    raise."""
    import quiverlab.resolutions_cs.diagonal as diagmod
    tc = TensorComplex(_kx2())
    assert tc.diagonal(0)                                  # base case: no solve, fine
    monkeypatch.setattr(diagmod, "solve", lambda M, rhs, dom: None)
    with pytest.raises(NotImplementedError, match="higher CS homotopy"):
        tc.diagonal(1)
