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
