"""QPA (GAP) as an EXTERNAL oracle for module Tor_n^A(M, N) (Plan 29 Part 4).

QPA ships NO native Tor: a live probe finds neither ``TorOverAlgebra`` nor ``Tor``
(``IsBoundGlobal`` both False). It DOES have ``HomOverAlgebra``, ``NthSyzygy`` and
``ExtOverAlgebra`` -- so we cross-check Tor through the classical duality

    dim Tor_n^A(M, N) = dim Ext_A^n(M, DN)   (D side-aware, Plan 24: DN of the LEFT
                                              module N is a RIGHT module),

computing the right-hand Ext in QPA with ITS OWN homological machinery and comparing to
our ``tor_dims``:  Ext^0 = dim Hom(M, DN); Ext^n = dim Ext^1(Omega^{n-1} M, DN) for
n >= 1 (dimension shifting via ``NthSyzygy``). QPA builds M and DN, resolves, and reads
the Ext dimensions independently of our engine, so a match is a genuine external
certificate for Tor. (Our side-aware D is itself QPA-validated elsewhere via
``DualOfModule`` in the tau / injective-resolution crosschecks.)

This test does NOT edit ``qpa/crosscheck.py`` or ``qpa/scripts.py`` -- it uses the
read-only script/module helpers and the live session.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import GF, QuantumCI, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.qpa_module import graded_form
from quiverlab.modules.tor import tor_dims
from quiverlab.qpa import scripts, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_has_no_native_tor():
    """Document the probe result the docstring relies on: QPA exposes no Tor entry
    point, so the duality-via-Ext route is the named coverage."""
    lg = session.libgap_handle()
    assert not bool(lg.IsBoundGlobal("TorOverAlgebra"))
    assert not bool(lg.IsBoundGlobal("Tor"))
    # the pieces we DO use exist
    for name in ("HomOverAlgebra", "NthSyzygy", "ExtOverAlgebra"):
        assert bool(lg.IsBoundGlobal(name)), name


def _qpa_ext_dims(A, M, Nright, top):
    """dim Ext_A^0..top(M, Nright) computed in QPA (M, Nright both RIGHT A-modules):
    Ext^0 = dim Hom; Ext^n = dim Ext^1(Omega^{n-1} M, N) via NthSyzygy dimension shift."""
    dvM, arrM = graded_form(M)
    dvN, arrN = graded_form(Nright)
    base = scripts.quiver_and_algebra_script(A)
    base += "\n" + scripts.module_decl(A, dvM, arrM, "MM")
    base += "\n" + scripts.module_decl(A, dvN, arrN, "NN")
    dims = [int(session.run(base + "\nhh := HomOverAlgebra(MM, NN);;\nLength(hh);"))]
    for n in range(1, top + 1):
        script = (base + f"\nsy := NthSyzygy(MM, {n - 1});;"
                  "\ne := ExtOverAlgebra(sy, NN);;\nLength(e[2]);")
        dims.append(int(session.run(script)))
    return dims


def _crosscheck(A, M, Nleft, top):
    ours = tor_dims(A, M, Nleft, top)
    qpa = _qpa_ext_dims(A, M, Nleft.dualize(), top)         # DN = D(N left) is right
    assert ours == qpa, (getattr(M, "name", "M"), getattr(Nleft, "name", "N"), ours, qpa)


@pytest.mark.parametrize("i,j", [(1, 1), (1, 2), (1, 3), (2, 3), (3, 3)])
def test_tor_kA3_simples_vs_qpa(i, j):
    """kA_3 simples: our Tor_n(S_i^right, S_j^left) vs QPA Ext^n(S_i, S_j^right)."""
    A = linear_path_algebra(3, field=QQ)
    _crosscheck(A, A.simple(i), A.simple(j, side="left"), 4)


def test_tor_kx2_infinite_pd_vs_qpa():
    """k[x]/(x^2): Tor_n(k, k) = k for all n, cross-checked against QPA to depth 5
    (exercises NthSyzygy on an infinite-pd module)."""
    B = truncated_polynomial(2, field=GF(7))
    _crosscheck(B, B.simple(1), B.simple(1, side="left"), 5)


def test_tor_quantum_ci_vs_qpa():
    """Quantum CI k<x,y>/(x^2, y^2, xy - 2 yx) over GF(7): Tor_n(k, k) = n + 1, checked
    against QPA's Ext of the trivial module against its dual."""
    Q = QuantumCI(2, field=GF(7))
    _crosscheck(Q, Q.simple(1), Q.simple(1, side="left"), 4)


def test_tor_nonsimple_pair_vs_qpa():
    """A non-simple pair on kA_3 (right interval [1,2] vs a left injective) vs QPA."""
    A = linear_path_algebra(3, field=QQ)
    M = A.module({1: 1, 2: 1, 3: 0}, {"a1": [[0, 0], [1, 0]], "a2": [[0, 0], [0, 0]]},
                 name="I12")
    _crosscheck(A, M, A.injective(3, side="left"), 3)
