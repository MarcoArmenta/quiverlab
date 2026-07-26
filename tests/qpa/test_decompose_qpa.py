"""QPA (GAP) as the independent oracle for Krull-Schmidt decomposition (Plan 30 Part A).

Crosschecks quiverlab's ``decompose`` / ``is_indecomposable`` against QPA's
``DecomposeModuleWithMultiplicities`` / ``IsIndecomposableModule`` across the zoo
(simples, projectives, injectives of kA_2 / kA_3 / commutative square / the Plan-18 line
algebra; the quantum complete intersection and a truncated polynomial ring) plus
constructed direct sums (S (+) S, P (+) S (+) P, mixed). QPA's decomposition primitives
require a FINITE field, and our trace-form certificate needs char > dim M, so every
crosscheck runs over GF(7). qpa-marked: skips locally, mandatory under
QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import GF, Quiver, QuantumCI, linear_path_algebra, truncated_polynomial
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

_F = GF(7)                                       # char 7 > dim of every module below


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=_F)


def _line_abc_cde():
    # Plan-18 multi-vertex record: 1->2->...->6 with relations a*b*c and c*d*e.
    Q = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 5), "e": (5, 6)})
    return Q.algebra(relations=["a*b*c", "c*d*e"], field=_F)


def _indec_modules(A):
    return ([A.simple(v) for v in A.quiver.vertices]
            + [A.projective(v) for v in A.quiver.vertices]
            + [A.injective(v) for v in A.quiver.vertices])


def _direct_sum(mods, name="DS"):
    A = mods[0].algebra
    dom = A.domain
    labels = set().union(*(m.action.keys() for m in mods))
    dims = [m.dim for m in mods]
    N = sum(dims)
    offs, o = [], 0
    for d in dims:
        offs.append(o)
        o += d
    action = {}
    for lab in labels:
        Mb = lm.zeros(N, N, dom)
        for m, off in zip(mods, offs):
            blk = m.action.get(lab)
            if blk is None:
                continue
            for i in range(m.dim):
                for j in range(m.dim):
                    Mb[off + i][off + j] = blk[i][j]
        action[lab] = Mb
    return Module(A, N, action, name=name, side=mods[0].side)


# --- indecomposables: our verdict + decomposition agree with QPA ------------
@pytest.mark.parametrize("A", [linear_path_algebra(2, field=_F),
                               linear_path_algebra(3, field=_F),
                               _square(), _line_abc_cde()])
def test_indecomposables_crosscheck(A):
    for M in _indec_modules(A):
        A.crosscheck("indecomposable", M).assert_agree()
        A.crosscheck("decompose", M).assert_agree()      # single summand, multiplicity 1


# --- constructed direct sums: full Krull-Schmidt data agrees with QPA -------
def test_direct_sums_crosscheck():
    A = linear_path_algebra(2, field=_F)
    P1, S1, S2, I2 = A.projective(1), A.simple(1), A.simple(2), A.injective(2)
    sums = [_direct_sum([S1, S1], "S1+S1"),
            _direct_sum([P1, S2, P1], "P1+S2+P1"),
            _direct_sum([P1, I2, S1, S1], "P1+I2+S1+S1")]
    for M in sums:
        assert M.check_module()[0]
        A.crosscheck("indecomposable", M).assert_agree()
        A.crosscheck("decompose", M).assert_agree()


def test_direct_sum_multivertex_crosscheck():
    A = _line_abc_cde()
    M = _direct_sum([A.projective(1), A.injective(6), A.projective(1)], "P1+I6+P1")
    A.crosscheck("decompose", M).assert_agree()
    A.crosscheck("indecomposable", M).assert_agree()


# --- single-vertex algebras with loops + relations (QCI, k[x]/(x^n)) --------
def test_quantum_ci_and_truncated_crosscheck():
    for A, M in [(QuantumCI(-1, field=_F), QuantumCI(-1, field=_F).projective(1)),
                 (truncated_polynomial(3, field=_F),
                  truncated_polynomial(3, field=_F).projective(1))]:
        A.crosscheck("indecomposable", M).assert_agree()
        A.crosscheck("decompose", M).assert_agree()
        # rad P is indecomposable too on these; QPA agrees.
        A.crosscheck("indecomposable", M.radical()).assert_agree()
