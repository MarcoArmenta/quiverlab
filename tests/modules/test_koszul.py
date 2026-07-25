"""Koszulity, quadraticity, the quadratic dual, and the Fröberg identity (Plan 27).

Theory/literature oracles (all run without [qpa]):
- quadraticity + G-quadratic certifier on the oracle battery (Priddy, `priddy`);
- the quadratic dual A^! and the arrow-reversal it lives on, both directions plus
  the double dual (Polishchuk-Positselski, `polishchuk_positselski`);
- the Fröberg numeric Koszulity test P(t)*C_A(-t) = I, hand-verified on kA_2 and
  failing at a computed degree on k[x]/x^3 (Fröberg, `froberg_koszul`);
- the E(A) = (A^!)^op dimension cross-check, corner transpose included.
"""
import pytest

from quiverlab import (
    CC, GF, Quiver, QuantumCI, linear_path_algebra, truncated_polynomial,
)
from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.modules.koszul import (
    QuadraticDual, dual_dims_crosscheck, froberg_obstruction,
    g_quadratic_certificate, is_quadratic, quadratic_dual,
)


# -- fixtures --------------------------------------------------------------
def _radsq_line(field=CC):
    """rad^2 = 0 on 1 -> 2 -> 3 (a: 1->2, b: 2->3), relation a*b."""
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=field)


def _hereditary_line(field=CC):
    """hereditary path algebra kQ on 1 -> 2 -> 3 (no relations)."""
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=field)


def _square(field=CC):
    """commutative square, relation a*b - c*d (quadratic, gl.dim 2)."""
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _cubic_relation(field=CC):
    """a genuinely cubic minimal relation a*b*c on 1->2->3->4 (NOT quadratic)."""
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=["a*b*c"], field=field)


def _quadratic_cubic_tip(field=GF(7)):
    """A finite-dimensional (dim 4) QUADRATIC algebra whose length-lex Gröbner basis
    acquires a CUBIC tip x^3: one vertex, loops x, y, relations x*y, y*y, y*x - x*x.
    Completion of {xy->0, yy->0, yx->x^2} resolves the overlap x*y*x to x^3, forcing
    x^3 -> 0. So it is quadratic (presentation) but NOT G-quadratic (a length-3 tip)."""
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*y", "y*y", "y*x - x*x"], field=field)


# -- is_quadratic ----------------------------------------------------------
def test_is_quadratic_battery():
    assert is_quadratic(truncated_polynomial(2)) is True          # k[x]/x^2
    assert is_quadratic(truncated_polynomial(3)) is False         # k[x]/x^3 (cubic relation)
    assert is_quadratic(_radsq_line()) is True                    # rad^2 = 0
    assert is_quadratic(_hereditary_line()) is True               # hereditary kQ (no relations)
    assert is_quadratic(QuantumCI(3)) is True                     # tips {x^2, y^2, yx}
    assert is_quadratic(_square()) is True                        # commutative square
    assert is_quadratic(_cubic_relation()) is False               # length-3 minimal relation


def test_is_quadratic_needs_presentation():
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(QuiverlabError):
        is_quadratic(A)


# -- G-quadratic certifier -------------------------------------------------
def test_g_quadratic_certifier_true_on_quadratic_battery():
    # Every quadratic battery member has a confluent length-2 Gröbner basis -> Koszul.
    assert g_quadratic_certificate(truncated_polynomial(2)) is True
    assert g_quadratic_certificate(_radsq_line()) is True
    assert g_quadratic_certificate(_hereditary_line()) is True    # empty rule set, vacuously
    assert g_quadratic_certificate(QuantumCI(3)) is True
    assert g_quadratic_certificate(_square()) is True


def test_g_quadratic_false_is_inconclusive_not_a_disproof():
    # A quadratic (presentation) algebra whose Gröbner basis has a cubic tip: the PBW
    # certificate FAILS, but this is inconclusive -- it disproves nothing.
    A = _quadratic_cubic_tip()
    assert A.dim == 4                          # finite-dimensional
    assert is_quadratic(A) is True             # quadratic by the presentation
    assert g_quadratic_certificate(A) is False  # cubic tip x^3 -> not G-quadratic


# -- quadratic dual --------------------------------------------------------
def test_dual_of_kx_x2_is_free_loop_ky():
    # dual(k[x]/x^2): R = span{x^2} is the whole length-2 space, so R^perp = 0 and
    # A^! = k[y] (free loop) -- infinite-dimensional -> QuadraticDual record, all
    # graded dims 1.
    D = quadratic_dual(truncated_polynomial(2))
    assert isinstance(D, QuadraticDual)
    assert D.relations == []                             # no relations on the free loop
    assert D.graded_dims(6) == [1, 1, 1, 1, 1, 1, 1]


def test_dual_of_radsq_is_free_hereditary_kQ_op():
    # dual(kQ/J^2): R = all length-2 paths, R^perp = 0 -> A^! = kQ^op FREE (no
    # relations). On the acyclic line it is finite-dimensional -> a genuine Algebra.
    A = _radsq_line()
    D = quadratic_dual(A)
    assert isinstance(D, Algebra)
    assert (D.relations or []) == []                     # zero relations: path algebra kQ^op
    # kQ^op on 1 <- 2 <- 3 keeps the length-2 path b*a (3 -> 1): dim = 3 + 2 + 1 = 6.
    assert D.dim == 6


def test_dual_of_hereditary_is_radical_square_zero():
    # dual(kQ): R = 0, so R^perp = all length-2 Q-paths -> A^! = kQ^op/J^2. On the
    # line 1->2->3 the only length-2 path a*b reverses to b*a, killed in the dual.
    A = _hereditary_line()
    D = quadratic_dual(A)
    assert isinstance(D, Algebra)
    assert len(D.relations or []) == 1                   # one relation (the reversed a*b)
    assert D.dim == 5                                    # 3 idempotents + 2 arrows


def test_double_dual_recovers_dimensions():
    # (A^!)^! recovers A (dims) on finite acyclic quadratic algebras, both directions.
    for A in (_radsq_line(GF(7)), _hereditary_line(GF(7))):
        D = quadratic_dual(A)
        assert isinstance(D, Algebra)
        DD = quadratic_dual(D)
        assert isinstance(DD, Algebra)
        assert DD.dim == A.dim


def test_dual_of_quantum_ci_is_the_quantum_plane():
    # dual(k<x,y>/(x^2, y^2, yx - q xy)) = the quantum plane k<x,y>/(yx - q' xy):
    # 2 generators, 1 relation -> infinite-dimensional, dim (A^!)_n = n + 1.
    D = quadratic_dual(QuantumCI(3))
    assert isinstance(D, QuadraticDual)
    assert len(D.relations) == 1
    assert D.graded_dims(6) == [1, 2, 3, 4, 5, 6, 7]


def test_quadratic_dual_refuses_non_quadratic():
    with pytest.raises(QuiverlabError):
        quadratic_dual(truncated_polynomial(3))          # k[x]/x^3 is not quadratic
    with pytest.raises(QuiverlabError):
        quadratic_dual(_cubic_relation())


# -- Fröberg numeric Koszulity test ---------------------------------------
def _hilbert_all_ones(d):
    """One-vertex Hilbert matrices with dim E^n = 1 for all n (the 1x1 record [[1]])."""
    return [[[1]] for _ in range(d + 1)]


def test_froberg_passes_on_koszul_kx_x2_through_6():
    # E(k[x]/x^2) = k[y], dim E^n = 1 all n; C_A = 1 + t. P*C_A(-t) = 1 -> no obstruction.
    A = truncated_polynomial(2)
    assert froberg_obstruction(A, _hilbert_all_ones(6), 6) is None


def test_froberg_passes_on_kA2_through_6():
    # kA_2: C = [[1, t], [0, 1]], P = [[1, t], [0, 1]] (Ext^1(S_1, S_2) = 1, E^{>=2} = 0).
    A = linear_path_algebra(2)
    zero = [[0, 0], [0, 0]]
    E0 = [[1, 0], [0, 1]]
    E1 = [[0, 1], [0, 0]]
    P = [E0, E1] + [zero] * 5
    assert froberg_obstruction(A, P, 6) is None


def test_froberg_fails_on_kx_x3_at_degree_2():
    # E(k[x]/x^3) = k[y,z]/(y^2), dim E^n = 1 all n; C_A = 1 + t + t^2.
    # M_2 = P_0 C_2 - P_1 C_1 + P_2 C_0 = 1 - 1 + 1 = 1 != 0 -> the identity FAILS at
    # degree 2 (certified NOT Koszul over this field).
    A = truncated_polynomial(3)
    assert froberg_obstruction(A, _hilbert_all_ones(6), 6) == 2


def test_froberg_needs_the_presentation():
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(QuiverlabError):
        froberg_obstruction(A, [[[1]]], 0)


# -- E(A) = (A^!)^op dimension cross-check --------------------------------
def _kQ_hilbert_line():
    """The E = kQ Hilbert matrices for the rad^2 = 0 line 1 -> 2 -> 3 (vertices in
    order 1, 2, 3). E^n[i][j] = #paths_n(Q) i -> j: E^0 = I, E^1 has a (1->2) and
    b (2->3), E^2 has the length-2 path a*b (1->3), E^{>=3} = 0."""
    E0 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    E1 = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
    E2 = [[0, 0, 1], [0, 0, 0], [0, 0, 0]]
    zero = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    return [E0, E1, E2, zero]


def test_ext_dims_of_radsq_line_are_path_counts():
    # rad^2 = 0 on 1->2->3 has E = kQ, so dim E^n = #paths_n(Q) = [3, 2, 1, 0, ...].
    hil = _kQ_hilbert_line()
    totals = [sum(sum(row) for row in mat) for mat in hil]
    assert totals == [3, 2, 1, 0]


def test_dual_dims_crosscheck_radsq_line_E_equals_kQ():
    # E(A) = (A^!)^op: A^! = kQ^op (free), and the corner transpose closes the
    # identification E^n[i][j] == D_n[j][i]. E is kQ, NOT kQ^op.
    A = _radsq_line()
    assert dual_dims_crosscheck(A, _kQ_hilbert_line(), 3) is True


def test_dual_dims_crosscheck_detects_a_wrong_table():
    # Corrupt one corner: the cross-check must reject it (not silently pass).
    A = _radsq_line()
    hil = _kQ_hilbert_line()
    hil[2][0][2] = 5                                     # the surviving a*b path had dim 1
    assert dual_dims_crosscheck(A, hil, 3) is False


def test_dual_dims_crosscheck_needs_the_transpose():
    # Feeding the TRANSPOSED (un-op'd) table must fail -- pins the corner convention.
    A = _radsq_line()
    hil = _kQ_hilbert_line()
    transposed = [[[mat[j][i] for j in range(3)] for i in range(3)] for mat in hil]
    assert dual_dims_crosscheck(A, transposed, 3) is False


def test_dual_dims_crosscheck_none_when_not_quadratic():
    assert dual_dims_crosscheck(truncated_polynomial(3), _hilbert_all_ones(6), 6) is None


# -- domain spread ---------------------------------------------------------
def test_quadratic_and_dual_over_gfp_and_char0_agree():
    # The dual presentation is field-agnostic on these battery members; graded dims
    # match across GF(2), GF(3), and char 0.
    for field in (GF(2), GF(3), CC):
        assert quadratic_dual(_radsq_line(field)).dim == 6
        assert quadratic_dual(truncated_polynomial(2, field=field)).graded_dims(5) == [1] * 6
