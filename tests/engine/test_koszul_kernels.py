"""The Koszul kernels K_n = intersection of V^i (x) R (x) V^j inside V^{(x)n}.

For a KOSZUL algebra, dim K_n is the n-th Koszul-dual Hilbert coefficient = the total
graded Betti number = the rank of the n-th term of the minimal bimodule resolution
(GHMS 2005). Plan 75 / R10.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.families import IncidenceAlgebra
from quiverlab.families.exterior import ExteriorAlgebra
from quiverlab.families.preprojective import PreprojectiveAlgebra
from quiverlab.fields import GF, QQ
from quiverlab.hochschild.koszul_ghms import koszul_betti, koszul_kernels

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("m, expect", [
    (2, [1, 2, 3, 4, 5, 6, 7]),               # C(n+1, n)   -- LIVE-VERIFIED
    (3, [1, 3, 6, 10, 15, 21, 28]),           # C(n+2, 2)   -- LIVE-VERIFIED
])
def test_exterior_kernels_are_koszul_dual_hilbert(m, expect):
    assert koszul_betti(ExteriorAlgebra(m, field=GF(32003)), 6) == expect


@lit
def test_exterior_kernels_are_domain_general():
    """The GHMS route runs over ANY exact Domain -- unlike the GF(p)-only syzygy engine."""
    for dom in (QQ, GF(2), GF(32003)):
        assert koszul_betti(ExteriorAlgebra(2, field=dom), 5) == [1, 2, 3, 4, 5, 6]


@xeng
def test_kernels_match_ext_algebra_betti():
    """dim K_n == the total graded Betti number from the INDEPENDENT Plan-27 Ext-algebra
    (Yoneda) route."""
    A = ExteriorAlgebra(2, field=QQ)
    hm = A.ext_algebra(6).hilbert_matrix_through(6)
    betti = [sum(sum(row) for row in hm[n]) for n in range(7)]
    assert koszul_betti(A, 6) == betti == [1, 2, 3, 4, 5, 6, 7]


@lit
def test_diamond_incidence_bridge():
    """The commutative-square incidence algebra: Koszul AND an incidence algebra -- the
    (a)<->(b) bridge of the plan. K_0 = |Q_0| = 4 (NOT 1 -- the multi-vertex point),
    K_1 = #arrows = 4, K_2 = the single commutativity relation, then 0 (gl.dim 2)."""
    D = IncidenceAlgebra([("0", "a"), ("0", "b"), ("a", "1"), ("b", "1")], field=QQ)
    assert D.dim == 9
    assert koszul_betti(D, 5) == [4, 4, 1, 0, 0, 0]


@lit
def test_cyclic_nakayama_kernels_are_periodic():
    """kZ_3/rad^2 -- multi-vertex, Koszul, self-injective: the kernels are PERIODIC of
    constant rank 3 (the resolution never terminates), unlike the diamond's gl.dim 2."""
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)})
    Z = Q.algebra(relations=["a*b", "b*c", "c*a"], field=QQ)
    assert Z.dim == 6
    assert koszul_betti(Z, 5) == [3, 3, 3, 3, 3, 3]


@selfcert
def test_preprojective_a3_is_the_negative_example():
    """Pi(A_3) is NOT Koszul -- ext_algebra reports the POSITIVE obstruction "a new
    Ext-algebra generator appears in degree 3" (Brenner-Butler-King almost-Koszul: Dynkin
    preprojectives break at h - 1 = 3 for A_3).

    And the kernels show CONCRETELY why the Koszul hypothesis is not decoration: they
    VANISH from degree 3 ([3,4,3,0,0,0]) while the true minimal bimodule ranks are
    [3,4,3,3,4,3,3] (syzygy engine). So a GHMS resolution built from these kernels would
    silently claim gl.dim 2 for a SELF-INJECTIVE algebra of infinite global dimension --
    which is exactly why `engine="ghms"` must REFUSE a non-Koszul input rather than
    compute something.
    """
    P = PreprojectiveAlgebra("A3", field=QQ)
    assert P.dim == 10
    ea = P.ext_algebra(4)
    assert ea.koszul is False
    assert ea.koszul_obstruction[0] == 3
    assert "generator" in ea.koszul_obstruction[1]
    assert koszul_betti(P, 5) == [3, 4, 3, 0, 0, 0]      # kernels DIE; the resolution does not


@selfcert
def test_kernels_are_nested_and_corner_typed():
    """K_n sits inside BOTH V (x) K_{n-1} and K_{n-1} (x) V (that is what makes the GHMS
    differential comultiplicative), so dim K_n <= dim V * dim K_{n-1} always."""
    for A in (ExteriorAlgebra(2, field=QQ),
              IncidenceAlgebra([("0", "a"), ("0", "b"), ("a", "1"), ("b", "1")], field=QQ)):
        b = koszul_betti(A, 5)
        n_arrows = len(A.quiver.arrows)
        for n in range(2, len(b)):
            assert b[n] <= n_arrows * b[n - 1]


@selfcert
def test_kernel_vectors_are_independent():
    K = koszul_kernels(ExteriorAlgebra(3, field=QQ), 4)
    from quiverlab.fields.linalg import rank
    for n, basis in enumerate(K):
        if basis:
            assert rank(basis, QQ) == len(basis), f"K_{n} basis is dependent"


@selfcert
def test_non_quadratic_refuses():
    from quiverlab.errors import QuiverlabError
    Q = Quiver([1], {"x": (1, 1)})
    A = Q.algebra(relations=["x*x*x"], field=QQ)          # cubic: NOT quadratic
    with pytest.raises(QuiverlabError, match="QUADRATIC|quadratic"):
        koszul_kernels(A, 3)
