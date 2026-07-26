"""Cross-invariant Hochschild identity oracles pinned against the literature
(Plan 29, Part 2 items 1-2). Two *independent* invariants must agree:

* **Happel's trace identity** -- Happel, "The trace of the Coxeter matrix and
  Hochschild cohomology," Linear Algebra Appl. 258 (1997) 169-177 (registry key
  ``happel_trace``). For a finite-global-dimension algebra,
  ``tr(Coxeter) == -sum_i (-1)^i dim HH^i``. The sign is fixed by quiverlab's
  Coxeter convention ``Phi = -C^{-T} C`` and is PINNED on A_3 first
  (tr = -1, Euler characteristic = +1).

* **Derived invariance** of ``HH^*`` / ``HH_*`` / ``HC_*`` -- Keller, J. Pure
  Appl. Algebra 123 (1998) 223-273 (``keller_cyclic_invariance``); Rickard,
  "Morita theory for derived categories," J. London Math. Soc. 39 (1989)
  436-456 (``rickard_derived``). Two orientations of one underlying graph are
  derived equivalent, so all three invariants must coincide.

The trace matrix is integer (field-independent) and these HH values are
characteristic-independent, so the fast GF(p) rank engine is used off A_3.
Provenance: ``docs/plans/2026-07-25-literature-oracles-deep-research.md``
(Happel / Keller / Rickard cluster; values re-verified live here).
"""
import pytest

from quiverlab import CC, GF, IncidenceAlgebra, Quiver

pytestmark = [pytest.mark.oracle_literature]

F = GF(32003)               # the fast GF(p) rank engine (identities are char-independent)
_DIAMOND = [("b", "x"), ("b", "y"), ("x", "t"), ("y", "t")]   # 2x2 grid = comm. square


def _trace(A):
    M = A.coxeter_matrix()
    return sum(M[i][i] for i in range(len(M)))


def _euler(dims):
    return sum((-1) ** i * d for i, d in enumerate(dims))


# ------------------------------------------------------------------ 1a. Happel

def test_happel_trace_sign_pinned_on_A3():
    """Fix the sign FIRST on linear A_3: quiverlab's Coxeter matrix has trace -1,
    HH^* = [1,0,0] so the Euler characteristic is +1, and the verified identity
    is tr == -Euler (Happel 1997, ``happel_trace``)."""
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=CC)
    hh = A.hochschild_cohomology(2).dims
    assert hh == [1, 0, 0]
    assert _trace(A) == -1
    assert _euler(hh) == 1
    assert _trace(A) == -_euler(hh)          # the pinned identity


def test_happel_trace_hereditary_battery():
    """Hereditary members (gl.dim 1, so HH^{>=2} = 0). top=2 fits the default
    guard for each and OBSERVES the vanishing tail HH^2 = 0, making the Euler
    characteristic complete. A_4 (dim 10) needs a raised guard at top=2 and lives
    in the deep companion ``tests/engine/test_identity_oracles_deep.py``; here it
    is checked at top=1 (its Euler characteristic is complete by gl.dim 1)."""
    top2 = {
        "2-Kronecker": (Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}), [1, 3, 0]),
        "3-Kronecker": (Quiver([1, 2], {"a": (1, 2), "b": (1, 2), "c": (1, 2)}), [1, 8, 0]),
        "D4": (Quiver([0, 1, 2, 3], {"a": (0, 1), "b": (0, 2), "c": (0, 3)}), [1, 0, 0]),
        "acyclic-triangle": (Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3)}), [1, 2, 0]),
    }
    for name, (Q, expect) in top2.items():
        A = Q.algebra(relations=[], field=F)
        hh = A.hochschild_cohomology(2).dims
        assert hh == expect, name
        assert hh[2] == 0, f"{name}: vanishing tail confirms the Euler sum is complete"
        assert _trace(A) == -_euler(hh), name

    # A_4: top=1 under the default guard; gl.dim 1 => HH^{>=2} = 0 completes Euler.
    A4 = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(relations=[], field=F)
    hh = A4.hochschild_cohomology(1).dims
    assert hh == [1, 0]
    assert _trace(A4) == -_euler(hh)         # Euler = HH^0 - HH^1 (higher HH vanish)


def test_happel_trace_finite_gldim_nonhereditary_comm_square():
    """Commutative square = incidence algebra of the 2x2 grid poset, gl.dim 2 and
    NOT hereditary. HH vanishes above the global dimension, so HH^2 = 0 is the
    vanishing tail and top=2 captures the whole Euler characteristic. (top=3 is
    unnecessary and infeasible: the bar differential d^3 is 37k x 4.6k = 170M
    cells, far past the guard.)"""
    A = IncidenceAlgebra(_DIAMOND, field=F)
    hh = A.hochschild_cohomology(2).dims
    assert hh == [1, 0, 0]
    assert hh[2] == 0                         # tail vanishes -> Euler complete (gl.dim 2)
    assert _trace(A) == -_euler(hh)


# --------------------------------------------------- 1b. Derived invariance

def _triple(A, top=2):
    return (A.hochschild_cohomology(top).dims,
            A.hochschild_homology(top).dims,
            A.cyclic_homology(top).dims)


def test_derived_invariance_atilde2_triangle_orientations():
    """The two acyclic orientations of the A~_2 triangle (extended Dynkin A_2)
    are derived equivalent, so HH^* / HH_* / HC_* coincide. Pinned from the
    Happel/Keller/Rickard cluster: [1,2,0] / [3,0,0] / [3,0,3] (``keller_cyclic_
    invariance`` / ``rickard_derived``)."""
    o1 = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3)}).algebra(relations=[], field=F)
    o2 = Quiver([1, 2, 3], {"a": (2, 1), "b": (3, 2), "c": (3, 1)}).algebra(relations=[], field=F)
    t1, t2 = _triple(o1), _triple(o2)
    assert t1 == t2                                          # derived invariance
    assert t1 == ([1, 2, 0], [3, 0, 0], [3, 0, 3])          # pinned values


def test_derived_invariance_D4_orientations():
    """Two orientations of D_4 (a tree): both give the trivial HH^* = [1,0,0] and
    equal HH_* / HC_*. Pins the derived-invariance scheme on a Dynkin graph."""
    o1 = Quiver([0, 1, 2, 3], {"a": (0, 1), "b": (0, 2), "c": (0, 3)}).algebra(relations=[], field=F)  # center = source
    o2 = Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}).algebra(relations=[], field=F)  # center = sink
    t1, t2 = _triple(o1), _triple(o2)
    assert t1 == t2                                          # derived invariance
    assert t1[0] == [1, 0, 0]                                # trivial HH^*
