"""Plan 74 -- the TWO BASES of the automorphism transport, and the non-vacuous
arbitration of its sign/variance conventions.

``_autom_action_on_classes`` builds the automorphism's action on the twisted bar
complex as ``M-slot (x) bar-block^{(x)n}``. The two slots live in DIFFERENT bases:

- the BAR slots are unit-adapted (``hbar`` is cut from ``h``'s unit-adapted
  matrix, dropping the unit row/column), but
- the COEFFICIENT slot is NOT: ``bar._coeff_in_unit_basis`` transports a
  coefficient bimodule with ``Bimodule.change_of_basis``, which by contract
  re-indexes only the A-element slot and leaves ``M``'s own basis alone, so ``M``
  stays in ``A``'s ORIGINAL basis.

The two matrices coincide whenever ``A`` is already unit-adapted (every local
algebra -- the dual numbers and ``k[x]/(x^3)`` batteries) or whenever ``h`` fixes
the unit-adapting vertex (a vertex-FIXING automorphism). They differ exactly when
a VERTEX-PERMUTING automorphism acts on a multi-vertex algebra, which is what this
file pins on both sides.

The convention arbiter below is non-vacuous BY CONSTRUCTION: on the 3-cycle
rotation ``hbar`` is neither symmetric nor self-inverse, so ``hbar``, ``hbar^T``,
``hbar^{-1}`` and ``(hbar^{-1})^T`` are four DISTINCT matrices, and the test
demands that exactly one of the eight (M-slot, bar-block) combinations is a chain
map on each side -- the shipped one. A test that merely asserted the shipped
combination works would pass under a symmetric involution too.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.fields import GF
from quiverlab.families.skew_group import (
    GroupAction, QuiverAutomorphism, skew_group_algebra)
from quiverlab.hochschild import bar
from quiverlab.hochschild import skew_group as sg
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.hochschild.skew_group import stefan_decomposition

MAX_CELLS = 8_000_000


def _rotation_3_cycle(field=GF(7)):
    """kZ_3/rad^2 with the Z/3 rotation: hbar's arrow block is a 3-cycle
    permutation -- non-symmetric AND non-self-inverse."""
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (1, 2), "b": (2, 3), "c": (3, 1)})
    A = Q.algebra(relations=["a*b", "b*c", "c*a"], field=field)
    rot = QuiverAutomorphism({1: 2, 2: 3, 3: 1},
                             {"a": "b", "b": "c", "c": "a"})
    return A, GroupAction.cyclic(rot, 3)


def _swap_two_cycle(field=GF(7)):
    """The two-cycle Nakayama with the vertex-SWAPPING involution."""
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2), "b": (2, 1)})
    A = Q.algebra(relations=["a*b", "b*a"], field=field)
    swap = QuiverAutomorphism({1: 2, 2: 1}, {"a": "b", "b": "a"})
    return A, GroupAction.cyclic(swap, 2)


def _eq(dom, X, Y):
    if len(X) != len(Y) or len(X[0]) != len(Y[0]):
        return False
    return all(dom.eq(X[i][j], Y[i][j])
               for i in range(len(X)) for j in range(len(X[0])))


# --------------------------------------------------------------- the two bases
@pytest.mark.oracle_selfcert
def test_coefficient_slot_and_bar_slot_bases_really_differ():
    """The premise of the arbiter: on a vertex-permuting action the two matrices
    are NOT equal (so choosing the wrong one is a real choice), while on a local
    algebra they ARE (so the whole distinction is invisible there)."""
    A, ga = _rotation_3_cycle()
    MhA = ga.matrix_of(1, A)
    MhB = sg._to_unit_basis(A, MhA)
    assert not _eq(A.domain, MhA, MhB)

    from quiverlab.families.basic import truncated_polynomial
    L = truncated_polynomial(3, field=GF(7))
    laut = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": 2})
    lga = GroupAction.cyclic(laut, 3)
    MlA = lga.matrix_of(1, L)
    assert L.is_unit_adapted
    assert _eq(L.domain, MlA, sg._to_unit_basis(L, MlA))


@pytest.mark.oracle_selfcert
def test_hbar_is_neither_symmetric_nor_self_inverse():
    """Guards the arbiter against silently going vacuous if the fixture changes."""
    A, ga = _rotation_3_cycle()
    dom = A.domain
    hbar = sg._hbar(dom, sg._to_unit_basis(A, ga.matrix_of(1, A)))
    assert not _eq(dom, hbar, sg._transpose(hbar))
    assert not _eq(dom, hbar, sg._matinv(dom, hbar))
    assert not _eq(dom, hbar, sg._transpose(sg._matinv(dom, hbar)))


# ------------------------------------------------- the chain-map convention arbiter
@pytest.mark.oracle_selfcert
def test_transport_convention_is_the_unique_chain_map():
    """EXACTLY ONE of the eight (M-slot, bar-block) combinations is a chain map on
    each side, and it is the shipped one: coefficient slot in A's own basis, with
    ``(hbar^{-1})^T`` for cohomology and ``hbar`` for homology."""
    A, ga = _rotation_3_cycle()
    dom = A.domain
    B = A.unit_adapted()
    MhA = ga.matrix_of(1, A)
    MhB = sg._to_unit_basis(A, MhA)
    hbar = sg._hbar(dom, MhB)

    blocks = {
        "hbar": hbar,
        "hbar^T": sg._transpose(hbar),
        "hbar^-1": sg._matinv(dom, hbar),
        "(hbar^-1)^T": sg._transpose(sg._matinv(dom, hbar)),
    }
    mslots = {"MhA": MhA, "MhB": MhB}

    for gi in (0, 1):                       # untwisted AND genuinely twisted
        Mb = bar._coeff_in_unit_basis(
            A, B, Bimodule.twisted(A, phi=ga.matrix_of(gi, A)))
        dco = {n: bar.coboundary_matrix(B, n, MAX_CELLS, coefficients=Mb)[0]
               for n in (0, 1)}
        dho = {n: bar.boundary_matrix(B, n, MAX_CELLS, coefficients=Mb)[0]
               for n in (1, 2)}

        coh_ok, hom_ok = set(), set()
        for mname, mslot in mslots.items():
            for bname, blk in blocks.items():
                H = {n: (sg._kron(dom, mslot, sg._kpow(dom, blk, n)) if n
                         else mslot) for n in (0, 1, 2)}
                mul = lambda X, Y: sg._mat_mul(dom, X, Y)  # noqa: E731
                if all(_eq(dom, mul(dco[n], H[n]), mul(H[n + 1], dco[n]))
                       for n in (0, 1)):
                    coh_ok.add((mname, bname))
                if all(_eq(dom, mul(dho[n], H[n]), mul(H[n - 1], dho[n]))
                       for n in (1, 2)):
                    hom_ok.add((mname, bname))

        assert coh_ok == {("MhA", "(hbar^-1)^T")}, coh_ok
        assert hom_ok == {("MhA", "hbar")}, hom_ok


# ----------------------------------------------- end-to-end vs the DIRECT engine
@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("side", ["coh", "hom"])
def test_vertex_permuting_action_agrees_with_direct(side):
    """The regression: with the unit-adapted matrix on the coefficient slot this
    transport left the cycle span and the whole decomposition refused."""
    A, ga = _swap_two_cycle()
    assert skew_group_algebra(A, ga).dim == 8
    rep = stefan_decomposition(A, ga, 2, side=side, verify_direct=True)
    assert rep.status == "complete"
    assert rep.agrees is True
    assert rep.dims == rep.direct_dims == [2, 1, 1]
