"""Module Tor_n^A(M, N) for a RIGHT A-module M and a LEFT A-module N (Plan 29 Part 4).

Self-certifying + literature oracles:
  (i)   DUALITY (the anchor, asserted on EVERY case): dim Tor_n^A(M, N) =
        dim Ext_A^n(M, DN), with D side-aware (Plan 24: DN of a left module is a RIGHT
        module). Every Tor value is cross-checked against the EXISTING, independently
        written Ext engine (modules/ext.py). Hom-tensor duality
        Hom_A(P_i, DN) ~= D(P_i (x)_A N) makes the two complexes k-dual, so their
        (co)homology dimensions agree degreewise.
  (ii)  BALANCE: resolving M vs resolving N over A^op agree (``resolve="second"``).
  (iii) hereditary => Tor_{>=2} = 0; semisimple => Tor_{>0} = 0; projective (flat) M =>
        Tor_{>0} = 0 with Tor_0 = e_v N; k[x]/(x^2) => Tor_n(k, k) = k for all n; the
        quantum complete intersection => Tor_n(k, k) = n + 1 (the quantum-Koszul Betti,
        independent of q); kA_2 / kA_3 simples worked BY HAND.

Sources: Cartan-Eilenberg, "Homological Algebra", Princeton Univ. Press (1956)
(citations key ``tensor_product``); Assem-Simson-Skowronski, "Elements of the
Representation Theory of Associative Algebras, Vol. 1", Cambridge Univ. Press (2006),
Ch. III (key ``assem_book``).
"""
import pytest

from quiverlab import CC, GF, Quiver, QuantumCI, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.tor import tor_dims


def _anchor(A, M, N, top):
    """Return tor_dims(A, M, N, top) after asserting the duality oracle
    dim Tor_n(M, N) = dim Ext^n(M, DN) degreewise against the Ext engine."""
    tors = tor_dims(A, M, N, top)
    exts = [A.ext(M, N.dualize(), n) for n in range(top + 1)]
    assert tors == exts, (getattr(M, "name", "M"), getattr(N, "name", "N"), tors, exts)
    return tors


# -- hand-worked kA_2 (1 --a--> 2) ------------------------------------------
def test_tor_kA2_hand_worked():
    """kA_2 simples worked by hand from the resolution 0 -> P_2 -> P_1 -> S_1 -> 0 and
    P_2 = S_2 projective:  Tor(S_1^r, S_1^l) = [1,0], Tor(S_1^r, S_2^l) = [0,1]
    (Ext^1(S_1, S_2) = 1, so Tor_1(S_1, S_2^left) = 1 by duality; D S_j^left = S_j^right).
    Tor_0 = M (x)_A N: S_1 (x) S_2 = 0 while S_1 (x) S_1 = k."""
    A = linear_path_algebra(2)
    S1r, S2r = A.simple(1), A.simple(2)
    S1l, S2l = A.simple(1, side="left"), A.simple(2, side="left")
    assert tor_dims(A, S1r, S1l, 2) == [1, 0, 0]
    assert tor_dims(A, S1r, S2l, 2) == [0, 1, 0]
    assert tor_dims(A, S2r, S2l, 2) == [1, 0, 0]        # S_2 = P_2 projective, Tor_0 = e_2 N
    assert tor_dims(A, S2r, S1l, 2) == [0, 0, 0]        # e_2(S_1^left) = 0
    # duality anchor + single-value delegator surface
    assert _anchor(A, S1r, S2l, 3) == [0, 1, 0, 0]
    assert A.tor(S1r, S2l, 1) == 1
    assert A.tor(S1r, S1l, 0) == 1


def test_tor_projective_is_flat():
    """A projective right module is flat: Tor_{>0}(P_v, N) = 0 and Tor_0 = e_v N =
    (dim-vector of N at v). Checked over every vertex x every left simple of kA_3."""
    A = linear_path_algebra(3)
    for v in (1, 2, 3):
        P = A.projective(v)
        for w in (1, 2, 3):
            Nl = A.simple(w, side="left")
            tors = _anchor(A, P, Nl, 3)
            assert tors[0] == Nl.dimension_vector()[v]     # Tor_0 = e_v N
            assert tors[1:] == [0, 0, 0]                    # flat


@pytest.mark.parametrize("field", [QQ, GF(2), GF(7), CC])
def test_tor_kx2_all_ones(field):
    """k[x]/(x^2) is self-injective with pd(k) = infinity: the minimal resolution of k
    is ... -> P -> P -> P -> k with every differential = mult-by-x, and x acts as 0 on
    k, so Tor_n(k, k) = k for all n (dims all 1). Matches Ext by self-duality. Spread
    over QQ / GF(2) / GF(7) / CC."""
    B = truncated_polynomial(2, field=field)
    kr, kl = B.simple(1), B.simple(1, side="left")
    assert _anchor(B, kr, kl, 5) == [1, 1, 1, 1, 1, 1]


@pytest.mark.parametrize("field", [QQ, GF(2), GF(7), CC])
def test_tor_duality_anchor_kA3_spread(field):
    """The duality anchor across kA_3 simples over a field spread; kA_3 is hereditary
    so Tor_{>=2} = 0, and the only nonvanishing higher Tor is Tor_1(S_i, S_{i+1}^left)
    (from Ext^1(S_i, S_{i+1}) = 1)."""
    A = linear_path_algebra(3, field=field)
    Sr = {i: A.simple(i) for i in (1, 2, 3)}
    Sl = {i: A.simple(i, side="left") for i in (1, 2, 3)}
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            tors = _anchor(A, Sr[i], Sl[j], 4)
            assert tors[0] == (1 if i == j else 0)         # Tor_0 = Hom side: S_i (x) S_j
            assert tors[2:] == [0, 0, 0]                    # hereditary
    assert tor_dims(A, Sr[1], Sl[2], 2) == [0, 1, 0]
    assert tor_dims(A, Sr[2], Sl[3], 2) == [0, 1, 0]
    assert tor_dims(A, Sr[1], Sl[3], 2) == [0, 0, 0]        # no arrow 1 -> 3


def test_tor_hereditary_vanishing():
    """Hereditary (gl.dim <= 1) => Tor_{>=2} = 0 for ALL module pairs, including
    non-simple ones."""
    A = linear_path_algebra(3)
    M = A.module({1: 1, 2: 1, 3: 0}, {"a1": [[0, 0], [1, 0]], "a2": [[0, 0], [0, 0]]},
                 name="I12")                                # interval [1,2], right
    for Nl in (A.simple(2, side="left"), A.projective(2, side="left"),
               A.injective(3, side="left")):
        tors = _anchor(A, M, Nl, 4)
        assert tors[2:] == [0, 0, 0]


def test_tor_semisimple_vanishing():
    """A semisimple algebra k x k (two vertices, no arrows): every module is projective,
    so Tor_{>0} = 0 and Tor_0 = e_v N is the ordinary tensor over the semisimple base."""
    SS = Quiver([1, 2], {}).algebra(relations=[])
    for v in (1, 2):
        Sr = SS.simple(v)
        for w in (1, 2):
            Sl = SS.simple(w, side="left")
            tors = _anchor(SS, Sr, Sl, 3)
            assert tors[0] == (1 if v == w else 0)
            assert tors[1:] == [0, 0, 0]


def test_tor_balance_resolve_second():
    """Balance of Tor: resolving M (default) vs resolving N over A^op agree degreewise.
    Tor_n^A(M, N) = Tor_n^{A^op}(N_as_right_over_Aop, M_as_left_over_Aop)."""
    A2 = linear_path_algebra(2)
    S1r, S2l = A2.simple(1), A2.simple(2, side="left")
    assert tor_dims(A2, S1r, S2l, 3) == tor_dims(A2, S1r, S2l, 3, resolve="second")

    B = truncated_polynomial(2, field=GF(5))
    kr, kl = B.simple(1), B.simple(1, side="left")
    assert tor_dims(B, kr, kl, 4) == tor_dims(B, kr, kl, 4, resolve="second")

    A3 = linear_path_algebra(3)
    M = A3.module({1: 1, 2: 1, 3: 0}, {"a1": [[0, 0], [1, 0]], "a2": [[0, 0], [0, 0]]})
    Nl = A3.injective(3, side="left")
    assert tor_dims(A3, M, Nl, 3) == tor_dims(A3, M, Nl, 3, resolve="second")


def test_tor_quantum_ci():
    """Quantum complete intersection A = k<x,y>/(x^2, y^2, xy + q yx): the trivial module
    k has the quantum-Koszul minimal resolution with Betti numbers n + 1, so
    Tor_n(k, k) = n + 1 independent of q (checked for q = 2, 3 over GF(7) and q = -1 over
    GF(3)); anchored against Ext via duality."""
    for q, field in [(2, GF(7)), (3, GF(7)), (-1, GF(3))]:
        Q = QuantumCI(q, field=field)
        kr, kl = Q.simple(1), Q.simple(1, side="left")
        assert _anchor(Q, kr, kl, 5) == [1, 2, 3, 4, 5, 6]


def test_tor_from_arrow_action_pair():
    """An arbitrary from_arrow_action pair on kA_3: a non-simple right module M (interval
    [1,2]) against a non-simple left module N built with side="left" (interval [2,3] on
    A^op). Anchored via duality only."""
    A = linear_path_algebra(3)
    M = A.module({1: 1, 2: 1, 3: 0}, {"a1": [[0, 0], [1, 0]], "a2": [[0, 0], [0, 0]]},
                 name="I12")
    N = A.module({1: 0, 2: 1, 3: 1}, {"a1": [[0, 0], [0, 0]], "a2": [[0, 1], [0, 0]]},
                 side="left", name="L23")
    assert M.side == "right" and N.side == "left"
    tors = _anchor(A, M, N, 3)
    assert tors == [0, 1, 0, 0]                             # nonvanishing, hereditary
    # and against the left projective/injective builders too
    _anchor(A, M, A.projective(2, side="left"), 3)


def test_tor_wrong_side_and_cross_algebra_refuse_loudly():
    """Tor pairs a RIGHT module with a LEFT module over ONE fixed base algebra; every
    other combination is a category error, raised loudly (mirrors hom._assert_comparable),
    never a silent 0."""
    A = linear_path_algebra(3)
    B = linear_path_algebra(2)
    with pytest.raises(QuiverlabError):                     # right, right
        A.tor(A.simple(1), A.simple(2), 0)
    with pytest.raises(QuiverlabError):                     # left, left
        A.tor(A.simple(1, side="left"), A.simple(2, side="left"), 0)
    with pytest.raises(QuiverlabError):                     # left as first arg
        A.tor(A.simple(1, side="left"), A.simple(2, side="left"), 0)
    with pytest.raises(QuiverlabError):                     # right M over the wrong algebra
        A.tor(B.simple(1), A.simple(1, side="left"), 0)
    with pytest.raises(QuiverlabError):                     # left N over the wrong algebra
        A.tor(A.simple(1), B.simple(1, side="left"), 0)
    with pytest.raises(QuiverlabError):                     # bad resolve kwarg
        tor_dims(A, A.simple(1), A.simple(1, side="left"), 2, resolve="both")
