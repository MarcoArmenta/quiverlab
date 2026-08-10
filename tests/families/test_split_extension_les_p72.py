"""Plan 72 part (a): the split-extension / trivial-extension Hochschild LES (R5).

Oracles (Plan-32 markers):
- oracle_literature: frozen HH^•(T(kA_n))/HH^•(T(kD4)) dims; the flank pins; the
  HH^1(T(B)) != 0 sweep + directed formula.
- oracle_crossengine: assembled HH^•(L) == direct HH^•(L) (the LES route vs. the
  standalone CS engine); the leading-piece identity + summand inequality.
- oracle_selfcert: LES exactness / SES dim identity; the snake adversarial-lift;
  Bimodule.inflate().check(); dim T(B) = 2 dim B.

**Engine note (declared deviation).** The snake is computed on the Chouhy-Solotar
L^e-projective resolution (block-partitioned by the coefficient split L = M (+) B),
NOT the bar coboundary matrices: bar is intractable at these depths (T(kD4) delta^4
lands in a >2M-row bar group). CS reproduces every live-verified pin exactly.
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.fields import QQ
from quiverlab.families import TrivialExtension
from quiverlab.hochschild.coefficients import Bimodule
from quiverlab.hochschild import split_extension as _se


# ------------------------------------------------------------------ builders
def kA2(field=QQ):
    return Quiver([1, 2], {"a1": (1, 2)}).algebra(field=field)


def kA3(field=QQ):
    return Quiver([1, 2, 3], {"a1": (1, 2), "a2": (2, 3)}).algebra(field=field)


def kD4(field=QQ):
    return Quiver([1, 2, 3, 4], {"a2": (2, 1), "a3": (3, 1), "a4": (4, 1)}).algebra(field=field)


def kron(field=QQ):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=field)


def kx2(field=QQ):
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x^2"], field=field)


def kx3(field=QQ):
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x^3"], field=field)


# =====================================================================
# Task A1 -- Bimodule.inflate + split_extension constructor + the flanks
# =====================================================================
@pytest.mark.oracle_selfcert
def test_inflate_check_and_dims_kA2():
    B = kA2()
    L = TrivialExtension(B)
    van = [a for a in L.quiver.arrows if a not in B.quiver.arrows]
    Md = Bimodule.inflate(L, B, Bimodule.dual(B), vanishing_arrows=van)
    Mr = Bimodule.inflate(L, B, Bimodule.regular(B), vanishing_arrows=van)
    assert Md.check() is True and Md.dim_M == 3
    assert Mr.check() is True and Mr.dim_M == 3


@pytest.mark.oracle_selfcert
def test_inflate_multivertex_kD4_check():
    B = kD4()
    L = TrivialExtension(B)
    van = [a for a in L.quiver.arrows if a not in B.quiver.arrows]
    assert len(van) == 3                         # three dual arrows
    Md = Bimodule.inflate(L, B, Bimodule.dual(B), vanishing_arrows=van)
    assert Md.check() is True and Md.dim_M == B.dim == 7


@pytest.mark.oracle_selfcert
def test_split_extension_is_trivial_extension():
    B = kA2()
    L = _se.split_extension(B)
    T = TrivialExtension(B)
    assert L.dim == 6 == 2 * B.dim
    assert set(L.quiver.arrows) == set(T.quiver.arrows) == {"a1", "te0"}


@pytest.mark.oracle_selfcert
def test_split_extension_general_M_refused():
    B = kA2()
    with pytest.raises(Exception):
        _se.split_extension(B, Bimodule.dual(B))     # general M is out of scope (DD-A1)


@pytest.mark.oracle_literature
def test_flank_pins_kA2():
    rep = _se.split_extension_cohomology(kA2(), 5)
    assert rep.flank_M == [2, 1, 0, 1, 2, 1, 0]      # 0..6 (M-flank to top+1)
    assert rep.flank_B == [1, 0, 1, 2, 1, 0]         # 0..5


@pytest.mark.oracle_literature
def test_flank_pins_kD4_boundary_rule():
    # The BOUNDARY-RULE pin: HH^5(L,M) = 4 != 0, HH^4(L,B) = 4 != 0.
    rep = _se.split_extension_cohomology(kD4(), 4)
    assert rep.flank_M == [4, 1, 0, 0, 1, 4]         # 0..5
    assert rep.flank_B == [1, 0, 0, 1, 4]            # 0..4


@pytest.mark.oracle_crossengine
def test_inflate_flank_equals_report_flank():
    """The Bimodule.inflate(dual) flank == the report's flank_M (the CS partition):
    inflate(D(B)) and the ideal-M are the SAME L-bimodule up to iso."""
    B = kA2()
    L = TrivialExtension(B)
    van = [a for a in L.quiver.arrows if a not in B.quiver.arrows]
    Md = Bimodule.inflate(L, B, Bimodule.dual(B), vanishing_arrows=van)
    flank_via_inflate = L.hochschild_cohomology(6, coefficients=Md, engine="cs").dims
    assert flank_via_inflate == _se.split_extension_cohomology(B, 5).flank_M


# =====================================================================
# Task A2 -- SES / snake delta / assembly / exactness
# =====================================================================
@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("B,top,direct,delta", [
    (kA2(), 5, [3, 1, 1, 1, 1, 1], [0, 0, 0, 2, 0, 0]),
    (kA3(), 4, [4, 1, 1, 1, 1], [0, 0, 0, 0, 0]),
    (kD4(), 4, [5, 1, 0, 0, 2], [0, 0, 0, 1, 2]),   # the boundary-rule case
])
def test_assembled_equals_direct(B, top, direct, delta):
    rep = _se.split_extension_cohomology(B, top)
    assert rep.direct == direct
    assert rep.assembled == direct
    assert rep.agrees is True
    assert rep.delta_ranks == delta


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("B,top", [(kA2(), 5), (kA3(), 4), (kD4(), 4)])
def test_les_exactness_and_ses_dim_identity(B, top):
    rep = _se.split_extension_cohomology(B, top)
    assert rep.exact is True


@pytest.mark.oracle_selfcert
def test_snake_adversarial_lift_well_defined():
    # shifting the lift by any M-cochain leaves delta[phi]'s class unchanged
    assert _se.snake_lift_well_defined(kA2(), 4) is True
    assert _se.snake_lift_well_defined(kD4(), 3) is True


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("B,top,direct", [
    (kA2(), 4, [3, 1, 1, 1, 1]),
    (kA3(), 4, [4, 1, 1, 1, 1]),
])
def test_homology_twin_assembled_equals_direct(B, top, direct):
    rep = _se.split_extension_homology(B, top)
    assert rep.direct == direct
    assert rep.assembled == direct
    assert rep.agrees is True
    assert rep.exact is True


# =====================================================================
# Task A3 -- HH^1 != 0 witness, directed formula, leading pieces, scale
# =====================================================================
@pytest.mark.oracle_literature
@pytest.mark.parametrize("B,directed,hh1L,formula", [
    (kA2(), True, 1, True),
    (kA3(), True, 1, True),
    (kD4(), True, 1, True),
    (kron(), True, 4, True),
    (kx2(), False, 4, None),
    (kx3(), False, 7, None),
])
def test_hh1_grading_witness_sweep(B, directed, hh1L, formula):
    w = _se.hh1_grading_witness(B)
    assert w["is_cocycle"] is True and w["is_outer"] is True
    assert w["nonzero"] is True                  # HH^1(L) != 0 always (M != 0)
    assert w["hh1_L"] == hh1L
    assert w["directed"] is directed
    assert w["formula"] is formula               # 1 + HH^1(B) on directed B; None else


@pytest.mark.oracle_crossengine
def test_leading_piece_identity_and_inequality():
    # p = 0 leading piece: the n=0 identity + the summand inequality (NOT equality).
    for B, top in [(kA2(), 5), (kA3(), 4), (kD4(), 4)]:
        rep = _se.split_extension_cohomology(B, top)
        assert rep.flank_B[0] == rep.leading_B[0]          # identity at n=0
        assert rep.flank_M[0] == rep.leading_M[0]
        for n in range(top + 1):
            assert rep.flank_B[n] >= rep.leading_B[n]       # summand inequality
            assert rep.flank_M[n] >= rep.leading_M[n]


@pytest.mark.oracle_crossengine
def test_leading_piece_tight_on_2kronecker():
    # T(2-Kronecker): HH^1(B) = 3 = flank_B[1] -- the inequality is TIGHT/non-vacuous.
    B = kron()
    rep = _se.split_extension_cohomology(B, 3)
    assert rep.leading_B == [1, 3, 0, 0]
    assert rep.flank_B == [1, 3, 3, 5]
    assert rep.flank_B[1] == rep.leading_B[1] == 3


@pytest.mark.oracle_literature
def test_trivial_extension_hh_dims_pins():
    # frozen HH^•(T(B)) via the assembled LES (== direct)
    assert _se.split_extension_cohomology(kA2(), 4).assembled == [3, 1, 1, 1, 1]
    assert _se.split_extension_cohomology(kA3(), 4).assembled == [4, 1, 1, 1, 1]
    assert _se.split_extension_cohomology(kD4(), 4).assembled == [5, 1, 0, 0, 2]


@pytest.mark.oracle_crossengine
def test_gf32003_parity_kA2():
    # char > dim: no TrivialExtension/socle raise; the LES route agrees with direct.
    rep = _se.split_extension_cohomology(kA2(GF(32003)), 4)
    assert rep.assembled == rep.direct == [3, 1, 1, 1, 1]
    assert rep.agrees is True
