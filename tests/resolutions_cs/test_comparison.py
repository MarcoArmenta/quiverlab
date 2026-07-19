"""Task 12: CS<->bar comparison morphisms and windowed transport of cup/cap/bracket.

Two disjoint engines are bridged: the Chouhy-Solotar cochain complex (CS side,
domain-generic) and the normalized bar complex with its Tamarkin-Tsygan calculus
(bar side, engine.tt_calculus over GF(p)). The comparison cochain map Phi# (bar ->
CS, honest -- carries the relation/beta correction, not only the leading block
term) is asserted to be a genuine chain map, and the induced iso on HH^* is used to
transport cocycles and their cup/bracket classes between the two worlds.

The four tests below are transcribed verbatim from p04-task-12-brief.md.
"""
import pytest
from quiverlab import Quiver, GF
from quiverlab.resolutions_cs.comparison import Comparison
pytest.importorskip("quiverlab.groebner")


def _kx2_gf():
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(32003))


def _square_gf():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b - c*d"], field=GF(32003))


def test_phi_is_chain_map():
    Comparison(_square_gf()).assert_chain_map(upto=2)                 # b_bar ∘ Φ = Φ ∘ d_cs


def test_transport_roundtrip_identity_on_cohomology():
    Comparison(_kx2_gf()).assert_transport_roundtrip_identity(upto=3)  # Ψ*Φ* = id on HH^n


def test_transported_cup_consistent():
    comp = Comparison(_kx2_gf())
    u = comp.hh_class_cs(1, 0)
    assert comp.same_cohomology_class(comp.cup_of_cs_classes(u, u),
                                      comp.transport_then_bar_cup(u, u), degree=2)


def test_operation_window_boundary():
    comp = Comparison(_kx2_gf())
    with pytest.raises(NotImplementedError):
        comp.cup_of_cs_classes(comp.hh_class_cs(20, 0), comp.hh_class_cs(20, 0))
