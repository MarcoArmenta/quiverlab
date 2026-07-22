import pytest
from quiverlab.errors import DepthLimitError
from quiverlab.resolutions_cs.ambiguities import SSequence
pytest.importorskip("quiverlab.groebner")


def test_kx2_single_chain_per_degree(kx2_rs):
    ss = SSequence(kx2_rs, max_degree=6)
    assert [len(ss.S(n)) for n in range(7)] == [1, 1, 1, 1, 1, 1, 1]
    assert ss.S(2)[0].word == ("x", "x")
    assert ss.S(3)[0].word == ("x", "x", "x")
    assert ss.S(3)[0].blocks == (("x",), ("x",), ("x",))          # x|x|x, three blocks


def test_commutative_square_length_two(square_rs):
    ss = SSequence(square_rs, max_degree=5)
    assert len(ss.S(0)) == 4 and len(ss.S(1)) == 4 and len(ss.S(2)) == 1
    assert ss.S(2)[0].word == ("c", "d") and ss.S(2)[0].o == 1 and ss.S(2)[0].t == 4
    assert ss.S(3) == [] and ss.S(4) == []


def test_quantum_ci_chain_counts(qci_rs):
    ss = SSequence(qci_rs(xi="2"), max_degree=6)
    assert [len(ss.S(n)) for n in range(7)] == [1, 2, 3, 4, 5, 6, 7]  # |S_n| = n+1
    assert {c.word for c in ss.S(2)} == {("x", "x"), ("y", "y"), ("y", "x")}
    assert {c.word for c in ss.S(3)} == {("x",)*3, ("y", "x", "x"), ("y", "y", "x"), ("y",)*3}


def test_depth_guard_loud(qci_rs):
    ss = SSequence(qci_rs(xi="2"), max_degree=60, max_cells=5)
    with pytest.raises(DepthLimitError) as e:
        ss.S(40)
    assert "certified" in str(e.value).lower()


def test_out_of_range_degree_raises_loud(kx2_rs):
    """GO-LOUD (open-item #1): a degree strictly above max_degree must raise
    DepthLimitError naming the certified range, not silently return []."""
    ss = SSequence(kx2_rs, max_degree=3)
    for n in (4, 5):
        with pytest.raises(DepthLimitError) as e:
            ss.S(n)
        msg = str(e.value)
        assert "certified" in msg.lower()          # names the certified range
        assert "0..3" in msg                        # the range 0..max_degree
        assert str(n) in msg                        # the offending degree


def test_in_range_empties_and_negative_unchanged(square_rs):
    """In-range genuinely-empty cochain spaces (S(3), S(4) at max_degree=5 for the
    commutative square) and the n<0 convention are untouched by the loud guard."""
    ss = SSequence(square_rs, max_degree=5)
    assert ss.S(3) == [] and ss.S(4) == []          # in-range, genuinely empty: still []
    assert ss.S(-1) == []                            # negative degree convention: still []


def _rs(rels, field=None):
    from quiverlab import Quiver, CC
    from quiverlab.groebner import build_reduction_system
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    return build_reduction_system(Q, rels, field or CC)


def test_ssequence_straddle_chains_present():
    """Plan 12: tips {xx,yy,xyx} — S_3 must contain the straddling chains xyxx and
    xyxyx (5 chains total, = the minimal-A^e Betti number; the exact-pair code had 3)."""
    ss = SSequence(_rs(["x*x", "y*y", "x*y*x"]), max_degree=4)
    assert {c.word for c in ss.S(3)} == {("x",) * 3, ("x", "x", "y", "x"),
                                         ("x", "y", "x", "x"), ("x", "y", "x", "y", "x"),
                                         ("y",) * 3}
    straddle = next(c for c in ss.S(3) if c.word == ("x", "y", "x", "x"))
    assert straddle.blocks == (("x",), ("y", "x"), ("x",))
    assert len(ss.S(4)) == 9                              # Betti sequence 2,3,5,9,17


def test_ssequence_qci32_matches_cs_phi_formula():
    """CS §7.2 (TeX ~2110–2140): for k<x,y>/(x^n, y^m, yx−ξxy), 𝒜_N =
    {y^{φ(s,m)} x^{φ(t,n)} : s+t = N+1}, φ(s,k) = (s/2)k if s even else ((s−1)/2)k + 1.
    quiverlab S_N = 𝒜_{N-1}. Pin n=3, m=2 through S_5."""
    def phi(s, k):
        return (s // 2) * k if s % 2 == 0 else ((s - 1) // 2) * k + 1
    ss = SSequence(_rs(["x*x*x", "y*y", "y*x - 2*x*y"]), max_degree=5)
    for N in range(2, 6):
        expect = {("y",) * phi(s, 2) + ("x",) * phi(t, 3)
                  for s in range(N + 1) for t in range(N + 1) if s + t == N}
        assert {c.word for c in ss.S(N)} == expect, f"S_{N}"
