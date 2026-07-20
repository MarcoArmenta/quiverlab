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
