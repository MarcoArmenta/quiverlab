"""CS §3 verbatim block combinatorics: the witness tip may STRADDLE a block boundary
(pair reducible with minimal extension), which the exact-pair condition missed.
Repro discovered 2026-07-22: tips {xx,yy,xyx} lose the chain xyxx, making Bardzell
HH wrong from degree 2 (vs bar and minimal-A^e oracles). Uniform-relation families
(truncated poly, Nakayama, radsq) have exact pairs only: pinned unchanged below."""
import pytest

from quiverlab.engine.resolutions_bardzell import MonomialPresentation


def _s1():
    # one vertex, loops x,y; tips {xx, yy, xyx} — the minimal straddle presentation
    return MonomialPresentation(["v"], [("x", "v", "v"), ("y", "v", "v")],
                                [("x", "x"), ("y", "y"), ("x", "y", "x")])


def _line():
    # 1->2->..->6, arrows a..e; tips {abc, cde} overlap in the single arrow c
    return MonomialPresentation(
        [1, 2, 3, 4, 5, 6],
        [("a", 1, 2), ("b", 2, 3), ("c", 3, 4), ("d", 4, 5), ("e", 5, 6)],
        [("a", "b", "c"), ("c", "d", "e")])


def test_straddle_ap3_complete():
    """AP^3 = 5 chains, matching the minimal-A^e Betti number rks[3] = 5 (the old
    exact-pair code returned 3, missing the straddles xyxx and xyxyx)."""
    ap3 = set(_s1().associated_paths(3, 30))
    assert ap3 == {("x", "x", "x"), ("x", "x", "y", "x"), ("x", "y", "x", "x"),
                   ("x", "y", "x", "y", "x"), ("y", "y", "y")}


def test_line_quiver_short_overlap_chain():
    pres = _line()
    assert pres.associated_paths(3, 20) == [("a", "b", "c", "d", "e")]
    assert pres.associated_paths(4, 24) == []          # no further overlap


def test_left_decomposition_straddle_blocks():
    pres = _s1()
    assert pres.left_decomposition(("x", "y", "x", "x"), 3) == \
        [("x",), ("y", "x"), ("x",)]                   # pair yx·x ⊇ xx straddling
    assert pres.left_decomposition(("x", "x", "y", "x"), 3) == \
        [("x",), ("x",), ("y", "x")]                   # pairs xx, xyx exact
    assert _line().left_decomposition(("a", "b", "c", "d", "e"), 3) == \
        [("a",), ("b", "c"), ("d", "e")]               # pair bc·de ⊇ cde straddling


def test_uniform_families_unchanged():
    tp = MonomialPresentation.truncated_polynomial(3)
    assert tp.associated_paths(3, 20) == [(0, 0, 0, 0)]
    assert tp.left_decomposition((0, 0, 0, 0), 3) == [(0,), (0, 0), (0,)]
    cn = MonomialPresentation.cyclic_nakayama(3, 2)
    assert sorted(cn.associated_paths(3, 20)) == \
        [(0, 1, 2), (1, 2, 0), (2, 0, 1)]              # ell=2: quadratic, single-arrow blocks
    cn3 = MonomialPresentation.cyclic_nakayama(2, 3)   # ell=3: overlapping relations
    assert sorted(cn3.associated_paths(3, 30)) == [(0, 1, 0, 1), (1, 0, 1, 0)]
    assert cn3.left_decomposition((0, 1, 0, 1), 3) == [(0,), (1, 0), (1,)]
    rs = MonomialPresentation.local_radsq(2)
    assert len(rs.associated_paths(3, 20)) == 8        # all words of length 3 (quadratic)


def test_deep_straddle_chain_growth():
    # AP^n matches the minimal-A^e Betti numbers 2,3,5,9,17 (session-verified via
    # minimal_resolution over GF(32003)) and stays finite to depth 8
    pres = _s1()
    counts = [len(pres.associated_paths(n, 60)) for n in range(1, 9)]
    assert counts[:5] == [2, 3, 5, 9, 17]
    assert all(c >= 2 for c in counts)                 # y^k chain and an x-side chain persist


def test_right_decomposition_mirror_and_straddle():
    pres = _s1()
    # xyxx: right blocks (xy)(x)(x) — witness xyx STARTS at the pair's start, straddling
    assert pres.right_decomposition(("x", "y", "x", "x"), 3) == \
        [("x", "y"), ("x",), ("x",)]
    assert pres.right_decomposition(("x", "x", "y", "x"), 3) == \
        [("x",), ("x", "y"), ("x",)]
    assert _line().right_decomposition(("a", "b", "c", "d", "e"), 3) == \
        [("a", "b"), ("c", "d"), ("e",)]


def test_right_decomposition_quadratic_reverses_left():
    # CS Prop. "cuadratico": for S ⊆ Q_2, v_i = α_{n-i} — blocks are single arrows,
    # so right blocks == left blocks elementwise.
    rs = MonomialPresentation.local_radsq(2)
    for n in (2, 3, 4):
        for p in rs.associated_paths(n, 20):
            assert rs.right_decomposition(p, n) == rs.left_decomposition(p, n)


def test_right_decomposition_palindromic_families_agree():
    tp = MonomialPresentation.truncated_polynomial(3)
    assert tp.right_decomposition((0, 0, 0, 0), 3) == [(0,), (0, 0), (0,)]
    cn = MonomialPresentation.cyclic_nakayama(3, 3)
    for n in (2, 3):
        for p in cn.associated_paths(n, 30):
            L, R = cn.left_decomposition(p, n), cn.right_decomposition(p, n)
            assert [len(b) for b in R] == [len(b) for b in L][::-1] or L == R


def test_qci32_right_blocks():
    # tips {yx, xxx, yy} (QCI(3,2) tip algebra): yxxx right = (y)(xx)(x), left = (y)(x)(xx)
    pres = MonomialPresentation(["v"], [("x", "v", "v"), ("y", "v", "v")],
                                [("y", "x"), ("x", "x", "x"), ("y", "y")])
    assert pres.left_decomposition(("y", "x", "x", "x"), 3) == [("y",), ("x",), ("x", "x")]
    assert pres.right_decomposition(("y", "x", "x", "x"), 3) == [("y",), ("x", "x"), ("x",)]


def test_left_right_ambiguity_sets_coincide():
    """CS §3 Proposition (citing Bardzell, Sköldberg): right n-ambiguities = left
    n-ambiguities as SETS. right_decomposition must succeed on every associated path."""
    for pres in (_s1(), _line(),
                 MonomialPresentation(["v"], [("x", "v", "v"), ("y", "v", "v")],
                                      [("y", "x"), ("x", "x", "x"), ("y", "y")])):
        for n in (2, 3, 4):
            for p in pres.associated_paths(n, 40):
                blocks = pres.right_decomposition(p, n)
                assert len(blocks) == n
                assert tuple(x for b in blocks for x in b) == tuple(p)


def test_straddle_hh_matches_bar_and_minimal():
    """The decisive repro: A = k<x,y>/(xx,yy,xyx), dim 6. Bar and minimal-A^e oracles
    agree; Bardzell must now too (was [4,4,7,6] vs [4,4,6,9] at 32003 before the fix)."""
    import quiverlab as ql
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    from quiverlab.engine.resolutions_bardzell import BardzellResolution
    from quiverlab.engine.resolutions_minimal import minimal_homology_dims

    Q = ql.Quiver(["v"], {"x": ("v", "v"), "y": ("v", "v")})
    A = Q.algebra(relations=["x*x", "y*y", "x*y*x"], field=ql.GF(32003))
    eng = to_engine(A)
    res = BardzellResolution(_s1())
    N, PRIMES = 3, (32003, 2, 3, 5)
    bar = hochschild_homology_dims(eng, N, primes=PRIMES)
    bard = hochschild_homology_dims(eng, N, primes=PRIMES, resolution=res)
    minh = minimal_homology_dims(eng, N, primes=PRIMES)
    for p in PRIMES:
        assert bard[p] == bar[p], f"p={p}: {bard[p]} != {bar[p]}"
        assert minh[p] == bar[p][:len(minh[p])], f"minimal p={p}"
    assert bar[32003] == [4, 4, 6, 9]                  # session-verified literals
    assert bar[2] == [4, 7, 9, 12]


def test_straddle_bardzell_depth_unlock():
    """Bardzell runs past the bar window on the straddle presentation (smoke, N=10)."""
    import quiverlab as ql
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    from quiverlab.engine.resolutions_bardzell import BardzellResolution

    Q = ql.Quiver(["v"], {"x": ("v", "v"), "y": ("v", "v")})
    A = Q.algebra(relations=["x*x", "y*y", "x*y*x"], field=ql.GF(32003))
    dims = hochschild_homology_dims(to_engine(A), 10, primes=(32003,),
                                    resolution=BardzellResolution(_s1()))[32003]
    assert len(dims) == 11 and all(d >= 0 for d in dims)
