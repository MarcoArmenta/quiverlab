"""Plan 35 wave 3c -- Yoneda exact sequences realizing Ext classes.

Marco: "for ext the bases ARE exact sequences, and we need to SHOW these exact
sequences." Every Ext^n(M, N) class is constructed as an explicit n-fold exact sequence
0 -> N -> Q -> ... -> M -> 0 with Q a genuine module, and its exactness is SELF-CERTIFIED.

The crux pins:
  * kA_2: the unique nonsplit class of Ext^1(S_1, S_2) is the projective cover --
    0 -> S_2 -> P_1 -> S_1 -> 0 -- so the constructed middle module E is P_1, verified
    by the library's OWN is_isomorphic / identify_standard certificate;
  * for a battery of algebras / degrees the constructed sequence is exact at every joint
    (compositions zero, ranks add up, injective/surjective ends), certified by the
    engine's exact linear algebra;
  * a NON-cocycle is refused loudly -- a wrong sequence is never shown.
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.complex_reps import _reconstruct_cocycle
from quiverlab.modules.ext import ext_dims
from quiverlab.modules.hom import identify_standard, is_isomorphic
from quiverlab.modules.resolution import minimal_resolution
from quiverlab.modules.yoneda import (baer_extension, spliced_sequence,
                                      yoneda_sequence)


def _kA2():
    return Quiver(vertices=[1, 2], arrows={"a": (1, 2)}).algebra(relations=[], field=GF(7))


def _loop3(p=7):
    return Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x*x"], field=GF(p))


def _cocycles(A, M, N, top):
    """Reconstruct, per degree, the actual cocycle matrices of the Ext basis classes --
    the SAME reconstruction the capture layer uses -- returning (terms, dmats, {n: [f]})."""
    from quiverlab.modules.complex_reps import (_hom_adjunction_basis,
                                                _reps_from_complex, _terms_info)
    from quiverlab.modules.ext import _delta_matrix
    dom = A.domain
    terms, dmats = minimal_resolution(M, top + 1)
    tinfo = _terms_info(M.algebra, terms)
    vb = {}
    hb = [_hom_adjunction_basis(N, ti, dom, vb) for ti in tinfo]
    homs = [h for _e, h in hb]
    deltas = []
    for nn in range(len(homs) - 1):
        dn1 = dmats[nn + 1]
        deltas.append(_delta_matrix(homs[nn], homs[nn + 1], dn1, dom)
                      if (dn1 and dn1[0]) else
                      lm.zeros(len(homs[nn + 1]), len(homs[nn]), dom))
    out = {}
    for n in range(1, top + 1):
        space = len(homs[n]) if n < len(homs) else 0
        if not space:
            out[n] = []
            continue
        here = deltas[n] if n < len(deltas) else None
        prev = deltas[n - 1] if 0 <= n - 1 < len(deltas) else None
        cols = _reps_from_complex(here, prev, space, dom)
        width = len(homs[n][0][0]) if homs[n] else 0
        out[n] = [_reconstruct_cocycle(c, homs[n], N.dim, width, dom) for c in cols]
    return terms, dmats, out


# --------------------------------------------------------------------------- #
# (1) The kA_2 Baer-extension pin: E is the projective P_1, by the library's own iso.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_literature
def test_kA2_baer_extension_is_projective_cover():
    """0 -> S_2 -> P_1 -> S_1 -> 0: the unique nonsplit extension of S_1 by S_2 is the
    projective cover P_1 (a textbook fact for kA_2). The CONSTRUCTED middle module is
    isomorphic to P_1, certified by is_isomorphic AND identify_standard."""
    A = _kA2()
    M, N = A.simple(1), A.simple(2)
    terms, dmats, coc = _cocycles(A, M, N, 3)
    assert len(coc[1]) == 1                              # dim Ext^1(S_1, S_2) = 1
    seq = baer_extension(M, N, coc[1][0], terms, dmats)
    ok, facts = seq.check_exact()
    assert ok, facts
    assert seq.roles == ["sub", "middle", "quotient"]
    E = seq.middle
    assert E.dim == 2 and E.dimension_vector() == {1: 1, 2: 1}
    # the two independent certificates the plan asks for
    assert is_isomorphic(E, A.projective(1))
    assert identify_standard(E) == ("projective", 1)
    # the ends: 0 -> S_2 -> E -> S_1 -> 0
    assert is_isomorphic(seq.modules[0], A.simple(2))
    assert is_isomorphic(seq.modules[-1], A.simple(1))


# --------------------------------------------------------------------------- #
# (2) Exactness self-certification across a battery of algebras / degrees.
# --------------------------------------------------------------------------- #
def _battery():
    A = _loop3()                                         # k[x]/(x^3), infinite pd
    Mloop = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    B = _loop3(2)                                         # char 2
    Mloop2 = B.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    C = Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x"], field=GF(7))                  # dual numbers, infinite pd
    kA3 = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=GF(5))                       # hereditary: only Ext^1
    comm = Quiver(vertices=[1, 2, 3, 4],
                  arrows={"a": (1, 2), "b": (1, 3), "c": (2, 4), "d": (3, 4)}).algebra(
        relations=["a*c-b*d"], field=GF(7))              # commutative square, pd 2
    return [
        (A, Mloop, A.simple(1), 5),
        (B, Mloop2, B.simple(1), 5),
        (C, C.simple(1), C.simple(1), 5),
        (kA3, kA3.simple(1), kA3.simple(2), 3),
        (comm, comm.simple(1), comm.simple(4), 4),
    ]


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("idx", range(5))
def test_constructed_sequences_are_exact(idx):
    A, M, N, top = _battery()[idx]
    terms, dmats, coc = _cocycles(A, M, N, top)
    built = 0
    for n in range(1, top + 1):
        for f in coc.get(n, []):
            seq = yoneda_sequence(M, N, f, n, terms, dmats)
            ok, info = seq.check_exact()
            assert ok, "n=%d: %s" % (n, info)
            # the sequence starts at N and ends at M
            assert seq.modules[0].dim == N.dim and seq.modules[-1].dim == M.dim
            # n-fold: length is n + 2 modules (0 -> N -> Q -> P_{n-2} -> ... -> M -> 0)
            assert len(seq.modules) == n + 2
            built += 1
    # every algebra in the battery has at least one class to realize
    assert built >= 1


@pytest.mark.oracle_selfcert
def test_capture_ships_certified_sequences():
    """Through the public capture path: every degree's Yoneda sequence is certified and
    its facts are self-consistent (rank identities hold)."""
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    _dims, pl = ext_dims(A, M, A.simple(1), 4, with_reps=True, interpret=True)
    seqs = pl["interpretation"]["sequences"]
    assert set(seqs) == {"1", "2", "3", "4"}
    for dkey, classes in seqs.items():
        for s in classes:
            assert s["certified"] is True and "error" not in s
            for f in s["facts"]:
                if f["fact"] == "im=ker":
                    assert f["rank_in"] + f["rank_out"] == f["dim"]
                else:
                    assert f["rank"] == f["dim"]


# --------------------------------------------------------------------------- #
# (3) Yoneda degree bookkeeping + aliases.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_crossengine
def test_spliced_middle_terms_are_resolution_terms():
    """For n >= 2 the middle terms P_{n-2}..P_0 are exactly the minimal resolution's
    terms, so their dimension vectors match the resolution's."""
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    terms, dmats, coc = _cocycles(A, M, A.simple(1), 4)
    seq = spliced_sequence(M, A.simple(1), coc[3][0], 3, terms, dmats)
    # roles: sub, middle, then resolution_term P_1, P_0, then quotient M
    rt = [i for i, r in enumerate(seq.roles) if r == "resolution_term"]
    assert len(rt) == 2                                  # P_1 and P_0
    # the last resolution term (P_0) has the dimension vector of terms[0]
    assert seq.modules[rt[-1]].dimension_vector() == terms[0].module.dimension_vector()


@pytest.mark.oracle_selfcert
def test_aliases_and_degree_guards():
    A = _loop3()
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")
    terms, dmats, coc = _cocycles(A, M, A.simple(1), 3)
    # baer_extension is the n=1 alias
    assert baer_extension(M, A.simple(1), coc[1][0], terms, dmats).degree == 1
    # spliced_sequence refuses n < 2
    with pytest.raises(QuiverlabError):
        spliced_sequence(M, A.simple(1), coc[1][0], 1, terms, dmats)
    # yoneda_sequence refuses n < 1
    with pytest.raises(QuiverlabError):
        yoneda_sequence(M, A.simple(1), coc[1][0], 0, terms, dmats)


# --------------------------------------------------------------------------- #
# (4) Honesty: a non-cocycle is refused loudly, never shown as a sequence.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_noncocycle_refused_loudly():
    A = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=GF(5))
    M, N = A.simple(1), A.simple(3)
    terms, dmats = minimal_resolution(M, 2)
    d1 = dmats[1]
    dom = A.domain
    bad = [[dom.one()] * len(d1[0]) for _ in range(N.dim)]   # not a cocycle
    with pytest.raises(QuiverlabError):
        yoneda_sequence(M, N, bad, 1, terms, dmats)


@pytest.mark.oracle_selfcert
def test_capture_guards_each_class():
    """A class whose sequence fails self-cert ships an honest error entry, never a
    certified-but-wrong sequence. (Simulated by feeding a corrupted cocycle to the
    per-class serializer.)"""
    from quiverlab.modules.complex_reps import _one_ext_sequence
    A = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=GF(5))
    M, N = A.simple(1), A.simple(3)
    terms, dmats = minimal_resolution(M, 2)
    from quiverlab.modules.complex_reps import (_hom_adjunction_basis, _terms_info)
    dom = A.domain
    tinfo = _terms_info(M.algebra, terms)
    _e, homs1 = _hom_adjunction_basis(N, tinfo[1], dom, {})
    width = len(homs1[0][0]) if homs1 else 1
    bad_col = [dom.one()] * len(homs1)                   # not a cocycle combination
    entry = _one_ext_sequence(A, M, N, 1, bad_col, homs1, width, terms, dmats,
                              "\\alpha^{1}_{1}")
    assert entry["certified"] is False and "error" in entry
