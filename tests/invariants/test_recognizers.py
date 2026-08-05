"""Recognizer batch pinned on textbook examples (ASS; Butler-Ringel;
Assem-Skowronski gentle papers)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.invariants import recognizers as rec

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra()


def _gentle_a3():
    # 1 --a--> 2 --b--> 3 with a*b = 0: the standard gentle example
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"])


def _sb_not_string():
    # special biserial with a non-monomial relation: commutative square
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_semisimple_and_radsq():
    S = Quiver([1, 2], {}).algebra()
    assert rec.is_semisimple(S) is True
    A = _kA2()
    assert rec.is_semisimple(A) is False
    assert rec.is_radical_square_zero(A) is True
    B = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x"], field=GF(5))
    assert rec.is_radical_square_zero(B) is False


def test_hereditary():
    assert rec.is_hereditary(_kA2()) is True
    assert rec.is_hereditary(_gentle_a3()) is False        # relations present
    loop = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(5))
    assert rec.is_hereditary(loop) is False                # cycle


def test_basic():
    assert rec.is_basic(_kA2()) is True
    assert rec.is_basic(_sb_not_string()) is True


def test_nakayama():
    assert rec.is_nakayama(_kA2()) is True
    cyc = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b*a", "b*a*b"], field=GF(5))
    assert rec.is_nakayama(cyc) is True
    assert rec.is_nakayama(_sb_not_string()) is False      # vertex 1 out-deg 2


def test_gentle_string_special_biserial_hierarchy():
    G = _gentle_a3()
    assert rec.is_gentle(G) and rec.is_string(G) and rec.is_special_biserial(G)
    SB = _sb_not_string()
    assert rec.is_special_biserial(SB) is True
    assert rec.is_string(SB) is False                      # binomial relation
    assert rec.is_gentle(SB) is False
    # NOT special biserial: three arrows out of one vertex
    W = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (1, 4)}).algebra()
    assert rec.is_special_biserial(W) is False


def test_hereditary_path_algebra_is_gentle_iff_a_n():
    # kA_n with no relations is gentle; the 3-star is not (condition 1)
    assert rec.is_gentle(_kA2()) is True


def test_algebra_method_surface():
    A = _gentle_a3()
    assert A.is_gentle() is True
    assert A.is_special_biserial() is True
    assert A.is_hereditary() is False
    assert _kA2().is_hereditary() is True
