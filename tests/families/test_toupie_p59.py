"""Toupie algebras (Plan 59 / R35, Artenstein-Lanzilotta-Solotar 2020). Self-cert:
dimension certificate; recognizer accepts every constructed toupie and rejects
non-toupie shapes. Literature: the a-Kronecker HH = [1, a^2-1, 0, ..]. sl_a: char 0."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.families.toupie import (ToupieAlgebra, is_toupie,
                                       toupie_branch_count,
                                       toupie_direct_arrow_count)

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
@pytest.mark.parametrize("branches, dim", [
    ([1, 1], 4),            # 2-Kronecker: 2 + 2
    ([1, 1, 1], 5),         # 3-Kronecker: 2 + 3
    ([2], 6),               # single length-2 branch = kA3
    ([2, 2], 10),           # two length-2 branches, no rel: 4 verts + 6 paths
    ([1, 2], 7),            # mixed lengths: 3 verts + 4 paths
])
def test_relation_free_dimension_certificate(branches, dim):
    A = ToupieAlgebra(branches, field=QQ)
    assert A.dim == dim
    assert is_toupie(A) is True
    assert toupie_branch_count(A) == len(branches)


@selfcert
def test_recognizer_rejects_non_toupie():
    kD4 = Quiver([0, 1, 2, 3], {"a": (1, 0), "b": (2, 0), "c": (3, 0)}).algebra(
        relations=[], field=QQ)                            # centre has in-deg 3 -> not a toupie
    assert is_toupie(kD4) is False
    line = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)
    assert is_toupie(line) is True                         # single branch = a toupie (a=1)


@selfcert
def test_recognizer_rejects_disconnected_or_cyclic_trap():
    # VERIFIED LIVE: (path U oriented-cycle) passes the DEGREE checks -- unique source
    # [0], unique sink [1], all other vertices in/out-deg 1 -- yet is disconnected AND has
    # a cycle, so kQ is infinite-dimensional. is_toupie MUST reject it (connected+acyclic).
    trap = Quiver([0, 1, "c0", "c1"],
                  {"p": (0, 1), "e": ("c0", "c1"), "f": ("c1", "c0")}).algebra(
        relations=["e*f", "f*e"], field=QQ)                # kill the cycle to stay f.d.
    assert is_toupie(trap) is False                        # disconnected + (pre-truncation) cyclic


@selfcert
def test_non_monomial_cross_branch_relation():
    # two length-2 branches p0 = x1*x2, p1 = y1*y2; the commutative toupie p0 - p1 = 0.
    A = ToupieAlgebra([2, 2], relations=["p0 - p1"], field=QQ)   # branch-path token names
    assert is_toupie(A) is True
    assert A.dim == 9                                      # 10 (free) - 1 (identifies the two 0->w paths)
    assert toupie_branch_count(A) == 2
    assert toupie_direct_arrow_count(A) == 0               # no direct source->sink arrow


@selfcert
def test_direct_vs_branch_count_on_mixed():
    # 2 direct arrows source->sink + one length-2 branch: 3 branches, 2 direct arrows.
    A = Quiver([0, "v", 1],
               {"d0": (0, 1), "d1": (0, 1), "b1": (0, "v"), "b2": ("v", 1)}
               ).algebra(relations=[], field=QQ)
    assert is_toupie(A) is True
    assert toupie_branch_count(A) == 3
    assert toupie_direct_arrow_count(A) == 2


@selfcert
def test_a_kronecker_shortcut():
    A = ToupieAlgebra([1] * 4, field=QQ)                   # a-Kronecker = ToupieAlgebra([1]*a)
    assert A.dim == 2 + 4
    assert toupie_branch_count(A) == 4
    assert toupie_direct_arrow_count(A) == 4               # all branches length 1 -> all direct


@selfcert
def test_expected_dim_certificate_and_mismatch():
    A = ToupieAlgebra([2, 2], relations=["p0 - p1"], field=QQ, expected_dim=9)
    assert A.dim == 9
    with pytest.raises(QuiverlabError):
        ToupieAlgebra([2, 2], relations=["p0 - p1"], field=QQ, expected_dim=10)


@selfcert
def test_loud_refusals():
    with pytest.raises(QuiverlabError):
        ToupieAlgebra([], field=QQ)                        # no branches
    with pytest.raises(QuiverlabError):
        ToupieAlgebra([1, 0], field=QQ)                    # branch length 0
    with pytest.raises(QuiverlabError):
        ToupieAlgebra([1, 1], relations=["p0 - p1"], field=QQ)   # length-1 cross relation (rad, not rad^2)
    # a structure-constant (presentation-less) algebra: is_toupie refuses loudly
    from quiverlab.core.algebra import Algebra
    one = QQ.one()
    sc = Algebra(QQ, [[[one]]], [one])                     # 1-dim k, quiver is None (escape hatch)
    assert sc.quiver is None
    with pytest.raises(QuiverlabError):
        is_toupie(sc)
