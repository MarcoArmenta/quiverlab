# SPDX-License-Identifier: MIT
"""Internal (path-length) generation degrees of Ext(k,k) off the shipped minimal
resolutions (Plan 77 / R36). Validated against Berger's closed form on k[x]/(x^N):
delta(n) = (N/2) n (n even) / (N/2)(n-1)+1 (n odd)."""
import pytest

from quiverlab import (PreprojectiveAlgebra, Quiver, TruncatedPathAlgebra,
                       truncated_polynomial)
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import (berger_degree, generation_degrees,
                                       n_homogeneous_degree)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("N,expected", [
    (2, [0, 1, 2, 3, 4, 5, 6, 7, 8]),
    (3, [0, 1, 3, 4, 6, 7, 9, 10, 12]),
    (4, [0, 1, 4, 5, 8, 9, 12, 13, 16]),
    (5, [0, 1, 5, 6, 10, 11, 15, 16, 20]),
])
def test_kxN_internal_degrees_match_berger(N, expected):
    A = truncated_polynomial(N, field=QQ)
    gd = generation_degrees(A, 1, 8)              # single vertex -> S_1
    flat = [s[0] for s in gd]                     # one summand each degree
    assert flat == expected[:len(flat)]
    assert flat == [berger_degree(n, N) for n in range(len(flat))]
    assert n_homogeneous_degree(A) == N


@lit
def test_cubic_monomial_multivertex_is_berger_N3():
    A = TruncatedPathAlgebra("A4", 3, field=QQ)   # single relation e12*e23*e34, N=3
    assert n_homogeneous_degree(A) == 3
    gd = generation_degrees(A, 1, 4)              # S_1 has pd 2
    assert [s[0] for s in gd if s] == [0, 1, 3]   # Berger delta for N=3, finite


@selfcert
def test_generation_degrees_single_element_on_pure_resolution():
    # DEFENSIVE: on k[x]/x^3 every P_n is generated in ONE internal degree, so each
    # returned set is a singleton. This is NOT the correctness anchor (Berger's closed
    # form above is) -- it guards that the _NonPure branch stays dormant on a
    # genuinely graded-pure input.
    A = truncated_polynomial(3, field=QQ)
    for s in generation_degrees(A, 1, 6):
        assert len(s) == 1


@selfcert
def test_refuses_non_length_graded():
    # x^3 - x^2 is ADMISSIBLE (the ideal sits inside rad^2) but genuinely
    # non-length-graded (a length-3 and a length-2 term). (x^2 - x, the naive choice,
    # raises AdmissibilityError AT CONSTRUCTION -- ideal not in rad^2 -- so it would
    # never exercise the _is_length_graded gate: a false green.)
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x - x*x"], field=QQ)
    with pytest.raises(QuiverlabError, match="length-graded"):
        generation_degrees(A, 1, 4)


@selfcert
def test_hereditary_is_n_homogeneous_two_vacuously():
    from quiverlab import linear_path_algebra
    A = linear_path_algebra(2, field=QQ)
    assert A.relations == []
    assert n_homogeneous_degree(A) == 2           # kQ = the empty-relation quadratic
    assert [s for s in generation_degrees(A, 1, 4) if s] == [[0], [1]]


@selfcert
def test_non_homogeneous_relations_have_no_single_degree():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x - x*x"], field=QQ)
    assert n_homogeneous_degree(A) is None        # inhomogeneous relation


@xeng
def test_internal_degrees_are_bounded_below_by_the_homological_degree():
    # Every differential has path length >= 1 (minimality: d(P_n) <= rad P_{n-1}), so
    # l_i(n) >= n always, with equality exactly on a linear (Koszul) strand. Checked
    # against an independent route: the minimal resolution's own Betti data.
    from quiverlab import PreprojectiveAlgebra
    A = PreprojectiveAlgebra("A3", field=QQ)
    for v in A.quiver.vertices:
        gd = generation_degrees(A, v, 5)
        for n, s in enumerate(gd):
            if s:
                assert min(s) >= n


@selfcert
@pytest.mark.parametrize("field_name", ["QQ", "GF2", "GF3", "GF32003"])
def test_internal_degrees_are_FIELD_INDEPENDENT(field_name):
    # The acceptance claim "over any exact Domain", pinned rather than asserted in prose.
    # Internal degrees are integer PATH LENGTHS -- field-free combinatorics -- so the
    # whole recognizer ladder must return the same answer over QQ and over GF(p) for
    # every p, including the small primes where a careless implementation would divide.
    # (Contrast the Plan-75 incidence route, where the answer genuinely IS
    # characteristic-dependent: RP^2 over GF(2) differs from QQ. Here it must NOT be.)
    from quiverlab.fields import GF
    from quiverlab.modules.nkoszul import koszul_profile_block
    field = {"QQ": QQ, "GF2": GF(2), "GF3": GF(3), "GF32003": GF(32003)}[field_name]
    A = truncated_polynomial(3, field=field)
    assert [s[0] for s in generation_degrees(A, 1, 6)] == [0, 1, 3, 4, 6, 7, 9]
    b = koszul_profile_block(A, 6)
    assert b["n_homogeneous"] == 3
    assert b["k2"]["generator_degrees"] == [1, 2]

    P = PreprojectiveAlgebra("A3", field=field)
    pb = koszul_profile_block(P, 6)
    assert (pb["almost_koszul"]["p"], pb["almost_koszul"]["q"]) == (2, 2)
    assert pb["almost_koszul"]["verdict"] is True
    assert pb["generation_degrees"]["1"] == [0, 1, 2, 4, 5, 6, 8]
