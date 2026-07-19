"""Radical / top / socle and the radical series (spec §3.6). Fixtures A & B."""
from quiverlab import Quiver, CC, linear_path_algebra


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_p1_of_a2_radical_series():
    A = linear_path_algebra(2)
    P1 = A.projective(1)
    top = P1.top()
    assert top.dim == 1 and top.dimension_vector() == {1: 1, 2: 0}   # top P_1 = S_1
    rad = P1.radical()
    assert rad.dim == 1 and rad.dimension_vector() == {1: 0, 2: 1}   # rad P_1 = S_2
    soc = P1.socle()
    assert soc.dim == 1 and soc.dimension_vector() == {1: 0, 2: 1}   # soc P_1 = S_2
    # rad^2 P_1 = rad(rad P_1) = 0
    assert rad.radical().dim == 0


def test_square_p1_radical_layers():
    A = _square()
    P1 = A.projective(1)               # dim 4, dimvec {1:1,2:1,3:1,4:1}
    assert P1.top().dimension_vector() == {1: 1, 2: 0, 3: 0, 4: 0}
    r1 = P1.radical()                  # rad P_1, dim 3
    assert r1.dim == 3
    assert r1.top().dimension_vector() == {1: 0, 2: 1, 3: 1, 4: 0}   # S_2 (+) S_3
    r2 = r1.radical()                 # rad^2 P_1 = span{a*b}, dim 1 ~ S_4
    assert r2.dim == 1 and r2.dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}
    assert r2.radical().dim == 0     # rad^3 P_1 = 0  -> Loewy length 3
    assert P1.socle().dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}   # soc P_1 = S_4


def test_simple_is_its_own_top_and_socle_zero_radical():
    A = linear_path_algebra(2)
    S1 = A.simple(1)
    assert S1.radical().dim == 0
    assert S1.top().dim == 1
    assert S1.socle().dim == 1
