"""Radical / top / socle and the radical series (spec §3.6). Fixtures A & B."""
from quiverlab import Quiver, GF, linear_path_algebra
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.radtopsoc import quotient, _rad_image_cols


def _square(field=None):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


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


def test_quotient_direct_on_genuine_submodule_solves_and_is_a_module():
    # Fix 1a: quotient() now asserts its coset solve is not None before indexing [0] /
    # slicing. On a genuine submodule (rad P_1) the solve must succeed, the assert must
    # NOT fire, and the result must be a genuine A-module (top P_1 = S_1 here).
    A = _square()
    P1 = A.projective(1)
    Q = quotient(P1, _rad_image_cols(P1), name="q")
    ok, why = Q.check_module()
    assert ok, why
    assert Q.dimension_vector() == {1: 1, 2: 0, 3: 0, 4: 0}   # top P_1 = S_1


def test_radtopsoc_outputs_are_genuine_modules():
    # T4 review (binding constraint): radical/top/socle must each RETURN a genuine
    # A-module -- check_module must PASS on every output, not merely on the input P_1.
    A = _square()
    P1 = A.projective(1)
    for M in (P1.radical(), P1.top(), P1.socle()):
        ok, why = M.check_module()
        assert ok, why


def test_square_p1_radtopsoc_over_gfp():
    # T4 review (binding constraint): domain-genericity over GF(p). Recompute the square
    # fixture's radical series over GF(5) (the reviewer's verified case): identical dimension
    # vectors to the default-field test, and every output still a genuine module.
    A = _square(field=GF(5))
    P1 = A.projective(1)               # dim 4, dimvec {1:1,2:1,3:1,4:1}
    assert P1.top().dimension_vector() == {1: 1, 2: 0, 3: 0, 4: 0}
    r1 = P1.radical()                  # rad P_1, dim 3
    assert r1.dim == 3
    assert r1.top().dimension_vector() == {1: 0, 2: 1, 3: 1, 4: 0}   # S_2 (+) S_3
    r2 = r1.radical()                 # rad^2 P_1 = span{a*b}, dim 1 ~ S_4
    assert r2.dim == 1 and r2.dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}
    assert r2.radical().dim == 0     # rad^3 P_1 = 0  -> Loewy length 3
    assert P1.socle().dimension_vector() == {1: 0, 2: 0, 3: 0, 4: 1}   # soc P_1 = S_4
    for M in (P1.radical(), P1.top(), P1.socle()):
        ok, why = M.check_module()
        assert ok, why
