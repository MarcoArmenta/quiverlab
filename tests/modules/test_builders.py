"""Simples/projectives/injectives from quiver provenance (spec §3.6). Fixtures A & B."""
import pytest
from quiverlab import Quiver, CC, GF, linear_path_algebra, Algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def test_a2_simples():
    A = linear_path_algebra(2)
    assert A.simple(1).dimension_vector() == {1: 1, 2: 0}
    assert A.simple(2).dimension_vector() == {1: 0, 2: 1}
    assert A.simple(1).dim == 1
    for v in (1, 2):
        ok, why = A.simple(v).check_module()   # each S_v is a genuine module
        assert ok, why


def test_a2_projectives():
    A = linear_path_algebra(2)
    P1, P2 = A.projective(1), A.projective(2)
    assert P1.dim == 2 and P1.dimension_vector() == {1: 1, 2: 1}
    assert P2.dim == 1 and P2.dimension_vector() == {1: 0, 2: 1}
    ok, why = P1.check_module()
    assert ok, why


def test_a2_injectives():
    A = linear_path_algebra(2)
    I1, I2 = A.injective(1), A.injective(2)
    assert I1.dim == 1 and I1.dimension_vector() == {1: 1, 2: 0}   # I_1 = S_1
    assert I2.dim == 2 and I2.dimension_vector() == {1: 1, 2: 1}
    ok, why = I2.check_module()
    assert ok, why
    # check_module alone cannot pin the genuine I_2: a ZERO arrow action has the
    # same dimvec {1:1, 2:1} AND passes check_module -- but that module is the
    # decomposable S_1 (+) S_2, not the indecomposable injective. Pin the arrow
    # action as NONZERO to kill that impostor.
    dom = I2.domain
    assert not all(dom.is_zero(x) for row in I2.action["a1"] for x in row)


def test_square_projectives_match_cartan_rows():
    A = _square()
    C = A.cartan_matrix()
    verts = [1, 2, 3, 4]
    for i, v in enumerate(verts):
        Pv = A.projective(v)
        assert Pv.dimension_vector() == {verts[j]: C[i][j] for j in range(4)}
    assert A.projective(1).dim == 4
    ok, _ = A.projective(1).check_module()
    assert ok


def test_injective_dimvec_is_cartan_column():
    A = _square()
    C = A.cartan_matrix()
    verts = [1, 2, 3, 4]
    for j, v in enumerate(verts):
        Iv = A.injective(v)
        assert Iv.dimension_vector() == {verts[i]: C[i][j] for i in range(4)}
        ok, why = Iv.check_module()   # relational (Groebner-route) injective is a genuine module
        assert ok, why


def test_builders_over_gfp():
    A = linear_path_algebra(2, field=GF(7))
    P1 = A.projective(1)
    assert P1.dimension_vector() == {1: 1, 2: 1}
    ok, why = P1.check_module()   # builder is a genuine module over GF(7) too
    assert ok, why


def test_builders_need_provenance():
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(Exception):
        A.simple(1)
