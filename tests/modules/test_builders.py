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
    ok, _ = I2.check_module()
    assert ok


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


def test_builders_over_gfp():
    A = linear_path_algebra(2, field=GF(7))
    assert A.projective(1).dimension_vector() == {1: 1, 2: 1}


def test_builders_need_provenance():
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(Exception):
        A.simple(1)
