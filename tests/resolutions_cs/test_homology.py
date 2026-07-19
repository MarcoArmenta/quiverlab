import pytest
from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
pytest.importorskip("quiverlab.groebner")


def _A(field=CC, rels=("x*x",), arrows=None, verts=(1,)):
    return Quiver(list(verts), arrows or {"x": (1, 1)}).algebra(relations=list(rels), field=field)


def test_kx2_dims_char0_and_char2():
    assert cs_cohomology_dims(_A(), 6).dims == [2, 1, 1, 1, 1, 1, 1]
    assert cs_homology_dims(_A(), 6).dims == [2, 1, 1, 1, 1, 1, 1]
    assert cs_cohomology_dims(_A(field=GF(2)), 5).dims == [2, 2, 2, 2, 2, 2]


def test_square_dims():
    A = _A(rels=["a*b - c*d"], verts=(1, 2, 3, 4),
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    assert cs_cohomology_dims(A, 4).dims == [1, 0, 0, 0, 0]
    assert cs_homology_dims(A, 4).dims == [4, 0, 0, 0, 0]


def test_qci_homology_matches_bank_vector():
    A = _A(rels=["x*x", "y*y", "y*x - 2*x*y"], arrows={"x": (1, 1), "y": (1, 1)})
    assert cs_homology_dims(A, 12).dims == [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]


def test_engine_facade_is_resolution_protocol():
    from quiverlab.resolutions_cs.engine_facade import CSResolution
    from quiverlab.engine.resolutions import Resolution
    from quiverlab.engine.adapter import to_engine
    Ap = _A(field=GF(32003), rels=["a*b - c*d"], verts=(1, 2, 3, 4),
            arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    R = CSResolution(Ap)
    assert isinstance(R, Resolution)
    E = to_engine(Ap.unit_adapted())
    b2 = R.term_basis(E, 2)
    assert R.differential_matrix(E, 2, b2, {g: i for i, g in enumerate(R.term_basis(E, 1))}).shape[1] == len(b2)
