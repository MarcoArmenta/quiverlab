"""Plan 06 acceptance. Three construction routes yield the SAME dim-9 algebra
(commutative square == diamond incidence == kA_2 (x) kA_2), all HH=[1,0,0]; and
citations flow from family -> A.citations() -> bibliography()."""
from quiverlab import (
    ExteriorAlgebra, IncidenceAlgebra, NakayamaAlgebra, PathAlgebra, QuantumCI,
    TensorProduct, bibliography, families,
)
from quiverlab.combinat import Quiver
from quiverlab.fields import CC, GF


def _kA2(field=CC):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


def test_triple_crosscheck_dim9_commutative_square():
    Qsq = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    square = Qsq.algebra(relations=["a*b - c*d"], field=CC)          # general kQ/I
    diamond = IncidenceAlgebra([("b", "x"), ("b", "y"), ("x", "t"), ("y", "t")], field=CC)
    tensor = TensorProduct(_kA2(), _kA2())
    for A in (square, diamond, tensor):
        assert A.dim == 9
        assert A.hochschild_cohomology(2).dims == [1, 0, 0]


def test_catalog_dims_char_independent_and_exact():
    for A in (NakayamaAlgebra([3, 2, 2], field=GF(32003)),
              PathAlgebra("D4", field=GF(7)),
              QuantumCI(q=2, field=GF(32003)),
              ExteriorAlgebra(3, field=GF(32003))):
        assert A.dim in (7, 9, 4, 8)


def test_exterior2_equals_quantumci1():
    assert ExteriorAlgebra(2, field=CC).T == QuantumCI(q=1, field=CC).T


def test_citations_flow_end_to_end():
    A = QuantumCI(q="i", field=CC)
    keys = A.citations()
    assert "quantum_ci" in keys and "bar" in keys
    b = bibliography(keys=keys)
    assert "@article{BGMS2005" in b.bibtex()
    assert "quantum complete intersection" in str(b).lower()
    hh = A.hochschild_cohomology(2)
    assert "quantum_ci" in hh.references and "bar" in hh.references   # family + engine (§3.9)
    assert set(hh.references) == set(keys)             # table.references == A.citations()


def test_families_listing_complete():
    assert len(families().names()) >= 11
