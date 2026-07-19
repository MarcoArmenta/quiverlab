"""Every family stamps registry-valid citation keys; A.citations() resolves them."""
import pytest

from quiverlab import citations, families
from quiverlab.families import NakayamaAlgebra, PathAlgebra, QuantumCI


@pytest.mark.parametrize("A,expected", [
    (NakayamaAlgebra([3, 2, 2]), "nakayama"),
    (PathAlgebra("A3"), "path_algebra"),
    (QuantumCI(q=2), "quantum_ci"),
])
def test_family_citations_resolve(A, expected):
    keys = A.citations()
    assert expected in keys
    for k in keys:
        citations.reference(k)                     # loud if any stamped key is bogus


def test_catalog_citations_are_all_registered():
    for info in families():
        for k in info.citations:
            citations.reference(k)
