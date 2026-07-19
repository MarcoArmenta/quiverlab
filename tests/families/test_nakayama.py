"""NakayamaAlgebra by Kupisch series and by (n, l, cyclic). Fixture N1."""
import pytest

from quiverlab.errors import AdmissibilityError
from quiverlab.families import NakayamaAlgebra
from quiverlab.fields import CC, GF


def test_kupisch_322_is_cyclic_dim7_cartan_det1():
    A = NakayamaAlgebra([3, 2, 2], field=CC)
    assert A.dim == 7
    assert A.cartan_matrix() == [[1, 1, 1], [0, 1, 1], [1, 0, 1]]     # Fixture N1
    assert A.hochschild_cohomology(0).dims == [1]                    # HH^0 = center = 1


def test_kupisch_322_char_independent_dim():
    assert NakayamaAlgebra([3, 2, 2], field=GF(32003)).dim == 7


def test_homogeneous_cyclic_n4_l3():
    A = NakayamaAlgebra(n=4, l=3, cyclic=True, field=CC)             # kZ_4 / rad^3
    assert A.dim == 12                                              # 4 * 3


def test_linear_form_b_n4_l3():
    A = NakayamaAlgebra(n=4, l=3, cyclic=False, field=CC)           # kA_4 / rad^3, Kupisch [3,3,2,1]
    assert A.dim == 9                                               # 3 + 3 + 2 + 1


def test_linear_series_ends_in_one():
    A = NakayamaAlgebra([3, 2, 1], field=CC)                        # linear A_3, no run-off relations
    assert A.dim == 6                                               # 3 + 2 + 1
    B = NakayamaAlgebra([2, 2, 1], field=CC)                        # kills a1*a2 only
    assert B.dim == 5                                               # e1,e2,e3,a1,a2


def test_bad_series_is_loud():
    with pytest.raises(AdmissibilityError):
        NakayamaAlgebra([2, 1, 2], field=CC)                        # interior 1: neither shape


def test_citations_include_nakayama():
    assert "nakayama" in NakayamaAlgebra([3, 2, 2]).citations()
