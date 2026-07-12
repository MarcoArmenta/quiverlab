import pytest
from quiverlab import GF, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError


def test_truncated_polynomial_dims():
    A = truncated_polynomial(4)
    assert A.dim == 4
    A2 = truncated_polynomial(2, field=GF(5))
    assert A2.hochschild_cohomology(2).dims == [2, 1, 1]  # char 5 behaves like char 0 here


def test_linear_path_algebra_dim():
    # kA_n has dimension n(n+1)/2
    assert linear_path_algebra(2).dim == 3
    assert linear_path_algebra(4).dim == 10


def test_bad_arguments_fail_loudly():
    with pytest.raises(QuiverlabError):
        truncated_polynomial(1)
    with pytest.raises(QuiverlabError):
        linear_path_algebra(0)
