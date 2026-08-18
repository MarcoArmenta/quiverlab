# SPDX-License-Identifier: MIT
"""K2 (Cassidy-Shelton): E(A)=Ext(k,k) generated in cohomological degrees 1,2, decided
through an explicit certified window (Plan 77 / R36)."""
import pytest

from quiverlab import (PreprojectiveAlgebra, TruncatedPathAlgebra, linear_path_algebra,
                       truncated_polynomial)
from quiverlab.fields import QQ
from quiverlab.modules.nkoszul import k2_certificate

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("N", [3, 4, 5])
def test_kxN_is_K2_through_window(N):
    c = k2_certificate(truncated_polynomial(N, field=QQ), top=8)
    assert c["generator_degrees"] == [1, 2]       # K2: gens in degrees 1,2
    assert c["verdict"] is None                   # self-injective: window-certified
    assert c["complete"] is False
    assert c["window"] >= 8
    assert "window" in c["reason"]


@lit
def test_finite_gldim_K2_is_complete_true():
    c = k2_certificate(TruncatedPathAlgebra("A4", 3, field=QQ), top=6)
    assert c["generator_degrees"] == [1, 2]
    assert c["verdict"] is True and c["complete"] is True


@lit
@pytest.mark.parametrize("typ", ["A3", "A4"])
def test_preprojective_is_NOT_K2(typ):
    # THE discriminator: Pi(A_{n>=3}) has a genuine Yoneda generator in cohomological
    # degree 3, so K2 fails DECISIVELY (no window caveat -- the obstruction is
    # exhibited inside the window).
    c = k2_certificate(PreprojectiveAlgebra(typ, field=QQ), top=6)
    assert 3 in c["generator_degrees"]
    assert c["verdict"] is False


@selfcert
def test_koszul_implies_K2():
    # Koszul (degree-1 gens only) => K2. kA2 is hereditary, gl.dim 1.
    c = k2_certificate(linear_path_algebra(2, field=QQ), top=6)
    assert c["generator_degrees"] == [1] and c["verdict"] is True


@selfcert
def test_preprojective_A2_is_koszul_hence_K2_through_window():
    # Pi(A2) = kZ2/rad^2 is genuinely Koszul but has infinite gl.dim, so K2 is
    # certified through the window and NEVER upgraded to an unconditional True.
    A = PreprojectiveAlgebra("A2", field=QQ)
    c = k2_certificate(A, top=8)
    assert c["generator_degrees"] == [1]
    assert c["verdict"] is None and c["complete"] is False
    assert A.ext_algebra(8).koszul is True


@xeng
def test_k2_generator_degrees_are_plan27_generators_verbatim():
    # P77 reports Plan-27's generators_by_degree, unmodified -- the named overlap.
    for A in (truncated_polynomial(4, field=QQ),
              PreprojectiveAlgebra("A3", field=QQ),
              TruncatedPathAlgebra("A4", 3, field=QQ)):
        Y = A.ext_algebra(6)
        c = k2_certificate(A, top=6, yoneda=Y)
        assert c["generator_degrees"] == sorted(Y.generators_by_degree)
        assert c["window"] == Y.certified_through_degree
