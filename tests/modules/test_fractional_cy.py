"""Fractional Calabi-Yau dimension of self-injective algebras in the stable category
(Plan 53 / R24). Serre functor S = Omega.nu, suspension Sigma = Omega^{-1}
(Ivanov-Volkov 1212.2619 sec 1.3; Erdmann-Skowronski). Criterion: least ell with
S^ell ~ Sigma^m. Over QQ (decompose char caveat).

Tier: weak-on-generators (object-wise on the simples + orbit reps) -- a NECESSARY
condition for weak/strong CY, never a functor iso; the payload `tier` field says so."""
import pytest

from quiverlab import PreprojectiveAlgebra, Quiver, linear_path_algebra, \
    truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.fractional_cy import (fractional_calabi_yau,
                                             is_fractionally_calabi_yau)
from quiverlab.modules.fractional_cy import _stable_nf as _nf

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
def test_truncated_polynomial_is_1_CY():
    # k[x]/(x^a), a>=3: symmetric (nu = id), Omega^2 ~ id => weakly 1-CY, (m,ell)=(1,1).
    A = truncated_polynomial(3, field=QQ)
    fcy = fractional_calabi_yau(A)
    assert fcy.status == "certified"
    assert (fcy.m, fcy.ell) == (1, 1) and fcy.cy_dimension == "1/1"
    assert fcy.weakly_n_cy == 1
    A5 = truncated_polynomial(5, field=QQ)
    assert (fractional_calabi_yau(A5).m, fractional_calabi_yau(A5).ell) == (1, 1)


@lit
def test_dual_numbers_shift_trivial():
    A = truncated_polynomial(2, field=QQ)                # k[x]/(x^2): Omega ~ id
    fcy = fractional_calabi_yau(A)
    assert fcy.status == "shift-trivial" and (fcy.m, fcy.ell) == (0, 1)


@lit
def test_preprojective_Dynkin_is_2_CY():
    # stab Pi(Delta) is 2-Calabi-Yau (Geiss-Leclerc-Schroer): (m,ell)=(2,1).
    A = PreprojectiveAlgebra("A3", field=QQ)             # rep-finite (n<=4) => periodic
    fcy = fractional_calabi_yau(A)
    assert fcy.status == "certified"
    assert (fcy.m, fcy.ell) == (2, 1) and fcy.weakly_n_cy == 2


@xeng
def test_ell1_equals_ivanov_volkov_integer_form():
    # ell=1 certificate S ~ Sigma^n  <=>  Omega^{n+1} ~ nu^{-1} on the simples.
    A = truncated_polynomial(3, field=QQ)
    fcy = fractional_calabi_yau(A)
    assert fcy.ell == 1
    n = fcy.m
    for v in A.quiver.vertices:
        S = A.simple(v)
        lhs = S
        for _ in range(n + 1):
            lhs = lhs.syzygy()                           # Omega^{n+1} S
        from quiverlab.modules.hom import is_isomorphic
        # Omega^{n+1} S ~ nu^{-1} S  (stably; both non-projective here)
        assert is_isomorphic(_nf(lhs), _nf(S.nakayama_minus()))


@selfcert
def test_not_self_injective_refused():
    A = linear_path_algebra(3, field=QQ)                 # hereditary: not self-injective
    with pytest.raises(QuiverlabError, match="self-injective"):
        fractional_calabi_yau(A)
