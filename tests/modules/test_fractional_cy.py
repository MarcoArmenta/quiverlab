"""Fractional Calabi-Yau dimension of self-injective algebras in the stable category
(Plan 53 / R24). Serre functor S = Omega.nu, suspension Sigma = Omega^{-1}
(Ivanov-Volkov 1212.2619 sec 1.3; Erdmann-Skowronski). Criterion: least ell with
S^ell ~ Sigma^m. Over QQ (decompose char caveat).

Tier: weak-on-generators (object-wise on the simples + orbit reps) -- a NECESSARY
condition for weak/strong CY, never a functor iso; the payload `tier` field says so."""
import pytest

from quiverlab import ExteriorAlgebra, NakayamaAlgebra, PreprojectiveAlgebra, Quiver, \
    linear_path_algebra, truncated_polynomial
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


@selfcert
def test_fractional_cy_sign_anchor_kz3():
    # SIGN ANCHOR (P53 critic top recommendation): NakayamaAlgebra(n=3, l=2, cyclic=True)
    # = kZ_3/J^2 is self-injective with a Nakayama permutation of ORDER 3, so nu is NOT
    # its own inverse on the simples (nu(S_i) NOT iso nu^{-1}(S_i), asserted below). The
    # Serre functor S = Omega.nu then satisfies S(S_i) = Omega(nu S_i) ~ S_i = Sigma^0 S_i,
    # giving (m, ell) = (0, 1). The WRONG composition order (nu^{-1}, or Omega on the wrong
    # side) yields (1, 1) instead -- so this test PINS the nu-direction of S = Omega.nu
    # forever (the discrimination the critic verified live: Omega(nu S_i) = S_i).
    from quiverlab.modules.hom import is_isomorphic
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ)
    assert A.is_selfinjective() and not A.is_symmetric()      # nu of order 3, not identity
    for v in A.quiver.vertices:                              # nu != nu^{-1} on the simples
        S = A.simple(v)
        assert not is_isomorphic(_nf(S.nakayama()), _nf(S.nakayama_minus()))
    fcy = fractional_calabi_yau(A)
    assert fcy.status == "certified"
    assert (fcy.m, fcy.ell) == (0, 1) and fcy.weakly_n_cy == 0
    assert fcy.cy_dimension == "0/1"
    assert len(fcy.certificate) == 3                         # one witness per simple


@selfcert
def test_fractional_cy_budget_honest():
    # BUDGET PATH (P53 critic): ExteriorAlgebra(2) is a representation-INFINITE
    # self-injective local algebra whose syzygies grow (dims 3, 5, 7, ...) -- NOT
    # Omega-periodic, so no (m, ell) exists. At a SMALL dim_budget the search trips
    # PROMPTLY (< 1 s) and returns an honest status="budget" with NO certified claim
    # (never a fabricated "not CY"). Runtime is kept tight per the critic's measurement.
    E = ExteriorAlgebra(2, field=QQ)
    assert E.is_selfinjective()
    fcy = fractional_calabi_yau(E, dim_budget=8, ell_max=2, m_window=4)
    assert fcy.status == "budget"
    assert fcy.m is None and fcy.ell is None and fcy.cy_dimension is None
    assert is_fractionally_calabi_yau(E, dim_budget=8, ell_max=2, m_window=4) is False
