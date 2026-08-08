"""Coefficient identities: HH^0=M^A, HH_0=M/[A,M], the exact duality
dim HH^n(A,DA)=dim HH_n(A,A), and the symmetric-algebra consistency web.

The non-symmetric discriminator is dim 9 (bar is exponential there), so its deep
identities are checked through the CS engine (deep, any Domain)."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import Bimodule


@pytest.mark.oracle_selfcert
def test_degree0_closed_forms():
    A = ql.truncated_polynomial(3, field=ql.CC)
    M = Bimodule.dual(A)
    c0 = A.hochschild_cohomology(0, coefficients=M).dims[0]
    h0 = A.hochschild_homology(0, coefficients=M).dims[0]
    assert c0 == M.invariants_dim()           # dim M^A = joint ker(Lact[j]-Ract[j])
    assert h0 == M.coinvariants_dim()          # dim M/[A,M]


@pytest.mark.oracle_selfcert     # W4: a theory identity on one engine's own output
def test_duality_identity_smoke_symmetric():
    # SMOKE ONLY: truncated_polynomial(3) is symmetric, so D(A) ~ A and the
    # identity is (nearly) tautological -- it exercises the plumbing, not content.
    A = ql.truncated_polynomial(3, field=ql.CC)
    assert A.is_symmetric() is True
    DA = Bimodule.dual(A)
    assert A.hochschild_cohomology(6, engine="bar", coefficients=DA).dims == \
        A.hochschild_homology(6, engine="bar").dims


@pytest.mark.oracle_selfcert     # W4
@pytest.mark.parametrize("field", [ql.CC, ql.GF(32003)])
def test_duality_identity_nonsymmetric_discriminator(field):
    # C3: the CONTENT case -- a NON-symmetric self-injective algebra where
    # D(A) not-iso A, so dim HH^n(A,DA)=dim HH_n(A,A) has real force (a wrong dual
    # or a wrong nu could not accidentally satisfy it). NakayamaAlgebra(3,3,cyclic),
    # dim 9 -- checked through CS (bar is exponential there).
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=field)
    assert A.is_symmetric() is False
    DA = Bimodule.dual(A)
    assert A.hochschild_cohomology(6, engine="cs", coefficients=DA).dims == \
        A.hochschild_homology(6, engine="cs").dims


@pytest.mark.oracle_selfcert
def test_symmetric_algebra_web():
    # A GENUINELY symmetric algebra (verified in-snippet): nu is inner => _1A_nu ~ A,
    # D(A) ~ A; all HH agree, and dim HH^n = dim HH_n follows.
    A = ql.truncated_polynomial(3, field=ql.CC)                # symmetric (verified below)
    assert A.is_symmetric() is True
    base = A.hochschild_cohomology(5).dims
    assert A.hochschild_cohomology(5, coefficients=Bimodule.dual(A)).dims == base
    assert A.hochschild_cohomology(5, coefficients=Bimodule.twisted_by_nakayama(A)).dims == base
    assert A.hochschild_homology(5).dims == base


@pytest.mark.oracle_selfcert
def test_twisted_equals_dual_nontrivial_nu():
    # C2/C3 discriminator: for self-injective A, D(A) ~ _1A_nu, so twisted and dual
    # coefficients give IDENTICAL HH dims -- non-vacuously on NON-symmetric A
    # (nu != id). Fails if twisted_by_nakayama built nu in the wrong basis (DD1b).
    # Deep (top=5) via CS (the default-engine bar route depth-falls-back to CS here).
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(32003))
    assert A.is_symmetric() is False
    assert A.hochschild_cohomology(5, engine="cs",
                                   coefficients=Bimodule.twisted_by_nakayama(A)).dims == \
        A.hochschild_cohomology(5, engine="cs", coefficients=Bimodule.dual(A)).dims
