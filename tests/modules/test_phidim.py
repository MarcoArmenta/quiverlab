"""phidim/psidim as algebra invariants (Plan 53 / R23a). Literature: gldim finite =>
findim = phidim = psidim = gldim (survey 2310.09283 chain collapse); self-injective =>
phidim = psidim = 0 (Plan 40 phi==0). Self-cert: the standing chain
findim <= phidim <= psidim <= gldim on computed terms; the rep-finite
phidim = phi(+ all indec) route == max phi over indecomposables. Over QQ (decompose
char caveat).

Cost note (W4): phi_dim on TruncatedPathAlgebra("A4", 2) (~10 indecomposables,
M0 dim ~20) must complete well under 30 s on the CI cell -- if it ever exceeds,
shrink the oracle to "A3". The knit budget is the guard against a runaway M0."""
import pytest

from quiverlab import GF, Quiver, TruncatedPathAlgebra, linear_path_algebra, \
    truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import phi_dim, psi_dim

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
def test_finite_gldim_collapse():
    # kA3/rad^2: gl.dim 2 => findim = phidim = psidim = gldim = 2 (chain collapse).
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    assert phi_dim(A) == 2 and phi_dim(A).exact is True
    assert psi_dim(A) == 2


@lit
def test_hereditary_kA3_phidim_is_gldim():
    A = linear_path_algebra(3, field=QQ)                # hereditary, gl.dim 1
    assert phi_dim(A) == 1 and psi_dim(A) == 1


@lit
def test_self_injective_phidim_zero():
    A = truncated_polynomial(4, field=QQ)               # k[x]/(x^4): self-injective
    pd = phi_dim(A)
    assert pd == 0 and pd.exact is True and pd.status == "self-injective"
    assert psi_dim(A) == 0


@xeng
def test_phidim_equals_max_phi_over_indecomposables():
    # the +-of-all route (add-monotonicity) == max phi over the knit indecomposables.
    from quiverlab.modules.ar import knit_ar_quiver
    from quiverlab.modules.homdims import igusa_todorov_phi
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    manual = max(igusa_todorov_phi(v["module"]) for v in ar.vertices)
    assert phi_dim(A).value == manual


@selfcert
def test_standing_chain_holds():
    from quiverlab.modules.homdims import _chain_selfcheck    # returns dict of checks
    for A in (linear_path_algebra(3, field=QQ),
              TruncatedPathAlgebra("A3", 2, field=QQ),
              truncated_polynomial(4, field=QQ)):
        chk = _chain_selfcheck(A)
        assert chk["ok"]                                # findim<=phidim<=psidim<=gldim


@selfcert
def test_char_caveat_refuses_loudly():
    from quiverlab.errors import QuiverlabError
    A = TruncatedPathAlgebra("A3", 2, field=GF(2))      # char 2 <= module dims
    with pytest.raises(QuiverlabError):
        phi_dim(A)
