"""Lat-Igusa-Todorov finitistic certificates (Plan 53 / R23c; Bravo-Lanzilotta-
Mendoza-Vivero 2002.07866). Each decidable family emits a proof-carrying finite
findim upper bound; the bound is never below the known exact findim. 'No known
decision procedure in general' -- NOT undecidable. Over QQ."""
import pytest

from quiverlab import Quiver, TruncatedPathAlgebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import (finitistic_dimension_bounds,
                                       lit_finitistic_certificate)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
def test_self_injective_findim_zero():
    A = truncated_polynomial(4, field=QQ)
    c = lit_finitistic_certificate(A)
    assert c.family == "self-injective" and c.findim_upper == 0


@lit
def test_gorenstein_nonselfinjective_findim_is_gorenstein_dim():
    # FAMILY 2 executing test (adversarial-review fix): kA3/J^2 is NOT self-injective
    # (family 1 does not preempt) and IS Iwanaga-Gorenstein because it has finite
    # global dimension (acyclic quiver => gl.dim < infinity => finite injective dim
    # both sides). For finite gl.dim, findim = gl.dim = id(A_A) = Gorenstein dim, so
    # the family-2 bound == int(A.global_dimension()). Asserting family == "gorenstein"
    # also confirms family 2 is the DECISIVE family (returned before family 3).
    A = TruncatedPathAlgebra("A3", 2, field=QQ)          # gl.dim 2, not self-injective
    c = lit_finitistic_certificate(A)
    assert c.family == "gorenstein"
    assert c.findim_upper == int(A.global_dimension())   # == 2 = Gorenstein dim


@selfcert
def test_lit_bound_respects_chain():
    # the LIT bound sits between the rigorous findim lower bound and gl.dim (family 3's
    # inequality findim <= phidim is the Task-A chain self-cert; here we tie the LIT
    # bound to [findim_lower, gl.dim]). Whichever family fires, the bound is valid.
    from quiverlab.modules.homdims import finitistic_dimension_bounds
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    c = lit_finitistic_certificate(A)
    lo = finitistic_dimension_bounds(A).lower
    assert c.findim_upper is not None
    assert lo <= c.findim_upper <= int(A.global_dimension())


@selfcert
def test_finitistic_upper_now_certified_where_it_was_none():
    # Plan 40 returned upper=None for self-injective (gl.dim infinite); the LIT
    # certificate now supplies the certified upper = 0 (family 1).
    A = truncated_polynomial(4, field=QQ)
    fb = finitistic_dimension_bounds(A)
    assert fb.lower == 0
    assert fb.upper == 0 and fb.note.startswith("LIT")       # was None pre-Plan-53


@selfcert
def test_lit_family4_machinery():
    # Family 4 (finite one-sided id(A_A)) is DEMOTED -- no end-to-end executing test.
    # It fires only when families 1-3 fail yet id(A_A) < infinity, i.e. id(A_A) finite
    # but id(_AA) infinite; the bounded engine never PROVES id = infinity, so this
    # precondition is not decidably reachable. The psi_D(V)+n+1 machinery is unit-covered
    # by a DIRECT call on constructed (D, n) data (never through
    # lit_finitistic_certificate).
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.homdims import _lit_family4_bound
    assert _lit_family4_bound(psi_D_of_V=2, n=1) == 4        # psi_D(V) + n + 1
    assert _lit_family4_bound(psi_D_of_V=0, n=0) == 1        # the +1 constant itself
    # arbitration safety gate: a bound below a KNOWN exact findim is a bug -> RAISE
    # (never clamp -- the Plan 40 honesty rule).
    with pytest.raises(QuiverlabError):
        _lit_family4_bound(psi_D_of_V=0, n=0, known_findim=5)


# NOTE (family 4 = finite one-sided id(A_A)): DEMOTED -- no end-to-end executing test.
# It fires only when families 1-3 fail yet id(A_A) < infinity, i.e. id(A_A) finite but
# id(_AA) infinite; the bounded engine never PROVES id = infinity, so this precondition
# is not decidably reachable. The psi_D(V)+n+1 machinery is unit-covered by a direct
# call on constructed (D, n) data (test_lit_family4_machinery), NOT through
# lit_finitistic_certificate.
