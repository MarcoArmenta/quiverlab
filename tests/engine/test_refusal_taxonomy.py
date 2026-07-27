"""Exception-taxonomy refusals: two loud refusals that used to raise a bare
`ValueError`, now raise `QuiverlabError` so `except quiverlab.QuiverlabError`
catches them uniformly across domains (spec 7 taxonomy).

  * nakayama_automorphism on a non-Frobenius algebra -- the GF(p) engine path
    (engine.coxeter) matched the QQ/generic path (invariants.frobenius) only
    after this fix; both now raise QuiverlabError.
  * a presentation-less (structure-constants) algebra routed to the
    Chouhy-Solotar engine (resolutions_cs.build.reduction_system_of).

These are contract tests (unmarked / not an oracle class)."""
import pytest

import quiverlab
from quiverlab import Quiver
from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import CC, GF, QQ


def _kA3(field):
    """kA_3 = 1 -> 2 -> 3, hereditary hence NOT Frobenius (no Nakayama auto)."""
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=field)


# --------------------------------------------------------------------------
# Fix 4: nakayama_automorphism refusal is QuiverlabError on BOTH domains
# --------------------------------------------------------------------------
def test_nakayama_engine_path_raises_quiverlab_error():
    """The GF(p) engine path (was a bare ValueError)."""
    from quiverlab.engine.coxeter import nakayama_automorphism
    from quiverlab.engine.scan2 import kA
    with pytest.raises(QuiverlabError):
        nakayama_automorphism(kA(3), 5)


@pytest.mark.parametrize("field", [GF(5), QQ], ids=["GF(5)", "QQ"])
def test_nakayama_public_refusal_same_type_across_domains(field):
    """kA_3 is not Frobenius; the public refusal is QuiverlabError over GF(5)
    (engine route) AND QQ (generic route) -- one taxonomy, both domains."""
    A = _kA3(field)
    assert A.is_frobenius() is False
    with pytest.raises(QuiverlabError):
        A.nakayama_automorphism()


# --------------------------------------------------------------------------
# Fix 5: CS on a presentation-less algebra is a QuiverlabError (was ValueError)
# --------------------------------------------------------------------------
def _dual_numbers_sc(field):
    """k[x]/(x^2) via raw structure constants: NO quiver, NO relations."""
    T = [[[1, 0], [0, 1]],
         [[0, 1], [0, 0]]]
    return Algebra.from_structure_constants(T, [1, 0], field=field)


@pytest.mark.parametrize("field", [GF(5), QQ, CC], ids=["GF(5)", "QQ", "CC"])
def test_cs_needs_presentation_raises_quiverlab_error(field):
    A = _dual_numbers_sc(field)
    with pytest.raises(QuiverlabError) as ei:
        A.hochschild_cohomology(2, engine="cs")
    assert "Quiver.algebra" in str(ei.value)


def test_cs_refusal_reduction_system_direct():
    """Same refusal at the source (resolutions_cs.build.reduction_system_of)."""
    from quiverlab.resolutions_cs.build import reduction_system_of
    A = _dual_numbers_sc(QQ)
    with pytest.raises(QuiverlabError):
        reduction_system_of(A)
