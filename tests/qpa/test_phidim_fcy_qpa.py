"""QPA (GAP) probe + prerequisite crosschecks for the Plan-53 phidim/psidim algebra
invariants and the stable-category fractional Calabi-Yau dimension.

The honest question first: does QPA 1.37 ship a phi-dimension / Igusa-Todorov
phi-dimension / stable-Calabi-Yau-dimension surface we could crosscheck the VALUES
against? It does NOT (verified via NamesGVars() -- zero names matching "IgusaTodorov",
"PhiDimension", "CalabiYau", "FractionalCalabiYau", "StableCalabiYau"). So the value
probe SKIPS honestly and FAILS only if QPA ever ships such a surface.

What QPA CAN crosscheck -- and this battery does -- are the PREREQUISITES the values
rest on:
  * IsSelfinjectiveAlgebra <-> A.is_selfinjective()  (the fractional-CY scope gate);
  * IsSymmetricAlgebra <-> A.is_symmetric()          (nu = id <=> symmetric, the
    k[x]/(x^a) fractional-CY derivation: symmetric => nu = id => nu^{-1} = id);
  * the syzygy period underpinning the CY value: QPA NthSyzygy of the simple on
    k[x]/(x^a) <-> Module.syzygy() orbit (the Omega^2 ~ id fact => (m,ell) = (1,1)).

qpa-marked: skips locally when [qpa] is absent, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab import PreprojectiveAlgebra, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.qpa import scripts, session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a phi-dimension / stable-CY-dimension surface.
_PHIDIM_CY_NAMES = (
    "IgusaTodorovFunction",
    "PhiDimension",
    "PsiDimension",
    "PhiDimensionOfAlgebra",
    "FiniteDimensionOfAlgebra",
    "CalabiYauDimension",
    "StableCalabiYauDimension",
    "FractionalCalabiYauDimension",
)


def test_qpa_exposes_no_phidim_or_stable_cy_surface():
    lg = session.libgap_handle()
    # Scan the global name table FIRST -- IsBoundGlobal REGISTERS a name into
    # NamesGVars(), so scanning after the queries would echo them back (the
    # test_products_qpa precedent, verified live 2026-08-01).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    surface_like = sorted(
        n for n in gvar_names
        if any(tok in n.lower() for tok in
               ("igusatodorov", "phidimension", "calabiyau")))
    # A name in NamesGVars() is NOT evidence on its own: `IsBoundGlobal("Foo")` REGISTERS
    # "Foo" into the table as a known-but-UNBOUND name, and the QPA session is SHARED across
    # every test in one `-m qpa` run -- so a SIBLING probe's queries land in this scan and the
    # first draft of these probes failed each other in file order (4 of 5 red on a full run,
    # pre-dating Plan 75). Keep only names that are actually BOUND: that makes the verdict
    # order-independent, and boundness is what "QPA ships a surface" means anyway. Re-querying
    # already-registered names adds nothing new to the table.
    surface_like = [n for n in surface_like if bool(lg.eval(f'IsBoundGlobal("{n}")'))]
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
             for name in _PHIDIM_CY_NAMES}
    present = [name for name, ok in bound.items() if ok]
    if surface_like or present:
        pytest.fail("QPA now ships a phidim / stable-CY surface -- wire a real value "
                    f"crosscheck: {surface_like or present}")
    pytest.skip("QPA 1.37 has NO phi-dimension / Igusa-Todorov-dimension / stable-"
                "Calabi-Yau-dimension surface; the covering oracle is the Plan-53 "
                "literature battery (k[x]/(x^a), Pi(Delta)) + the internal-consistency "
                "self-cert. This SKIP FAILS if QPA ever ships such a surface.")


def _selfinjective_qpa(A):
    base = scripts.quiver_and_algebra_script(A)
    return bool(session.run(base + "\nIsSelfinjectiveAlgebra(A);"))


def _symmetric_qpa(A):
    base = scripts.quiver_and_algebra_script(A)
    return bool(session.run(base + "\nIsSymmetricAlgebra(A);"))


@pytest.mark.parametrize("A,expected", [
    (truncated_polynomial(3, field=QQ), True),     # k[x]/(x^3): self-injective
    (truncated_polynomial(4, field=QQ), True),     # k[x]/(x^4): self-injective
    (PreprojectiveAlgebra("A3", field=QQ), True),  # Pi(A3): self-injective (periodic)
    (linear_path_algebra(3, field=QQ), False),     # kA3 (hereditary): NOT self-injective
])
def test_is_selfinjective_agrees_with_qpa(A, expected):
    # the fractional-CY SCOPE gate: our is_selfinjective must agree with QPA's
    # IsSelfinjectiveAlgebra (a non-self-injective input is refused loudly).
    assert bool(A.is_selfinjective()) == expected
    assert _selfinjective_qpa(A) == expected


@pytest.mark.parametrize("A,expected", [
    (truncated_polynomial(3, field=QQ), True),     # symmetric => nu = id (the (1,1) derivation)
    (PreprojectiveAlgebra("A3", field=QQ), False), # self-injective but NOT symmetric
    (linear_path_algebra(3, field=QQ), False),     # hereditary: not symmetric
])
def test_is_symmetric_agrees_with_qpa(A, expected):
    # nu = id <=> symmetric is exactly the k[x]/(x^a) fractional-CY derivation anchor.
    assert bool(A.is_symmetric()) == expected
    assert _symmetric_qpa(A) == expected


def test_syzygy_period_of_simple_agrees_with_qpa():
    # The Omega^2 ~ id fact on k[x]/(x^3) underpins the stable CY dimension (m,ell)=(1,1):
    # QPA's NthSyzygy(S, 2) has the same dimension as our Omega^2(S), and both equal S's.
    A = truncated_polynomial(3, field=QQ)
    S = A.simple(next(iter(A.quiver.vertices)))
    ours = S.syzygy().syzygy().dim
    base = scripts.quiver_and_algebra_script(A)
    qpa = int(session.run(base + "\nS := SimpleModules(A)[1];;\n"
                          "Dimension(NthSyzygy(S, 2));"))
    assert ours == qpa == S.dim        # Omega^2(S) ~ S (dim 1)
