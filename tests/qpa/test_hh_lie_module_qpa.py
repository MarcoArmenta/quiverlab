"""QPA/GAP cross-check of HH^* as a Lie module over HH^1 (Plan 71, Task 5).

QPA 1.37 has NO Hochschild-Lie-module surface (there is no HH^*-module verb at all),
so the covering oracle is GAP's OWN MeatAxe: our action matrices rho_n(D) are fed to
GModuleByMats + MTX.CompositionFactors for an independent count of the irreducible
constituents. The SUBSTANTIVE case is the MULTI-summand k[x,y]/(x,y)^2 (gl2 acting on
HH^n with dim End = 2, 2, 6 -> >= 2 indecomposables); the sl2 = kK2 case is a trivial
irreducible sentinel. The first test is an HONEST SKIP that FAILS if QPA ever grows a
Hochschild-Lie-module surface (the Plan-35 precedent)."""
import pytest

pytestmark = pytest.mark.qpa

from quiverlab.qpa import session
from quiverlab.qpa.session import should_skip_qpa

# Every name QPA might plausibly use for a Hochschild-Lie-module surface.
_HH_LIE_MODULE_NAMES = [
    "HochschildLieModule",
    "HochschildCohomologyLieModule",
    "HochschildCohomologyAsLieModule",
    "LieModuleStructureOfHochschildCohomology",
]


@pytest.mark.skipif(should_skip_qpa(), reason="[qpa] libgap backend not installed")
def test_qpa_exposes_no_hh_lie_module_surface():
    """Honest skip: QPA has no Hochschild-Lie-module surface. FAILS if one ever appears
    so the honest-scope entry is revisited and a real crosscheck wired in."""
    lg = session.libgap_handle()
    # (1) Scan NamesGVars FIRST -- IsBoundGlobal("Foo") REGISTERS "Foo" as a known name,
    #     so a probe-before-scan would echo our probes back (Plan-35 precedent).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    hh_lie_like = sorted(
        n for n in gvar_names
        if "hochschild" in n.lower() and ("lie" in n.lower() or "module" in n.lower()))
    # A name in NamesGVars() is NOT evidence on its own: `IsBoundGlobal("Foo")` REGISTERS
    # "Foo" into the table as a known-but-UNBOUND name, and the QPA session is SHARED across
    # every test in one `-m qpa` run -- so a SIBLING probe's queries land in this scan and the
    # first draft of these probes failed each other in file order (4 of 5 red on a full run,
    # pre-dating Plan 75). Keep only names that are actually BOUND: that makes the verdict
    # order-independent, and boundness is what "QPA ships a surface" means anyway. Re-querying
    # already-registered names adds nothing new to the table.
    hh_lie_like = [n for n in hh_lie_like if bool(lg.eval(f'IsBoundGlobal("{n}")'))]
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _HH_LIE_MODULE_NAMES}
    present = [name for name, ok in bound.items() if ok]
    if not present and not hh_lie_like:
        pytest.skip(
            "QPA 1.37 exposes no Hochschild-Lie-module surface (NamesGVars has no "
            "Hochschild+Lie/Module name). Covering oracle: GAP MeatAxe on the fed "
            "rho_n matrices (the multi-summand gl2 case) + the internal module-axiom / "
            "inner-zero / weight self-certs and the ALS-toupie / truncated-Witt lit pins.")
    pytest.fail(
        "QPA now exposes a Hochschild-Lie-module surface -- wire a real crosscheck "
        f"against the Plan-71 lie_module_action. Found: present={present}, "
        f"NamesGVars-hits={hh_lie_like}")


@pytest.mark.skipif(should_skip_qpa(), reason="[qpa] libgap backend not installed")
def test_hh_lie_module_gap_multisummand_decomposition():
    """The SUBSTANTIVE GAP oracle: k[x,y]/(x,y)^2, gl2 acting on HH^n with dim End =
    2, 2, 6 -> >= 2 indecomposables. GAP's MTX.CompositionFactors on the fed
    representation is a genuine independent count (== our summands for these semisimple
    reductive modules), per degree. (Over GF(p), p > dim HH^n keeps them semisimple.)"""
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.fields import GF
    from quiverlab.qpa.crosscheck import crosscheck_hh_lie_module
    p = 101
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "x*y", "y*x", "y*y"], field=GF(p))
    r = crosscheck_hh_lie_module(A, top=2)               # gl2 modules, multi-summand
    assert any(len(o) >= 2 for _, o, _, _ in r.per_degree)   # a genuine >= 2-summand degree
    r.assert_agree()                                     # GAP factor dims == ours, per degree
    # trivial sentinel: sl2 = HH^1(kK2) acts irreducibly (1 summand) -- kept, proves little
    kron = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=GF(p))
    crosscheck_hh_lie_module(kron, top=1).assert_agree()
