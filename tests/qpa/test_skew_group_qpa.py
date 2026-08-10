"""QPA (GAP) battery for the skew group algebra A⋊G (Plan 74, acceptance #6).

QPA 1.37 has NO skew-group-algebra HH-DECOMPOSITION surface (no isotypic /
conjugacy-class / Z(g)-invariants verb) -- the honest question, asked with a
`NamesGVars()` fail-if-appears guard (the Plan-35 products precedent). What QPA
CAN do is build a skew group algebra BY HAND (as its presented kQ'/I') and compute
its Hochschild cohomology dims; we use that as an INPUT-LEVEL crosscheck of the
DIRECT HH^•(A⋊G) on the small Z/2 / Z/3 examples (a construction anchor, NOT a
decomposition compare -- there is nothing on the QPA side to compare the
decomposition against).

qpa-marked: skips locally without GAP, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.families.basic import truncated_polynomial
from quiverlab.families.skew_group import GroupAction, QuiverAutomorphism
from quiverlab.fields import GF, QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Names QPA might plausibly use for a skew-group / smash-product HH surface.
_SKEW_HH_NAMES = (
    "SkewGroupAlgebra",
    "SmashProduct",
    "HochschildCohomologyOfSkewGroupAlgebra",
    "ConjugacyClassDecomposition",
    "GroupActionOnHochschild",
)


def test_qpa_exposes_no_skew_group_hh_decomposition_surface():
    lg = session.libgap_handle()
    # (1) Scan the global name table FIRST (IsBoundGlobal REGISTERS a queried name,
    #     so scanning after would echo our probes back -- the P35 lesson).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    skew_like = sorted(n for n in gvar_names
                       if "smash" in n.lower()
                       or ("skew" in n.lower() and "group" in n.lower()))
    assert not skew_like, (
        "QPA appears to ship a skew-group/smash surface now: " + str(skew_like)
        + " -- revisit the honest-scope entry and add a real decomposition compare.")
    # (2) None of the named skew-group-HH entry points are bound (callable).
    present = [name for name in _SKEW_HH_NAMES
               if bool(lg.eval(f'IsBoundGlobal("{name}")'))]
    assert not present, (
        "QPA now binds a skew-group-HH surface: " + str(present)
        + " -- add a decomposition crosscheck.")


def test_direct_hh_of_smash_matches_qpa_z2_dual_numbers():
    # Z/2 on k[x]/(x^2): A⋊G is the orbit 2-cycle Nakayama kQ'/I' (presented_form),
    # HH^• = [1,1,1,1] (direct). Feed QPA the presented algebra and compare its HH.
    A = truncated_polynomial(2, field=QQ)
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": -1})
    S = A.skew_group(GroupAction.cyclic(sig, 2))
    pf = S.presented_form()               # 1⇄2 / rad^2, ≅ A⋊G (HH-faithful)
    assert pf.hochschild_cohomology(3).dims == [1, 1, 1, 1]
    pf.crosscheck("hochschild", 3).assert_agree()


def test_direct_hh_of_smash_matches_qpa_z3_cubic():
    # Z/3 on k[x]/(x^3): A⋊G is dim 9; presented_form needs char > dim (finding 3c),
    # so this non-involution example uses GF(13) (13 = 1 mod 3; 3^3 = 27 = 1 mod 13;
    # dim 9 < 13). Feed QPA the presented form and crosscheck the DIRECT HH^•.
    A = truncated_polynomial(3, field=GF(13))
    sig = QuiverAutomorphism({1: 1}, {"x": "x"}, {"x": 3})    # 3^3 = 1 mod 13
    S = A.skew_group(GroupAction.cyclic(sig, 3))
    pf = S.presented_form()
    assert pf.hochschild_cohomology(2).dims == S.hochschild_cohomology(2).dims == [1, 1, 1]
    pf.crosscheck("hochschild", 2).assert_agree()
