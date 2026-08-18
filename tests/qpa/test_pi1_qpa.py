"""QPA (GAP) probe for a fundamental-group / simple-connectivity surface (Plan 56).

Plan 56 exposes pi1(Q, I), its abelianization, the three-valued is_simply_connected
verdict, and the strongly-simply-connected / separation-condition recognizer. This
file asks the honest question: does QPA 1.37 ship any surface we could crosscheck
those against?

It does NOT. QPA 1.37 has no fundamental-group surface: no FundamentalGroup, no
IsSimplyConnected, no SeparationCondition, and NamesGVars() contains no name
matching "fundamental" / "simplyconnected" / "separation". So there is nothing to
compare against; the covering oracles are the Plan-56 literature/theory pins (the
corrected commutative square, Le Meur Example 1, the Zito simply-but-not-strongly
example, the star/tree separation discriminators) and the Hom(pi1,k+) <= dim HH^1
Hurewicz cross-check over the triangular zoo.

Consequently this probe SKIPS honestly and FAILS if a future QPA ever ships one of
those verbs, so the honest-scope claim on the verification page cannot silently rot.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a pi1 / simple-connectivity surface.
_PI1_NAMES = (
    "FundamentalGroup",
    "FundamentalGroupOfQuiverAlgebra",
    "IsSimplyConnected",
    "IsSimplyConnectedAlgebra",
    "IsStronglySimplyConnected",
    "SeparationCondition",
    "SeparationProperty",
)


def test_qpa_exposes_no_fundamental_group_surface():
    lg = session.libgap_handle()

    # (1) Scan the global name table FIRST (IsBoundGlobal REGISTERS names into
    #     NamesGVars(), so scanning after the queries would echo them back).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    pi1_like = sorted(
        n for n in gvar_names
        if any(tok in n.lower()
               for tok in ("fundamental", "simplyconnected", "separation")))

    # (2) None of the named entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _PI1_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not pi1_like:
        pytest.skip(
            "QPA 1.37 exposes no fundamental-group / simple-connectivity surface "
            "(no FundamentalGroup / IsSimplyConnected / SeparationCondition; "
            "NamesGVars has no fundamental/simplyconnected/separation name). "
            "Covering oracles: the Plan-56 literature/theory pins + the "
            "Hom(pi1,k+) <= dim HH^1 Hurewicz cross-check over the triangular zoo."
        )

    # If a future QPA grows such a surface, FAIL loudly so this honest-scope skip is
    # revisited and a real crosscheck is wired in.
    pytest.fail(
        "QPA now exposes a fundamental-group / simple-connectivity surface -- wire a "
        f"real crosscheck against the Plan-56 pi1 / separation results. Found: "
        f"present={present}, pi1_like_names={pi1_like}."
    )
