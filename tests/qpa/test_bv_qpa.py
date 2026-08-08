"""QPA (GAP) probe for a Batalin-Vilkovisky / BV operator surface (Plan 54, R2).

Plan 54 exposes the BV operator Delta: HH^n(A) -> HH^{n-1}(A) on the Hochschild
cohomology of a Frobenius / symmetric algebra with semisimple Nakayama automorphism.
This file asks the honest question: does QPA 1.37 ship any BV / Delta surface we
could crosscheck it against?

It does NOT. Plan 35 already established that QPA 1.37 has no Hochschild-cohomology
PRODUCT surface at all (no CupProduct / HochschildCohomologyRing*), and the BV
operator is a further structure ON that ring, so a fortiori QPA exposes no BV
operator: NamesGVars() contains no name matching "Batalin"/"Vilkovisky"/"BVOperator",
and none of the plausible entry points are bound.

Consequently this probe SKIPS with an honest message (the expected outcome). The
covering oracles for the BV operator are internal: Delta^2 = 0, the seven-term
relation, the perfect-pairing certificate, and -- decisively -- the derived bracket
equalling the INDEPENDENT Gerstenhaber bracket in-window (a cross-engine arbiter with
no Frobenius input), plus the Tradler k[x,y]/(x^2,y^2) and BIKLZ k[x]/(x^N) literature
pins. The verification page records this honest-scope entry.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a BV / Delta surface on HH^*.
_BV_NAMES = (
    "BVOperator",
    "BatalinVilkoviskyOperator",
    "BVDifferential",
    "HochschildBVOperator",
    "GerstenhaberBracket",         # not a BV operator, but the calculus it lives in
    "ConnesBOperator",
)


def test_qpa_exposes_no_bv_operator_surface():
    lg = session.libgap_handle()

    # (1) Scan the global name table FIRST (IsBoundGlobal REGISTERS a queried name,
    #     so the scan must precede the queries -- the Plan-35 precedent).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    bv_like = sorted(n for n in gvar_names
                     if "batalin" in n.lower() or "vilkovisky" in n.lower()
                     or "bvoperator" in n.lower())

    # (2) None of the named BV entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _BV_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not bv_like:
        pytest.skip(
            "QPA 1.37 exposes no BV / Delta surface on Hochschild cohomology "
            "(no BVOperator / BatalinVilkovisky*; NamesGVars has no "
            "Batalin/Vilkovisky/BVOperator name; QPA has no HH product ring at all, "
            "Plan 35). Covering oracle: Delta^2=0 + the seven-term relation + the "
            "derived-bracket == independent-Gerstenhaber arbiter + the Tradler / "
            "BIKLZ literature pins (tests/hochschild/test_bv_*.py)."
        )

    # If a future QPA ever grows a BV surface, FAIL loudly so this honest-scope skip
    # is revisited and a real crosscheck is wired in.
    pytest.fail(
        "QPA now exposes a BV/Delta surface -- wire a real crosscheck against the "
        f"Plan-54 bv_operator matrices. Found: present={present}, bv_like={bv_like}."
    )
