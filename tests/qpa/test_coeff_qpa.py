"""QPA (GAP) probe for Hochschild-with-COEFFICIENTS / RELATIVE Hochschild (Plan 52).

Plan 52 computes HH^*(A, M) / HH_*(A, M) for an arbitrary bimodule M and relative
HH_*(A|kQ_0, M). This file asks the honest question: does QPA 1.37 ship any
surface we could crosscheck those against?

It does NOT. QPA 1.37 exposes ordinary Hochschild homology/cohomology DIMENSIONS
of an algebra with itself (regular coefficients) only -- there is no name in
NamesGVars() for a with-coefficients or relative Hochschild surface. So this probe
SKIPS honestly and FAILS if a matching verb ever appears (then a real crosscheck
is owed). The covering oracles are the internal identities (test_coeff_identities),
the ported GF(p) bank (engine/bimodule.py), the CS/minimal cross-engine web, and
the literature pins (BLOCKED-until-transcribed).

qpa-marked (by directory): skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Names QPA might plausibly use for a with-coefficients / relative Hochschild surface.
_COEFF_HH_NAMES = (
    "HochschildCohomologyWithCoefficients",
    "HochschildHomologyWithCoefficients",
    "RelativeHochschildCohomology",
    "RelativeHochschildHomology",
    "HochschildCohomologyOfBimodule",
    "HochschildHomologyOfBimodule",
)


def test_qpa_exposes_no_coefficient_or_relative_hochschild_surface():
    lg = session.libgap_handle()
    # (1) Scan the global name table FIRST (IsBoundGlobal REGISTERS the queried
    #     name into NamesGVars(), so scanning after the queries would echo them).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    coeff_like = sorted(
        n for n in gvar_names
        if ("hochschild" in n.lower()
            and ("coeff" in n.lower() or "bimodule" in n.lower() or "relative" in n.lower())))
    # (2) None of the named entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _COEFF_HH_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not coeff_like:
        pytest.skip(
            "QPA 1.37 exposes no Hochschild-with-coefficients / relative-HH surface "
            "(NamesGVars has no matching name; only ordinary HH dimensions of A with "
            "itself). Covering oracles: the internal identities + the ported GF(p) "
            "bank + the CS/minimal cross-engine web + the literature pins.")

    pytest.fail(
        "QPA now exposes a with-coefficients / relative Hochschild surface -- wire a "
        f"real crosscheck. Found: present={present}, coeff_like_names={coeff_like}.")
