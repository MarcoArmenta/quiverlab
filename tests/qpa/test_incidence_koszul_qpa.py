"""QPA (GAP) probe for the Plan-75 surfaces: the incidence-vs-nerve theorem (R9) and the
GHMS comultiplicative Koszul resolution (R10).

The honest question this file asks: does QPA 1.37 ship anything we could crosscheck them
against? It does NOT, on three separate counts.

  1. **No Hochschild surface at all** -- QPA has no ``HochschildCohomology*``, so there is
     nothing to compare ``incidence_cohomology`` (or any other HH route) against.
  2. **No simplicial/order-complex surface** -- no ``OrderComplex`` / ``NerveOfPoset`` /
     ``SimplicialCohomology``: QPA is a representation-theory package, not a topology one,
     so the identity ``HH^*(kP) = H^*(Delta(P))`` has no QPA-side leg.
  3. **No Koszul surface and no GHMS resolution** -- QPA has ``IsQuadraticIdeal`` but NO
     ``IsKoszul`` (already recorded as honest scope by Plan 27), and no comultiplicative
     bimodule resolution.

So this probe SKIPS with a message naming what IS there, and FAILS loudly if QPA ever
grows any of the three -- at which point a real crosscheck must be wired in.

**The one live QPA leg that DOES exist** is the Koszul-BETTI direction, and it is already
wired: the Plan-27 ``ExtAlgebraGenerators`` crosschecks in ``tests/qpa/test_ext_algebra_qpa.py``
pin the graded Ext-algebra generator counts that ``koszul_kernels``' ``dim K_n`` must
match. That is the covering QPA oracle; the rest of P75 is covered by the theory pins
(B_3 / crown / RP^2 char split) and the cross-engine agreements (order complex == the
general CS engine; GHMS == the minimal syzygy engine == bar/CS), per the verification page.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for the three surfaces.
_PROBE_NAMES = (
    # (1) a Hochschild surface of any kind
    "HochschildCohomology", "HochschildCohomologyRing", "HochschildHomology",
    # (2) an order complex / nerve / simplicial (co)homology
    "OrderComplex", "NerveOfPoset", "SimplicialCohomology", "SimplicialHomology",
    "IncidenceAlgebra",
    # (3) Koszulity and a comultiplicative resolution
    "IsKoszul", "KoszulDual", "KoszulComplex", "IsKoszulAlgebra",
    "BimoduleResolution", "MinimalBimoduleResolution",
)


def test_qpa_exposes_no_incidence_nerve_or_koszul_surface():
    lg = session.libgap_handle()

    # Scan the global name table FIRST: in GAP, IsBoundGlobal("Foo") REGISTERS "Foo"
    # into NamesGVars() as a known-but-unbound name, so scanning AFTER the queries
    # below would echo them back and falsely "find" a surface (the Plan-35 probe's
    # lesson, verified live 2026-08-01).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    lowered = [n.lower() for n in gvar_names]
    suspicious = sorted({n for n, low in zip(gvar_names, lowered)
                         if "hochschild" in low or "koszul" in low
                         or "simplicial" in low or "ordercomplex" in low
                         or "nerve" in low})

    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _PROBE_NAMES}
    present = sorted(name for name, ok in bound.items() if ok)

    # What QPA DOES offer, so the skip says what the covering oracle rides on.
    available = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
                 for name in ("ExtAlgebraGenerators", "IsQuadraticIdeal")}

    if not present and not suspicious:
        pytest.skip(
            "QPA 1.37 exposes no Hochschild surface, no order-complex/nerve surface and "
            "no IsKoszul / comultiplicative bimodule resolution (NamesGVars has no "
            "Hochschild / Koszul / simplicial / nerve name). What it does offer: "
            f"{available} -- the Ext-algebra generator counts, already crosschecked in "
            "tests/qpa/test_ext_algebra_qpa.py, which is the live QPA leg for the "
            "Koszul-Betti direction. Covering oracles for the rest: the theory pins "
            "(B_3, crown, RP^2 char split) and the cross-engine agreements "
            "(order complex == CS; GHMS == minimal syzygy == bar/CS)."
        )

    pytest.fail(
        "QPA now exposes one of the Plan-75 surfaces -- wire a REAL crosscheck against "
        "incidence_cohomology / koszul_kernels / the GHMS resolution instead of this "
        f"honest skip. Found: present={present}, suspicious_names={suspicious}."
    )
