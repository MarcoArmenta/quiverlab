"""QPA (GAP) probe for a WALL-AND-CHAMBER / stability-fan surface (Plan 63, Task 4).

Plan 63 exposes the wall-and-chamber structure of a finite-dimensional algebra: the wall
``D(B)`` of every brick as an exact rational inequality system, the chambers as g-vector
cones, and the honest brick-finite <=> tau-tilting-finite completeness gate. This file asks
the honest question: does QPA 1.37 ship any surface we could crosscheck those walls /
chambers against?

It does NOT. QPA 1.37 has NO wall-and-chamber / stability-fan surface: no WallAndChamber,
no StabilitySpace, no BrickWall, no ChamberOfAlgebra, no StabilityFunction verb, and
NamesGVars() contains zero names matching "chamber", "stability", "semistab", "brickwall",
"wallandchamber", "tautilt", or "supporttau". More fundamentally QPA has no tau-tilting /
support-tau-tilting surface at all (the Plan-45 finding) -- what it DOES ship is CLASSICAL
tilting/cotilting (``IsTiltingModule``, ``AllComplementsOfAlmostCompleteTiltingModule``,
...), a different theory. So ``#chambers = #support tau-tilting`` itself inherits NO live QPA
cross-check either -- the external cross-checks named on the verification page are the
FD-Applet / Demonet-Iyama-Jasso tables, not a QPA call. (Bare "wall" is deliberately NOT a
scan token: GAP core ships ``WallForm`` -- a reflection-group wall form, unrelated -- and bare
"tilt" is CLASSICAL tilting, which QPA legitimately has and which is NOT the tau-tilting fan.)

Consequently this probe SKIPS with an honest message (the expected outcome named in the
Task-4 brief). The covering oracles are the Plan-63 literature pins (kA2 = 5 chambers / 3
walls with D(P1) a ray; the cyclic rad^2-Nakayama N3^2 Ex. 15 = 14 chambers; the
non-thin discriminator kZ2/rad^2 = 4 walls) and the internal self-cert / cross-engine
batteries. The verification page records this honest-scope entry.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a wall-and-chamber / stability-fan surface.
_WALL_CHAMBER_NAMES = (
    "WallAndChamber",
    "WallAndChamberStructure",
    "StabilitySpace",
    "StabilityFunction",
    "BrickWall",
    "ChamberOfAlgebra",
    "GVectorFan",
    "SupportTauTiltingModules",
    "TauTiltingModules",
)


def test_qpa_exposes_no_wall_and_chamber_surface():
    lg = session.libgap_handle()

    # (1) Scan the global name table FIRST -- for no name matching a PRECISE wall-and-chamber
    #     / tau-tilting token (bare "wall" and bare "tilt" are excluded: "WallForm" is GAP
    #     core, classical "tilting" is a different theory QPA does have). This must precede any
    #     IsBoundGlobal query below: in GAP, IsBoundGlobal("Foo") REGISTERS "Foo" into
    #     NamesGVars() (as an unbound known name), so scanning after the queries would echo
    #     them back and falsely "find" a surface.
    _TOKENS = ("chamber", "stability", "semistab", "brickwall", "wallandchamber",
               "wallchamber", "tautilt", "supporttau")
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    wall_like = sorted(n for n in gvar_names
                       if any(t in n.lower() for t in _TOKENS))

    # (2) None of the named wall-and-chamber / tau-tilting entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
             for name in _WALL_CHAMBER_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not wall_like:
        pytest.skip(
            "QPA 1.37 exposes no wall-and-chamber / stability-fan surface "
            "(no WallAndChamber / StabilitySpace / BrickWall / ChamberOfAlgebra; "
            "NamesGVars has no wall/chamber/stability/brick name), and no support "
            "tau-tilting surface at all -- so #chambers = #support tau-tilting has no "
            "live QPA cross-check either (FD-Applet / DIJ tables are the named external "
            "checks). Covering oracles: the Plan-63 literature pins (kA2 5/3; cyclic "
            "rad^2-Nakayama N3^2 = 14; kZ2/rad^2 = 4 walls) + the self-cert / cross-engine "
            "batteries in tests/modules/test_wall_chamber_*.py."
        )

    # If a future QPA ever grows a wall-and-chamber / tau-tilting surface, FAIL loudly so this
    # honest-scope skip is revisited and a real crosscheck is wired in.
    pytest.fail(
        "QPA now exposes a wall-and-chamber / tau-tilting surface -- wire a real crosscheck "
        f"against the Plan-63 wall/chamber payload. Found: present={present}, "
        f"wall_like_names={wall_like}."
    )
