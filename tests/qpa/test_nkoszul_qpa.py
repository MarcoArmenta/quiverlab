"""QPA (GAP) probe for the Plan-77 surfaces: N-Koszul (Berger), K2 (Cassidy-Shelton),
(p,q)-almost-Koszul (Brenner-Butler-King) and multi-Koszul (Herscovich).

The honest question: does QPA 1.37 ship any GENERALIZED Koszulity surface we could
crosscheck against? It does NOT. Plan 27 already recorded that QPA has ``IsQuadraticIdeal``
but no ``IsKoszul``; the generalized ladder is even further out -- there is no
``IsNKoszul`` / ``IsK2Algebra`` / ``IsAlmostKoszul`` / ``MultiKoszul``, and no surface that
reports the INTERNAL (path-length) generation degrees of Ext(k,k), which is Plan 77's
headline primitive.

So this probe SKIPS with a message naming what IS there, and FAILS loudly if QPA ever
grows any of them -- at which point a real crosscheck must be wired in.

**The live QPA legs that DO exist** are one step below: ``ExtAlgebraGenerators`` (the
HOMOLOGICAL Yoneda generator degrees, already crosschecked in
``tests/qpa/test_ext_algebra_qpa.py``) and ``IsQuadraticIdeal``. Those cover the K2 INPUT
data -- Plan 77's ``k2_certificate`` reports exactly Plan-27's ``generators_by_degree``
verbatim, so a QPA-verified generator table is a QPA-verified K2 verdict. What QPA cannot
reach at all: the internal degrees, Berger's alternation, and the (p,q) label. Those ride
on the literature pins (Berger's delta(n) on k[x]/x^N; BBK's (h-2,2) on the Dynkin
preprojectives) and the cross-engine agreement of the two N-Koszul certificates, per the
verification page.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a generalized Koszulity surface.
_PROBE_NAMES = (
    # Koszulity itself (Plan-27 finding, re-probed here so the claim cannot rot)
    "IsKoszul", "IsKoszulAlgebra", "KoszulDual",
    # N-Koszul / D-Koszul (Berger, Green-Marcos-Martinez-Villa-Zhang)
    "IsNKoszul", "IsNKoszulAlgebra", "IsDKoszul", "IsDKoszulAlgebra", "NKoszulDegree",
    # K2 (Cassidy-Shelton)
    "IsK2", "IsK2Algebra",
    # almost Koszul (Brenner-Butler-King)
    "IsAlmostKoszul", "AlmostKoszulPair", "IsPQKoszul",
    # multi-Koszul (Herscovich)
    "IsMultiKoszul", "MultiKoszul",
    # the internal-degree primitive
    "GenerationDegrees", "InternalDegrees", "GradedBettiNumbers",
)


def test_qpa_exposes_no_generalized_koszul_surface():
    lg = session.libgap_handle()

    # Scan the global name table FIRST: in GAP, IsBoundGlobal("Foo") REGISTERS "Foo"
    # into NamesGVars() as a known-but-unbound name, so scanning AFTER the queries below
    # would echo them back and falsely "find" a surface.
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    suspicious = sorted({n for n in gvar_names
                         if "koszul" in n.lower() or "internaldegree" in n.lower()})
    # A name in NamesGVars() is NOT evidence on its own, and the QPA session is SHARED
    # across a whole `-m qpa` run, so a SIBLING probe's queries land in this scan (the
    # bug Plan 75 fixed: five probes failing EACH OTHER in file order). Keep only names
    # that are actually BOUND -- that is what "QPA ships a surface" means, and it makes
    # the verdict order-independent.
    suspicious = [n for n in suspicious if bool(lg.eval(f'IsBoundGlobal("{n}")'))]

    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _PROBE_NAMES}
    present = sorted(name for name, ok in bound.items() if ok)

    # What QPA DOES offer, so the skip names what the covering oracle rides on.
    available = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
                 for name in ("ExtAlgebraGenerators", "IsQuadraticIdeal")}

    if not present and not suspicious:
        pytest.skip(
            "QPA 1.37 exposes no generalized Koszulity surface: no IsKoszul (the Plan-27 "
            "finding, re-probed), no IsNKoszul/IsDKoszul, no K2, no almost-Koszul, no "
            "multi-Koszul, and nothing reporting the INTERNAL generation degrees of "
            f"Ext(k,k). What it does offer: {available} -- ExtAlgebraGenerators is the "
            "live leg, already crosschecked in tests/qpa/test_ext_algebra_qpa.py, and it "
            "covers K2's input data (k2_certificate reports Plan-27's "
            "generators_by_degree verbatim). Covering oracles for the rest: Berger's "
            "delta(n) on k[x]/x^N, BBK's (h-2,2) on the Dynkin preprojectives, and the "
            "cross-engine agreement of the two N-Koszul certificates."
        )

    pytest.fail(
        "QPA now exposes a generalized Koszulity surface -- wire a REAL crosscheck "
        "against n_koszul_certificate / k2_certificate / almost_koszul_certificate "
        f"instead of this honest skip. Found: present={present}, "
        f"suspicious_names={suspicious}."
    )
