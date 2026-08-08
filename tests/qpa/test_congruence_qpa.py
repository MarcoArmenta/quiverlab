"""QPA (GAP) probe for a torsion-lattice / congruence / wide-subcategory surface (Plan 64).

Plan 64 exposes the lattice theory of torsion classes: the finite lattice tors A, its
congruence lattice Con(tors A), the forcing order on bricks, the canonical join
representations, and the wide-subcategory poset (Enomoto's core label order). This file asks
the honest question: does QPA 1.37 ship any surface we could crosscheck those against?

It does NOT. QPA 1.37 has NO tau-tilting surface at all (the Plan-45 finding), hence no
torsion-class lattice, no congruence lattice, no brick-forcing order, no core-label /
wide-subcategory poset. NamesGVars() contains zero REPRESENTATION-THEORY torsion/wide names
(the compound tokens "tautilt", "torsionclass", "torsionpair", "widesubcat", "corelabel",
"forcingorder" -- NOT the generic "torsion"/"congruence" substrings, which GAP's Semigroups
and Polycyclic packages use heavily for semigroup congruences and group torsion subgroups,
objects unrelated to the lattice of torsion CLASSES of a finite-dimensional algebra). So
there is no lattice / congruence / wide surface to compare against -- and even #torsion /
#bricks themselves have no QPA cross-check (Plan 45 stated it), so the P64 lattice/congruence
invariants inherit no QPA oracle.

Consequently this probe SKIPS with an honest message (the expected outcome). The covering
oracles for P64 are the DIRRT / BCZ / Enomoto worked examples (the kA2 pentagon N5,
|Con(tors kA3)| = 14, #wide(kA2) = 5 = M3, #wide(kA3) = 14 = NC(A3)) + the internal
self-certs (is_lattice, Con distributive, #J(Con) == #bricks, kappa order == core label
order), NOT a live QPA call. The verification page records this honest-scope entry. If a
future QPA ever ships any of these verbs, this test FAILS loudly so the honest-scope claim
cannot silently rot.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a torsion-lattice / congruence / wide surface.
_LATTICE_NAMES = (
    "TorsionClasses",
    "TorsionClassLattice",
    "LatticeOfTorsionClasses",
    "CongruenceLattice",
    "LatticeCongruences",
    "WideSubcategories",
    "WideSubcategoryPoset",
    "Bricks",
    "BrickLabelling",
    "ForcingOrder",
    "CanonicalJoinRepresentation",
    "CoreLabelOrder",
)


def test_qpa_exposes_no_torsion_lattice_surface():
    lg = session.libgap_handle()

    # (1) Scan the global name table FIRST -- for no name matching a REP-THEORY torsion-class
    #     / tau-tilting / wide-subcategory surface. The tokens are COMPOUND on purpose: bare
    #     "torsion" / "congruence" would false-positive on GAP's Semigroups + Polycyclic
    #     packages (SemigroupCongruence, TorsionSubgroup, ...), which are UNRELATED to the
    #     lattice of torsion CLASSES of a finite-dimensional algebra. This must precede any
    #     IsBoundGlobal query below: in GAP, IsBoundGlobal("Foo") REGISTERS "Foo" into
    #     NamesGVars() (as an unbound known name), so scanning after the queries would echo
    #     them back and falsely "find" a surface.
    _TOKENS = ("tautilt", "tau_tilt", "torsionclass", "torsionpair", "widesubcat",
               "corelabel", "forcingorder", "bricklabel")
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    lattice_like = sorted(n for n in gvar_names
                          if any(t in n.lower() for t in _TOKENS))

    # (2) None of the named entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _LATTICE_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not lattice_like:
        pytest.skip(
            "QPA 1.37 exposes no torsion-lattice / congruence / wide-subcategory surface "
            "(no tau-tilting surface at all; NamesGVars has no Torsion/Congruence/Wide/"
            "Forcing name). Covering oracles: the DIRRT/BCZ/Enomoto worked examples (kA2 "
            "pentagon, |Con(tors kA3)|=14, #wide=M3/NC) + the P64 self-certs, NOT a live "
            "QPA call."
        )

    # If a future QPA ever grows one of these surfaces, FAIL loudly so this honest-scope skip
    # is revisited and a real crosscheck is wired in.
    pytest.fail(
        "QPA now exposes a torsion-lattice / congruence / wide surface -- wire a real "
        f"crosscheck against the Plan-64 lattice/congruence/wide invariants. Found: "
        f"present={present}, lattice_like_names={lattice_like}."
    )
