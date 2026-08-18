"""QPA (GAP) probe for a tau-cluster-morphism-category / picture-group surface (Plan 66).

Plan 66 exposes the tau-cluster morphism category W(A) (Buan-Marsh IMRN 2021): its objects
(= the tau-perpendicular wide subcategories), its rank-graded morphisms, the Hanson-Igusa
cube-complex classifying space (face vector + Euler characteristic + K(pi,1) verdict), and the
picture-group presentation (generators = bricks, relations per rank-2 wide, abelianization).
This file asks the honest question: does QPA 1.37 ship any surface we could crosscheck those
against?

It does NOT. QPA 1.37 has NO tau-tilting surface at all (the Plan-45/63/64 finding), hence no
tau-perpendicular reduction, no cluster-morphism category, no cube-complex classifying space,
no picture group. NamesGVars() contains zero names for any of these (the compound tokens
"clustermorphism", "picturegroup", "tauperpendicular", "cubecomplex", "widesubcat" -- NOT the
bare "cluster"/"picture"/"group" substrings, which GAP uses heavily for cluster algebras and
group theory, objects unrelated to the tau-cluster MORPHISM category of a finite-dimensional
algebra). So there is no category / classifying-space / picture-group surface to compare
against.

Consequently this probe SKIPS with an honest message (the expected outcome). The covering
oracles for P66 are the Igusa-Todorov-Weyman / Hanson-Igusa / Buan-Marsh worked examples (the
kA2 pentagon = 3 generators + 1 atom relation, the Nakayama K(pi,1), the kA3 picture group =
6 generators + 4 atom + 2 commutation, object count = #wide) + the internal cross-engine /
self-cert oracles (object count == wide_subcategories.size, #generators == #bricks, the
closed-star = #faces(C_W) Jasso link identity, face_vector[0] == object_count,
sum(face_vector) == morphism_count), NOT a live QPA call. The verification page records this
honest-scope entry. If a future QPA ever ships any of these verbs, this test FAILS loudly so
the honest-scope claim cannot silently rot.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a cluster-morphism / picture-group / cube-complex
# surface.
_CLUSTER_NAMES = (
    "TauClusterMorphismCategory",
    "ClusterMorphismCategory",
    "PictureGroup",
    "PictureGroupPresentation",
    "TauPerpendicular",
    "TauPerpendicularCategory",
    "WideSubcategories",
    "WideSubcategoryPoset",
    "CubeComplex",
    "ClassifyingSpace",
)


def test_qpa_exposes_no_cluster_morphism_surface():
    lg = session.libgap_handle()

    # (1) Scan the global name table FIRST -- for no name matching a tau-cluster-morphism /
    #     picture-group / cube-complex / tau-perpendicular surface. The tokens are COMPOUND on
    #     purpose: bare "cluster" / "picture" / "group" would false-positive on GAP's cluster-
    #     algebra + group-theory machinery, UNRELATED to the tau-cluster MORPHISM category of a
    #     finite-dimensional algebra. This must precede any IsBoundGlobal query below: in GAP,
    #     IsBoundGlobal("Foo") REGISTERS "Foo" into NamesGVars() (as an unbound known name), so
    #     scanning after the queries would echo them back and falsely "find" a surface.
    _TOKENS = ("clustermorphism", "picturegroup", "tauperpendicular", "tau_perpendicular",
               "cubecomplex", "widesubcat", "tauclustermorphism")
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    cluster_like = sorted(n for n in gvar_names
                          if any(t in n.lower() for t in _TOKENS))
    # A name in NamesGVars() is NOT evidence on its own: `IsBoundGlobal("Foo")` REGISTERS
    # "Foo" into the table as a known-but-UNBOUND name, and the QPA session is SHARED across
    # every test in one `-m qpa` run -- so a SIBLING probe's queries land in this scan and the
    # first draft of these probes failed each other in file order (4 of 5 red on a full run,
    # pre-dating Plan 75). Keep only names that are actually BOUND: that makes the verdict
    # order-independent, and boundness is what "QPA ships a surface" means anyway. Re-querying
    # already-registered names adds nothing new to the table.
    cluster_like = [n for n in cluster_like if bool(lg.eval(f'IsBoundGlobal("{n}")'))]

    # (2) None of the named entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _CLUSTER_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not cluster_like:
        pytest.skip(
            "QPA 1.37 exposes no tau-cluster-morphism-category / picture-group / cube-complex "
            "surface (no tau-tilting surface at all; NamesGVars has no ClusterMorphism/"
            "PictureGroup/TauPerpendicular/CubeComplex/WideSubcat name). Covering oracles: the "
            "ITW/HI/BM worked examples (kA2 pentagon, Nakayama K(pi,1), kA3 = 4 atom + 2 comm, "
            "object count = #wide) + the P66 cross-engine / self-cert oracles, NOT a live QPA "
            "call."
        )

    # If a future QPA ever grows one of these surfaces, FAIL loudly so this honest-scope skip
    # is revisited and a real crosscheck is wired in.
    pytest.fail(
        "QPA now exposes a tau-cluster-morphism / picture-group / cube-complex surface -- wire "
        "a real crosscheck against the Plan-66 category / classifying-space / picture-group "
        f"invariants. Found: present={present}, cluster_like_names={cluster_like}."
    )
