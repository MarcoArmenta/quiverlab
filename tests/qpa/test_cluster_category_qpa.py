"""QPA (GAP) probe for the Plan-79 cluster-category surface.

The honest question: does QPA 1.37 ship a cluster category, a cluster-tilting-object
construction, or a cluster-tilted-algebra builder we could crosscheck against? It does
NOT -- QPA has tilting-theory machinery but no CLUSTER category, no orbit category
`D^b/tau^-1[1]`, no cluster-tilting enumeration and no Jacobian-algebra-of-a-potential
constructor. Macaulay2 has none of it either.

So this probe SKIPS with a message naming what IS there, and FAILS loudly if QPA ever
grows any of them -- at which point a real crosscheck must be wired in.

**The covering oracles instead** are all internal or literature: the cluster numbers
themselves (Catalan(n+1) in type A, Buan-Marsh 50 for D4, 182 for D5) as frozen
literature values; the AIR bijection making the shipped tau-tilting exchange graph an
INDEPENDENT route to the same count; Keller-Reiten's Gorenstein-dimension-<=-1 theorem on
every cluster-tilted algebra built; and the 0-mismatch Auslander-Reiten sweep behind the
2-CY certificate. See docs/verification.md.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

_PROBE_NAMES = (
    # a cluster category / orbit category of any kind
    "ClusterCategory", "ClusterCategoryOfQuiver", "OrbitCategory",
    # cluster-tilting objects / the exchange graph
    "ClusterTiltingObject", "ClusterTiltingObjects", "IsClusterTilting",
    "ClusterExchangeGraph", "ClusterMutation",
    # cluster-tilted algebras and quivers with potential
    "ClusterTiltedAlgebra", "IsClusterTiltedAlgebra",
    "JacobianAlgebra", "QuiverWithPotential", "PotentialOfQuiver",
    # the 2-Calabi-Yau property
    "Is2CalabiYau", "CalabiYauDimension",
)


def test_qpa_exposes_no_cluster_category_surface():
    lg = session.libgap_handle()

    # Scan the global name table FIRST: in GAP IsBoundGlobal("Foo") REGISTERS "Foo" into
    # NamesGVars() as a known-but-unbound name, and the QPA session is SHARED across the
    # whole -m qpa run, so scanning after our own queries -- or trusting a sibling probe's
    # leftovers -- would "find" surfaces that do not exist (the bug Plan 75 fixed).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]

    def _looks_like_a_cluster_category_name(name):
        """A bare "cluster" substring is NOT evidence -- GAP's own group theory ships
        ``ClusterConjugacyPermgroups`` / ``RefineClusterConjugacyPermgroups``, which are
        about conjugacy classes of permutation groups and made the first draft of this
        probe FAIL on a live run. Require a second token that places the name in
        representation theory, so the probe stays loud about a real cluster-category
        surface without crying wolf over an unrelated namespace collision.
        """
        low = name.lower()
        if "calabiyau" in low:
            return True
        if "potential" in low and "quiver" in low:
            return True
        return "cluster" in low and any(
            tok in low for tok in ("categ", "tilt", "mutat", "exchange", "quiver",
                                   "algebra"))

    suspicious = sorted(n for n in gvar_names
                        if _looks_like_a_cluster_category_name(n))
    # Keep only names that are actually BOUND -- that is what "QPA ships a surface" means,
    # and it makes the verdict order-independent.
    suspicious = [n for n in suspicious if bool(lg.eval(f'IsBoundGlobal("{n}")'))]

    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _PROBE_NAMES}
    present = sorted(name for name, ok in bound.items() if ok)

    available = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
                 for name in ("TiltingModule", "DTr", "IsTiltingModule")}

    if not present and not suspicious:
        pytest.skip(
            "QPA 1.37 exposes no cluster category, no cluster-tilting-object "
            "construction, no cluster-tilted-algebra builder and no quiver-with-potential "
            f"Jacobian constructor (NamesGVars carries no such name). What it does offer: "
            f"{available} -- classical tilting theory, one level below the cluster "
            "category. Covering oracles for Plan 79: the cluster numbers as literature "
            "values (Catalan(n+1); Buan-Marsh 50 for D4, 182 for D5), the AIR bijection "
            "making the tau-tilting exchange graph an independent route to the same "
            "count, Keller-Reiten Gorenstein-dim <= 1 on every cluster-tilted algebra "
            "built, and the 0-mismatch Auslander-Reiten sweep behind the 2-CY certificate."
        )

    pytest.fail(
        "QPA now exposes a cluster-category surface -- wire a REAL crosscheck against "
        "ClusterCategory instead of this honest skip. Found: "
        f"present={present}, suspicious_names={suspicious}."
    )
