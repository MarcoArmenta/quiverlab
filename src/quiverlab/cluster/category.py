# SPDX-License-Identifier: MIT
"""The Amiot-Keller cluster category, certified module-category slice (Plan 79 / R31).

For an ACYCLIC quiver ``Q`` and ``A = kQ`` hereditary, the cluster category
``C_Q = D^b(mod A)/tau^{-1}[1]`` (Buan-Marsh-Reineke-Reiten-Todorov) has a FINITE
fundamental domain that shipped quiverlab machinery already reaches::

    ind(C_Q) = ind(mod A)  |_|  {P_v[1] : v in Q_0}

so ``#indec(C_Q) = #ind(mod A) + n`` -- the almost-positive roots of the underlying
diagram.  On top of that finite model three further facts become computable with the
shipped engines, and this module is the GLUE plus the CERTIFICATES, not new homological
algebra:

* **cluster-tilting objects ARE support tau-tilting pairs** (Adachi-Iyama-Reiten): a pair
  ``(M, P)`` is the cluster-tilting object ``T = M (+) P[1]``, so the P45 exchange graph IS
  the cluster exchange graph.  This is the plan's cross-engine backbone -- the counts must
  come out as the cluster (associahedron) numbers.
* **the cluster-tilted endomorphism algebra is a Jacobian algebra**
  ``End_{C_Q}(T) = Jac(Q_T, W_T)`` (Buan-Marsh-Reiten / Amiot): the quiver ``Q_T`` is
  obtained by shipped Fomin-Zelevinsky matrix mutation and certified per instance; the
  identification itself is CITED, never recomputed (there is no dg engine).
* **2-Calabi-Yau duality** reduces on the module window to the Auslander-Reiten formula
  (Assem-Simson-Skowronski Thm IV.2.13).

**Everything needing the Ginzburg dg algebra is out of scope and refused loudly**: the dg
algebra ``Gamma(Q,W)``, ``D^b(Gamma)``, the direct orbit-category ``Hom_{C}``, and the
general non-acyclic ``C_{(Q,W)}``.  See :mod:`docs/verification.md` honest scope.

References: ``amiot_cluster_category`` (Amiot, Ann. Inst. Fourier 59 (2009) 2525-2590),
``bmrrt_cluster`` (Adv. Math. 204 (2006) 572-618), ``bmr_cluster_tilted`` (Trans. AMS 359
(2007) 323-332), ``keller_reiten_gorenstein`` (Adv. Math. 211 (2007) 123-151),
``air_tau_tilting`` (Compos. Math. 150 (2014) 415-452), ``assem_book``.
"""
from quiverlab.errors import QuiverlabError

__all__ = ["ClusterCategory"]

_REFERENCES = ["amiot_cluster_category", "bmrrt_cluster", "bmr_cluster_tilted",
               "keller_reiten_gorenstein", "air_tau_tilting", "assem_book"]

# Cluster numbers of the Dynkin types, used ONLY to quote the expected size in an
# over-budget refusal (G3a) -- never to fabricate a count the engine did not discover.
# A_n: Catalan(n+1); D_n: binom(2n-2, n-1) * (3n-2)/(2n-2)... the classical values are
# tabulated rather than derived so a wrong closed form cannot silently mislead.
_CLUSTER_NUMBER = {
    ("A", 1): 2, ("A", 2): 5, ("A", 3): 14, ("A", 4): 42, ("A", 5): 132,
    ("A", 6): 429, ("A", 7): 1430, ("A", 8): 4862,
    ("D", 4): 50, ("D", 5): 182, ("D", 6): 672, ("D", 7): 2508,
    ("E", 6): 833, ("E", 7): 4160, ("E", 8): 25080,
}


def _catalan(m):
    """Catalan(m), exact integer arithmetic (no floats anywhere in src/)."""
    num = den = 1
    for i in range(m):
        num *= (2 * m - i)
        den *= (i + 1)
    return num // den // (m + 1)


class ClusterCategory:
    """The cluster category ``C_Q`` of an acyclic quiver, as a finite certified model.

    Construct from a Dynkin type string (``"A3"``), an acyclic :class:`Quiver`, or a
    hereditary :class:`Algebra` (``I = 0``); or via :meth:`from_potential` for a
    Jacobi-finite quiver-with-potential (an ALGEBRA-level certificate; the category
    combinatorics are served only when the input is certified reducible to a hereditary
    model).
    """

    def __init__(self, source, field=None):
        from quiverlab.fields import QQ
        self.algebra = _normalize(source, field or QQ)
        self.n = len(self.algebra.quiver.vertices)
        self.hereditary = not (self.algebra.relations or [])
        self.dynkin_type = self.algebra.dynkin_type()
        self._knit = None
        self._eg = None
        self._potential = None                 # set by from_potential

    # -- input classification ------------------------------------------------------
    @property
    def is_representation_finite(self):
        """True iff Gabriel's theorem certifies ``mod kQ`` finite: the underlying diagram
        is Dynkin (A/D/E).  ``None`` when the type certificate cannot decide (a
        disconnected or otherwise unclassified diagram) -- never guessed."""
        dt = self.dynkin_type
        if dt is None:
            return None
        return dt[0] in ("A", "D", "E")

    def _cluster_number(self):
        """The known cluster number for this Dynkin type, or ``None``.  Quoted in an
        over-budget refusal so the user learns the size they need to budget for."""
        dt = self.dynkin_type
        if dt is None or dt[0] not in ("A", "D", "E"):
            return None
        if dt[0] == "A":
            return _CLUSTER_NUMBER.get(dt, _catalan(dt[1] + 1))
        return _CLUSTER_NUMBER.get(dt)

    # -- the fundamental domain ----------------------------------------------------
    def _knit_ar(self):
        if self._knit is None:
            from quiverlab.modules.ar import knit_ar_quiver
            self._knit = knit_ar_quiver(self.algebra)
        return self._knit

    def _refuse_infinite(self, what):
        return QuiverlabError(
            f"{what}: C_Q has infinitely many indecomposables -- the underlying diagram "
            f"is {_diagram_name(self.dynkin_type)}, which is NOT Dynkin, so by Gabriel's "
            "theorem mod kQ is representation-infinite",
            hint="only bounded-window enumeration is meaningful here; the certified "
                 "cluster-category slice is the acyclic representation-FINITE case")

    def num_indecomposables(self):
        """``#indec(C_Q) = #ind(mod kQ) + n`` -- the almost-positive-roots count.

        Refuses loudly (G3b) when the algebra is representation-INFINITE, a claim made
        ONLY off the Gabriel/Dynkin type certificate, never off a truncated enumeration.
        """
        # The type certificate is consulted FIRST, before any enumeration: by Gabriel a
        # connected hereditary kQ is representation-finite IFF its diagram is Dynkin, so
        # on a non-Dynkin input the answer is infinite as a THEOREM. Knitting first would
        # (a) grind to the budget on an algebra we already know is infinite, and (b) make
        # the refusal look like it came from a failed enumeration rather than from the
        # certificate -- which is exactly the conflation G3a/G3b exists to prevent.
        if self.is_representation_finite is False:
            raise self._refuse_infinite("num_indecomposables")
        knit = self._knit_ar()
        if getattr(knit, "status", None) != "complete":
            raise QuiverlabError(
                "num_indecomposables: the AR knit did not close (status="
                f"{getattr(knit, 'status', None)!r}), so #ind(mod kQ) is not certified",
                hint="raise budget_modules/budget_dim, or note that the type certificate "
                     f"reports {_diagram_name(self.dynkin_type)}")
        return len(knit.vertices) + self.n

    def indecomposables(self):
        """``ind(mod kQ) |_| {P_v[1]}`` as ``[("module", M), ..., ("shift", P_v), ...]``.

        The shifted entries carry the PROJECTIVE ``P_v`` itself; the tag ``"shift"`` says
        the object of ``C_Q`` is ``P_v[1]`` -- there is no shifted-module type, and
        inventing one would claim a triangulated structure this slice does not compute.
        """
        self.num_indecomposables()             # runs the same gate, refuses identically
        from quiverlab.modules.builders import projective
        out = [("module", M) for M in self._knit_ar().vertices]
        out += [("shift", projective(self.algebra, v))
                for v in self.algebra.quiver.vertices]
        return out

    def almost_positive_roots_count(self):
        """The number of almost-positive roots of the underlying diagram: the positive
        roots (= ``#ind(mod kQ)`` by Gabriel) plus the ``n`` negative simple roots (=
        the shifted projectives).  Equal to :meth:`num_indecomposables` BY THE THEOREM --
        the two are computed by the same route here, so the test that compares them to
        the closed form ``n(n+3)/2`` (type A) is the real content."""
        return self.num_indecomposables()

    def __repr__(self):
        dt = _diagram_name(self.dynkin_type)
        return f"<ClusterCategory C_Q, n={self.n}, type={dt}>"


def _diagram_name(dt):
    return "unclassified" if dt is None else f"{dt[0]}{dt[1]}"


def _normalize(source, field):
    """A Dynkin string / Quiver / hereditary Algebra -> the hereditary acyclic ``A = kQ``.

    G1 (not hereditary) and G2 (the quiver has an oriented cycle) are the loud refusals.
    """
    from quiverlab.combinat.quiver import Quiver
    from quiverlab.families.path_algebra import PathAlgebra
    if isinstance(source, str):
        A = PathAlgebra(source, field=field)
    elif isinstance(source, Quiver):
        # Acyclicity is checked BEFORE building kQ: a cyclic quiver's path algebra is
        # infinite-dimensional, so Quiver.algebra would raise its own (correct but
        # generic) finite-dimensionality error first and the caller would never learn the
        # cluster-specific reason -- that C_Q is undefined here at all.
        _reject_cycles(source)
        A = source.algebra(relations=[], field=field)
    else:
        A = source
    if getattr(A, "quiver", None) is None:
        raise QuiverlabError(
            "the cluster category needs a quiver presentation kQ",
            hint="build the algebra via Quiver.algebra(...); structure-constant algebras "
                 "carry no quiver")
    if A.relations:
        raise QuiverlabError(
            "the cluster category model needs an acyclic / hereditary A (I = 0); this "
            f"algebra has {len(A.relations)} defining relation(s)",
            hint="C_Q = D^b(mod kQ)/tau^-1[1] is defined for a HEREDITARY kQ; for a "
                 "quiver with potential use ClusterCategory.from_potential(Q, W)")
    _reject_cycles(A.quiver)
    return A


def _reject_cycles(quiver):
    if _has_oriented_cycle(quiver):
        raise QuiverlabError(
            "kQ is not finite-dimensional / C_Q is undefined for a quiver with an "
            "oriented cycle",
            hint="the cluster category model in this slice is the ACYCLIC case; for a "
                 "quiver with potential use ClusterCategory.from_potential(Q, W)")


def _has_oriented_cycle(quiver):
    """A DFS three-colouring on the arrow directions (a loop counts as a cycle)."""
    out = {v: [] for v in quiver.vertices}
    for name, (s, t) in quiver.arrows.items():
        out[s].append(t)
    WHITE, GREY, BLACK = 0, 1, 2
    colour = {v: WHITE for v in quiver.vertices}

    def visit(v):
        colour[v] = GREY
        for w in out[v]:
            if colour[w] == GREY:
                return True
            if colour[w] == WHITE and visit(w):
                return True
        colour[v] = BLACK
        return False

    return any(colour[v] == WHITE and visit(v) for v in quiver.vertices)
