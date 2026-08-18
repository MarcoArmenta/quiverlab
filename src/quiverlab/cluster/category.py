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
        self._require_category_model("num_indecomposables")
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
        # ARQuiver.vertices is a list of RECORD DICTS {"name", "dimvec", "module"};
        # the module itself is what belongs in the fundamental domain.
        out = [("module", rec["module"]) for rec in self._knit_ar().vertices]
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


# ---------------------------------------------------------------------------------
# Cluster-tilting objects = support tau-tilting pairs (Adachi-Iyama-Reiten)
# ---------------------------------------------------------------------------------
def _exchange_graph_of(A, budget_pairs):
    from quiverlab.tautilting.mutation import exchange_graph
    return exchange_graph(A, budget_pairs=budget_pairs)


def _cluster_tilting_methods():
    """Attached to :class:`ClusterCategory` below; kept in one block so the AIR
    bijection and its certification gate read together."""


def _num_cluster_tilting(self, budget_pairs=512):
    """The number of cluster-tilting objects of ``C_Q``, with its provenance.

    **The AIR bijection is the backbone:** a support tau-tilting pair ``(M, P)`` IS the
    cluster-tilting object ``T = M (+) P[1]`` (Adachi-Iyama-Reiten), and every
    cluster-tilting object of ``C_Q`` decomposes uniquely so.  Hence the shipped P45
    exchange graph IS the cluster exchange graph and its vertex count IS the cluster
    (associahedron) number -- which is what makes this a genuine CROSS-ENGINE oracle
    rather than a definition.

    **Certification is a RUNTIME property, not a type list** -- the count is certified
    iff, at compute time, ``status == "complete"`` OR (``status == "error"`` AND the
    Plan-63 n-regularity recovery closes the graph).  Anything else is refused loudly and
    the discovered vertex count is NEVER emitted as if certified.  The three refusals are
    kept distinct because conflating them would slander a perfectly finite algebra:

    * **G3a over-budget** -- a ``"budget"`` stop on a known-finite Dynkin input means
      "rep-finite, but the cluster number exceeds the budget"; the known cluster number is
      quoted in ``expected_count`` so the caller can size the budget.  NEVER an
      infiniteness claim.
    * **G3b infinite** -- claimed ONLY off the Gabriel/Dynkin type certificate.
    * **G4 non-recoverable error** -- the BFS hit a ``mutate`` boundary and the graph is
      not n-regular, so completeness cannot be certified either way.

    ``count`` is what the ENGINE certified and is ``None`` in every refusal; the known
    literature value never leaks into it (it lives in ``expected_count``).
    """
    self._require_category_model("num_cluster_tilting")
    dt = self.dynkin_type
    dt_name = None if dt is None else _diagram_name(dt)
    out = {"count": None, "certified": False, "status": None,
           "dynkin_type": dt_name, "expected_count": self._cluster_number(),
           "budget_pairs": budget_pairs, "note": None}
    eg = _exchange_graph_of(self.algebra, budget_pairs)
    self._eg = eg
    out["status"] = eg.status
    discovered = len(eg.vertices)

    if eg.status == "complete":
        out.update(count=discovered, certified=True,
                   note="the exchange-graph BFS closed natively (status='complete')")
        return out

    if eg.status == "error":
        from quiverlab.tautilting.wallchamber import _closed_by_n_regularity
        if _closed_by_n_regularity(eg, self.n):
            out.update(count=discovered, certified=True,
                       note=("the BFS reported status='error' (a mutate boundary), but "
                             f"the discovered graph is {self.n}-regular and was not "
                             "budget-capped, so it IS closed under mutation -- "
                             "completeness recovered by the Plan-63 n-regularity "
                             "certificate"))
            return out
        out["note"] = (
            f"the exchange graph did not close (status='error') and the {self.n}-"
            "regularity recovery FAILED, so completeness is not certified; the "
            f"{discovered} discovered pairs are NOT reported as the cluster number. "
            "This is a mutate-boundary limitation (Plan 65 owns the root cause), NOT a "
            "claim that C_Q is infinite")
        return out

    # status == "budget": rep-finite-but-over-budget vs genuinely infinite.
    if self.is_representation_finite:
        expected = self._cluster_number()
        quoted = ("" if expected is None
                  else f"; its cluster number is {expected}")
        out["note"] = (
            f"rep-finite (Dynkin {dt_name}){quoted}, which exceeds budget_pairs="
            f"{budget_pairs} -- raise the budget. This is an OVER-BUDGET stop, NOT an "
            "infiniteness claim")
        return out
    if self.is_representation_finite is False:
        out["note"] = (
            f"C_Q has infinitely many cluster-tilting objects: the diagram {dt_name} is "
            "NOT Dynkin, so by Gabriel mod kQ is representation-infinite and the "
            "exchange graph does not close")
        return out
    out["note"] = (
        f"the BFS stopped at budget_pairs={budget_pairs} and the type certificate could "
        "not classify the diagram, so neither finiteness nor infiniteness is claimed")
    return out


def _cluster_tilting_objects(self, budget_pairs=512):
    """The cluster-tilting objects as ``T = M (+) P[1]`` (the AIR reading of each
    support tau-tilting pair).  Refuses whenever :meth:`num_cluster_tilting` refuses --
    an uncertified enumeration must not be handed out as if it were the whole set."""
    info = self.num_cluster_tilting(budget_pairs=budget_pairs)
    if not info["certified"]:
        raise QuiverlabError(
            "cluster_tilting_objects: " + str(info["note"]),
            hint="the enumeration is not certified complete, so the list would be a "
                 "partial set presented as the whole -- refused")
    out = []
    for rec in self._eg.vertices:
        pair = rec["pair"]
        out.append({"pair": pair, "summands": list(pair.summands),
                    "support": sorted(pair.support), "label": rec["label"],
                    "as_object": f"{rec['label']}: T = M (+) P[1]",
                    "rank": len(pair.summands) + len(pair.support)})
    return out


ClusterCategory.num_cluster_tilting = _num_cluster_tilting
ClusterCategory.cluster_tilting_objects = _cluster_tilting_objects


# ---------------------------------------------------------------------------------
# The cluster-tilted End-algebra End_{C_Q}(T) = Jac(Q_T, W_T)  (BMR / Amiot)
# ---------------------------------------------------------------------------------
def _quiver_from_exchange_matrix(B, verts):
    """The quiver of a skew-symmetric exchange matrix: ``B[i][j] > 0`` means ``B[i][j]``
    arrows ``verts[i] -> verts[j]`` (and the negative entry is the same information read
    from the other side, so only ``i < j`` is walked)."""
    from quiverlab.combinat.quiver import Quiver
    arrows, idx = {}, 0
    for i in range(len(B)):
        for j in range(i + 1, len(B)):
            m = B[i][j]
            if m == 0:
                continue
            s, t, mult = ((verts[i], verts[j], m) if m > 0
                          else (verts[j], verts[i], -m))
            for _ in range(mult):
                idx += 1
                arrows[f"a{idx}"] = (s, t)
    return Quiver(list(verts), arrows)


def _oriented_3_cycles(Q):
    """EVERY oriented 3-cycle of ``Q``, deduplicated up to cyclic rotation.

    Enumerating all of them (not just the first) is what makes the canonical type-A
    potential ``W = sum of oriented 3-cycles`` right on a cluster-tilted quiver carrying
    several triangles -- a single-triangle flagship would never exercise this.
    """
    names = sorted(Q.arrows)
    seen, out = set(), []
    for a in names:
        for b in names:
            if b == a or Q.target(a) != Q.source(b):
                continue
            for c in names:
                if c in (a, b) or Q.target(b) != Q.source(c):
                    continue
                if Q.target(c) != Q.source(a):
                    continue
                cyc = (a, b, c)
                # canonicalize up to rotation: start at the smallest arrow name
                k = min(range(3), key=lambda r: cyc[r])
                key = cyc[k:] + cyc[:k]
                if key in seen:
                    continue
                seen.add(key)
                out.append(key)
    return out


def _cluster_tilted_algebra(self, mutation_seq, field=None):
    """``End_{C_Q}(T) = Jac(Q_T, W_T)`` for ``T`` reached from the initial cluster-tilting
    object by the given Fomin-Zelevinsky mutation sequence (BMR / Amiot).

    **The identification is CITED, not computed.** There is no dg engine here, so
    ``Hom_{C_Q}`` is never formed directly; what IS computed and certified per instance is
    (a) the cluster-tilted QUIVER ``Q_T`` by shipped FZ matrix mutation, and (b) the
    Jacobian algebra of ``(Q_T, W_T)`` with its presentation verified (dimension, Cartan
    matrix, self-injectivity).

    ``mutation_seq`` entries are VERTEX POSITIONS in ``Q.vertices`` order (the
    ``matrix_mutation`` convention).

    **G5:** the canonical potential is the sum of oriented 3-cycles, which is the right
    potential in type A.  When the mutated quiver carries no 3-cycle at all, or the input
    is not type A, the QUIVER is still returned (FZ-certified) with ``algebra=None`` and a
    note -- inferring a general potential needs DWZ mutation / right-equivalence, deferred
    since P48.1.  A potential is never guessed.
    """
    self._require_category_model("cluster_tilted_algebra")
    from quiverlab.surfaces.flip import exchange_matrix, matrix_mutation
    Q0 = self.algebra.quiver
    verts = list(Q0.vertices)
    B = exchange_matrix(Q0)
    for k in mutation_seq:
        if not (0 <= k < len(verts)):
            raise QuiverlabError(
                f"cluster_tilted_algebra: mutation index {k} is out of range for "
                f"{len(verts)} vertices",
                hint="mutation_seq entries are VERTEX POSITIONS in Q.vertices order")
        B = matrix_mutation(B, k)
    QT = _quiver_from_exchange_matrix(B, verts)
    out = {"quiver": {"vertices": list(QT.vertices),
                      "arrows": {a: [QT.source(a), QT.target(a)] for a in QT.arrows}},
           "exchange_matrix": [row[:] for row in B],
           "algebra": None, "dim": None, "self_injective": None,
           "three_cycles": 0, "verified": [], "note": None}

    cycles = _oriented_3_cycles(QT)
    out["three_cycles"] = len(cycles)
    # The FZ certificate: rebuilding the quiver from B and re-reading its exchange matrix
    # must return B itself -- the per-instance quiver-level check (surfaces.flip precedent).
    out["fz_certified"] = (exchange_matrix(QT) == B)
    if not out["fz_certified"]:
        out["note"] = ("the rebuilt quiver's exchange matrix does not match the mutated "
                       "matrix -- the quiver step is NOT certified for this instance")
        return out
    out["verified"].append("fz_quiver_mutation")

    if not cycles:
        out["note"] = ("the cluster-tilted quiver carries no oriented 3-cycle, so the "
                       "canonical type-A potential is empty; inferring a general "
                       "potential needs DWZ mutation / right-equivalence, which is out "
                       "of scope (deferred since P48.1) -- the QUIVER is certified, the "
                       "algebra is not built")
        return out
    if self.dynkin_type is None or self.dynkin_type[0] != "A":
        out["note"] = (
            "the canonical 'sum of oriented 3-cycles' potential is the type-A potential; "
            f"this input is type {_diagram_name(self.dynkin_type)}, where the potential "
            "needs DWZ mutation / right-equivalence (out of scope, P48.1) -- the QUIVER "
            f"is FZ-certified and carries {len(cycles)} oriented 3-cycle(s), but no "
            "potential is guessed")
        return out

    from quiverlab.families.jacobian import JacobianAlgebra, Potential
    from quiverlab.fields import QQ
    W = Potential(QT, [(1, cyc) for cyc in cycles])
    JA = JacobianAlgebra(QT, W, field=field or self.algebra.domain or QQ)
    out["algebra"] = JA
    out["dim"] = JA.dim
    out["verified"].append("jacobian_built")
    from quiverlab.modules.ext import is_selfinjective
    out["self_injective"] = bool(is_selfinjective(JA))
    out["verified"].append("self_injective")
    out["note"] = (
        f"W = the sum of all {len(cycles)} oriented 3-cycle(s) of Q_T (the canonical "
        "type-A potential). End_C(T) = Jac(Q_T, W_T) is Buan-Marsh-Reiten / Amiot -- "
        "CITED, not recomputed: no dg engine forms Hom_C directly. What is verified here "
        "is the FZ quiver mutation and the Jacobian algebra's own presentation")
    return out


ClusterCategory.cluster_tilted_algebra = _cluster_tilted_algebra


# ---------------------------------------------------------------------------------
# The 2-Calabi-Yau certificate, scoped to the module window
# ---------------------------------------------------------------------------------
def _two_cy_certificate(self):
    """The 2-CY duality on the module block, via the Auslander-Reiten formula.

    In ``C_Q`` one has (BMRRT) ``Ext^1_C(X,Y) = Ext^1_A(X,Y) (+) D Ext^1_A(Y,X)``, so
    ``dim Ext^1_C(X,Y) = e_XY + e_YX`` is symmetric BY CONSTRUCTION -- that half is
    bookkeeping, and this payload says so rather than dressing ``a+b = b+a`` up as a
    verified theorem.

    **The computable content is the Auslander-Reiten formula** (Assem-Simson-Skowronski
    Thm IV.2.13): ``Ext^1_A(X,Y) = D \\overline{Hom}_A(Y, tau X)``.  The certificate uses
    the ORDINARY ``Hom(Y, tau X)``; the theorem is about the injectively-stable
    ``\\overline{Hom}``.  On a representation-finite HEREDITARY algebra the two coincide
    across the tested window -- that is the theorem-backed proxy, and the evidence is the
    0-mismatch sweep this function performs and reports (``mismatches``), not an
    assumption.  A mismatch is REPORTED, never smoothed away.
    """
    from quiverlab.modules.ext import ext_dims
    from quiverlab.modules.hom import hom_dim
    self.num_indecomposables()                       # the same G3b/knit gate
    mods = [rec["module"] for rec in self._knit_ar().vertices]
    A = self.algebra
    taus = [M.tau() for M in mods]
    ext1 = {}
    mismatches = []
    for i, X in enumerate(mods):
        for j, Y in enumerate(mods):
            e = ext_dims(A, X, Y, 1)[1]
            ar = hom_dim(Y, taus[i]) if taus[i] is not None else 0
            if e != ar:
                mismatches.append({"X": i, "Y": j, "ext1_A": e, "hom_Y_tauX": ar})
            ext1[(i, j)] = e
    n_mods = len(mods)
    ext1_C = {(i, j): ext1[(i, j)] + ext1[(j, i)]
              for i in range(n_mods) for j in range(n_mods)}
    symmetric = all(ext1_C[(i, j)] == ext1_C[(j, i)]
                    for i in range(n_mods) for j in range(n_mods))
    return {
        "ar_formula_holds": not mismatches,
        "mismatches": mismatches,
        "pairs_checked": n_mods * n_mods,
        "modules_checked": n_mods,
        "ext1_A": {f"{i},{j}": v for (i, j), v in sorted(ext1.items())},
        "ext1_C": {f"{i},{j}": v for (i, j), v in sorted(ext1_C.items())},
        "ext1_C_symmetric": symmetric,
        "nonzero_ext1_A_pairs": sorted(k for k, v in ext1.items() if v),
        "scope": "module_window",
        "shifted_by_citation": True,
        "note": ("the AR formula Ext^1_A(X,Y) = D Hom-bar(Y, tau X) (ASS Thm IV.2.13) is "
                 "the computed content, checked on every ORDERED pair of ind(mod kQ) with "
                 "ordinary Hom as the theorem-backed proxy for the injectively-stable "
                 "Hom-bar (they coincide on a representation-finite hereditary algebra "
                 "across this window -- the 0-mismatch sweep is the evidence). The "
                 "symmetry of dim Ext^1_C = e_XY + e_YX is BY CONSTRUCTION, not a "
                 "verified theorem. The shifted P_v[1] pairs are NOT computed: they hold "
                 "by BMRRT, cited"),
    }


def _is_2_calabi_yau(self):
    """A SCOPED verdict, never a bald categorical ``True``.

    The module block is verified by the AR formula; the pairs involving the shifted
    projectives ``P_v[1]`` are asserted by the BMRRT theorem and are NOT computed here
    (that would need the triangulated structure of the orbit category).  The payload
    always says which is which.
    """
    cert = self.two_cy_certificate()
    return {
        "verdict": bool(cert["ar_formula_holds"] and cert["ext1_C_symmetric"]),
        "scope": "module_window",
        "shifted_by_citation": True,
        "pairs_checked": cert["pairs_checked"],
        "note": ("2-CY verified on the module block via the AR formula (ASS IV.2.13) "
                 f"across all {cert['pairs_checked']} ordered pairs of ind(mod kQ); the "
                 "shifted P_v[1] pairs hold by BMRRT (cited, not computed)"),
    }


ClusterCategory.two_cy_certificate = _two_cy_certificate
ClusterCategory.is_2_calabi_yau = _is_2_calabi_yau


# ---------------------------------------------------------------------------------
# from_potential: the Jacobi-finite (Q, W) ALGEBRA-level certificate
# ---------------------------------------------------------------------------------
def _is_single_oriented_cycle(Q, length):
    """True iff ``Q`` is exactly one oriented cycle through ``length`` vertices."""
    verts = list(Q.vertices)
    if len(verts) != length or len(Q.arrows) != length:
        return False
    outdeg = {v: 0 for v in verts}
    indeg = {v: 0 for v in verts}
    for a in Q.arrows:
        s, t = Q.source(a), Q.target(a)
        if s == t:
            return False
        outdeg[s] += 1
        indeg[t] += 1
    return all(outdeg[v] == 1 and indeg[v] == 1 for v in verts)


def _reducible_hereditary_model(Q, JA):
    """The NARROW, DOCUMENTED reducibility recognition.

    ``(3-cycle, alpha beta gamma)`` has Jacobian algebra ``kZ_3/J^2`` (dim 6), which is
    the cluster-tilted algebra of ``A_3`` -- so its cluster category IS ``C_{A_3}`` and
    the category combinatorics may be served.  Recognition is per-instance CERTIFIED (the
    shape of ``Q`` and ``dim Jac == 6``), never inferred.

    A GENERAL mutation-equivalence-to-acyclic test is deliberately NOT attempted: deciding
    it needs quiver-mutation-class enumeration, and guessing would let the category
    surface answer for an algebra whose cluster category is not hereditary at all.
    Everything else returns ``None`` and the category invariants are refused (G7).
    """
    if _is_single_oriented_cycle(Q, 3) and JA.dim == 6:
        return ClusterCategory("A3", field=JA.domain)
    return None


@classmethod
def _from_potential(cls, Q, W, field=None):
    """A Jacobi-finite quiver-with-potential: the ALGEBRA-level Amiot certificate.

    Builds ``Jac(Q, W)``; a Jacobi-INFINITE input propagates the shipped
    ``NotFiniteDimensionalError`` (G6).  Certifies what Amiot / Keller-Reiten give at the
    algebra level -- finite-dimensional, and Gorenstein of dimension at most one.

    **The CATEGORY combinatorics are served only when ``(Q, W)`` is certified reducible to
    a hereditary model** (the documented 3-cycle instance).  Otherwise ``C_{(Q,W)}`` is
    not a hereditary cluster category, there is no finite mod-A model, and every category
    invariant refuses loudly (G7) -- the Jacobian algebra's own module theory is still
    served by the rest of quiverlab.
    """
    from quiverlab.families.jacobian import JacobianAlgebra
    from quiverlab.fields import QQ
    self = object.__new__(cls)
    JA = JacobianAlgebra(Q, W, field=field or QQ)          # G6 propagates from here
    self.jacobian = JA
    self._potential = (Q, W)
    self._knit = None
    self._eg = None
    model = _reducible_hereditary_model(Q, JA)
    self._model = model
    if model is None:
        self.algebra = None
        self.n = len(list(Q.vertices))
        self.hereditary = False
        self.dynkin_type = None
    else:
        self.algebra = model.algebra
        self.n = model.n
        self.hereditary = model.hereditary
        self.dynkin_type = model.dynkin_type
    return self


def _jacobian_certificate(self):
    """The algebra-level Amiot / Keller-Reiten certificate for a ``from_potential`` input:
    finite-dimensional, and Gorenstein of dimension at most one."""
    JA = getattr(self, "jacobian", None)
    if JA is None:
        raise QuiverlabError(
            "jacobian_certificate is only defined for a from_potential input",
            hint="this ClusterCategory was built from an acyclic quiver, not from (Q, W)")
    g = JA.gorenstein_dimension()
    reducible = self._model is not None
    return {
        "dim": JA.dim,
        "finite_dimensional": True,                # a non-finite Jac raised at build time
        "gorenstein": bool(g.is_gorenstein),
        "gorenstein_dimension": max(g.left_id, g.right_id) if g.is_gorenstein else None,
        "gorenstein_at_most_one": bool(g.is_gorenstein
                                       and max(g.left_id, g.right_id) <= 1),
        "self_injective": _self_injective(JA),
        "reducible_to_hereditary": reducible,
        "model": (None if not reducible else _diagram_name(self._model.dynkin_type)),
        "note": ("Jac(Q,W) is finite-dimensional and Gorenstein of dimension <= 1 "
                 "(Keller-Reiten); by Amiot, when (Q,W) is Jacobi-finite the generalized "
                 "cluster category C_{(Q,W)} carries a cluster-tilting object whose "
                 "endomorphism algebra IS this Jacobian algebra -- CITED, not computed. "
                 + ("This (Q,W) is certified reducible to a hereditary model, so the "
                    "category invariants are served through it."
                    if reducible else
                    "This (Q,W) is NOT certified reducible to a hereditary model, so the "
                    "CATEGORY invariants need D^b(Gamma) and are refused (out of scope)")),
    }


def _self_injective(A):
    from quiverlab.modules.ext import is_selfinjective
    return bool(is_selfinjective(A))


def _require_category_model(self, what):
    if getattr(self, "algebra", None) is None:
        raise QuiverlabError(
            f"{what}: C_(Q,W) category invariants need D^b(Gamma) -- the generalized "
            "cluster category of a quiver with potential that is not certified reducible "
            "to a hereditary model is out of scope; only the Jacobian-algebra Amiot "
            "certificate is available",
            hint="call jacobian_certificate() for the algebra-level statement")


ClusterCategory.from_potential = _from_potential
ClusterCategory.jacobian_certificate = _jacobian_certificate
ClusterCategory._require_category_model = _require_category_model
