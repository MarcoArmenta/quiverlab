# Plan 56: π₁(Q,I) + the Strongly-Simply-Connected Recognizer (R14 + R16) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The combinatorial fundamental group of a presentation, exact and honest.
Two public surfaces on `Algebra`: (1) **π₁(Q,I)** — a finitely presented group
(generators = non-tree arrows, relators from the relation-homotopy along a spanning
tree) with its **abelianization π₁^{ab}** computed by exact Smith normal form (the
always-computable headline), the `Hom(π₁,k⁺) ↪ HH¹` Hurewicz cross-check, and a
loud typed refusal for the *intrinsic* π₁; (2) **`is_simply_connected`** — a
three-valued verdict (True only via decidable sufficient criteria; False via a
decidable witness; else `None`/inconclusive, honest per Adian–Rabin) that carries
the **strongly-simply-connected certificate** (R16 separation condition for every
full convex subcategory, witness on failure). The strongly-simply-connected verdict
is a **clean consumable certificate that GATES P62** (Tits-form tame/wild). Two
no-code compute kinds — `fundamental_group` and `simply_connected` — expose the
whole surface end-to-end.

**Architecture:** A new `src/quiverlab/invariants/coverings.py` holds the group
construction (spanning-forest walk group → homotopy quotient → SNF abelianization →
Tietze-lite triviality search), the bypass/double-bypass detection, the three-valued
`is_simply_connected`, and the separation-condition / `is_strongly_simply_connected`
recognizer. Pure-combinatorics graph primitives (spanning forest, predecessors,
induced subquiver, undirected components) are added to
`src/quiverlab/combinat/quiver.py` as `Quiver` methods (they must not import
`invariants`, to avoid the `combinat → invariants` cycle; the existing `is_connected`
delegate at `combinat/quiver.py:75` and the `_undirected` helper at
`invariants/dynkin_type.py:19` are the shape to mirror, but the new primitives live
in `combinat`). Every result is a frozen dataclass following the `GlobalDimension`
honesty pattern (`modules/ext.py`): a certified value, or a labelled marker
(`inconclusive`/`refused`), never a fabricated verdict. Thin lazy-import wrapper
methods on `Algebra` sit in the `# -- invariants` block (`core/algebra.py:566-618`,
the `dynkin_type` wrapper at `600-611` is the template for the presentation guard).

Two mathematical engines are reused, not rebuilt: the **exact ℤ SNF** via the same
sympy helper `invariant_factors` (`sympy.matrices.normalforms`, imported at
`derived/fingerprint.py:17` and used at `:38` *without* a domain), here called with
an explicit `domain=ZZ` (verified valid: `invariant_factors(Matrix([[-1]]),
domain=ZZ) == [1]`); and the **Krull–Schmidt `decompose`** (`modules/decompose.py:361`,
returning `(module, multiplicity)` pairs) applied to `rad P_x` (`P.radical()`) for
the separation supports. The **minimal relations** that generate the homotopy
relation are computed by **exact linear algebra on the finite-dimensional space
I/(rad·I + I·rad)** (Task 2, H1) — NOT by any greedy reduction of the stored
generators (which is unsound). `combinat/relations.py::Relation` already stores each
relation *uniform* (its `parse` enforces a common source/target across terms,
`relations.py:187-193`), and `resolutions_cs.build.reduction_system_of(A)` — a
Gröbner *completion* — is used ONLY for the normal-form/ideal-membership map (a path
combination lies in I iff it reduces to 0), never to "minimize" (a completion is the
opposite of a minimal generating set).

**Tech Stack:** sympy exact integer linear algebra (SNF via `invariant_factors(...,
domain=ZZ)`, rank) — already a dep; `quiverlab.fields.linalg` for the exact
field-linear-algebra kernels (the path→A map's kernel = `I ∩ V_{x,y}`, and the
`rad·I + I·rad` subspace); the reduction-system normal form
(`ReductionSystem.normal_form`, `groebner/system.py:35`) as the path→A /
ideal-membership map only. No floats in `src/` (the group theory is over ℤ;
`Hom(π₁,k⁺)` counts free rank + char-divisible torsion — an integer/char
computation).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- **Bucket decision (stated, not guessed):** new tests live in a **new
  `tests/coverings/` directory**, which `tests/conftest.py::_bucket` (`_DEEP_DIRS`
  at line 50) does **not** list ⇒ it collects in the **fast** bucket. Keep every
  covering test fast: the zoo `Hom(π₁,k⁺) ≤ HH¹` sweep is bounded to the
  **triangular** zoo members up to a small `dim_max` so it fits the fast budget;
  the separation batteries use hand-sized oracles (kA_n, squares, the Zito
  example). QPA honest-skips go in `tests/qpa/` (bucket = the class); cross-runner
  webapp/gui tests go in `tests/webapp/` (unmarked — extras-gated dir, Plan-32 rule).
- New scalar kinds act on the **existing algebra block** ⇒ **schema stays v1**
  (Plan-26/38 precedent: request-schema bumps only for new request *blocks*).
- **Canonical keys are pinned to the frozen `_V` constant by design** (metaplan §3):
  adding `fundamental_group` / `simply_connected` to the compute list must not
  re-key any existing request. Verify by running `test_runner_delegation.py`
  BEFORE and AFTER adding goldens; existing golden entries stay byte-identical.
- **Every citation key is BibTeX-verified before use** (Plan-29 rule); `_r(...)`
  registry entry + `references.bib` entry both land, gated by
  `tests/citations/test_bib_structure.py`.
- All refusals are loud (`QuiverlabError(message, hint=...)`, the two-space
  `[hint: ...]` convention, `errors.py:4-7`); presentation-less (structure-constant)
  algebras refuse every kind here — the quiver is required — mirroring the
  `dynkin_type` guard (`core/algebra.py:604-609`).
- **Merge-train discipline (metaplan §3, binding):** stage conflict resolutions by
  EXPLICIT path (never `git add -A`); merge `_runner_goldens.json` SEMANTICALLY
  (dict-union), never textually; grep for conflict markers after every merge;
  compile-check every touched `.py` before staging; run
  `tests/citations/test_bib_structure.py` + `tests/webapp/test_js_parses.py` +
  `test_runner_delegation.py` after every merge.
- Conventional commits; green at every commit; branch `plan-56-pi1` off `dev`.
  P51/P52/P53 run in locked worktrees (`plan-51-bracket`, `plan-52-hh-coefficients`,
  `plan-53-invariants`); P56 is **greenfield** (no `pi1`/`fundamental_group`/
  `simply_connected`/`separation` code exists in `src/`). Shared-file overlap with
  those plans is confined to the GUI runner / i18n / golden / citation files —
  the metaplan merge-train hazard; follow the discipline above.

---

## Records (verbatim — the mandate for this plan)

> **R14 — Fundamental group π₁(Q,I) of a presentation + Hurewicz to HH¹.**
> [merged B-P5 + C-P1; adjudicated corrections applied] Object: the COMBINATORIAL
> presentation group π₁(Q,I) (walks mod relation-homotopy along a spanning tree) —
> finite presentation always computable; abelianization = universal abelian grading
> group via exact Smith normal form; Hom(π₁, k⁺) ↪ HH¹ cross-check. HONEST
> DECIDABILITY (adjudicated): triviality of an f.p. group is undecidable
> (Adian–Rabin) — emit π₁^{ab} always, 'simply connected' only via decidable
> sufficient criteria (separation condition R16; Le Meur's privileged presentation:
> char 0 + no double bypasses). The INTRINSIC π₁ (inverse limit over connected
> gradings) is NOT bounded-computable in general — refuse loudly; the char-p oracle
> π₁(k[x]/(x^p)) = ℤ × C_p (verified verbatim in arXiv:0906.3069) is the INTRINSIC
> group, not the presentation group (a monomial loop relation gives presentation
> π₁ = ℤ) — use it only if the grading route is implemented. CORRECTED ORACLE:
> commutative square WITH the commutativity relation ⇒ π₁ = 1 (the relation
> identifies the parallel paths); ℤ only WITHOUT it. Refs: Assem–de la Peña Comm.
> Alg. 24 (1996) 187–208; Cibils–Redondo–Solotar arXiv:0706.2491, arXiv:1010.6296,
> arXiv:0906.3069; Le Meur arXiv:math/0503302; Briggs–Rubio y Degrassi
> arXiv:2109.03704 (IMRN 2023: every maximal torus of HH¹ is dual to some π₁ —
> verified). Size M.

> **R16 — Strongly-simply-connected recognizer (separation condition).**
> [C-scout P9; keep] Object: for triangular A: every convex subcategory satisfies
> separation (rad P_x decomposes with supports in distinct components of the
> non-predecessor subquiver) — pure finite combinatorics + Krull–Schmidt; witness
> on failure. Gate for R19. Refs: Skowroński CMS Conf. Proc. 14 (1993); ASS book;
> arXiv:1905.06028. Size M.

---

## Reference re-verification (mandatory; findings recorded)

All statements below were re-read from the sources during authoring (WebFetch of
the arXiv PDFs; the definitions are quoted where load-bearing). Keep this section —
the implementer transcription-checks the two `# PIN`s against primary sources.

- **Le Meur, arXiv:math/0503302 — VERIFIED VERBATIM (PDF §2, Thm 1.1).** The
  homotopy relation and minimal relation are defined exactly as we implement them:
  *"Let r = t₁u₁+…+t_nu_n ∈ I where t_i ∈ k\* and the u_i's are distinct paths;
  then r is called a **minimal relation** if n ≥ 2 and if for any non-empty proper
  subset E of {1,…,n} the term Σ_{i∈E} t_i·u_i does not lie in I. The **homotopy
  relation ~_I** is the smallest equivalence relation on walks compatible with
  concatenation such that (i) for any arrow α: x→y, αα⁻¹ ~_I e_y and α⁻¹α ~_I e_x;
  (ii) u₁ ~_I u₂ as soon as t₁u₁+…+t_nu_n is a minimal relation."* π₁(Q,I,x₀) =
  ~_I-classes of closed walks at x₀. **Example 1 (VERIFIED):** I = ⟨da⟩ gives
  π₁ ≅ ℤ; J = ⟨da − dcb⟩ gives π₁ = 0 — a **monomial single-path relation creates
  NO identification** (n = 1, not minimal), a binomial does. **Bypass** = a pair
  (α, u), α an arrow and u an oriented path parallel to α, u ≠ α. **Double bypass**
  = a quadruple (α, u, β, v) where (α,u),(β,v) are bypasses and β appears in u.
  **Theorem 1.1 (VERIFIED):** if char k = 0, A triangular, Q has no double
  bypasses, there is a presentation kQ/I₀ ≅ A such that for any other kQ/J ≅ A
  there is a *surjective* group morphism π₁(Q,I₀) → π₁(Q,J). *"Our assumption on
  double bypasses implies Q has no double arrows."* **Remark 6 (VERIFIED):** if
  p = char k ≠ 0 and Q has fewer than p bypasses, Theorem 1.1 still holds.
  **Prop 4.2/Remark 4 (VERIFIED):** if Q has m bypasses the comparison graph Γ has
  a unique source and ≤ 1+m+…+mᵐ vertices — in particular **m = 0 ⇒ π₁ is
  presentation-independent** (one homotopy relation).
- **Cibils–Redondo–Solotar, arXiv:1010.6296 (Hurewicz) — VERIFIED VERBATIM
  (abstract + §1).** *"The Hurewicz morphism from the vector space of abelian
  characters of π₁(C) to the first Hochschild–Mitchell cohomology vector space of C
  is an isomorphism"* for a Schurian category (hom-spaces dim ≤ 1). For matrix
  algebras M_p(k), p prime, char 0: π₁ = F_{p-1} × C_p (free × cyclic) — the sibling
  of the k[x]/(x^p) oracle.
- **Cibils–Redondo–Solotar, arXiv:0906.3069 "Connected gradings and the fundamental
  group" (Algebra Number Theory 4 (2010) 625–648) — # PIN (adjudicated; not
  byte-extractable from the PDF).** The INTRINSIC group π₁(k[x]/(x^p)) = ℤ × C_p in
  char p. Corroborated by the M_p(k) = F_{p-1} × C_p sibling above. We DO NOT compute
  this (it is the intrinsic/grading group); it is recorded only as documentation of
  *why* the presentation group (ℤ, which we compute for the monomial loop) is not the
  algebra invariant. Transcription-check at implementation; if the PDF cannot be
  read, cite the record as adjudicated and keep the k[x]/(x^p) presentation-π₁ = ℤ
  test (that one IS ours to verify).
- **Briggs–Rubio y Degrassi, arXiv:2109.03704 (IMRN 2023) — VERIFIED (abstract).**
  *"Every maximal torus in HH¹(A) arises as the dual of some fundamental group of
  A."* Documentation of the maximal-torus refinement of the Hurewicz cross-check;
  not a v1 computable.
- **Kasjan–Skowroński, arXiv:1905.06028 "On tame strongly simply connected algebras"
  — VERIFIED (PDF Thm 1 + §2).** For strongly simply connected A: *A is tame ⟺ the
  Tits form q_A is weakly nonnegative* (Brüstle–de la Peña–Skowroński). Tits form
  q_A(z) = Σ_i z_i² − Σ_{i→j} z_iz_j + Σ_{i,j} r_{ij} z_iz_j (r_{ij} = number of
  minimal relations from i to j). This is the theorem P62 consumes downstream of our
  strongly-simply-connected certificate.
- **Zito, arXiv:2006.15729 "Several Results Concerning Convex Subcategories" —
  VERIFIED VERBATIM (PDF §2, §3.2).** *Full* subquiver: for any two vertices v,w in
  C, all arrows of Q from v to w are in C. *Convex*: for any v,w ∈ C and any path p
  from v to w, every vertex of p is in C. *"An algebra Λ is **simply connected**
  provided its ordinary quiver Q has no oriented cycles and, for any presentation
  Λ ≅ KQ/I, the fundamental group π₁(Q,I) is trivial."* *"Λ is **strongly simply
  connected** if every full, convex subcategory is simply connected."* **The
  discriminating example (VERIFIED):** Q = (1→2, 2→3→5, 2→4→5), I = ⟨αβγ − αδε⟩,
  gives π₁(Q,I) = 0 for every presentation (so Λ is simply connected), but the full
  convex subcategory on {2,3,4,5} is hereditary with π₁ = ℤ (NOT simply connected) —
  hence Λ is simply connected but NOT strongly simply connected. This is our key
  strongly-s.c.-vs-s.c. oracle.
- **Separation condition — # PIN (secondary sources; transcription-check against
  ASS2006 + Skowroński 1993 at implementation).** An indecomposable projective P_a
  has a **separated radical** if the pairwise non-isomorphic indecomposable direct
  summands of rad P_a have supports in *distinct connected components* of the full
  subquiver **Q_a = the full subquiver on the vertices from which a is UNREACHABLE**,
  i.e. delete `a` together with its entire **transitive predecessor closure**
  `{ x : there is a directed path x ⇝ a }` (which includes `a` itself, via the trivial
  path). A satisfies the **separation condition** if every P_a does. **Skowroński
  (CMS 14, 1993):** a triangular algebra is strongly simply connected ⟺ every full
  convex subcategory satisfies the separation condition. This is R16.
  **Why the closure (and `a`) MUST be removed — hand-derivation on the star**
  `1→2, 1→3` (a tree, hence strongly simply connected, so it MUST come out
  separated): the WRONG convention (delete only immediate proper predecessors, keep
  `a`) leaves `Q_1 = {1,2,3}` connected, so rad P₁ = S₂ ⊕ S₃ has supports {2},{3} in
  ONE component (joined through vertex 1) → "separation fails" → trees mislabelled
  not strongly simply connected → the P62 gate corrupted. Under the CORRECT
  convention `Q_1 = {2,3}` (delete `1` and its closure `{1}`) has TWO components
  {2},{3} → separated ✓. rad P_a is always supported on proper successors of `a`,
  which survive in Q_a (a successor cannot also be a predecessor — triangular).

---

### Task 0: citations (do first — later tasks reference the keys)

**Files:**
- Modify: `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`
- Test: `tests/citations/test_bib_structure.py` (existing gate)

**Interfaces:** the `_r(key, bibtex_key, kind, title, annotation, *tags)` registry
helper (`registry.py:24`) + the `@article`/`@incollection` bib patterns
(`references.bib`). `assem_book` (→ `ASS2006`) and `skowronski_yamagata` already
exist; do NOT re-add. `_citation_pairs` (`hpc/spec.py:1557`) resolves keys →
`[key, formatted]` pairs for the payloads.

- [ ] **Step 1: add the keys** (BibTeX-verify each before committing — Plan-29 rule):
  - `assem_delapena` → Assem, de la Peña, *The fundamental groups of a triangular
    algebra*, Comm. Algebra **24** (1996) 187–208, doi 10.1080/00927879608825561.
    (foundation; π₁ of a presentation + the HH¹ link)
  - `martinez_villa_delapena` → Martínez-Villa, de la Peña, *The universal cover of
    a quiver with relations*, J. Pure Appl. Algebra **30** (1983) 277–292.
    (foundation; the ORIGINAL homotopy relation ~_I)
  - `le_meur_pi1` → Le Meur, *The fundamental group of a triangular algebra without
    double bypasses*, C. R. Acad. Sci. Paris Sér. I **341** (2005), arXiv:math/0503302.
    (foundation; Thm 1.1 privileged presentation) — # verify exact pages.
  - `crs_hurewicz` → Cibils, Redondo, Solotar, *Fundamental group of Schurian
    categories and the Hurewicz isomorphism*, arXiv:1010.6296. (foundation)
  - `crs_gradings` → Cibils, Redondo, Solotar, *Connected gradings and the
    fundamental group*, Algebra Number Theory **4** (2010) 625–648, arXiv:0906.3069.
    (foundation; intrinsic π₁, universal abelian grading, k[x]/(x^p) oracle)
  - `crs_intrinsic` → Cibils, Redondo, Solotar, *The intrinsic fundamental group of
    a linear category*, Algebr. Represent. Theory, arXiv:0706.2491. (foundation)
  - `briggs_ryd_tori` → Briggs, Rubio y Degrassi, *Maximal tori of the first
    Hochschild cohomology and fundamental groups*, IMRN (2023), arXiv:2109.03704.
    (foundation) — # verify exact title.
  - `skowronski_ssc` → Skowroński, *Simply connected algebras and Hochschild
    cohomologies*, in Representations of Algebras (Ottawa, 1992), CMS Conf. Proc.
    **14**, Amer. Math. Soc. (1993) 431–447. (foundation; separation condition,
    strongly simply connected) — use the `@incollection` CMS-Conf-Proc pattern
    already present at `references.bib` (the `cibils1998radsq` shape).
- [ ] **Step 2:** `... -m pytest tests/citations/test_bib_structure.py -q` green
  (registry ↔ bib in sync; every `_r` key resolvable via `bibtex(key)`).
- [ ] **Step 3: Commit** — `docs(citations): π₁ / simply-connected references (Assem–de la Peña, Le Meur, CRS ×3, Briggs–RyD, Skowroński 1993)`

---

### Task 1: quiver graph primitives (combinat, no invariants import)

**Files:**
- Modify: `src/quiverlab/combinat/quiver.py`
- Test: `tests/coverings/test_quiver_graph.py`

**Interfaces:** add five pure-combinatorics `Quiver` methods (no `invariants`
import — they must be usable from `combinat`):
```python
def predecessors(self, v) -> set          # {s for (s,t) in arrows.values() if t==v} (proper: excludes v)
def successors(self, v) -> set            # {t for (s,t) in arrows.values() if s==v}
def undirected_components(self, vertices=None) -> list[frozenset]
                                          # connected components of the underlying
                                          # undirected graph induced on `vertices`
                                          # (default all); loops ignored for adjacency
def spanning_forest(self, root=None) -> tuple[set[str], dict]
                                          # (tree_arrow_names, parent_map). One
                                          # spanning tree per component (BFS). A
                                          # tree arrow is chosen once per newly-
                                          # reached vertex; parallel/back arrows are
                                          # non-tree. Deterministic (sorted order).
def induced_subquiver(self, vertices) -> "Quiver"
                                          # full subquiver on `vertices` (arrows with
                                          # both endpoints in the set). Reuses Quiver.__init__.
```
Mirror `_undirected` (`invariants/dynkin_type.py:19`) for the adjacency build but
keep it inside `combinat`. `spanning_forest` returns the SET of tree-arrow names
(a chosen tree per component) so π₁ generators = `set(arrows) - tree_arrows`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/coverings/test_quiver_graph.py
"""Pure-combinatorics graph primitives underpinning pi1 (spanning forest,
predecessors, induced subquiver, undirected components). Self-certifying."""
import pytest
from quiverlab import Quiver

pytestmark = pytest.mark.oracle_selfcert


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_predecessors_successors():
    Q = _square()
    assert Q.predecessors(4) == {2, 3}
    assert Q.successors(1) == {2, 3}
    assert Q.predecessors(1) == set()


def test_spanning_forest_betti_number():
    Q = _square()
    tree, _ = Q.spanning_forest()
    # connected: |tree| = |V| - 1; non-tree count = first Betti number = |E|-|V|+c
    assert len(tree) == 3
    non_tree = set(Q.arrows) - tree
    assert len(non_tree) == 4 - 4 + 1                  # Betti = 1


def test_undirected_components_disconnected():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)})
    comps = Q.undirected_components()
    assert sorted(map(sorted, comps)) == [[1, 2], [3, 4]]


def test_induced_subquiver_drops_boundary_arrows():
    Q = Quiver([1, 2, 3, 4, 5],
               {"al": (1, 2), "be": (2, 3), "ga": (3, 5), "de": (2, 4), "ep": (4, 5)})
    sub = Q.induced_subquiver({2, 3, 4, 5})            # drop vertex 1 and arrow "al"
    assert set(sub.vertices) == {2, 3, 4, 5}
    assert set(sub.arrows) == {"be", "ga", "de", "ep"}


def test_spanning_forest_one_tree_per_component():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)})
    tree, _ = Q.spanning_forest()
    assert len(tree) == 2                              # one edge per component, c=2
```

- [ ] **Step 2: Run to verify failure** (`AttributeError: predecessors`).
- [ ] **Step 3: Implement** the five methods (all BFS/set comprehensions over
  `self.arrows`; `spanning_forest` iterates components in sorted vertex order,
  sorted arrow order, so the tree is deterministic — important for reproducible
  π₁ presentations and byte-stable payloads).
- [ ] **Step 4: Run** — `... -m pytest tests/coverings/test_quiver_graph.py -q` green.
- [ ] **Step 5: Commit** — `feat(combinat): quiver graph primitives (spanning forest, predecessors, induced subquiver, components) for pi1`

---

### Task 2: π₁(Q,I) presentation + π₁^{ab} via exact SNF (the headline)

**Files:**
- Create: `src/quiverlab/invariants/coverings.py`
- Modify: `src/quiverlab/invariants/__init__.py`, `src/quiverlab/core/algebra.py`
  (`Algebra.fundamental_group`, `Algebra.intrinsic_fundamental_group`)
- Test: `tests/coverings/test_pi1.py`

**Interfaces:**
```python
@dataclass(frozen=True)
class FundamentalGroup:
    generators: tuple          # non-tree arrow labels (base-point loop names)
    relators: tuple            # words in generators (strings), from minimal relations
    free_rank: int             # rank of the free part of pi1^ab (= |generators| - rank(R))
    invariant_factors: tuple   # the SNF factors d_i > 1 (torsion of pi1^ab)
    components: int            # connected components of Q (pi1^ab = direct sum over them)
    base_vertex: object
    presentation_note: str     # "presentation invariant of (Q, I); NOT of the algebra A"
    def abelianization_repr(self) -> str      # "Z^r (+) Z/d1 (+) ... " or "0"
    def hom_to_additive_dim(self, char: int) -> int
        # dim_k Hom(pi1^ab, (k,+)) = free_rank + #{d_i : char>0 and char | d_i}

def fundamental_group(A, base=None) -> FundamentalGroup
def intrinsic_fundamental_group(A)              # ALWAYS raises QuiverlabError (see below)
```

**The construction (exact linear algebra; ALWAYS terminates, NEVER refuses):**
1. Guard: `A.quiver is None` → loud refusal (presentation-less). Admissibility is
   guaranteed by the build routes (`AdmissibilityError` at build time); no re-check.
   There is **no** minimality-refusal path — the algorithm below always succeeds.
2. **Minimal relations = a basis of I/(rad·I + I·rad) (exact linear algebra).** ~_I is
   generated by the *minimal relations* of I, which are the minimal generators of I as
   a two-sided ideal — a basis of the finite-dimensional space `I/(rad·I + I·rad)`
   (Martínez-Villa–de la Peña; Bardzell–Marcos). **This is NOT a greedy reduction of
   the stored generators** — a stored generator can be irredundant yet not a minimal
   relation. *(Concrete unsoundness the old design hit: with `g₁ = u₁−u₂` and
   `g₂ = (u₁−u₂)+(u₃−u₄)` both stored, neither is greedily redundant, but g₂ has the
   proper sub-relation g₁ ∈ I; a naive "use g₂'s paths" glues u₁~u₂~u₃~u₄ — WRONG,
   the truth is u₁~u₂ and u₃~u₄ separately — and a "refuse non-minimal generators"
   rule would reject a perfectly computable input.)*
   Since I is a two-sided ideal, `I = ⊕_{(x,y)} (I ∩ e_y·kQ·e_x)` decomposes over
   source/target blocks; compute block-by-block. Let `N = A.loewy_length()` so
   `rad(kQ)^N ⊆ I` (exact truncation; `core/algebra.py:729`). For each ordered pair
   (x,y):
   - `V_{x,y}` = all paths x→y of length in `[2, N]` (finite BFS; I ⊆ rad² excludes
     length-1 arrows; length > N paths lie in rad^N ⊆ I automatically).
   - `W_{x,y} = I ∩ V_{x,y}` = **kernel of the path→A map**: the matrix with rows
     `normal_form(p)` (coords in A's basis of `e_y A e_x`) for `p ∈ V_{x,y}`
     (`reduction_system_of(A).normal_form`); a path-combination is in I iff it maps to
     0 in A. Kernel exact over k (`fields.linalg`).
   - `R_{x,y} = (rad·I + I·rad) ∩ V_{x,y}` = span of
     `{ α∗w : α arrow x→x', w ∈ W_{x',y} } ∪ { w∗α : α arrow y'→y, w ∈ W_{x,y'} }`
     (each extends a shorter relation by one arrow; drop length > N). `R_{x,y} ⊆ W_{x,y}`.
   - **minimal relations of (x,y)** = a complement basis of `W_{x,y}/R_{x,y}` (extend
     a basis of R to a basis of W by reduced row echelon; the added, reduced-support
     rows are the minimal relations — each a k-combination of **parallel** paths x→y,
     automatically uniform, every path nonzero in A).
   - **Self-cert:** total minimal-relation count `== dim I/(rad·I + I·rad)
     == Σ_{x,y}(dim W_{x,y} − dim R_{x,y})` (gated by a test).
   A single-path (monomial, support size 1) minimal relation glues nothing (Le Meur
   Example 1). *Cite `martinez_villa_delapena` + `assem_delapena`.*
3. **Spanning forest** `tree, _ = A.quiver.spanning_forest(root=base)`; generators =
   `sorted(set(arrows) - tree)`.
4. **π₁^{ab} (SNF).** The abelianized cycle space of the underlying graph is
   ℤ^{generators} (rank = first Betti number = |Q₁| − |Q₀| + c, handling
   disconnected automatically). For each minimal relation with parallel paths
   u₀,…,u_{n-1} (support ≥ 2), each pair contributes the integer row
   `vec(u₀) − vec(u_i)`, where `vec(u)` ∈ ℤ^{generators} counts each non-tree arrow in
   u (tree arrows drop out; tree-path corrections cancel since u₀,u_i share endpoints).
   Stack rows → matrix R. Factors = `invariant_factors(sp.Matrix(R), domain=ZZ)` (the
   same helper `fingerprint.py:38` uses, here with explicit `domain=ZZ`);
   `free_rank = len(generators) − rank(R)` (`from quiverlab.fields.linalg import rank`);
   `invariant_factors` field = the factors > 1. Empty R (monomial or hereditary) ⇒
   `free_rank = |generators|`, no torsion. The abelianization is a function of ~_I
   alone (determined by I — MVdlP), hence **basis-independent** (the complement-basis
   choice in Step 2 does not affect it).
5. **Presentation (relators as words).** For each minimal relation and each pair,
   emit `word(u₀)·word(u_i)⁻¹` as a reduced string in the generators (tree arrows →
   ε via the parent map). This is the finite presentation payload.

`intrinsic_fundamental_group(A)` **ALWAYS** raises:
```python
raise QuiverlabError(
    "the intrinsic fundamental group (inverse limit over connected gradings) is "
    "not bounded-computable in general; quiverlab emits the PRESENTATION group "
    "pi1(Q, I) and its abelianization instead",
    hint="e.g. the presentation pi1(k[x]/(x^p)) = Z (monomial loop) differs from "
         "the intrinsic pi1 = Z x C_p in characteristic p (Cibils-Redondo-Solotar, "
         "arXiv:0906.3069) -- the C_p torsion is invisible to any presentation")
```

`Algebra.fundamental_group(self, base=None)` / `Algebra.intrinsic_fundamental_group(self)`
wrappers in `core/algebra.py:566-618` with the `dynkin_type`-style quiver guard;
re-export `fundamental_group`, `intrinsic_fundamental_group`, `FundamentalGroup`
from `invariants/__init__.py`.

- [ ] **Step 1: Write the failing tests** (the CORRECTED oracle is the centrepiece)

```python
# tests/coverings/test_pi1.py
"""pi1(Q,I): presentation + abelianization by exact SNF. Oracles: the CORRECTED
commutative square (WITH relation => 1, WITHOUT => Z), Le Meur Example 1, trees,
single loop, multi-loop free ranks, the annulus. Refs: Le Meur math/0503302;
Martinez-Villa-de la Pena 1983; Assem-de la Pena 1996."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.coverings import fundamental_group, intrinsic_fundamental_group

pytestmark = pytest.mark.oracle_literature


def _square(rel):
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=rel, field=GF(7))


def test_commutative_square_with_relation_is_trivial():
    A = _square(["a*b - c*d"])                 # commutativity: identifies the parallel paths
    g = fundamental_group(A)
    assert g.free_rank == 0 and g.invariant_factors == ()   # pi1^ab = 0
    assert g.abelianization_repr() == "0"


def test_commutative_square_without_relation_is_Z():
    A = _square([])                            # hereditary square: a hole survives
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()   # pi1^ab = Z
    assert g.abelianization_repr() == "Z"


def test_le_meur_example_1_monomial_vs_binomial():
    # Le Meur Ex.1 concretely + ADMISSIBLY (all relations in rad^2). d*a (1->2->3) and
    # d*c*b (1->2->4->3) are parallel 1->3 paths. I = <d*a> monomial (n=1) => NO
    # identification => pi1^ab = Z; J = <d*a - d*c*b> binomial minimal => identifies
    # the parallel paths => pi1^ab = 0.
    Q = Quiver([1, 2, 3, 4], {"d": (1, 2), "a": (2, 3), "c": (2, 4), "b": (4, 3)})
    Amono = Q.algebra(relations=["d*a"], field=GF(7))
    Abino = Q.algebra(relations=["d*a - d*c*b"], field=GF(7))
    assert fundamental_group(Amono).free_rank == 1        # I = <da>     -> Z
    assert fundamental_group(Abino).free_rank == 0        # J = <da-dcb> -> 0


def test_two_independent_pairs_not_over_glued():
    # DISCRIMINATING ORACLE (adversarial H1). Four parallel routes 1->3 and a
    # NON-minimal generating set I = <u1-u2, (u1-u2)+(u3-u4)>. The exact
    # I/(rad.I + I.rad) linear algebra recovers the TWO minimal relations u1-u2 and
    # u3-u4 => glue u1~u2 and u3~u4 ONLY => pi1^ab = Z (one route-class comparison
    # survives). The OLD greedy/minimality-refusal design would REFUSE g2 (proper
    # sub-relation g1 in I) or naively glue all four routes => pi1^ab = 0 (WRONG).
    Q = Quiver([1, 3, "A", "B", "C", "D"],
               {"a1": (1, "A"), "a2": ("A", 3), "b1": (1, "B"), "b2": ("B", 3),
                "c1": (1, "C"), "c2": ("C", 3), "d1": (1, "D"), "d2": ("D", 3)})
    A = Q.algebra(relations=["a1*a2 - b1*b2", "a1*a2 - b1*b2 + c1*c2 - d1*d2"],
                  field=GF(7))
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()     # Z, NOT 0


def test_tree_is_trivial():
    A = linear_path_algebra(5, field=GF(7))    # kA5: underlying graph a tree
    g = fundamental_group(A)
    assert g.free_rank == 0 and g.invariant_factors == () and g.components == 1


def test_single_loop_monomial_is_Z():
    A = truncated_polynomial(4, field=GF(7))   # k[x]/(x^4): one loop, monomial relation
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()   # presentation pi1 = Z (NOT ZxC_p)


def test_multi_loop_free_rank():
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y", "y*x"], field=GF(5))   # all monomial
    g = fundamental_group(A)
    assert g.free_rank == 2 and g.invariant_factors == ()      # free on 2 loops


def test_disconnected_is_direct_sum():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)}).algebra(field=GF(5))
    g = fundamental_group(A)
    assert g.components == 2 and g.free_rank == 0               # two trees


def test_intrinsic_refused_loudly():
    A = truncated_polynomial(3, field=GF(7))
    with pytest.raises(QuiverlabError, match="intrinsic"):
        intrinsic_fundamental_group(A)
```

```python
# tests/coverings/test_pi1_annulus.py
"""The annulus has one hole: pi1^ab = Z. VERIFIED LIVE during authoring:
jacobian_of(annulus_triangulation(2,2)) is 4 vertices, 4 arrows, EMPTY relations
(hereditary) -- so free_rank = 1 is the first Betti number of a 4-cycle (a DERIVED
fact, not a published value): oracle_selfcert."""
import pytest
from quiverlab import annulus_triangulation, jacobian_of
from quiverlab.invariants.coverings import fundamental_group

pytestmark = pytest.mark.oracle_selfcert


def test_annulus_C22_pi1ab_is_Z():
    A = jacobian_of(annulus_triangulation(2, 2))   # 4 vertices, 4 arrows, no relations
    assert A.relations == []                        # hereditary (re-verify the shape)
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()   # Betti(4-cycle) = 1 => Z
```

- [ ] **Step 2: Run to verify failure** (`ModuleNotFoundError: coverings`).
- [ ] **Step 3: Implement `coverings.py`** per the construction (Step 2 = the
  I/(rad·I + I·rad) block linear algebra; NO refusal path). The Le Meur Example-1 and
  two-pair quivers are concrete above. The annulus is confirmed hereditary
  (4 vertices, 4 arrows, `relations == []`) so `free_rank = 1` is Betti(4-cycle) — no
  binomial-torsion subtlety.
- [ ] **Step 4: Run** — `... -m pytest tests/coverings/test_pi1.py tests/coverings/test_pi1_annulus.py -q` green.
- [ ] **Step 5: Commit** — `feat(invariants): pi1(Q,I) presentation + abelianization via exact SNF; intrinsic pi1 refused loudly`

---

### Task 3: Hom(π₁,k⁺) ↪ HH¹ Hurewicz cross-check

**Files:**
- Modify: `src/quiverlab/invariants/coverings.py` (nothing new to expose;
  `hom_to_additive_dim` lives on `FundamentalGroup`)
- Test: `tests/coverings/test_hurewicz_bound.py`

**Interfaces:** `dim_k Hom(π₁^{ab}, (k,+)) = free_rank + #{d_i : char > 0 and
char | d_i}` (over char 0 only the free rank counts). **Assertion:** for a
**triangular** algebra, `Hom(π₁,k⁺) ↪ HH¹(A)` (Assem–de la Peña, `assem_delapena`;
iso for Schurian, `crs_hurewicz`), so `hom_to_additive_dim(char) ≤ dim HH¹(A)`.
`dim HH¹ = A.hochschild_cohomology(1)[1]` (`hochschild/table.py`, `HHTable.dims[1]`).

**Honest scope:** the embedding is a theorem for TRIANGULAR A only. The sweep runs
over the triangular members of the zoo; non-triangular members (loops, cyclic
Nakayama) are recorded but excluded from this particular assertion (their π₁^{ab}
is still emitted and pinned in Task 2).

- [ ] **Step 1: Write the failing tests**

```python
# tests/coverings/test_hurewicz_bound.py
"""Hom(pi1, k+) ↪ HH^1 for triangular algebras (Assem-de la Pena; CRS Hurewicz iso
for Schurian). Cross-engine: pi1 via SNF vs HH via the shipped engines."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra, zoo
from quiverlab.invariants.coverings import fundamental_group

pytestmark = pytest.mark.oracle_crossengine

P = 32003


def _dimhh1(A):
    return A.hochschild_cohomology(1)[1]


def test_square_equality_schurian():
    # commutative square (Schurian, triangular): Hom(pi1,k+) = dim HH^1 (Hurewicz iso)
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=["a*b - c*d"], field=GF(P))
    g = fundamental_group(A)
    assert g.hom_to_additive_dim(P) == 0
    assert _dimhh1(A) == 0                                 # square with relation: HH^1 = 0


def test_hereditary_square_Z_bound():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=[], field=GF(P))
    g = fundamental_group(A)
    assert g.hom_to_additive_dim(P) == 1 <= _dimhh1(A)     # Z: Hom-dim 1 <= dim HH^1


def test_bound_holds_on_triangular_zoo():
    for A in zoo(dim_max=9, field=GF(P)):
        if A.quiver is None or not A.quiver.is_acyclic():  # embedding is a theorem for triangular only
            continue
        g = fundamental_group(A)
        assert g.hom_to_additive_dim(P) <= _dimhh1(A)
```

- [ ] **Step 2: Run to verify failure**, then confirm PASS after Task 2 lands.
- [ ] **Step 3: Implement** — nothing beyond `hom_to_additive_dim` (Task 2). If any
  triangular zoo member violates `≤`, STOP: it is either a non-triangular record
  mis-tagged or an engine disagreement to escalate — do not weaken the assertion.
- [ ] **Step 4: Run** — `... -m pytest tests/coverings/test_hurewicz_bound.py -q` green.
- [ ] **Step 5: Commit** — `test(coverings): Hom(pi1,k+) ↪ HH^1 Hurewicz bound over the triangular zoo`

---

### Task 4: bypass detection + the three-valued `is_simply_connected`

**Files:**
- Modify: `src/quiverlab/invariants/coverings.py`, `src/quiverlab/core/algebra.py`
  (`Algebra.is_simply_connected`)
- Test: `tests/coverings/test_simply_connected.py`

**Interfaces:**
```python
def bypasses(A) -> list            # [(arrow, path_word), ...]: alpha and an oriented
                                   # path parallel to alpha, distinct (bounded path
                                   # enumeration -- triangular => finite; loud budget)
def has_double_bypass(A) -> bool   # a (alpha,u,beta,v) with beta appearing in u

@dataclass(frozen=True)
class SimpleConnectivity:
    verdict: bool | None           # True / False / None (inconclusive -- Adian-Rabin)
    reason: str                    # the criteria trail (which route fired / why None)
    witness: dict | None           # on False: {kind: "disconnected"|"oriented_cycle"|
                                   #   "nontrivial_pi1ab", ...}
    abelianization: FundamentalGroup
    strongly: "StrongSimpleConnectivity | None"   # the R16 certificate (Task 5);
                                   # None when not computed (cheap route decided and
                                   # strong not requested); GATES P62 when present

def is_simply_connected(A, strong="auto", convex_budget=20000) -> SimpleConnectivity
    # strong: "auto" (compute the R16 certificate only if R1/R3 do not decide) |
    #         True (always compute it -- the simply_connected kind + P62 gate) |
    #         False (cheap routes only -- NEVER raises/refuses on R1/R3 inputs)
```

**Decision logic (SOUND; True NEVER from a failed search):**
- Presentation-less → loud refusal (needs the quiver).
- **False guards (decidable witnesses, checked first):**
  - **Not connected** (`components > 1`) → False, witness `disconnected` (s.c.
    requires a connected quiver — Zito §3.2).
  - **Not triangular** (`not quiver.is_acyclic()`) → False, witness `oriented_cycle`
    (s.c. requires no oriented cycles — Zito §3.2).
  - **π₁(Q,I)^{ab} ≠ 0** (`free_rank > 0 or invariant_factors`) → False, witness
    `nontrivial_pi1ab` (the stored presentation has π₁ ≠ 1; s.c. requires *every*
    presentation trivial — one nontrivial presentation is a decidable False).
- **True routes, cheapest first (the first that fires wins — B2 reorder):**
  - **R1 tree (O(1)):** `len(generators) == 0` (underlying graph a tree, Betti 0) ⇒
    π₁(Q,J) = 1 for *every* presentation ⇒ simply connected. (Also strongly s.c.)
  - **R3 Le Meur no-bypass:** triangular and `len(bypasses(A)) == 0` (⇒ no double
    bypasses either, and Γ has a unique vertex, so π₁ is presentation-independent —
    Le Meur Prop 4.2/Rem 4). π₁^{ab} = 0 was already checked, so run **Tietze-lite**:
    reduces to the empty group ⇒ True; else fall through. *(The char-p
    "`< char` bypasses" refinement is UNREACHABLE here — with 0 bypasses it is
    vacuous — so it is deleted, not a separate branch.)*
  - **R2 separation (LAST — the expensive convex sweep):** compute the R16 certificate
    `strong_cert = is_strongly_simply_connected(A, convex_budget)` (Task 5); if
    `strong_cert.verdict is True` ⇒ strongly simply connected ⇒ simply connected.
- **Else → `verdict = None`** (inconclusive; honest per Adian–Rabin): `reason` names
  the obstruction — bypasses present (π₁ is presentation-dependent, no privileged-
  presentation route is built), a perfect-π₁ corner (π₁^{ab} = 0 but Tietze-lite did
  not trivialise), or `strong_cert.verdict is None` (budget/char, B4). NEVER a
  search-derived True.

**Laziness (B2).** `strong=False` skips R2 entirely (cheap routes only ⇒
`is_simply_connected` can NEVER raise/refuse on an R1/R3-decidable input);
`strong="auto"` computes the certificate only if R1/R3 did not decide; `strong=True`
always computes it (the `simply_connected` compute kind + the P62 gate pass this).
`SimpleConnectivity.strongly` is `None` when not computed. Because
`is_strongly_simply_connected` maps every internal char/budget event to a `None`
verdict (B4), R2 never raises mid-sweep — so `is_simply_connected` is total (only the
presentation-less guard raises).

**Tietze-lite** (decidable *sufficient* triviality test): drop trivial relators;
repeatedly eliminate a generator that appears with total exponent ±1 in some
relator (Nielsen/Tietze substitution), simplify; if all generators are eliminated,
the group is trivial. It only ever certifies triviality — it never returns False.

- [ ] **Step 1: Write the failing tests**

```python
# tests/coverings/test_simply_connected.py
"""Three-valued is_simply_connected. The CORRECTED square oracle drives the
True/False split; trees True; disconnected/cyclic/nontrivial-pi1ab False;
bypass-bearing perfect corner => None (Adian-Rabin honesty)."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.invariants.coverings import (
    bypasses, fundamental_group, is_simply_connected)

pytestmark = pytest.mark.oracle_literature


def test_square_with_relation_true():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=["a*b - c*d"], field=GF(7))
    assert is_simply_connected(A).verdict is True          # no bypasses + pi1 trivial (R3)


def test_square_without_relation_false():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
               ).algebra(relations=[], field=GF(7))
    r = is_simply_connected(A)
    assert r.verdict is False and r.witness["kind"] == "nontrivial_pi1ab"


def test_tree_true():
    assert is_simply_connected(linear_path_algebra(5, field=GF(7))).verdict is True


def test_loop_not_triangular_false():
    A = truncated_polynomial(3, field=GF(7))               # loop: not triangular
    r = is_simply_connected(A)
    assert r.verdict is False and r.witness["kind"] == "oriented_cycle"


def test_disconnected_false():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)}).algebra(field=GF(5))
    r = is_simply_connected(A)
    assert r.verdict is False and r.witness["kind"] == "disconnected"


def test_never_true_from_failed_search():
    # ADMISSIBLE bypass example (B3): alpha:1->3 bypasses the path a*b:1->2->3; the
    # suffix arrow s:3->4 lets the ADMISSIBLE (rad^2) relation alpha*s - a*b*s kill the
    # triangle cycle => pi1^ab = 0 (NOT caught by the nontrivial-pi1ab False guard).
    # With only the cheap routes (strong=False): NOT a tree (R1 off), a bypass is
    # present (R3 off) => the honest verdict is None -- the cheap routes can NEVER
    # fabricate True from a failed search.
    Q = Quiver([1, 2, 3, 4], {"al": (1, 3), "a": (1, 2), "b": (2, 3), "s": (3, 4)})
    A = Q.algebra(relations=["al*s - a*b*s"], field=GF(7))
    assert fundamental_group(A).free_rank == 0             # pi1^ab = 0 (cycle killed)
    assert bypasses(A)                                     # (al, a*b) is a bypass
    assert is_simply_connected(A, strong=False).verdict is None   # never a search-True
```

- [ ] **Step 2: Run to verify failure**, implement, iterate.
- [ ] **Step 3: Implement** `bypasses`/`has_double_bypass` (bounded oriented-path
  enumeration between each arrow's endpoints — finite because triangular; loud
  budget refusal if a pathological quiver blows the cap), the `SimpleConnectivity`
  dataclass, `is_simply_connected` with the **cheapest-first route order R1→R3→R2**
  and the **lazy `strong`** parameter (R2 / `is_strongly_simply_connected` invoked
  only when `strong` requires it), and Tietze-lite (a sufficient triviality test:
  drop trivial relators; eliminate a generator occurring with exponent ±1 in some
  relator by substitution; iterate; empty ⇒ trivial — it NEVER returns False). The
  never-True test is admissible above; the load-bearing assertion is that the cheap
  routes (`strong=False`) yield `None`, never a search-derived `True`.
- [ ] **Step 4: Run** — `... -m pytest tests/coverings/test_simply_connected.py -q` green.
- [ ] **Step 5: Commit** — `feat(invariants): three-valued is_simply_connected (tree/no-bypass True routes, decidable False witnesses, Adian-Rabin None)`

---

### Task 5: the separation condition + `is_strongly_simply_connected` (R16, the P62 gate)

**Files:**
- Modify: `src/quiverlab/invariants/coverings.py`, `src/quiverlab/core/algebra.py`
  (`Algebra.separation_condition`, `Algebra.is_strongly_simply_connected`)
- Test: `tests/coverings/test_separation.py`

**Interfaces:**
```python
@dataclass(frozen=True)
class Separation:
    holds: bool | None         # True / False / None (undecided_char: decompose refused)
    witness: dict | None       # on False: {vertex, summand_supports: [set, set]}
    reason: str                # "separated" | "vertex <a>: summands share a component"
                               #   | "undecided_char (decompose refused over char<=dim)"
def separation_condition(A) -> Separation
    # For triangular A, at each vertex a: decompose rad P_a into (module, mult) pairs
    #   for M, mult in decompose(A.projective(a).radical()):  ...
    # the DISTINCT indecomposable summands' supports (nonzero entries of
    # M.dimension_vector()) must lie in DISTINCT connected components of
    #   Q_a = induced_subquiver({ v : a is UNREACHABLE from v })
    #       = Q minus a and its TRANSITIVE predecessor closure (incl. a itself).
    # RAISES loudly ONLY on non-triangular input. A decompose char refusal
    # (char <= dim M, Plan-30) is CAUGHT -> holds=None, reason="undecided_char"
    # (never a raise).

@dataclass(frozen=True)
class StrongSimpleConnectivity:
    verdict: bool | None       # True / False / None (budget_exceeded OR undecided_char)
    witness: dict | None       # on False: {convex_subset, vertex, summand_supports}
    reason: str                # "all convex separated" | "<subset>@<vertex> fails"
                               #   | "budget_exceeded" | "undecided_char@<subset>"
    checked_convex: int        # number of convex subcategories tested
def is_strongly_simply_connected(A, convex_budget=20000) -> StrongSimpleConnectivity
    # Skowronski (CMS 14, 1993): triangular A is strongly simply connected iff every
    # full convex subcategory satisfies the separation condition. Enumerate convex
    # subsets of Q_0 with pruning; test separation on each induced subquiver algebra.
    # First HARD failure (Separation.holds is False) -> verdict False + witness.
    # A per-subcat undecided_char (Separation.holds is None) OR exceeding
    # convex_budget -> verdict None with the reason recorded (B4): NEVER a raise from
    # inside the sweep, NEVER a silent True. Non-triangular TOP input raises loudly.
```

**Definition pinned** (# PIN — transcription-check ASS2006 + `skowronski_ssc`; the
convention was CORRECTED after adversarial review — see the Reference section
hand-derivation): the subquiver Q_a at vertex a = the full subquiver on
`{ v : a is UNREACHABLE from v }` — i.e. delete `a` together with its **transitive
predecessor closure** `{ x : ∃ path x ⇝ a }` (which includes `a`). Implement via
`A.quiver.induced_subquiver(set(vertices) - transitive_predecessors(a) - {a})`. rad P_a
is supported on proper successors of a ⊆ Q_a (a successor cannot be a predecessor —
triangular). **Convex subcategory** = full subquiver on a convex vertex set (Zito §2).
The algebra of a convex subset S is `A.quiver.induced_subquiver(S).algebra(relations=<I
restricted to S>, field=...)` — restrict the stored relations to those whose paths stay
in S (drop any touching a removed vertex); admissibility is inherited (a convex subset
of an acyclic quiver stays acyclic).

**Honest scope (record on the verification page):** (a) requires TRIANGULAR — raise
loudly on non-triangular TOP input; (b) `decompose` is rigorous only over char 0 or
char > dim (Plan-30) — a per-vertex refusal is CAUGHT as `Separation.holds=None`
(undecided_char) and mapped by the sweep to `verdict=None` (B4), NEVER a raise; run
oracles over QQ / GF(32003) where decompose decides; (c) convex-subset enumeration is
worst-case exponential — a `convex_budget` bounds it and a `verdict=None`
(budget_exceeded) is returned honestly, never a silent True.

**The P62 handshake:** `StrongSimpleConnectivity` is the clean consumable certificate
P62 consumes — its `verdict is True` scopes P62's Tits-form tame/wild verdict layer
(Kasjan–Skowroński: strongly s.c. ⇒ tame ⟺ Tits form weakly nonnegative). Keep the
dataclass stable; P62 reads `verdict` + `witness` only.

- [ ] **Step 1: Write the failing tests**

```python
# tests/coverings/test_separation.py
"""Separation condition + strongly-simply-connected recognizer (R16, gate for P62).
Oracles: STAR/tree separated (the corrected-convention discriminator); a branching
source separated inside a NON-tree; the Zito example simply connected but NOT strongly;
the hereditary square genuinely NOT separated. Refs: Skowronski CMS 14 (1993); ASS2006;
arXiv:1905.06028.

Oracle discrimination (H2): `test_star_separated` + `test_deeper_tree_separated` catch
the B1 subquiver-convention bug (the WRONG Q_a mislabels trees as non-separated);
`test_branching_source_separated_in_nontree` shows the fix is not tree-only (a
branching source with 2 summands, separated inside a Betti-1 quiver);
`test_hereditary_square_not_separated` + `test_zito_...` confirm genuine
non-separation is still detected (no false True). The H1 minimal-relations bug is
caught by `test_two_independent_pairs_not_over_glued` in tests/coverings/test_pi1.py."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.coverings import (
    is_strongly_simply_connected, separation_condition, is_simply_connected)

pytestmark = pytest.mark.oracle_literature


def test_linear_An_separated_and_strongly():
    A = linear_path_algebra(4, field=QQ)                 # uniserial: rad P_a has <=1 summand
    assert separation_condition(A).holds is True
    assert is_strongly_simply_connected(A).verdict is True


def test_star_separated():
    # THE convention discriminator (B1): the star 1->2, 1->3 is a tree (strongly s.c.).
    # rad P_1 = S_2 (+) S_3 (two summands). CORRECT Q_1 = {2,3} (delete 1 and its
    # closure) has TWO components {2},{3} -> separated. The OLD (wrong) Q_1 = {1,2,3}
    # keeps 1, joins {2},{3} into one component -> would spuriously FAIL.
    A = Quiver([1, 2, 3], {"a": (1, 2), "c": (1, 3)}).algebra(field=QQ)
    assert separation_condition(A).holds is True
    assert is_strongly_simply_connected(A).verdict is True
    assert is_simply_connected(A).verdict is True         # tree => R1


def test_deeper_tree_separated():
    # 1->2->4, 1->3 : Q_1 = {2,3,4} has components {2,4},{3}; rad P_1 = a.A (+) c.A
    # with supports {2,4},{3} -> separated. (A deeper branching tree.)
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3)}).algebra(field=QQ)
    assert separation_condition(A).holds is True
    assert is_strongly_simply_connected(A).verdict is True


def test_branching_source_separated_in_nontree():
    # H2: a NON-tree (Betti 1: the commutative square 2-4-6-5) with a branching source
    # at 1. rad P_1 = P_2 (+) P_3 (two NON-iso summands), supports {2,4,5,6} and {3} in
    # DISTINCT components of Q_1 = {2,3,4,5,6} -> vertex 1 separated. The relation p*r-q*s
    # keeps vertex 2 separated too. Shows the corrected convention is not tree-only.
    Q = Quiver([1, 2, 3, 4, 5, 6],
               {"a": (1, 2), "c": (1, 3), "p": (2, 4), "q": (2, 5),
                "r": (4, 6), "s": (5, 6)})
    A = Q.algebra(relations=["p*r - q*s"], field=QQ)
    assert A.quiver.is_acyclic()                          # triangular; underlying graph non-tree (Betti 1)
    assert separation_condition(A).holds is True          # branching source 1 separated


def test_zito_example_simply_but_not_strongly():
    # Q: 1->2, 2->3->5, 2->4->5 ; I = <al*be*ga - al*de*ep>. pi1(Q,I)=0 (simply
    # connected) but separation FAILS (at vertex 2 already: rad P_2 = be.A (+) de.A,
    # supports {3,5},{4,5} share 5 in the connected Q_2 = {3,4,5}) -> NOT strongly s.c.
    Q = Quiver([1, 2, 3, 4, 5],
               {"al": (1, 2), "be": (2, 3), "ga": (3, 5), "de": (2, 4), "ep": (4, 5)})
    A = Q.algebra(relations=["al*be*ga - al*de*ep"], field=QQ)
    assert is_simply_connected(A, strong=True).verdict is True   # pi1(Q,I) = 0 (Zito, R3)
    ssc = is_strongly_simply_connected(A)
    assert ssc.verdict is False and ssc.witness["vertex"] == 2


def test_hereditary_square_not_separated():
    # pure commutative square (no relation): rad P_2 = be.A (+) de.A, supports {3,5}
    # and {4,5} share 5; Q_2 = {3,4,5} connected -> separation fails at vertex 2.
    Q = Quiver([2, 3, 4, 5], {"be": (2, 3), "ga": (3, 5), "de": (2, 4), "ep": (4, 5)})
    A = Q.algebra(relations=[], field=QQ)
    sep = separation_condition(A)
    assert sep.holds is False and sep.witness["vertex"] == 2


def test_non_triangular_raises():
    from quiverlab import truncated_polynomial
    with pytest.raises(QuiverlabError, match="triangular"):
        separation_condition(truncated_polynomial(3, field=QQ))


@pytest.mark.oracle_selfcert
def test_char_caveat_is_undecided_not_a_raise():
    # over a small GF(p) with char <= dim rad P_a, decompose refuses; separation must
    # CATCH it as undecided_char (B4), and the sweep maps it to verdict None -- NEVER a
    # raise from inside the sweep, NEVER a silent verdict. (Pick p and a quiver where
    # some rad P_a has dim >= p at implementation.)
    A = linear_path_algebra(6, field=GF(2))
    sep = separation_condition(A)
    assert sep.holds in (True, False, None)               # decided OR undecided_char
    ssc = is_strongly_simply_connected(A)
    assert ssc.verdict in (True, False, None)             # never raises; None if char-blocked
```

- [ ] **Step 2: Run to verify failure**, implement, iterate.
- [ ] **Step 3: Implement** `separation_condition` (raise loudly if not triangular;
  per vertex a: `for M, mult in decompose(A.projective(a).radical()):` collect each
  summand's support = nonzero entries of `M.dimension_vector()`; distinct-component
  test against the CORRECTED `Q_a = induced_subquiver(vertices − transitive_predecessors(a) − {a})`'s
  `undirected_components`; witness = two distinct summands sharing a component. Wrap
  the `decompose` call in `try/except QuiverlabError` → `holds=None,
  reason="undecided_char"` — never re-raise). Then `is_strongly_simply_connected`
  (enumerate convex subsets with pruning — grow convex sets, prune non-convex; bound
  by `convex_budget`; build each convex subalgebra via `induced_subquiver(S)` +
  restricted relations; call `separation_condition`: `holds is False` → verdict
  False + witness; `holds is None` (undecided_char) → verdict None; exceeding the
  budget → verdict None; else all separated → verdict True). NEVER raise from inside
  the sweep. The char-caveat test asserts the undecided path (None), NOT a raise.
- [ ] **Step 4: Run** — `... -m pytest tests/coverings/test_separation.py -q` green;
  wire `is_simply_connected`'s R2 route to call `is_strongly_simply_connected` and
  re-run Task 4.
- [ ] **Step 5: Commit** — `feat(invariants): separation condition + strongly-simply-connected recognizer (R16, the P62 gate), witness on failure`

---

### Task 6: QPA honest-scope skip

**Files:**
- Create: `tests/qpa/test_pi1_qpa.py`

**Interfaces:** QPA 1.37 has **no** fundamental-group / simply-connected surface.
Probe live via `NamesGVars()` (the `tests/qpa/test_products_qpa.py` precedent): the
test SKIPS honestly and **FAILS if QPA ever ships** `FundamentalGroup` /
`IsSimplyConnected` / `SeparationCondition` (so our honest-scope claim on the
verification page cannot silently rot).

- [ ] **Step 1: Write the probe test** (header
  `pytestmark = pytest.mark.skipif(session.should_skip_qpa(), reason=...)`); sweep
  `NamesGVars()` for the three verb names; `pytest.skip(...)` if absent, `assert
  False` (loud) if any appears.
- [ ] **Step 2: Run** — `... -m pytest tests/qpa/test_pi1_qpa.py -q -m qpa` (venv has
  `[qpa]`): expect a clean skip.
- [ ] **Step 3: Commit** — `test(qpa): honest no-pi1-surface probe (fails if QPA ships one)`

---

### Task 7: no-code exposure — `fundamental_group` + `simply_connected` compute kinds

The organizing invariant is **byte-identity across the two runners**, achieved by
routing BOTH through a single **shared library builder** (the `strings_block` /
`quasi_hereditary_block` / `recognizers_block` precedent — each runner only adds
`citations` locally). Anchors below are current-tree (verified during authoring);
still "adjust to reality" if a merge shifts them.

**Files (the seven-touchpoint checklist as it exists in the tree):**
1. **Shared builders** — Create `src/quiverlab/invariants/coverings_block.py`:
   `fundamental_group_block(A) -> dict` and `simply_connected_block(A) -> dict`,
   each returning a dict carrying a `"references"` list (both runners import these,
   so the block is byte-identical by construction).
2. **Server/HPC** — Modify `src/quiverlab/hpc/spec.py`: two `if kind == ...:`
   branches in `_dispatch` (~line 1516 `strings`/`recognizers` pattern) calling the
   builder + `block["citations"] = _citation_pairs(block["references"])`
   (`_citation_pairs` at `spec.py:1557`); two `_snip` reproduce entries in the map
   at `spec.py:2268-2334`. Do **NOT** add to `MODULE_KINDS` (`spec.py:72-77`) — these
   are algebra-level scalars routed through `_dispatch`, not `_dispatch_module`.
   There is NO compute-kind enum/whitelist: the `_dispatch` elif IS the registration
   (the only unknown-kind guard is `spec.py:1554`).
3. **Pyodide twin** — Modify `docs/gui/runner.py`: two byte-identical
   `elif name == ...:` branches in `compute_one` (~line 842, mirroring the server);
   two entries in the `calls` reproduce map (`runner.py:1160-1195` — **hard-indexed,
   KeyErrors if missing**); two `ETA_MODEL["scalars"]` costs (`runner.py:1255-1291`;
   default is a silent `0.1` if absent — set `fundamental_group` ~0.5,
   `simply_connected` ~2.0).
4. **GUI JS** — Modify `docs/gui/gui.js` AND `webapp/static/gui/gui.js` (they are
   `diff`-identical — apply every edit to both): (a) checkbox
   `<input id="qlgui-fundamental_group">` / `qlgui-simply_connected` in the
   invariants block (~lines 99-118); (b) the bare id in the id-registry array
   (~201-240); (c) the scalar push-list array (~857-861); (d) `KIND_CTRL` entries
   (~3894-3931, `{cb: "<kind>"}`); (e) a `renderBlock` `else if (name === ...)`
   branch before the citations mount (~line 3249); (f) the `scheduleProbe`
   change-listener list (~3340-3348); (g) the `THEMES` `kinds` array (invariants
   theme ~line 3427) inside the `QLGUI-THEMES-BEGIN/END` sentinels. Optional:
   structured renderers in `webapp/static/app.js` (~605-631) — otherwise the kinds
   surface via the JSON dump.
5. **i18n (FOUR locales) + picker template** — `webapp/server/i18n/` has
   `LANGS = ("en","es","fr","zh")` (`i18n/__init__.py:12`). Add to **all four**
   `{en,es,fr,zh}.json`: `pick.kind.fundamental_group` / `pick.kind.simply_connected`
   (REQUIRED — gated), plus `inv.fundamental_group` / `inv.simply_connected`
   (checkbox labels) and `block.fundamental_group.title` /
   `block.simply_connected.title` (+ per-row labels: abelianization, intrinsic_note,
   verdict_{yes,no,undecided}, strongly). Add one
   `data-pick-kind-<kind>="{{ t('pick.kind.<kind>') }}"` line per kind in
   `webapp/templates/draw.html` (~lines 46-81).
6. **Report renderer** — Modify `src/quiverlab/trace/results_html.py`: two
   `if kind == ...:` branches in `_block_html` (~718-731 scalar pattern; render the
   witness + the honest-decidability sentence + the intrinsic-refusal note as
   documentation) and two `_HEADINGS` entries (~lines 42-67).
7. **Tests / gates** — add both kinds to `ALL_KINDS` in
   `tests/webapp/test_layout_picker.py:26-37` (must equal the THEMES kind set,
   each once — `test_themes_cover_every_kind_exactly_once` + the four-locale
   `pick.kind.*` translation gate); add two goldens to
   `tests/webapp/_runner_goldens.json` with a dated change-log bullet in
   `tests/webapp/test_runner_delegation.py` (verify existing entries byte-identical
   FIRST); add a twin-parity test (copy `tests/webapp/test_quasi_hereditary_p47.py::
   test_twin_parity`).
- Test files: `tests/webapp/test_coverings_exposure_p56.py`,
  `tests/hpc/test_p56_kinds.py`.

**Schema stays v1** (`ComputeRequest.schema_version`, `schema.py:236`; v2 is required
only for `module`/`ext_target`/`tor_target` blocks). **Canonical keys are
request-derived** (`cache.py::canonical_key`, sha256 over the sorted-keys blob +
`library_version()`): a scalar-only request that does not use the new kinds keys
byte-identically to before, so existing goldens are untouched.

Two new blocks:
- `fundamental_group`: `{"generators": [...], "relators": [...], "abelianization":
  {"free_rank": int, "invariant_factors": [ints]}, "hom_to_additive_dim": int,
  "components": int, "presentation_note": str, "latex": "\\pi_1^{ab} = ...",
  "intrinsic_note": "<the loud-refusal sentence, as DOCUMENTATION — the intrinsic
  group is not computed>", "references": ["assem_delapena","martinez_villa_delapena",
  "crs_hurewicz","crs_gradings","briggs_ryd_tori"], "citations": ...}`. The HH¹ bound
  is NOT computed in the payload (keeps the kind cheap; the `≤ HH¹` cross-check is
  verification-side).
- `simply_connected`: `{"verdict": true|false|null, "reason": str, "witness":
  {...}|null, "strongly": {"verdict": true|false|null, "witness": {...}|null,
  "reason": str, "checked_convex": int}, "references": ["assem_book",
  "skowronski_ssc","le_meur_pi1"], "citations": ...}` (the compute kind passes
  `strong=True` so the `strongly` certificate — the P62 consumable — is always
  present; a char/budget event surfaces as `strongly.verdict = null` with the reason,
  never an error). The three-valued verdict
  renders as **"simply connected" / "not simply connected (witness: …)" / "undecided
  (triviality of a finitely presented group is undecidable — Adian–Rabin)"** in GUI +
  report + all four locales. A refusal on this input (presentation-less, or the
  separation char-caveat / budget) is reported as `{"error": "<the loud message>"}`
  per-verdict, never a silent default (the Plan-30 honest-per-entry precedent).

- [ ] **Step 1: Write the failing cross-runner test** (copy
  `tests/webapp/test_module_blocks_m0729.py` / `test_quasi_hereditary_p47.py`
  runner-pair fixtures; unmarked file — extras-gated dir):

```python
# tests/webapp/test_coverings_exposure_p56.py
"""fundamental_group + simply_connected kinds: served by hpc.spec, mirrored by the
Pyodide twin, three-valued verdict rendered honestly, both runners byte-identical."""
import json


def test_fundamental_group_block_shape(tmp_path):
    # request: commutative square WITH relation over GF(7), compute=["fundamental_group"]
    #   block["abelianization"]["free_rank"] == 0
    #   block["abelianization"]["invariant_factors"] == []
    #   "assem_delapena" in block["references"]
    #   "intrinsic_note" in block            # the honest-decidability documentation
    ...


def test_simply_connected_three_valued(tmp_path):
    # square WITHOUT relation: verdict False, witness kind "nontrivial_pi1ab";
    # the Zito example: verdict True, strongly.verdict False.
    ...


def test_twin_parity(tmp_path):
    # run both requests through docs/gui/runner.py::compute_one and hpc.spec.run;
    # assert json.dumps(block, sort_keys=True) equality on both blocks.
    ...
```

- [ ] **Step 2: Implement** the shared builders (touchpoint 1), then the server +
  twin dispatch (2, 3), the GUI JS in both `gui.js` files (4), the four-locale i18n
  + `draw.html` (5), and the `results_html.py` branches (6). `_snip` recipes:
  `"fundamental_group": lambda it: "A.fundamental_group()"`,
  `"simply_connected": lambda it: "A.is_simply_connected()"`.
- [ ] **Step 3: Gate the layout/i18n first** —
  `... -m pytest tests/webapp/test_layout_picker.py tests/webapp/test_js_parses.py -q`
  (THEMES == `ALL_KINDS`; `pick.kind.*` present in all four locales).
- [ ] **Step 4: Add two golden fixtures** (`fundamental_group_square_gf7`,
  `simply_connected_zito`) to `_runner_goldens.json`; note them in the
  `test_runner_delegation.py` docstring change-log. Run the delegation test BEFORE
  adding to confirm existing entries stay byte-identical; on a rebase, merge the
  JSON **semantically** (dict-union), never textually.
- [ ] **Step 5: Run the gates** —
  `... -m pytest tests/webapp/test_coverings_exposure_p56.py tests/webapp/test_runner_delegation.py tests/hpc -q`.
  Expected PASS; both runners byte-identical.
- [ ] **Step 6: Commit** — `feat(gui,webapp,hpc): fundamental_group + simply_connected compute kinds -- pi1 and strong simple connectivity clickable end-to-end (4 locales, v1 schema)`

---

### Task 8: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Test: existing release gates

- [ ] **Step 1: Verification page.** Add the Plan-56 subsystem row
  (`invariants/coverings.py` | oracles: `oracle_literature` — the corrected square
  (WITH ⇒ 1 / WITHOUT ⇒ ℤ), Le Meur Example 1, trees, single loop, multi-loop, the
  annulus, the Zito simply-but-not-strongly example, kA_n separated;
  `oracle_crossengine` — Hom(π₁,k⁺) ≤ dim HH¹ over the triangular zoo;
  `oracle_selfcert` — the graph primitives, SNF-abelianization consistency, the
  separation witness validity, three-valued-verdict-never-True-from-search;
  `qpa` — the honest no-surface probe). Add the **honest-scope entries**:
  (a) the INTRINSIC π₁ is refused loudly — not bounded-computable; the k[x]/(x^p) =
  ℤ × C_p intrinsic oracle documents *why* the presentation group (ℤ) is not the
  algebra invariant; (b) `is_simply_connected` is three-valued True/False/None —
  `None` is honest per Adian–Rabin (triviality of an f.p. group is undecidable);
  True only via decidable sufficient criteria (tree / no-bypass Le Meur / separation);
  (c) the separation recognizer requires triangular + char 0 or char > dim (Plan-30
  `decompose` caveat) and bounds the convex sweep with a loud budget refusal;
  (d) the `Hom(π₁,k⁺) ↪ HH¹` bound is asserted for TRIANGULAR algebras only
  (Assem–de la Peña embedding); (e) QPA has no π₁ surface. Recount the class table
  (`tests/release/test_oracle_classes.py` drives the numbers — collect, paste,
  re-run to green).
- [ ] **Step 2: README.** One features line: "the fundamental group π₁(Q,I) with
  exact abelianization, the Hurewicz `Hom(π₁,k⁺) ↪ HH¹` check, and a
  strongly-simply-connected recognizer (separation condition) — clickable via
  `fundamental_group` / `simply_connected`."
- [ ] **Step 3: Full gate** — `... -m pytest -q -m fast` green;
  `... -m pytest tests/coverings -q` green; `... -m pytest tests/qpa -q -m qpa`
  green; `... -m pytest tests/release -q` green.
- [ ] **Step 4: Commit** — `docs(verification): Plan-56 oracle rows + honest scope (intrinsic pi1 refused, Adian-Rabin None, separation char/budget caveats) + recounted classes`

---

## Acceptance (Plan-56 definition of done)

1. `A.fundamental_group()` public — finite presentation (generators + relators) +
   `π₁^{ab}` (free rank + invariant factors) by exact ℤ SNF. **Genuinely always
   computable**: minimal relations = exact linear algebra on I/(rad·I + I·rad)
   (no greedy-reduction unsoundness, NO refusal path), self-cert
   count == dim I/(rad·I + I·rad). `A.intrinsic_fundamental_group()` ALWAYS raises the
   loud honest refusal (documenting the k[x]/(x^p) = ℤ × C_p intrinsic-vs-presentation
   gap). Discriminating oracle: the two-independent-pairs non-minimal generating set
   is NOT over-glued (π₁^{ab} = ℤ, not 0).
2. The CORRECTED oracle green: commutative square WITH the commutativity relation
   ⇒ π₁ = 1; WITHOUT ⇒ ℤ. Trees ⇒ trivial; single monomial loop ⇒ ℤ (presentation
   group, not the intrinsic ℤ × C_p); multi-loop free ranks; Le Meur Example 1;
   the annulus ⇒ π₁^{ab} = ℤ.
3. `Hom(π₁,k⁺) ↪ HH¹` holds degreewise over the triangular zoo (`oracle_crossengine`),
   equality on the Schurian square (Hurewicz iso).
4. `A.is_simply_connected()` three-valued (True only via tree / no-bypass-Le Meur /
   separation; False via disconnected / oriented-cycle / nontrivial-π₁^{ab};
   else `None`), NEVER True from a failed search.
5. `A.is_strongly_simply_connected()` (R16) — separation (CORRECTED convention:
   `Q_a` deletes `a` and its transitive predecessor closure) for every full convex
   subcategory; raises only on non-triangular top input; maps per-subcat char refusal
   (`undecided_char`) and budget to `verdict = None` (never a mid-sweep raise);
   witness on failure. Discriminating oracles green: the star + deeper tree +
   branching-source-in-a-non-tree come out separated (the convention fix); the Zito
   example is simply connected but NOT strongly (split green); the hereditary square
   genuinely fails. The `StrongSimpleConnectivity` dataclass is the clean consumable
   certificate P62 reads.
6. `fundamental_group` + `simply_connected` clickable end-to-end (GUI canvas → block
   → report) in all four locales (en/es/fr/zh), both runners byte-identical via a
   shared block builder, two goldens added with a documented change-log entry,
   `test_layout_picker.py` `ALL_KINDS`/THEMES gate green, schema still v1, canonical
   keys unchanged.

---

## Change log

- **2026-08-07 adversarial review: 3 blocking + 2 majors + 7 minors applied** —
  minimal relations = exact I/(rad·I + I·rad) linear algebra with lifting (always
  computable, refusal path deleted; two-independent-pairs oracle added); separation
  subquiver = transitive-predecessor-closure deletion incl. `a` (star/tree
  counterexample fixed; star + deeper-tree + branching-in-non-tree oracles added);
  lazy `strongly` + True-route reorder R1→R3→R2; typed `undecided_char` (char/budget
  → `verdict None`, never a mid-sweep raise); annulus re-marked oracle_selfcert
  (hereditary, Betti-of-4-cycle); fingerprint provenance corrected (same helper,
  explicit `domain=ZZ`); unreachable R3 char-p branch deleted; explicit `field=` on
  every oracle snippet; `decompose` unpacked as `(module, mult)`; Le Meur Ex.1 and
  never-True-from-failed-search rewritten as admissible rad² constructions.
7. QPA honest-scope probe skips (fails if QPA ever ships a π₁ surface);
   `docs/verification.md` recounted with the five honest-scope entries; README line
   added; fast + coverings + qpa + release suites green.

---

## Addendum (2026-08-07, P62 contract): expose the per-pair minimal-relation counts

**Why.** P62 (Tits-form tame/wild, R19) needs the Tits form
`q_A(x) = Σ x_i² − Σ_{arrows i→j} x_i x_j + Σ_{i,j} r_{ij} x_i x_j`, where
`r_{ij}` = the number of **minimal relations** from `i` to `j` =
`dim_k e_j·(I / (rad·I + I·rad))·e_i`. **This is exactly the per-`(x,y)`-block quantity
Task 2's `fundamental_group` already computes** (`dim W_{x,y} − dim R_{x,y}`, the
`I/(rad·I + I·rad)` complement-basis count). Rather than P62 forking that unsound-if-done-
naively linear algebra, P56 retains and exposes it.

**Contract (P56 delivers; P62 consumes read-only).** Add to
`src/quiverlab/invariants/coverings.py`:
```python
def minimal_relation_counts(A) -> dict:
    """{(src, tgt): count} for every ordered vertex pair with >= 1 minimal relation
    (= dim_k e_tgt (I/(rad.I + I.rad)) e_src). Presentation-less A: loud QuiverlabError.
    Empty dict for a hereditary (relation-free) algebra."""
```
**Honest scope of the change (NOT a five-line extraction).** Task 2's block loop currently
folds each block's minimal-relation basis straight into the ℤ SNF rows of π₁^{ab} without
retaining the per-`(x,y)` *count*. So this is a **small refactor**: retain the per-pair
count `dim W_{x,y} − dim R_{x,y}` in a dict while the loop runs, return it from a shared
helper, and have BOTH `fundamental_group` (unchanged output) and the new public
`minimal_relation_counts` read that single computation. **P56's existing outputs and
goldens stay byte-identical** — this only ADDS a public accessor over data already formed.
Gate it with a one-line self-cert test in `tests/coverings/` (the total count ==
`dim I/(rad·I + I·rad)` == the count `fundamental_group`'s self-cert already asserts), and
P62's `test_r_ij_matches_ext2_of_simples` cross-checks it against `dim Ext²(S_i, S_j)`.
