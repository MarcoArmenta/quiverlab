# 15 — The τ-tilting engine

## What this computes

`tautilting/` (Plan 45, Adachi–Iyama–Reiten, *Compos. Math.* 150 (2014)) is the C4
engine: support τ-tilting pairs of an algebra `A`, their g-vectors, the mutation
**exchange graph**, the torsion-class lattice with brick / semibrick labels, 2-term
silting, King θ-stability with the wall-and-chamber fan, and maximal green sequences.
Two disciplines run through all of it. First, every enumeration is **budget-capped with
the honest complete-iff-τ-tilting-finite contract** — a run either closes as
`status="complete"` or stops loudly at `status="budget"`, never a silently truncated
answer. Second, every geometric payload is **exact rational** (no floats in `src/`); the
only fraction→pixel conversion happens in the browser.

## g-vectors and τ-rigidity (`rigid.py`)

`is_tau_rigid(M)` is the definition `Hom_A(M, τM) = 0` (AIR); since `τ(projective) = 0`,
every projective is τ-rigid. `g_vector(M)` reads the **minimal projective presentation**
`P_1 → P_0` (`minimal_resolution(M, 1)`) and returns the vertex-keyed dict
`g^M(v) = (multiplicity of P_v in P_0) − (multiplicity of P_v in P_1)` — the convention
"P₀ positive, P₁ negative", so `g^{P_v} = e_v` and g is additive on direct sums, valued in
`K_0(proj A)`. For a support τ-tilting pair, `g_columns(pair)` is the ordered list of
columns — the module summands' g-vectors, then `−e_v` for each killed projective `P_v[1]`
in the support — and `g_matrix(pair)` assembles them with **rows = vertices, columns =
g-columns**. The set of columns determines the pair (`SupportTauTiltingPair.g_key`, a
frozenset of column tuples); the matrix is unimodular (`det = ±1`, AIR Thm 5.1).

## The mutation arbiter — a uniform 2-term silting engine (`mutation.py`, `_twoterm.py`)

The naive picture — mutate a pair by an add-approximation at the module level — is
**incomplete for the module ↔ support crossing**, and the engine does not use it: over
kA₂, un-killing vertex 1 of the pair `(S₂, {1})` brings back `P₁`, which no approximation
of `P₁` against `S₂` produces. Only the cone/cocone in `K^b(proj A)` plus a
minimal-complex reduction recovers it. So mutation is done uniformly on **2-term complexes
of projectives** (`_twoterm.py::PComplex`, cohomological, degrees `−1, 0`):

- `mutate(pair, k)` forms the two candidate summands at position `k` — the **left**
  mutation `cone(f)` of the minimal left add(U)-approximation and the **right** mutation
  `cocone(g) = cone(g)[−1]` of the minimal right add(U)-approximation — where `U` is the
  direct sum of the other summands. Approximations use `hom_kb` (chain maps modulo
  null-homotopies, pure Domain linear algebra); minimality is by retract pruning.
- `reduce_complex` / `_cancel` drive each candidate to its **minimal complex**: repeatedly
  find an isomorphism block `P_v → P_v` inside a differential and cancel it by a Schur
  complement, until every differential's entries lie in the radical (the delooping lemma;
  homotopy-equivalent, so the summand iso-type is preserved). `summand_class` then reads
  the reduced complex as a `"module"` (the cokernel of `P_1 → P_0`) or a `"support"` shift
  `P_v[1]`.

**What certifies a mutation.** `mutate` builds each candidate pair through
`make_pair(..., check=True)` (the full four-axiom validation) and accepts the first for
which the g-key change is *exactly one column*: the k-th g-column removed and exactly one
new column added. Mutation is thereby an **involution swapping exactly one g-column**
(`mutate(mutate(pair, k), k') = pair`); if neither branch validates it refuses loudly
(naming the char caveat). Unimodularity is a property of the resulting g-matrix, asserted
in the docstrings and pinned by the oracle tests, not recomputed inside `mutate`.

## The exchange graph and its budget honesty (`mutation.py`)

`exchange_graph(A, budget_pairs=512)` is a breadth-first search from `initial_pair(A) =
(A_A, ∅)`, deduping discovered pairs by `g_key()`. It mutates every pair at every position
`k`. The `ExchangeGraph` result carries the pair records, the undirected `arrows`
(edges `(i, j)`, `i < j`, each labelled by the wall brick; orientation is the separate
`hasse_orientation` step), the adjacency, `is_complete`, `status`, and `n_regular`. The
contract is exact:

- It closes as `status="complete"`, `is_complete=True`, and `n_regular` (every pair has
  exactly `n` neighbours) **iff `A` is τ-tilting-finite** — the exchange graph is then
  connected (AIR Cor 2.38).
- On a τ-tilting-infinite algebra (e.g. the 2-Kronecker) it hits `budget_pairs` and returns
  immediately with `status="budget"`, `is_complete=False` — never a silently truncated
  graph.
- A `mutate` that refuses loudly sets `status="error"`, `is_complete=False`; the failure is
  surfaced, never swallowed.

Every downstream enumeration (`bricks`, `semibricks`, `maximal_green_sequences`, the fan,
the four-way counts) inherits this honesty and omits its value when the BFS did not close.

## King θ-stability and the wall-and-chamber fan (`stability.py`)

`is_theta_semistable(M, θ)` / `is_theta_stable(M, θ)` are King's (1994) conditions —
`θ·dim M = 0` and `θ·dim N ≤ 0` (resp. `< 0`) for every (proper nonzero) submodule `N` —
with `θ` an exact-rational vertex-ordered vector and submodules enumerated exactly by
`_submodule_dimvecs` (BFS with A-closure, deduped by RREF key, loud past its budget). No
floats.

`wall_and_chamber_fan(A, budget=512)` returns, when the exchange graph closes, one
**chamber per pair** — its `g_matrix` and the exact-`Fraction` `rays` (the g-columns) —
and one **wall per exchange edge**, whose normal is the shared g-facet and whose label is
the wall brick's dimension vector (King). For `n = 3` each chamber additionally carries the
L1/octahedron projection `rays_l1`, the `faces`, and the 2-D `net2d` unfolding (`_l1_project`,
exact Fractions). The `n ∈ {2, 3}` gate itself lives in the caller
`block.py::tau_tilting_block` (the fan is `None` for other `n`). The *completeness*
certificates are **oracle tests**, not payload fields (see the verification page): for
`n = 2` an exact angular sweep proves the cones tile ℝ²; for `n = 3` the L1 unfolding is a
rendering certified only per chamber (each chamber g-matrix unimodular) plus a net-sanity
check — there is no 3-D fan-tiling certificate. The payload stays exact-rational and the
browser does the only fraction→pixel conversion.

## Bricks and iso-class labelling (`torsion.py`)

`_torsion_universe(A)` collects every indecomposable summand across all pairs of the
exchange graph, deduped by the exact `is_isomorphic` (with a dimension-vector prefilter),
cached by algebra object in a `WeakKeyDictionary`. Labelling bricks by **iso-class, not
dimension vector**, is load-bearing: over the symmetric Nakayama algebra kZ₂/rad² the two
projective-injectives `P₁` and `P₂` are non-isomorphic bricks **sharing** the dimension
vector `(1, 1)`, sitting on opposite King-stability rays of the same hyperplane. A
first-dimension-vector-match label would collapse them, undercounting `bricks()` (3 instead
of 4) and merging the distinct semibricks `{P₁}`, `{P₂}` (5 instead of 6) — silently
breaking the AIR four-way identity.

`_edge_brick` therefore uses the DIRRT labelling when several non-isomorphic bricks share a
wall normal: on the strict cover `T > T'`, the wall brick is the unique candidate whose
`Gen`-membership flips across the wall, and it **refuses to guess** (loud `QuiverlabError`)
if membership fails to single out exactly one. `bricks(A)` are the edge bricks deduped by
iso-class (`end_dim == 1` by construction); `semibricks(A)` are the Asai down-labels of each
pair, deduped by iso-class multiset (`#semibricks = #pairs` requires exactly this
dedup). `hasse_orientation(eg)` orients each edge by torsion-class inclusion (via the
`Gen(M)` size from `torsion_class_data`), with `(A, 0)` the unique source and `(0, A)` the
unique sink. A brick is `end_dim == 1`, which reads "`End_A(B) = k`" only over an
algebraically closed / char-0 base; the GF(pⁿ) proper-division-ring caveat is honest scope,
and the batteries run over ℚ.

## The four-way identity, char scope, and the payload (`block.py`, `silting.py`, `green.py`)

`tau_tilting_block(A, budget=512)` is the algebra-level `tau_tilting` compute kind. On a
closed BFS it reports the pairs (each with g-matrix, label, summand dimension vectors,
support), the oriented Hasse edges, the fan (for `n ∈ {2, 3}`), the maximal-green-sequence
count, and the empirically computed **four-way counts**
`#s-τ-tilt = #f.f. torsion = #2-term silting = #semibricks`. On an incomplete BFS it sets
`hasse`, `fan`, `green_count`, and `counts` to the honest empty / `None` — "a partial fan
would be a lie". There is **no proactive characteristic guard**: the engine runs over ℚ and
inherits the loud refusals of `decompose` / `is_isomorphic` / the trace-form radical, which
are rigorous over char 0 or char > dim (Dickson/CIW); where every module involved is a brick
or splits it computes correctly even at char ≤ dim (GF(2) kZ₂/rad² agrees with the certified
ℚ counts).

## The oracles

- **Literature** — `#s-τ-tilt(kA_n) = Catalan(n+1)` (2, 5, 14), exchange-graph
  n-regularity, the AIR four-way count identity on kA₂/kA₃, hereditary `τ-rigid ⇔ rigid`,
  kA₂ = 2 maximal green sequences, and the kZ₂/rad² four-brick / six-semibrick non-thin gate.
- **Self-cert** — `g^{P_v} = e_v` and additivity, the four-axiom pair certification,
  mutation as a one-column-swap involution, g-matrix unimodularity (det ±1), the n=2 fan
  tiling and the n=3 per-chamber-unimodular net-sanity, and King θ-stability on the worked
  kA₂ example.
- **Cross-engine** — pair ↔ `Gen(M)` torsion-class injectivity, the fan wall normals ⊥ the
  shared g-facet, and the GF(2) no-proactive-char-guard agreement.
- **QPA** — **none**: QPA 1.37 exposes no support-τ-tilting / mutation / g-vector /
  stability surface, so there is no `qpa` battery here (stated in honest scope); the
  Demonet–Iyama–Jasso tables and Iyama's `fd-applet` are named, not wired live.

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| g-vectors, τ-rigidity, g-matrix | `tautilting/rigid.py` | `is_tau_rigid`, `g_vector`, `g_columns`, `g_matrix` |
| support τ-tilting pairs + axioms | `tautilting/pairs.py` | `SupportTauTiltingPair`, `make_pair`, `initial_pair`, `terminal_pair` |
| mutation + exchange-graph BFS | `tautilting/mutation.py` | `mutate`, `exchange_graph`, `ExchangeGraph` |
| the uniform 2-term silting engine | `tautilting/_twoterm.py` | `PComplex`, `cone`, `reduce_complex`, `_cancel`, `min_left_approx`, `min_right_approx` |
| torsion lattice + brick / semibrick labels | `tautilting/torsion.py` | `torsion_class_data`, `hasse_orientation`, `bricks`, `semibricks`, `_edge_brick` |
| King θ-stability + wall-and-chamber fan | `tautilting/stability.py` | `is_theta_semistable`, `is_theta_stable`, `wall_and_chamber_fan` |
| maximal green sequences | `tautilting/green.py` | `maximal_green_sequences` |
| 2-term silting bridge | `tautilting/silting.py` | `two_term_silting`, `silting_count` |
| the `tau_tilting` block | `tautilting/block.py` | `tau_tilting_block` |
