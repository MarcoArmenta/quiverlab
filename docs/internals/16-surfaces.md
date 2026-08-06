# 16 — Marked surfaces to gentle algebras

## What this computes

`surfaces/` (Plan 48) is a no-code **input** method: pick or build a marked bordered
surface, triangulate it, and read off the associated gentle algebra
`Jac(Q(T), W(T))` — the quiver-with-potential of Fomin–Shapiro–Thurston (2008) and
Labardini-Fragoso (2009), whose Jacobian is gentle in the unpunctured case (Assem–Brüstle–
Charbonneau-Jodoin–Plamondon 2010). The produced algebra then flows through every existing
compute kind (Hochschild, resolutions, modules, products, …). All surface data is exact
integers and tuples — float-free.

**v1 scope is binding: unpunctured surfaces with non-empty boundary only.** That is the
ABCP/LFS regime where every arc-adjacency is clean, there are no self-folded triangles, and
the Jacobian is certifiably gentle. Punctures, closed surfaces, and self-folded triangles
refuse loudly, each naming the successor **P48.1**.

## The marked-surface datum (`marked.py`)

`MarkedSurface` is a frozen dataclass with `genus`, `boundary_marked` (the tuple
`(k₁, …, k_b)`, each `k_i ≥ 1` marked points on boundary component `i`), and `punctures`.
Derived accessors:

- `euler_characteristic() = 2 − 2·genus − b` (with `b = num_boundary_components`).
- `arc_count()` — the **FST arc count** `n = 6·genus − 6 + 3(b + p) + c`, where `p = punctures`
  and `c = Σ kᵢ`. This is the number of interior arcs in any ideal triangulation.
- `triangle_count() = (2n + c) / 3`.
- `in_v1_scope() = (punctures == 0 and b ≥ 1)`.

Construction is self-validating (each failure a `QuiverlabError`): non-integer / negative /
sub-1 fields, and — the FST-admissibility gate — the degenerate `n < 1` cases (monogon,
digon, triangle) and the explicit exclusion set `_FST_EXCLUDED` (sphere with 1, 2, or 3
punctures; the once-punctured monogon). Booleans are rejected explicitly so `True`/`False`
cannot smuggle in as 1/0. A `MarkedSurface` may hold punctured/closed data (so it can be
refused downstream), but the gentle pipeline only runs `in_v1_scope` surfaces.

## Triangulations (`triangulation.py`)

`Triangulation` is a frozen dataclass of a `MarkedSurface` plus `triangles` — a tuple of
anticlockwise 3-tuples of *sides*, where a side is either an interior-arc integer or a
boundary-segment string. Construction self-certifies every combinatorial constraint:

1. each triangle has exactly 3 distinct sides (a repeated side is a self-folded triangle —
   refused, naming P48.1);
2. **arc-adjacency** — each interior arc borders exactly 2 triangles, each boundary segment
   exactly 1 (a `Counter`);
3. the arc count equals `surface.arc_count()` and the triangle count equals
   `surface.triangle_count()`;
4. a boundary-component guard — per-component segment counts match the `kᵢ` (so a
   disc-shaped triangle list cannot pass on an annulus surface).

Four constructors ship: `fan_triangulation(marked)` (the fan of a `marked`-gon from one
vertex — the disc), `annulus_triangulation(n, m)` (the balanced band triangulation of the
annulus `C(n, m)`), `hexagon_with_internal_triangle()` (one central internal triangle plus
three boundary "ears"), and `once_punctured_torus()` — the **pinned loud-refusal oracle**:
it constructs as a valid FST-admissible `Triangulation` (genus 1, one puncture, three
loop-arcs) but `quiver_of` / `jacobian_of` refuse it as out of v1 scope.

## Triangulation → quiver → potential → Jacobian (`qp.py`)

`quiver_of(T)` takes the interior arcs as vertices and reads the arrows off each
anticlockwise triangle. The angle→arrow **orientation is arbitrated, not assumed**:
`_raw_arrows` walks the cyclically-consecutive ordered pairs `(a, b)` of arc-sides in each
triangle and emits the arrow `b → a` ("later anticlockwise side → earlier"). The disc oracle
fixes this: the fan of the `(n+3)`-gon must give the *linear* `Aₙ` quiver `1 → 2 → … → n`.
The naive anticlockwise reading `sᵢ → sᵢ₊₁` gave the reversed chain `n → … → 1`, so v1 ships
the flipped convention `sᵢ₊₁ → sᵢ`, documented in `qp.py` and pinned by the disc arbiter
test. `_reduce_two_cycles` cancels oriented 2-cycles — a defensive no-op in v1, since
2-cycles provably need a self-folded triangle (`C(1,1)` yields *parallel* arrows, not a
2-cycle), kept for the P48.1 punctured successor.

`potential_of(T)` adds one `+1` **3-cycle term** for each internal triangle (a triangle all
of whose sides are arcs — `_internal_triangles`), ordering its arrows into a composable cycle
with `_cycle_word` (which refuses a degenerate internal triangle). A fan has no internal
triangle, so `W = 0`. `jacobian_of(T)` reads the quiver back off the potential and returns
`JacobianAlgebra(W.quiver, W, …)` — the P44 constructor, whose success **is** the finiteness
certificate. The Jacobian is gentle (ABCP/Labardini), self-certified by P46's `is_gentle`;
the single internal triangle of the hexagon yields the dimension-6 triangle algebra.

Loud refusals live in `_require_v1_scope` (punctures out of scope; closed / no-boundary out
of scope) and `_cycle_word` (degenerate internal triangle); self-folded triangles are already
refused at `Triangulation` construction. Every message names successor P48.1.

## Flip ≡ Fomin–Zelevinsky matrix mutation (`flip.py`)

`flip(T, arc)` is the combinatorial quadrilateral flip: it reads the two anticlockwise
triangles sharing the interior `arc`, forms the surrounding quadrilateral, and replaces the
diagonal — **reusing `arc`'s id** so the exchange matrices stay index-aligned. It refuses a
boundary segment ("not flippable") or a non-2-bordered arc ("self-folded configuration, out
of scope").

`certify_flip_mutation(T, arc)` is the per-instance certificate that this combinatorial flip
agrees with cluster-algebra mutation **at the quiver level**. It builds the skew-symmetric
exchange matrix `B` of `quiver_of(T)` (the entry at row i, column j is
`#(i→j) − #(j→i)`), applies the Fomin–Zelevinsky mutation `μ_k` (`matrix_mutation`, with `k`
the position of `arc`), builds the exchange matrix of `quiver_of(flip(T, arc))` realigned on
the same vertex order, and returns whether the two matrices are equal. Full DWZ **potential**
right-equivalence under mutation is *not attempted* — it is deferred to P48.1; for the gentle
v1 scope the quiver-level certificate plus `is_gentle` on both sides is the shipped guarantee.

## The no-code surface (`block.py`)

`surface_block(T)` returns the surface invariants (`genus`, `boundary_marked`, `punctures`,
`arc_count`, `triangle_count`, `euler_characteristic`), the produced quiver (`vertices`,
`arrows`), the Jacobian's `relations` and `dim`, and `is_gentle`. The AG derived invariant is
attached **only when the Jacobian is gentle** (`ag_invariant`, Avella-Alaminos–Geiss, from
P46) — it is a derived invariant, provably *not* complete, so completeness is never claimed.
`surface_block` is not itself a webapp compute kind: surfaces are an input method, and the
produced algebra runs through the existing kinds. The free-form draw-a-surface canvas is a
named post-release successor (see the
[v0.2.0 GUI-deferral ledger](../verification.md#v020-gui-deferral-ledger)); v1 ships three
build-time presets (disc fan `A₃`, annulus `C(2,2)`, hexagon).

## The oracles

- **Literature** — the disc fan of the `(n+3)`-gon gives linear `Aₙ` (the orientation
  arbiter); the annulus `C(n, m)` gives the *acyclically* oriented affine `Ã_{n+m-1}` (no
  oriented cycle, so the P44 Jacobian is finite); the once-punctured torus is the loud-refusal
  oracle. These are Labardini/FST worked examples.
- **Self-cert** — the arc-adjacency, arc-count, triangle-count, and boundary-component
  certificates at `Triangulation` construction, and the flip ↔ FZ-matrix-mutation identity.
- **Cross-engine** — `is_gentle(jacobian_of(T))` holds across the disc/annulus/hexagon zoo
  (tying P44's Jacobian, P46's `is_gentle`, and this plan).
- **QPA** — QPA has no surface / triangulation constructor (an `IsBoundGlobal` sweep), so the
  crosschecks are at the resulting **gentle-algebra** level (`IsGentleAlgebra` /
  `IsSpecialBiserialAlgebra`), with a standing guard that fails if QPA ever ships one — a
  no-code surface *input* method is white space even in QPA.

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| marked-surface datum, FST arc count | `surfaces/marked.py` | `MarkedSurface`, `arc_count`, `triangle_count`, `in_v1_scope`, `_FST_EXCLUDED` |
| triangulations + constructors | `surfaces/triangulation.py` | `Triangulation`, `fan_triangulation`, `annulus_triangulation`, `hexagon_with_internal_triangle`, `once_punctured_torus` |
| quiver / potential / Jacobian | `surfaces/qp.py` | `quiver_of`, `potential_of`, `jacobian_of`, `_raw_arrows`, `_require_v1_scope` |
| flip ≡ FZ matrix mutation | `surfaces/flip.py` | `flip`, `certify_flip_mutation`, `exchange_matrix`, `matrix_mutation` |
| the `surface` input block | `surfaces/block.py` | `surface_block` |
