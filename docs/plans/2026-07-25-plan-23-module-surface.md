# Plan 23 — the module-theoretic surface (Tier 1b)

Date: 2026-07-25. Branch: `plan-23-module-surface`. Backlog: Tier 1b
(`docs/plans/DEEPER-ENGINES-BACKLOG.md`, Marco 2026-07-24).

**Vision (Marco):** every representation theorist can use quiverlab — specify a
module in the GUI without writing code, and read off the classical module-level
invariants. Tier 1b has **four** items that share ONE engine — the opposite
algebra `A^op` and the duality functor `D` — so they are planned together and
delivered in slices.

## Scope of this plan

All four Tier 1b items are **specified** here. This branch **implements** the
engine slice — items 1, 2, 4:

1. **AR translates τ / τ⁻** for f.d. right modules. *(delivered)*
2. **Injective resolutions + injective dimension.** *(delivered)*
3. **No-code module input (GUI + webapp).** *(specified only — the NEXT slice;
   checkbox stays unchecked.)*
4. **QPA (GAP) as the oracle** — `-m qpa` crosschecks for τ/τ⁻, projective /
   injective resolution terms & Betti numbers, injective dimension. *(delivered)*

Plus, per the 2026-07-25 verification-transparency standing constraint,
**theory/literature oracles** wherever a theorem or a worked example pins a
value (these run without the `[qpa]` extra).

## Mathematical design

### The shared engine: `A^op` and `D`

`A = kQ/I` over an exact `Domain`. Basis = trivial paths `e_v` and irreducible
paths (label `a*b*c`, read left-to-right, `target(a)=source(b)`, Assem–Simson–
Skowroński). Structure constants `T[i][j]` = coords of `b_i · b_j`.

**Opposite algebra** (`modules/opposite.py::opposite_algebra`). `A^op` is the
SAME `k`-space with reversed product `x ·^op y := y · x`. Realized as a
first-class `Algebra` with

- reversed quiver `Q^op`: arrow `a: s→t` becomes `a: t→s`;
- **transposed structure constants** `T^op[i][j] = T[j][i]` (same coordinate
  pages, index set unchanged);
- reversed basis labels: `e_v ↦ e_v`, `a*b*c ↦ c*b*a` (a valid path in `Q^op`);
- reversed relations (each word reversed) so `A^op`-modules validate against
  `I^op`.

Key: the index set is preserved (index `i` in `A` ↔ index `i` in `A^op`, same
vector, reversed label), which makes `D` and the corner transpose coordinate-
clean and avoids re-running Gröbner completion. `(A^op)^op` is `A` (cached
cross-link; label/`T`/quiver reversals are involutive). Correctness of
`T^op[i][j]=T[j][i]` with reversed labels: if `p,q` are paths, `reverse(q·p) =
reverse(p)*reverse(q)` in `Q^op`, matching `p ·^op q = q · p`. Verified by an
associativity/round-trip test and by module `check_module` passing over `A^op`.

**Duality `D = Hom_k(−,k)`** (`modules/duality.py::dualize`). A right `A`-module
`M` (column convention `m·b = action[b] @ m`, `action[x*y]=action[y]@action[x]`)
dualizes to a right `A^op`-module `DM = M*` with

    action_{DM}[reverse(ℓ)] = action_M[ℓ]^T   for every basis label ℓ.

Consistency with the anti-homomorphism convention: `x ·^op y = y · x`, so
`action_{DM}[y·x]=action_{DM}[y]@action_{DM}[x]` reduces to `R_{y·x}^T =
(R_x R_y)^T`, which holds since `R_{y·x}=R_x R_y`. `D` **preserves dimension
vectors** (`dim_v DM = rank action_M[e_v]^T = dim_v M`) but transposes/reverses
arrow actions. `D∘D = id` and `D` is contravariant (a map `M→N` dualizes to
`DN→DM`). `builders.injective(A,v) = I_v = D(A e_v)` was the only prior implicit
use of `D`; we test the explicit `D(projective(A^op,v)) ≅ injective(A,v)` (the
implicit-vs-explicit agreement — a bug if it fails).

### Transpose `Tr` and the AR translates

Minimal projective presentation `P_1 --d1--> P_0 --> M --> 0` (Plan 05
`modules/resolution.py::minimal_resolution`, minimal by construction — top-
generator covers, so `Tr` carries no spurious projective summands). Apply
`(−)* = Hom_A(−,A)` (contravariant, right `A`-mod → right `A^op`-mod):
`(P_v)* ≅ A e_v = e_v A^op = projective(A^op,v)`. For `φ_y: P_w→P_v` (left-mult
by `y ∈ e_v A e_w`), the transpose is `φ_y^* = ` left-mult by `ȳ` in `A^op`
(`ȳ ∈ e_w A^op e_v`) — right-mult by `y` on `A e_v → A e_w`.

    Tr M = coker( d1^* : ⊕_i P_{v_i}^{op} → ⊕_j P_{w_j}^{op} )   (right A^op-module)

Concretely (`modules/duality.py::transpose_module`): read the corner elements
`y_{ij} ∈ e_{v_i} A e_{w_j}` off the canonical-generator columns of `d1`; because
the `A ↔ A^op` index set is preserved, `y_{ij}`'s global coordinate vector IS
`ȳ_{ij}`'s. Build `d1^*` as the right-`A^op` map sending each source generator
`g_i` to `Σ_j ȳ_{ij}` (mirrors `projective_cover`'s cover assembly:
`basis·^op p ↦ action[p] @ (generator image)`), then take the cokernel
(`radtopsoc.quotient`).

    τ M  = D (Tr M)          τ⁻ M = Tr (D M)

Variance: `Tr: R-mod → R^op-mod`, `D: R-mod → R^op-mod`. `τM = D(Tr M)` is over
`(A^op)^op = A`; `τ⁻M = Tr(DM)` is over `(A^op)^op = A`. Public surface:
`M.tau()`, `M.tau_minus()` (plus `M.dualize()`, `M.transpose()`,
`Algebra.opposite()`).

**Honest gates (tested):**
- `τ(projective) = 0` and `τ⁻(injective) = 0` — the minimal presentation of a
  projective has `P_1 = 0` so `Tr = 0`; `D(injective)=projective over A^op` so
  `Tr(D I)=0`.
- `τ⁻ τ M ≅ M` **only for M with no projective summands** (indecomposable
  non-projective is the clean statement). Iso is certified by an exact
  invertible-hom witness (`modules/hom.py::is_isomorphic`, positive certificate,
  never a silent wrong answer).
- **Decomposability caveat (loud in docs):** `τ` works summand-wise; on a
  decomposable module `τ⁻τM` recovers only the non-projective part. The iso gate
  is asserted only for indecomposables. `is_isomorphic` decides over the base
  field (not the algebraic closure).
- Minimality is load-bearing: a non-minimal presentation gives `Tr ⊕
  projectives`; `minimal_resolution` guarantees minimality.

### Injective resolutions + injective dimension

Dual of Plan 05's projective machinery on the same op+D engine
(`modules/injective.py`):

    E(M) = D(projective cover of DM over A^op)          (injective envelope)
    E^n  = D(P_n)   where P_• = min. proj. resolution of DM over A^op
    inj.dim_A(M) = pd_{A^op}(DM)

`injective_resolution(M, length)` dualizes the projective resolution of `DM`
term-by-term (`E^n` is a sum of `I_u`, `betti(n)` = # injective summands, term
dim-vectors = those of `P_n` since `D` preserves dim vectors);
`injective_dimension(M, bound)` = `DM.projective_resolution(bound).pd()` (int, or
`None` when unresolved within `bound` = infinite, e.g. self-injective
non-projective). Public: `M.injective_resolution(length)`,
`M.injective_dimension(bound=32)`.

### Item 3 — no-code module input (SPECIFIED, not built here)

The constructor exists (`modules/module.py::Module.from_arrow_action`: dim vector
+ one exact matrix per arrow, loud `RelationError` on relation violation). The
NEXT slice exposes it with zero code:
- **webapp** (Plan 09 server): request schema **v2** gains a `module` block
  `{dims: {v: n_v}, maps: {arrow: [[…]]}}` + module compute kinds (rad/top/soc
  from `radtopsoc.py`, Ext from `modules/ext.py`, and τ/τ⁻ + both resolutions
  from this plan); versioned bump, served by BOTH tiers.
- **GUI** (Pyodide): per-vertex dimension picker + per-arrow matrix grid, and
  zero-typing pick-lists for `S(v)/P(v)/I(v)` (builders exist).
The engine slice here is exactly what those surfaces call; the graded-form
emitter (`modules/qpa_module.py::graded_form`) added for the QPA oracle is also
the natural serialization for the schema.

## Item 4 — QPA / theory oracles

**QPA idioms confirmed live** (`[qpa]` extra, libgap single-statement, per-line
eval): `DTr`, `TrD`, `DualOfModule`, `RightModuleOverPathAlgebra(A, dimvec,
[[arrow, matrix], …])` (arrow matrix is `dim_s × dim_t`, ROW convention),
`DimensionVector`, `IsomorphicModules`, `ProjectiveResolution` +
`ObjectOfComplex`, `InjDimensionOfModule(M, bound)` (returns `false` when inj.dim
> bound), `GlobalDimensionOfAlgebra(A, bound)`, `IsSelfinjectiveAlgebra`. New
crosscheck kinds in `qpa/crosscheck.py` + `qpa/scripts.py`:

- `A.crosscheck("tau", M, minus=False)` / `("tau_minus", M)` — τ/τ⁻ dimension
  vectors vs `DTr`/`TrD`, AND iso class via translating our module to QPA
  (`graded_form` → `RightModuleOverPathAlgebra`, transpose to the row
  convention) then `IsomorphicModules`.
- `A.crosscheck("proj_resolution", M, top)` — resolution term dim-vectors +
  Betti (top-of-syzygy) vs QPA `ProjectiveResolution`.
- `A.crosscheck("inj_dimension", M, bound)` — vs `InjDimensionOfModule`
  (`None ↔ false`).
- `A.crosscheck("inj_resolution", M, top)` — our `E^n` dim-vectors vs QPA proj
  resolution of `DualOfModule(M)` over the opposite algebra (the defining
  identity `inj.res_A M ↔ proj.res_{A^op} DM`).

Zoo: kA_2, kA_3, commutative square, cyclic Nakayama `kZ_3/rad²`
(self-injective), and the Plan-18 multi-vertex `line_abc_cde`. QPA field scope =
QQ / prime GF(p) (inherited from `qpa/scripts.py`).

**Theory / literature oracles (plain `-m deep`, no extra):**
- **Hereditary Coxeter transformation** (ASS Cor. IV/VIII): for kA_n and `M`
  non-projective indecomposable, `dim τM = Φ^{-T}·dim M`; non-injective,
  `dim τ⁻M = Φ^{T}·dim M`, where `Φ = A.coxeter_matrix()`. (quiverlab's Coxeter
  uses the `e_iAe_j` Cartan convention, so the transform is `Φ^{-T}` — an
  independent path-counting cross-check of the op+D+Tr τ.)
- kA_2 / kA_3 explicit τ/τ⁻ tables (worked AR quiver).
- **Self-injective ⟺ inj.dim 0 for all modules; inj.dim 0 ⟺ projective**
  (k[x]/(x²)): non-projective ⇒ `injective_dimension = None` (infinite);
  projective ⇒ 0.
- **Nakayama τ-orbits**: `kZ_n/rad²` — τ permutes the simples cyclically.
- **inj.dim vs gl.dim on the commutative square** (gl.dim 2): `max_v inj.dim
  S_v = gl.dim`.
- `(A^op)^op ≅ A`, `D∘D ≅ id`, `D(projective(A^op,v)) ≅ injective(A,v)`
  (implicit-vs-explicit `D`).

## Acceptance

- New src: `modules/opposite.py`, `modules/duality.py`, `modules/injective.py`,
  `modules/qpa_module.py`; extensions to `modules/hom.py` (`is_isomorphic`),
  `modules/module.py` (public methods), `core/algebra.py` (`opposite()`),
  `qpa/scripts.py` + `qpa/crosscheck.py`. No floats in `src/` (AST gate). Domain-
  generic (GF(p), QQ, CC).
- Tests: `tests/modules/test_opposite.py`, `test_duality_tau.py`,
  `test_injective.py`, `test_module_iso.py` (deep bucket);
  `tests/qpa/test_module_ar_crosscheck.py` (qpa bucket).
- Docs: `docs/internals/` module chapter updated; `CLAUDE.md`, ROADMAP, backlog
  status lines (τ/τ⁻ + injective + QPA oracle delivered; GUI/webapp item stays
  unchecked).
- **Verification page (2026-07-25 constraint):** `docs/verification.md` is being
  created on a concurrent branch; Plan 23's oracles (this section) MUST be added
  to it at merge time — the module-surface subsystem row maps to the test files
  above.

## Deferred / out of scope

- Item 3 (no-code GUI/webapp) — next slice.
- Native AR quiver enumeration (Tier 2 non-goal; `[qpa]` covers it).
- Distinct-module Ext(M,N), M≠N (Plan-05-flagged QPA extension) — untouched.
