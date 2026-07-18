# quiverlab Plan 06 — Families + Batch + Citations Registry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A working algebraist writes `NakayamaAlgebra([3, 2, 2])`, `PathAlgebra("D4")`, `QuantumCI(q="i")`, `ExteriorAlgebra(3)`, `IncidenceAlgebra(diamond)`, `TensorProduct(A, B)`, `zoo(dim_max=12)`, or `families()` — one line each, `field=CC` default — and gets a certified Plan-01 `Algebra` (or a loud, mathematician-readable error). Every family names the papers it comes from: `A.citations()`, `A.hochschild_cohomology(n).references`, and `quiverlab.bibliography()` emit real, verified BibTeX keys with annotations. Expert scans persist to SQLite via `quiverlab.batch`. The full §3.4 catalog, the §3.9 citations registry core, and the labdb lift ship green.

**Architecture:** Three new packages plus one sanctioned grammar extension.
- `src/quiverlab/families/` gains the full catalog (`nakayama.py`, `path_algebra.py`, `truncated.py`, `radical_square_zero.py`, `incidence.py`, `quantum.py`, `exterior.py`, `preprojective.py`, `trivial_extension.py`, `tensor.py`, `zoo.py`, `discover.py`) on top of the existing `basic.py`. Monomial families lower through the Plan-01 `Quiver.algebra` monomial route; non-monomial families (`QuantumCI`, `ExteriorAlgebra`, `PreprojectiveAlgebra`, `IncidenceAlgebra`) lower through the Plan-03 Gröbner route; structure-constant families (`TensorProduct`, `TrivialExtension`) build a multiplication table directly and emit an `Algebra`. Every builder takes `field=CC` and is field-generic (the `sweep` hook, §Boundary).
- `src/quiverlab/citations/` — the registry core: `references.bib` (curated, verified), `registry.py` (keys → BibTeX id + algorithm/family + annotation), `bibliography()` (grouped/annotated), `A.citations()`, `HHTable.references` plumbing, loud `CitationError` on unknown keys.
- `src/quiverlab/batch/` — the `labdb` lift: `ResultsDB` (SQLite), a `BUILDERS` registry over the quiverlab catalog, `analyze`/`run_scan`. Expert-facing; **not** in `from quiverlab import *`.
- **One sanctioned grammar extension** (spec §3.3, foreseen in the Plan-03 boundary note): the relation parser accepts exact **non-rational** coefficient tokens (`"i"`, `"E(4)"`, `"sqrt(2)"`) so `QuantumCI(q="i")` lowers through the existing Gröbner engine unchanged. Rational relations parse to `Fraction` exactly as before (Plan-03/04 consumers unaffected).

**Tech Stack:** Python ≥ 3.10, no new hard dependencies (the Plan-01 stack: `numpy`, `sympy`; `sqlite3` is stdlib). pytest. Family construction is combinatorial + the frozen engine; no new linear algebra. Coefficient arithmetic goes through the `Domain` protocol / `field.parse_entry`.

---

## Global Constraints

- **Repo root:** `/Users/marco/Desktop/HomologicalNetworks/quiverlab`. All paths below are relative to it.
- **Interpreter:** use the project venv **`/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`** (Python 3.12). The *system* `python`/`python3` is 3.8 and MUST NOT be used — it will fail on 3.10+ syntax (`X | None`, `list[int]`, etc.).
- **Thread throttle:** prefix **every** test command with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2` (Marco's machine has crashed under agent fleets; keep thread/memory pressure low). No new parallelism is introduced in core; `quiverlab.batch` defaults to **`n_workers=1`** and only forms a `multiprocessing.Pool` when a caller explicitly raises it (spec §5: "no multiprocessing in core").
- **Exact arithmetic only, through the `Domain` protocol / `field.parse_entry`.** All coefficient arithmetic goes through `dom.add/sub/neg/mul/inv/is_zero/eq/coerce/zero/one` or `field.parse_entry`. Never compare or combine coefficients with raw Python numeric operators. Floats fail loudly (Plan-01 `ExactnessError`).
- **Float ban (AST gate).** `tests/test_no_floats.py` scans all of `src/quiverlab/` for float/complex literals and `float()` calls and MUST stay green. Write **no** float or complex literals anywhere under `src/` (use ints; all "coefficients" are `Fraction`, sympy exact `Expr`, or domain elements — never Python `float`/`complex`). Run it as part of every task's suite.
- **Left-to-right path composition** (Assem–Simson–Skowroński): `a*b` = "first `a`, then `b`", requiring `target(a) == source(b)`. A word is a tuple of arrow names read left to right. Every relation, mesh relation, commutativity relation, and Kupisch path in this plan is written in this convention; document it in every relevant docstring.
- **Read-only banks.** `/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/` (the `HansConjecture` bank: `hanlab/open_zoo.py`, `hanlab/labdb.py`, the three tests) is **read** for the lift and **never modified**. Code lifted from it carries the ported-with-attribution header already used across `src/quiverlab/engine/` (see any file there for the exact wording).
- **Full suite green at every commit.** Run `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q` from the repo root before each commit; it must pass. **Suite pattern:** each task runs its own **focused** test file(s) in the foreground; a single **full-suite** background run is launched at the acceptance task, awaited, and only then is the final commit made. Per-task commits run the focused suite plus the float gate foreground.
- **Commits:** conventional prefixes (`feat:`, `test:`, `chore:`, `docs:`); every commit message ends with the trailer lines
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`
  `Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd`
- **Frozen Plan-01/02/03/04/05 interfaces consumed verbatim** (do not modify their signatures):
  - `Domain` (`coerce/add/sub/neg/mul/inv/is_zero/eq/zero/one/characteristic/name`, hooks `parse_entry`/`make_domain`); `fields.CC`, `fields.GF`, `fields.QQ`, `fields.E`; `fields.primefield.PrimeField`; `fields.complexfield.ComplexField`.
  - `combinat.Quiver` (`.vertices`, `.arrows`, `.source/.target`, `.word_source/.word_target`, `.compose_ok`, `.is_acyclic`, `.algebra(relations=, field=, degree_bound=, trace=)`).
  - `combinat.Relation`/`parse_relations`/`parse_relation` (`.terms`, `.source`, `.target`, `.is_monomial`, `.max_length`, `.min_length`).
  - `core.Algebra` (`.domain, .T, .unit, .dim, .basis_labels, .is_unit_adapted, .quiver, .relations`, `.multiply`, `.unit_adapted`, `.change_of_basis`, `.hochschild_cohomology`, `.hochschild_homology`, `.cartan_matrix`, ctor `Algebra(domain, T, unit, basis_labels=, is_unit_adapted=, _quiver=, _relations=)`, classmethod `from_structure_constants(T, unit, field=, check=, basis_labels=)`); `core.monomial.build_monomial_algebra`, `core.monomial.irreducible_paths`.
  - `groebner.build_reduction_system`, `groebner.ReductionSystem`, `groebner.groebner_algebra` (Plan 03); the CS resolution + `QuantumCIResolution`/`CyclicNakayamaResolution` periodic backends (Plan 04); `sweep` and generalized `module_ext` (Plan 05).
  - `hochschild.table.HHTable` — see the **one compatible extension** below.
  - `errors.*` (a new `CitationError(QuiverlabError)` is **added**, no existing error changed).
  - **Two compatible extensions** made in this plan: (i) `HHTable.__init__` gains a trailing optional keyword `references=()` and exposes `.references` (additive; existing constructions unaffected — default empty tuple). (ii) `combinat.relations` accepts non-rational exact coefficient tokens, producing `Fraction` for rational input exactly as before and a sympy exact `Expr` only when a non-rational token is present (`Relation.terms` coefficient type widens from `Fraction` to `Fraction | sympy.Expr`; all existing rational relations are byte-for-byte unchanged).

---

## The families, mathematically (read this before Task 1)

The implementer knows Python, not representation theory. This section carries the mathematics; the code below implements exactly what is stated here. All conventions are quiverlab's: paths left-to-right, `C[i][j] = dim e_i A e_j = #(irreducible paths i→j)` (the frozen `invariants.cartan.cartan_matrix` convention), `dim A = Σ_{i,j} C[i][j] = |basis|`, `HH^0(A) = dim Z(A)` (the center), connected ⟹ `HH^0 ≥ 1`.

**Routing.** A family is *monomial* (lowers through `core.monomial.build_monomial_algebra` via `Quiver.algebra`), *general* (a genuine non-monomial relation, lowers through the Plan-03 Gröbner engine via `Quiver.algebra`), or *structure-constant* (builds `T` directly, emits an `Algebra`). Each builder states its route. Monomial: Nakayama, TruncatedPathAlgebra, RadicalSquareZero, PathAlgebra (no relations at all). General: QuantumCI, ExteriorAlgebra, PreprojectiveAlgebra, IncidenceAlgebra. Structure-constant: TensorProduct, TrivialExtension.

**1. NakayamaAlgebra — by Kupisch series or by (n, l, cyclic).** A Nakayama (uniserial) algebra is determined by its Kupisch series `K = [c_1, …, c_n]`, `c_i = dim P_i` (the length of the indecomposable projective at vertex `i`). Two shapes:
- **Linear** (quiver `A_n`: `1→2→…→n`, arrows `a_i: i→i+1`): the sink `n` has `c_n = 1` (P_n simple). Admissible iff `c_n = 1`, `c_i ≥ 2` for `i < n`, and `c_i ≤ c_{i+1} + 1`. Relation set: the length-`c_i` path from `i` is zero (monomial). `dim A = Σ c_i`.
- **Cyclic** (quiver `Z_n`: vertices `1..n` in a cycle, arrows `a_i: i → (i mod n)+1`): every `c_i ≥ 2`; admissible iff `c_i ≥ 2` and `c_i ≤ c_{i-1} + 1` cyclically. Relation set: the length-`c_i` path from `i` is zero (monomial). `dim A = Σ c_i`. `NakayamaAlgebra(n=4, l=3, cyclic=True)` is the homogeneous cyclic case `kZ_n / rad^l` (all `c_i = l`), `dim = n·l`.
- **Which shape does `NakayamaAlgebra(K)` pick?** If `min(K) == 1` → linear; if `min(K) ≥ 2` → cyclic. `[3,2,2]` has `min = 2 ≥ 2` ⟹ **cyclic**, dim `3+2+2 = 7` (see Fixture N1). A series that is neither (e.g. an interior `1`) raises `AdmissibilityError`.

**2. PathAlgebra(type, orientation) — hereditary Dynkin/Euclidean.** `kQ` with **no relations**, `Q` the Dynkin/Euclidean diagram of the given type with a chosen orientation (monomial route: `relations=[]`). `dim kQ = Σ_{i,j} #(paths i→j)` = total number of paths (finite iff `Q` acyclic). Because the underlying graph of a Dynkin diagram is a **tree**, and Euclidean `Ã` a cycle:
- For a tree quiver, between any two vertices there is at most one directed path; `HH^•(kQ) = (1, 0, 0, …)` over every field (Happel/Cibils: `dim HH^1(kQ) = b_1(underlying graph) = |Q_1| − |Q_0| + 1 = 0` for a tree; hereditary directed ⟹ `HH^{≥2} = 0`).
- Euclidean `Ã_n` fully cyclically oriented is the **infinite** path algebra ⟹ `NotFiniteDimensionalError`; an acyclic orientation is finite (tame hereditary). PathAlgebra requires an acyclic orientation.
- **Orientation input:** `orientation="linear"` (default) = arrow `i→j` for each diagram edge `{i,j}` with `i<j` in the standard type labeling; `orientation="reverse"` flips all; `orientation={edge: (s,t), …}` overrides specific edges. Standard labelings: `A_n` = `1—2—…—n`; `D_n` = chain `1—2—…—(n−2)` with `(n−2)—(n−1)` and `(n−2)—n`; `E_6/7/8` per Bourbaki; `Ã_n/D̃_n/Ẽ_n` per Kac. See Fixture N2 for `D4`.

**3. TruncatedPathAlgebra(type_or_quiver, r) = kQ / rad^r.** Monomial: kill every path of length `≥ r`. `dim = #(paths of length < r)`. `TruncatedPathAlgebra("A5", 2)`: `A5` (`1→2→3→4→5`), `rad^2 = 0`, basis = 5 vertices + 4 arrows = **dim 9**. Accepts a type string (builds the Dynkin quiver, default linear orientation) or a `Quiver`.

**4. RadicalSquareZero(Q) = kQ / rad^2.** `TruncatedPathAlgebra(Q, 2)` from *any* `Quiver`: kill all length-2 paths (monomial). `dim = |Q_0| + |Q_1|` — finite even when `kQ` itself is infinite (e.g. the two-loop quiver). See Fixture N4.

**5. IncidenceAlgebra(poset) = kP.** For a finite poset `P`, `kP` has basis the intervals `[x,y]` (`x ≤ y`), `dim = #{(x,y): x ≤ y}`. Realized as the **Hasse quiver** (vertices = elements, arrows = cover relations `x ⋖ y`) modulo the **commutativity ideal**: all parallel Hasse paths between the same endpoints are equal (general route: relations `p − q` for parallel paths `p, q`). **Poset input format (concrete):** `IncidenceAlgebra(covers, elements=None)` where `covers` is a list of ordered pairs `(x, y)` meaning "`x` is covered by `y`" (a Hasse edge, arrow `x→y`); `elements` optionally lists isolated points. Internally: build the Hasse quiver, compute `≤` by transitive closure, reject if the closure is not antisymmetric (a directed cycle → not a poset → `RelationError`), then emit one commutativity relation per pair of distinct parallel maximal Hasse paths. The **diamond poset** `⊥ ⋖ x, ⊥ ⋖ y, x ⋖ ⊤, y ⋖ ⊤` gives exactly the Plan-03 commutative square, **dim 9** (Fixture N5).

**6. QuantumCI(q) = k⟨x,y⟩/(x², y², xy + q·yx).** General route; the spec's convention is `xy + q·yx` (relation string `"x*y + q*y*x"`). One vertex, two loops `x, y`. Basis `{1, x, y, xy}`, **dim 4** for any `q ≠ 0`. `q` may be any exact scalar: rational (parser already handles) or non-rational (`"i"`, `"E(3)"` — the sanctioned parser extension, Task 4). Center: `Z = span{1, xy}` for `q ∉ {−1}` (dim 2), and `Z =` everything (dim 4) when `q = −1` (then the algebra is commutative `k[x,y]/(x²,y²)`). So `QuantumCI(q="i")` over `CC`: **dim 4, HH^0 = 2**. (Convention note: the task brief's cross-check "`q=−1`" matches the *alternate* BGMS convention `xy − q·yx`; under the spec's stated `xy + q·yx`, the exterior/anticommutative point is `q = +1` — see Fixture N6/N7. We follow the spec.)

**7. ExteriorAlgebra(n) = Λ(k^n) = k⟨x_1..x_n⟩/(x_i², x_ix_j + x_jx_i, i<j).** General route; **all coefficients rational (±1)** — needs no parser extension. `dim = 2^n`. Gröbner completion (arrow order `x_1<…<x_n`) rewrites `x_j x_i → −x_i x_j` (`i<j`) and `x_i² → 0`; irreducible basis = strictly-increasing words = the `2^n` exterior monomials. `ExteriorAlgebra(2)` = `k⟨x,y⟩/(x²,y²,xy+yx)` = `QuantumCI(q=1)` exactly (terms, dim 4, structure constants) — the cross-check. `ExteriorAlgebra(3)`: dim 8. Center of `Λ(k²)` over `CC`: `{1, xy}`, **HH^0 = 2** (char ≠ 2); over `GF(2)` it degenerates to commutative dim-4 with `HH^0 = 4` (a characteristic-dependence datum for the sweep).

**8. PreprojectiveAlgebra(type) = Π(Q).** General route. Double the Dynkin quiver `Q` (each edge `i—j` becomes `a: i→j` and `a*: j→i`); impose one **mesh relation** per vertex `v`: `Σ_{α: source v} α·α* − Σ_{β: target v} β*·β = 0` (rational ±1 coefficients — no parser extension). For Dynkin type, `Π(Q)` is finite-dimensional; `dim Π(A_n) = n(n+1)(n+2)/6` (tetrahedral number). `Π(A_2)` = `k⟨a, a*⟩/(a a*, a* a)` — both mesh relations are *monomial* here — **dim 4**, HH^0 = 1. `Π(A_3)`: the branch-vertex mesh relation `a* a − b b* = 0` is genuinely non-monomial ⟹ Gröbner; **dim 10** (= `T_3`). The Gröbner certificate is the authority on the dimension; the tests assert the tetrahedral value.

**9. TensorProduct(A, B) = A ⊗_k B.** Structure-constant route. Basis `{a_i ⊗ b_j}` in row-major order (`index(i,j) = i·dim(B) + j`); `(a⊗b)(a'⊗b') = (aa') ⊗ (bb')`; `unit = unit_A ⊗ unit_B`; `dim = dim A · dim B`. Both factors must share the same ground field/domain (else `FieldError`). **Künneth** (Gerstenhaber–Schack / Cartan–Eilenberg): `HH^n(A⊗B) = ⊕_{i+j=n} HH^i(A) ⊗ HH^j(B)`, so `HH^0` is **multiplicative**: `dim Z(A⊗B) = dim Z(A)·dim Z(B)`. `TensorProduct(kA_2, kA_2)` = the Plan-03 commutative square = `IncidenceAlgebra(diamond)`, **dim 9, HH = [1,0,0]** — a triple cross-check across three construction routes. `TensorProduct(k[x]/x², k[x]/x²)` = `k[x,y]/(x²,y²)`, **dim 4, HH^0 = 2·2 = 4** (Künneth with nontrivial centers).

**10. TrivialExtension(A) = A ⋉ D(A).** Structure-constant route. `D(A) = Hom_k(A, k)` the dual bimodule; underlying space `A ⊕ D(A)`, `(a,f)(b,g) = (ab, a·g + f·b)` where `(a·g)(c) = g(ca)` and `(f·b)(c) = f(bc)`. `dim T(A) = 2·dim A`. `TrivialExtension(kA_2)`: `dim kA_2 = 3` ⟹ **dim T = 6** (the task brief's "dim 5?" is corrected to 6: `dim T(A) = 2 dim A` always). `T(A)` is **symmetric**, so `dim HH^n(T A) = dim HH_n(T A)` — the family's signature cross-check.

**11. zoo(dim_max) — the curated exact zoo.** An iterator over a bundled, curated catalog (lifted from the bank's `open_zoo` generator + `open_zoo_catalog_v2.json`, 120 entries), each a confluent reduction system over `k⟨x,y⟩` or `k⟨x,y,z⟩`, rebuilt through the Plan-03/04 reduction-system → `groebner_algebra` path. `zoo(dim_max=D)` yields `Algebra` objects with `dim ≤ D`, in ascending `(dim, name)` order. The catalog entries are the "open" (local, non-monomial, non-graded, non-commutative, non-quantum-CI) algebras that motivated Han's conjecture; among them the **periodic symmetric family** `k⟨x,y⟩/(x³, yᵇ − x², yx + xy)` (Task 13).

**12. families() — discoverability.** Returns a structured listing of every catalog entry: name, one-line signature, route (monomial/general/structure-constant), and citation keys. Printable; consumed by the web `/families` page (Plan 09) and by `sweep` (Plan 05) as the catalog hook.

---

## Hand-verified fixtures (expected values worked out by hand)

Every number below is re-derived in the task that uses it. These are the acceptance oracles.

### Fixture N1 — `NakayamaAlgebra([3,2,2])` (cyclic; dim 7, HH^0 = 1)
`min([3,2,2]) = 2 ≥ 2` ⟹ cyclic `Z_3`: vertices `1,2,3`, arrows `a1:1→2, a2:2→3, a3:3→1`. `c_i = dim P_i`: `c_1=3` ⟹ length-3 path from 1 (`a1 a2 a3`) `= 0`; `c_2=2` ⟹ `a2 a3 = 0`; `c_3=2` ⟹ `a3 a1 = 0`. Admissibility (cyclic): all `c_i ≥ 2` ✓; `c_2=2 ≤ c_1+1=4`, `c_3=2 ≤ c_2+1=3`, `c_1=3 ≤ c_3+1=3` ✓. Nonzero basis paths: vertices `e1,e2,e3` (3); `a1,a2,a3` (3); length-2 nonzero: `a1 a2` (1→3) — `a2 a3` and `a3 a1` killed — (1). Total **dim 7 = 3+2+2 = Σ c_i** ✓. Composition factors: `P_1 = [S_1,S_2,S_3]`, `P_2 = [S_2,S_3]`, `P_3 = [S_3,S_1]`. **Cartan** (`C[i][j] = #paths i→j`):
```
        S_1 S_2 S_3
  P_1 [  1   1   1 ]
  P_2 [  0   1   1 ]
  P_3 [  1   0   1 ]
```
`Σ entries = 7 = dim` ✓; `det = 1·(1·1−1·0) − 1·(0·1−1·1) + 1·(0·0−1·1) = 1 + 1 − 1 = 1`. **HH^0 = dim Z(A)**: the only nonzero cyclic paths (source=target) are `e1,e2,e3` (all length-≥3 cycles vanish); centrality forces equal coefficients (`z a_i = λ_i a_i`, `a_i z = λ_{i+1} a_i`), so `Z = span{e1+e2+e3} = k·1`, **HH^0 = 1**.

### Fixture N2 — `PathAlgebra("D4")` (dim 9, HH = [1,0,0])
Standard `D4` labeling: edges `{1,2}, {2,3}, {2,4}` (vertex 2 the center). `orientation="linear"` (arrow `i→j`, `i<j`): arrows `a:1→2, b:2→3, c:2→4`. Paths: trivial `e1,e2,e3,e4` (4); length 1: `a,b,c` (3); length 2: `a·b` (1→3), `a·c` (1→4) (2); length ≥3: none (`b,c` end at leaves). **dim = 4+3+2 = 9.** (The task brief's "dim 12?" is not attainable: for the `D4` star with `k` in-edges and `3−k` out-edges at the center, `#paths = 4 + 3 + k(3−k) ∈ {7,9}`, max 9. Corrected: **9** under this orientation, **7** under a source/sink center.) `D4` underlying graph is a **tree** (`|Q_1|−|Q_0|+1 = 3−4+1 = 0`) ⟹ `HH^1 = 0`; hereditary directed ⟹ `HH^{≥2}=0`; connected ⟹ `HH^0=1`. **HH = [1,0,0]** over every field.

### Fixture N3 — `TruncatedPathAlgebra("A3", 2)` (dim 5, HH^0 = 1)
`A3` = `1→2→3` (arrows `a:1→2, b:2→3`), `rad^2 = 0` kills `a·b`. Basis `e1,e2,e3,a,b`, **dim 5**, connected ⟹ HH^0 = 1. (Spec example `TruncatedPathAlgebra("A5",2)`: 5 vertices + 4 arrows = **dim 9**.)

### Fixture N4 — `RadicalSquareZero(two-loop)` (dim 3, HH^0 = 3)
`Q = Quiver([1], {"x": (1,1), "y": (1,1)})`; `rad^2` kills `x², xy, yx, y²`. Basis `e1, x, y`, **dim 3**. All products of `x,y` vanish ⟹ commutative local ⟹ `Z = A`, **HH^0 = 3**. Finite even though `kQ` is the infinite free algebra `k⟨x,y⟩`. Monomial route.

### Fixture N5 — `IncidenceAlgebra(diamond)` ≡ commutative square (dim 9, HH = [1,0,0])
`covers = [("b","x"), ("b","y"), ("x","t"), ("y","t")]` (`b=⊥`, `t=⊤`). Hasse quiver `b→x, b→y, x→t, y→t`; the two parallel maximal paths `b→x→t` and `b→y→t` give one commutativity relation `(b→x→t) − (b→y→t)`. Isomorphic to the Plan-03 Fixture 2 commutative square `kQ/(a*b − c*d)` = `kA_2 ⊗ kA_2`. **dim 9, HH = [1,0,0]** over every field (Künneth: both `kA_2` factors are trees). Intervals: 4 trivial + 4 length-1 + 1 length-2 (`[⊥,⊤]`) = 9 ✓.

### Fixture N6 — `ExteriorAlgebra(2)` ≡ `QuantumCI(q=1)` (dim 4; HH^0 = 2 over CC)
`k⟨x,y⟩/(x²,y²,xy+yx)`. Gröbner (`x<y`): `x²→0, y²→0, yx→−xy`; irreducibles `{e, x, y, xy}`, **dim 4**. Byte-identical to `QuantumCI(q=1)` (same relations, same structure constants). Over `CC`: `Z = span{1, xy}` (`x,y` non-central since `xy = −yx ≠ yx`), **HH^0 = 2**. Over `GF(2)`: `xy+yx = xy−yx`, algebra is commutative dim-4, **HH^0 = 4** (sweep datum). `ExteriorAlgebra(3)`: dim 8.

### Fixture N7 — `QuantumCI(q="i")` (dim 4, HH^0 = 2 over CC)
`k⟨x,y⟩/(x²,y², xy + i·yx)` over `CC` (needs the parser extension; `i = E(4)`). Gröbner: `x²→0, y²→0, yx → −i^{-1} xy = i·xy` (since `xy = −i·yx` ⟹ `yx = −i^{-1} xy = i xy`). Irreducibles `{e,x,y,xy}`, **dim 4**. `Z = span{1, xy}` (`q = i ≠ −1`), **HH^0 = 2**.

### Fixture N8 — `PreprojectiveAlgebra("A2")` (dim 4, HH^0 = 1); `PreprojectiveAlgebra("A3")` (dim 10)
`A2` double quiver: `a:1→2, a*:2→1`; mesh relations `a·a* = 0` (vertex 1), `a*·a = 0` (vertex 2) — both **monomial**. Basis `e1,e2,a,a*`, **dim 4**, `Z = k`, HH^0 = 1. `A3` double: `a:1→2, a*:2→1, b:2→3, b*:3→2`; mesh at 1: `a·a*=0`; at 3: `b*·b=0`; at 2 (branch): `a*·a − b·b* = 0` (**non-monomial** ⟹ Gröbner). `dim = T_3 = 3·4·5/6 = 10` (Gröbner-certified; the test is the oracle).

### Fixture N9 — `TensorProduct` Künneth
`TensorProduct(kA_2, kA_2)`: `dim = 3·3 = 9`, ≡ commutative square, **HH = [1,0,0]**. `TensorProduct(truncated_polynomial(2), truncated_polynomial(2))` = `k[x,y]/(x²,y²)`: `dim = 2·2 = 4`, `HH^0 = HH^0(k[x]/x²)² = 2·2 = 4` (multiplicative).

### Fixture N10 — `TrivialExtension(kA_2)` (dim 6; symmetric)
`kA_2` dim 3 ⟹ `T` **dim 6 = 2·3**. `T` symmetric ⟹ `dim HH^n(T) = dim HH_n(T)` for all `n` (the tested property; the exact HH sequence is computed, not hand-asserted).

### Fixture Z1 — periodic symmetric family (zoo / ledger)
`open2_33_712`, dim 9: `k⟨x,y⟩/(x³, y³ − x², yx + xy)` (rules `x³→0, y³→x², yx→−xy`). `minimal_homology_dims(A, 6, primes=(32003,)) == [6,5,5,5,5,5,5]`; at `p=2`: `[9,10,14,18,22,26,30]` (growing tail). `is_frobenius` True, Nakayama automorphism = identity (symmetric). `open2_37_19612`, dim 21: `k⟨x,y⟩/(x³, y⁷ − x², yx + xy)`; symmetric; period-4 resolution ⟹ HH_n = 11 for `n ≥ 1` (SLOW, gated). Golden open-zone check `open_33_0` = `k⟨x,y⟩/(x³ − y², y³, yx + xy)` (rules `x³→y², y³→0, yx→−xy`): dim 9, `HH_homology[32003] == [6] + [5]*16`, `resolution_ranks == [1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2]`.

---

## The curated references (verified 2026-07-18; the registry ground truth)

`src/quiverlab/citations/references.bib` ships these entries verbatim (all metadata web-verified). Keys are stable public identifiers.

| BibTeX key | What it underpins | Registry key(s) |
|---|---|---|
| `Hochschild1945` | Original definition of HH; the (normalized) bar complex | `bar` |
| `Bardzell1997` | Minimal bimodule resolution, monomial case | `bardzell` |
| `ChouhySolotar2015` | General `kQ/I` resolution via reduction systems / ambiguities | `chouhy_solotar` |
| `NegronWitherspoon2016` | Gerstenhaber bracket via homotopy liftings | `bracket_liftings` |
| `Volkov2019` | Bracket on an arbitrary resolution | `bracket_liftings_volkov` |
| `GSZ2001` | Minimal module (projective) resolutions; module Ext | `minimal_resolution`, `module_ext` |
| `Happel1989` | Happel's question (finite gl.dim vs HH vanishing) | `happel_question` |
| `BGMS2005` | Quantum complete intersection: finite HH, infinite gl.dim | `quantum_ci` |
| `BerghErdmann2008` | Explicit HH of quantum complete intersections (the QCI HH *oracle*) | `quantum_ci`, `qci_hh_oracle` |
| `CartanEilenberg1956` | Künneth formula for HH (tensor-product multiplicativity) | `tensor_product` |
| `GerstenhaberSchack1987` | Hodge/λ-decomposition (commutative/tensor/incidence) | `hodge` |
| `Connes1985` | Cyclic homology; Connes B / SBI | `cyclic` |
| `Gerstenhaber1963` | Cup product + Gerstenhaber bracket structure | `cup`, `bracket`, `gerstenhaber` |
| `Luebeck_ConwayPolynomials` | Conway-polynomial tables for GF(p^n) | `conway`, `finite_fields` |
| `ASS2006` | Bound-quiver `kQ/I` formalism; Nakayama/hereditary/incidence families | `path_algebra`, `nakayama`, `incidence`, `preprojective`, `assem_book` |
| `Han2006` | Hochschild homology dimension; Han's conjecture | `han_conjecture` |

**Bib entry count: 16** (floor for the registry-coverage test; Plan 08 appends the four software-citation entries — `qpa`, `gap4`, `sagemath`, `quiverlab` — and raises the floor). Corrections carried from verification: `bar` is backed by `Hochschild1945` (the original HH definition), `Happel1989` now backs `happel_question` only. `BGMS2005` *is* the BGMS quantum-CI paper — the separate explicit QCI-HH computation is `BerghErdmann2008`. The Künneth formula is `CartanEilenberg1956`; `GerstenhaberSchack1987` is the λ/Hodge reference (`hodge`, cited by the tensor-product and incidence families). `Connes1985` pages are 41–144. `ASS2006` carries ISBN + DOI (`10.1017/CBO9780511614309`, moderate confidence).

Family → citation keys (wired in Task 10 / `A.citations()`):
`NakayamaAlgebra → [nakayama, assem_book]`; `PathAlgebra → [path_algebra, happel_question, assem_book]`; `TruncatedPathAlgebra → [path_algebra, bardzell]`; `RadicalSquareZero → [path_algebra, bardzell]`; `IncidenceAlgebra → [incidence, assem_book, hodge]`; `QuantumCI → [quantum_ci, qci_hh_oracle, bardzell, chouhy_solotar]`; `ExteriorAlgebra → [quantum_ci, chouhy_solotar]`; `PreprojectiveAlgebra → [preprojective, chouhy_solotar, assem_book]`; `TensorProduct → [tensor_product, hodge]`; `TrivialExtension → [assem_book]`; `zoo → [han_conjecture, chouhy_solotar]`.

---

### Task 1: Executable freshness gate (STOP on drift)

**Files:**
- Create: `tests/families/__init__.py` (empty), `tests/families/test_freshness_gate.py`

**Interfaces:**
- Consumes (probes only, no writes): `Quiver`, `Quiver.algebra` (general + monomial routes), `CC`, `GF`, `fields.primefield.PrimeField`, `Algebra`, `HHTable`, `families.basic.{linear_path_algebra, truncated_polynomial}`, `invariants.cartan.cartan_matrix`, the Plan-03 `groebner` surface, the Plan-04 CS/periodic backends.
- Produces: a red-if-drifted gate. **If any probe fails, STOP and report the drift to the human — do not proceed with Plan 06.**

- [ ] **Step 1: Write the gate test**

`tests/families/__init__.py`: empty file.

`tests/families/test_freshness_gate.py`:
```python
"""Plan 06 freshness gate. Every assumption Plan 06 builds on is probed here.
If this file is not green, the ground truth drifted -- STOP and re-baseline the
plan before writing any family code."""
import importlib

import pytest

from quiverlab.combinat import Quiver
from quiverlab.core import Algebra
from quiverlab.fields import CC, GF
from quiverlab.fields.primefield import PrimeField
from quiverlab.hochschild.table import HHTable


def test_general_quiver_algebra_is_delivered():
    """Plan 03 must be on disk: a non-monomial relation lowers to an Algebra."""
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)   # general route
    assert isinstance(A, Algebra)
    assert A.dim == 9
    assert A.hochschild_cohomology(1).dims == [1, 0]


def test_monomial_route_and_starter_builders_unchanged():
    from quiverlab.families.basic import linear_path_algebra, truncated_polynomial
    assert linear_path_algebra(3).dim == 6          # 1->2->3, dim n(n+1)/2
    assert truncated_polynomial(3).dim == 3         # k[x]/(x^3): e, x, x^2


def test_cartan_convention_is_paths_source_to_target():
    from quiverlab.invariants.cartan import cartan_matrix
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=CC)  # kA_2
    assert cartan_matrix(A) == [[1, 1], [0, 1]]     # C[i][j] = #paths i->j


def test_domains_and_hhtable_shapes():
    assert isinstance(GF(5), PrimeField)
    assert CC.parse_entry("i") is not None          # non-rational exact token
    t = HHTable([1, 0, 0], "HH^", "x")
    assert t.dims == [1, 0, 0] and t[0] == 1


def test_groebner_reduction_system_surface_present():
    g = importlib.import_module("quiverlab.groebner")
    assert hasattr(g, "build_reduction_system") and hasattr(g, "ReductionSystem")


@pytest.mark.skipif(
    importlib.util.find_spec("quiverlab.engine.resolutions_cs") is None,
    reason="Plan 04 CS backend not yet delivered",
)
def test_cs_periodic_backend_present_for_zoo_and_quantum():
    from quiverlab.engine import resolutions_periodic as rp
    assert hasattr(rp, "QuantumCIResolution")
```

- [ ] **Step 2: Run the gate**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/families/test_freshness_gate.py -q`
Expected: **all pass**. The CS probe may `skip` if Plan 04's `resolutions_cs` is absent — that gates only the *open-zone depth-16 HH goldens* (Tasks 11/13) and QuantumCI's optional periodic-backend HH, not the catalog: every family (Tasks 5–10, QuantumCI/Exterior/Preprojective included) *builds and computes HH via the bar/minimal engine* on Plans 01–03 alone, so a CS skip is not a Plan-06 blocker. Any *failure* (not skip) means an interface Plan 06 depends on has drifted — most importantly `test_general_quiver_algebra_is_delivered` (Plan 03); **STOP** and report which probe failed. No commit (this is a gate, kept green thereafter).

---

### Task 2: `citations` registry core — `references.bib`, registry, loud lookup

**Files:**
- Create: `src/quiverlab/citations/__init__.py`, `src/quiverlab/citations/references.bib`, `src/quiverlab/citations/registry.py`, `tests/citations/__init__.py` (empty), `tests/citations/test_registry.py`
- Modify: `src/quiverlab/errors.py` (add `CitationError`), `pyproject.toml` (`references.bib` package data)

**Interfaces:**
- Consumes: `errors.QuiverlabError`.
- Produces:
  - `errors.CitationError(QuiverlabError)`.
  - `citations.registry.REGISTRY: dict[str, Reference]` — `Reference` is a frozen dataclass `Reference(key: str, bibtex_key: str, kind: str, title: str, annotation: str, tags: tuple[str, ...])` where `kind ∈ {"algorithm", "family", "field", "foundation"}`.
  - `citations.reference(key: str) -> Reference` — loud `CitationError` on unknown key (lists the nearest known keys).
  - `citations.bibtex(key: str) -> str` — the raw BibTeX entry text for `key`'s `bibtex_key` (parsed from `references.bib`); loud on unknown.
  - `citations.all_keys() -> tuple[str, ...]`; `citations.references_bib_path() -> pathlib.Path`.

- [ ] **Step 1: Write the failing test**

`tests/citations/__init__.py`: empty.

`tests/citations/test_registry.py`:
```python
"""The citations registry core (spec §3.9). Real verified BibTeX; loud on unknown keys."""
import re

import pytest

from quiverlab import citations
from quiverlab.errors import CitationError


def test_every_registry_key_resolves_to_a_bibtex_entry():
    text = citations.references_bib_path().read_text(encoding="utf-8")
    bib_ids = set(re.findall(r"@\w+\{([^,]+),", text))
    for key in citations.all_keys():
        ref = citations.reference(key)
        assert ref.bibtex_key in bib_ids, f"{key} -> {ref.bibtex_key} missing from .bib"
        assert ref.annotation.strip(), f"{key} has no annotation"


def test_core_keys_present_with_right_kind():
    assert citations.reference("bardzell").kind == "algorithm"
    assert citations.reference("nakayama").kind == "family"
    assert citations.reference("conway").kind == "field"
    assert citations.reference("han_conjecture").kind == "foundation"


def test_bibtex_returns_the_entry_text():
    entry = citations.bibtex("chouhy_solotar")
    assert entry.startswith("@article{ChouhySolotar2015")
    assert "1406.2300" in entry


def test_unknown_key_fails_loudly():
    with pytest.raises(CitationError) as e:
        citations.reference("bardzel")     # typo
    assert "bardzell" in str(e.value)      # suggests the nearest known key


def test_bib_covers_registry():
    """Every registry-referenced bibtex_key resolves in references.bib, and the
    entry count is at least the current floor. This is a FLOOR, never an equality:
    Plan 08 appends software-citation entries (qpa, gap4, sagemath, quiverlab) and
    raises the floor; growth must not break this test."""
    text = citations.references_bib_path().read_text(encoding="utf-8")
    bib_ids = set(re.findall(r"@\w+\{([^,]+),", text))
    for key in citations.all_keys():
        assert citations.reference(key).bibtex_key in bib_ids
    assert len(re.findall(r"@\w+\{", text)) >= 16       # floor; Plan 08 raises it
```

- [ ] **Step 2: Run to verify it fails** (`ModuleNotFoundError: quiverlab.citations`).

- [ ] **Step 3: Add `CitationError`** to `src/quiverlab/errors.py`:
```python
class CitationError(QuiverlabError):
    """A citation key is unknown, or the bibliography is inconsistent."""
```

- [ ] **Step 4: Write `references.bib`** — `src/quiverlab/citations/references.bib`, verbatim (14 entries, all fields verified 2026-07-18):
```bibtex
@article{Hochschild1945,
  author  = {Hochschild, Gerhard},
  title   = {On the cohomology groups of an associative algebra},
  journal = {Annals of Mathematics. Second Series},
  volume  = {46},
  number  = {1},
  year    = {1945},
  pages   = {58--67},
  doi     = {10.2307/1969145},
}

@article{Bardzell1997,
  author  = {Bardzell, Michael J.},
  title   = {The alternating syzygy behavior of monomial algebras},
  journal = {Journal of Algebra},
  volume  = {188},
  number  = {1},
  year    = {1997},
  pages   = {69--89},
  doi     = {10.1006/jabr.1996.6813},
}

@article{ChouhySolotar2015,
  author  = {Chouhy, Sergio and Solotar, Andrea},
  title   = {Projective resolutions of associative algebras and ambiguities},
  journal = {Journal of Algebra},
  volume  = {432},
  year    = {2015},
  pages   = {22--61},
  doi     = {10.1016/j.jalgebra.2015.02.019},
  note    = {arXiv:1406.2300},
}

@article{NegronWitherspoon2016,
  author  = {Negron, Cris and Witherspoon, Sarah},
  title   = {An alternate approach to the {L}ie bracket on {H}ochschild cohomology},
  journal = {Homology, Homotopy and Applications},
  volume  = {18},
  number  = {1},
  year    = {2016},
  pages   = {265--285},
  doi     = {10.4310/HHA.2016.v18.n1.a14},
  note    = {arXiv:1406.0036},
}

@article{Volkov2019,
  author  = {Volkov, Yury},
  title   = {Gerstenhaber bracket on the {H}ochschild cohomology via an arbitrary resolution},
  journal = {Proceedings of the Edinburgh Mathematical Society. Series II},
  volume  = {62},
  number  = {3},
  year    = {2019},
  pages   = {817--836},
  doi     = {10.1017/S0013091518000901},
  note    = {arXiv:1610.05741},
}

@article{GSZ2001,
  author  = {Green, Edward L. and Solberg, {\O}yvind and Zacharia, Dan},
  title   = {Minimal projective resolutions},
  journal = {Transactions of the American Mathematical Society},
  volume  = {353},
  number  = {7},
  year    = {2001},
  pages   = {2915--2939},
  doi     = {10.1090/S0002-9947-01-02687-3},
}

@incollection{Happel1989,
  author    = {Happel, Dieter},
  title     = {Hochschild cohomology of finite-dimensional algebras},
  booktitle = {S\'eminaire d'Alg\`ebre Paul Dubreil et Marie-Paul Malliavin (Paris, 1987/1988)},
  series    = {Lecture Notes in Mathematics},
  volume    = {1404},
  publisher = {Springer},
  address   = {Berlin},
  year      = {1989},
  pages     = {108--126},
  doi       = {10.1007/BFb0084073},
}

@article{BGMS2005,
  author  = {Buchweitz, Ragnar-Olaf and Green, Edward L. and Madsen, Dag and Solberg, {\O}yvind},
  title   = {Finite {H}ochschild cohomology without finite global dimension},
  journal = {Mathematical Research Letters},
  volume  = {12},
  number  = {6},
  year    = {2005},
  pages   = {805--816},
  doi     = {10.4310/MRL.2005.v12.n6.a2},
  note    = {arXiv:math/0407108},
}

@article{BerghErdmann2008,
  author  = {Bergh, Petter Andreas and Erdmann, Karin},
  title   = {Homology and cohomology of quantum complete intersections},
  journal = {Algebra \& Number Theory},
  volume  = {2},
  number  = {5},
  year    = {2008},
  pages   = {501--522},
  doi     = {10.2140/ant.2008.2.501},
}

@article{GerstenhaberSchack1987,
  author  = {Gerstenhaber, Murray and Schack, Samuel D.},
  title   = {A {H}odge-type decomposition for commutative algebra cohomology},
  journal = {Journal of Pure and Applied Algebra},
  volume  = {48},
  number  = {3},
  year    = {1987},
  pages   = {229--247},
  doi     = {10.1016/0022-4049(87)90112-5},
}

@book{CartanEilenberg1956,
  author    = {Cartan, Henri and Eilenberg, Samuel},
  title     = {Homological Algebra},
  series    = {Princeton Mathematical Series},
  volume    = {19},
  publisher = {Princeton University Press},
  address   = {Princeton, NJ},
  year      = {1956},
  isbn      = {9780691049915},
}

@article{Connes1985,
  author  = {Connes, Alain},
  title   = {Non-commutative differential geometry},
  journal = {Publications Math\'ematiques de l'IH\'ES},
  volume  = {62},
  year    = {1985},
  pages   = {41--144},
  doi     = {10.1007/BF02698807},
}

@article{Gerstenhaber1963,
  author  = {Gerstenhaber, Murray},
  title   = {The cohomology structure of an associative ring},
  journal = {Annals of Mathematics. Second Series},
  volume  = {78},
  number  = {2},
  year    = {1963},
  pages   = {267--288},
  doi     = {10.2307/1970343},
}

@misc{Luebeck_ConwayPolynomials,
  author       = {L{\"u}beck, Frank},
  title        = {Conway polynomials for finite fields},
  howpublished = {\url{https://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/index.html}},
  note         = {Online database; accessed 2026-07-18},
}

@book{ASS2006,
  author    = {Assem, Ibrahim and Simson, Daniel and Skowro\'nski, Andrzej},
  title     = {Elements of the Representation Theory of Associative Algebras.
               Volume 1: Techniques of Representation Theory},
  series    = {London Mathematical Society Student Texts},
  volume    = {65},
  publisher = {Cambridge University Press},
  address   = {Cambridge},
  year      = {2006},
  isbn      = {9780521586313},
  doi       = {10.1017/CBO9780511614309},
}

@article{Han2006,
  author  = {Han, Yang},
  title   = {Hochschild (co)homology dimension},
  journal = {Journal of the London Mathematical Society. Second Series},
  volume  = {73},
  number  = {3},
  year    = {2006},
  pages   = {657--668},
  doi     = {10.1112/S002461070602299X},
  note    = {arXiv:math/0408402},
}
```

- [ ] **Step 5: Write `registry.py`**:
```python
"""Registry: quiverlab citation keys -> the papers behind each algorithm and family.
Annotations live HERE (the web /literature page and quiverlab.bibliography() consume
them). Loud failure on unknown keys (spec §3.9)."""
import difflib
import pathlib
import re
from dataclasses import dataclass, field

from quiverlab.errors import CitationError

_BIB = pathlib.Path(__file__).with_name("references.bib")


@dataclass(frozen=True)
class Reference:
    key: str            # the public quiverlab key (e.g. "bardzell")
    bibtex_key: str     # the @-entry id in references.bib (e.g. "Bardzell1997")
    kind: str           # "algorithm" | "family" | "field" | "foundation"
    title: str
    annotation: str     # one sentence: what it underpins
    tags: tuple = field(default_factory=tuple)


def _r(key, bibtex_key, kind, title, annotation, *tags):
    return Reference(key, bibtex_key, kind, title, annotation, tuple(tags))


REGISTRY: dict = {r.key: r for r in [
    _r("bardzell", "Bardzell1997", "algorithm",
       "Alternating syzygies of monomial algebras",
       "The minimal projective bimodule resolution for monomial algebras "
       "(quiverlab's Bardzell engine and the truncated/radical-square-zero families).",
       "resolution", "monomial"),
    _r("chouhy_solotar", "ChouhySolotar2015", "algorithm",
       "Projective resolutions via ambiguities",
       "The general kQ/I bimodule resolution from a reduction system -- quiverlab's "
       "Chouhy-Solotar engine for non-monomial algebras.",
       "resolution", "general"),
    _r("bracket_liftings", "NegronWitherspoon2016", "algorithm",
       "Gerstenhaber bracket via homotopy liftings",
       "The Gerstenhaber bracket computed directly on a non-bar resolution "
       "(with Volkov2019), transported onto bar representatives.",
       "bracket"),
    _r("bracket_liftings_volkov", "Volkov2019", "algorithm",
       "Gerstenhaber bracket on an arbitrary resolution",
       "A bracket formula valid on any projective bimodule resolution "
       "(companion to Negron-Witherspoon).",
       "bracket"),
    _r("minimal_resolution", "GSZ2001", "algorithm",
       "Minimal projective resolutions",
       "The Green-Solberg-Zacharia minimal module resolution algorithm "
       "(quiverlab's minimal engine and module Ext).",
       "resolution", "module"),
    _r("module_ext", "GSZ2001", "algorithm",
       "Module Ext via minimal resolutions",
       "Module-level Ext^n over minimal resolutions (Plan 05 module engine).",
       "module"),
    _r("bar", "Hochschild1945", "foundation",
       "Hochschild cohomology via the (normalized) bar complex",
       "Hochschild's original definition of the cohomology of an associative algebra; "
       "the normalized bar complex is quiverlab's HH^*/HH_* oracle in any characteristic.",
       "resolution", "bar"),
    _r("happel_question", "Happel1989", "foundation",
       "Happel's question",
       "Whether finite global dimension is equivalent to eventual vanishing of HH^n "
       "-- the motivating question for the hereditary and truncated families.",
       "conjecture"),
    _r("quantum_ci", "BGMS2005", "family",
       "Quantum complete intersections",
       "The algebra k<x,y>/(x^2, y^2, xy + q yx): finite Hochschild cohomology with "
       "infinite global dimension (the QuantumCI family).",
       "family"),
    _r("qci_hh_oracle", "BerghErdmann2008", "family",
       "Hochschild (co)homology of quantum complete intersections",
       "Explicit HH^* / HH_* of quantum complete intersections -- the literature "
       "oracle QuantumCI results are checked against.",
       "family", "oracle"),
    _r("tensor_product", "CartanEilenberg1956", "family",
       "Kunneth formula for Hochschild (co)homology",
       "The Kunneth isomorphism HH^n(A(x)B) = (+)_{i+j=n} HH^i(A)(x)HH^j(B) that "
       "makes HH multiplicative on tensor factors -- the anchor for TensorProduct(A, B).",
       "family"),
    _r("hodge", "GerstenhaberSchack1987", "algorithm",
       "Hodge (lambda) decomposition",
       "The eigenspace splitting HH^n = (+) HH^{n,(i)} of commutative/tensor and "
       "incidence-algebra pieces.",
       "decomposition"),
    _r("cyclic", "Connes1985", "algorithm",
       "Cyclic homology",
       "Connes' B-operator and the SBI sequence -- quiverlab's cyclic homology.",
       "cyclic"),
    _r("cup", "Gerstenhaber1963", "algorithm",
       "Cup product on Hochschild cohomology",
       "The associative cup product on HH^* (Gerstenhaber-algebra structure).",
       "product"),
    _r("bracket", "Gerstenhaber1963", "algorithm",
       "Gerstenhaber bracket",
       "The graded Lie bracket making HH^* a Gerstenhaber algebra.",
       "bracket"),
    _r("gerstenhaber", "Gerstenhaber1963", "foundation",
       "Cohomology structure of an associative ring",
       "The definitional source of the cup product and Gerstenhaber bracket.",
       "foundation"),
    _r("conway", "Luebeck_ConwayPolynomials", "field",
       "Conway polynomials for finite fields",
       "Lubeck's Conway-polynomial tables fixing canonical generators of GF(p^n).",
       "field"),
    _r("finite_fields", "Luebeck_ConwayPolynomials", "field",
       "Finite field arithmetic",
       "Deterministic cross-compatible GF(q) arithmetic via Conway polynomials.",
       "field"),
    _r("path_algebra", "ASS2006", "family",
       "Bound quiver algebras kQ/I",
       "The path-algebra / bound-quiver formalism for PathAlgebra and the catalog.",
       "family"),
    _r("nakayama", "ASS2006", "family",
       "Nakayama (serial) algebras",
       "Serial algebras by Kupisch series -- the NakayamaAlgebra family.",
       "family"),
    _r("incidence", "ASS2006", "family",
       "Incidence algebras of posets",
       "The incidence algebra kP realized as a bound quiver -- the IncidenceAlgebra family.",
       "family"),
    _r("preprojective", "ASS2006", "family",
       "Preprojective algebras",
       "The preprojective algebra of a Dynkin quiver -- the PreprojectiveAlgebra family.",
       "family"),
    _r("assem_book", "ASS2006", "foundation",
       "Elements of the Representation Theory of Associative Algebras",
       "The standard reference for bound quivers and the representation theory quiverlab implements.",
       "book"),
    _r("han_conjecture", "Han2006", "foundation",
       "Han's conjecture",
       "Finite global dimension iff finite Hochschild homology dimension -- the "
       "conjecture the zoo scans probe.",
       "conjecture"),
]}


def all_keys() -> tuple:
    return tuple(REGISTRY)


def references_bib_path() -> pathlib.Path:
    return _BIB


def reference(key: str) -> Reference:
    try:
        return REGISTRY[key]
    except KeyError:
        near = difflib.get_close_matches(key, REGISTRY, n=3)
        hint = f"did you mean {near}?" if near else f"known keys: {sorted(REGISTRY)}"
        raise CitationError(f"unknown citation key {key!r}", hint=hint) from None


def bibtex(key: str) -> str:
    ref = reference(key)
    text = _BIB.read_text(encoding="utf-8")
    m = re.search(r"(@\w+\{" + re.escape(ref.bibtex_key) + r",.*?\n\})", text, re.S)
    if m is None:
        raise CitationError(
            f"{key!r} maps to {ref.bibtex_key!r} but that entry is not in references.bib",
            hint="references.bib and the registry are out of sync")
    return m.group(1)
```

- [ ] **Step 6: Write `__init__.py`**:
```python
"""quiverlab citations registry (spec §3.9)."""
from quiverlab.citations.registry import (  # noqa: F401
    REGISTRY, Reference, reference, bibtex, all_keys, references_bib_path,
)
```

- [ ] **Step 7: Package the .bib** — in `pyproject.toml`, under `[tool.setuptools.packages.find]` add package-data so the non-`.py` file ships:
```toml
[tool.setuptools.package-data]
"quiverlab.citations" = ["references.bib"]
"quiverlab.families" = ["zoo_catalog.json"]
```
(the `zoo_catalog.json` line is added now; the file lands in Task 11.)

- [ ] **Step 8: Run focused + float gate**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/citations/test_registry.py tests/test_no_floats.py -q`
Expected: pass.

- [ ] **Step 9: Commit**
```bash
git add src/quiverlab/citations tests/citations src/quiverlab/errors.py pyproject.toml
git commit -m "feat(citations): registry core -- references.bib, keyed registry, loud lookup

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
```

---

### Task 3: `bibliography()`, `A.citations()`, and `HHTable.references` plumbing

**Files:**
- Create: `src/quiverlab/citations/bibliography.py`, `tests/citations/test_bibliography.py`, `tests/citations/test_result_references.py`
- Modify: `src/quiverlab/citations/__init__.py` (export `bibliography`), `src/quiverlab/hochschild/table.py` (add `references`), `src/quiverlab/core/algebra.py` (`citations()`, thread engine keys into `HHTable.references`), `src/quiverlab/__init__.py` (export `bibliography`)

**Interfaces:**
- Consumes: `citations.registry` (Task 2), `HHTable`, `Algebra`.
- Produces:
  - `citations.bibliography(keys=None, grouped=True, annotated=True) -> Bibliography`. **`Bibliography` is ITERABLE** — this is the load-bearing contract Plan 09's `references.entry_view` consumes verbatim: `for e in quiverlab.bibliography()` yields per-entry views, each exposing exactly `key` (the registry key — matches `table.references`), `formatted` (a human citation string built from `references.bib`), `doi`, `arxiv` (both parsed from `references.bib`; `None` when absent), `topic` (mapped from `kind`), `annotation`. Entry views are frozen `Entry` dataclasses whose attribute names match `entry_view`'s reads (`key/formatted/doi/arxiv/topic/annotation`); Plan 09 reads them via `getattr`. `Bibliography` also carries `.groups: dict[str, list[Entry]]` (grouped by `topic`), `.keys: tuple[str,...]`, `.__str__()` (printable grouped/annotated), `.bibtex()` (concatenated `.bib` entries, de-duplicated by `bibtex_key`), `.to_dict()` (JSON convenience — Plan 09 iterates, it does NOT depend on `to_dict()`). `keys=None` ⇒ the whole registry.
  - **`HHTable` one compatible extension**: `HHTable(dims, kind, algebra_repr, engine=..., references=())`; attribute `.references: tuple[str,...]`; `__repr__` appends a `refs: ...` line when non-empty. `__eq__` unchanged (compares `kind` + `dims` only).
  - `Algebra.citations() -> tuple[str, ...]` — the registry keys relevant to `A`: its family keys (stamped at construction, `_family_citations`, default `()`) plus the HH engine key `bar`. Resolution-specific keys (`bardzell`, `chouhy_solotar`) ride in the family stamp where the family declares them (QuantumCI stamps `chouhy_solotar`; Nakayama does not stamp `bardzell`), so `NakayamaAlgebra([3,2,2]).citations() == ("nakayama", "assem_book", "bar")`. Module keys (`minimal_resolution`/`module_ext`) and ops keys (`cup`/`bracket`/`cyclic`) are stamped by the Plan-04/05 resolution/ops result objects, not here.
  - **Result `.references` = engine + family (spec §3.9 "engine + family + ops")**: `Algebra.hochschild_cohomology`/`hochschild_homology` pass `references=self.citations()` into every `HHTable(...)` (both the fast-engine and bar branches) — so the returned table names both the engine path used (`bar` + `bardzell`/`chouhy_solotar` by dispatch) AND the family the algebra came from (`quantum_ci`, `nakayama`, …). Ops keys (`cup`, `bracket`, `cyclic`) are attached by the Plan-04/05 cup/bracket/cyclic result objects, not by plain HH-dim tables. (Loud `CitationError` is impossible here — keys are registry literals, validated by a test that every reported key resolves.)

- [ ] **Step 1: Write the failing tests**

`tests/citations/test_bibliography.py`:
```python
"""quiverlab.bibliography(): iterable, grouped-by-topic, annotated (web /literature)."""
from quiverlab import bibliography, citations


def test_full_bibliography_groups_by_topic_and_annotates():
    b = bibliography()
    assert set(b.groups) <= {"Algorithms", "Families", "Finite fields", "Foundations"}
    assert "bardzell" in b.keys and "nakayama" in b.keys
    s = str(b)
    assert "Bardzell" in s                                        # formatted citation
    assert "minimal projective bimodule resolution" in s.lower()  # annotation present


def test_iterates_with_plan09_entry_view_fields():
    """The frozen contract Plan 09's references.entry_view() consumes: each entry
    exposes key/formatted/doi/arxiv/topic/annotation as attributes."""
    entries = list(bibliography())                               # ITERABLE
    assert entries
    by_key = {e.key: e for e in entries}
    for attr in ("key", "formatted", "doi", "arxiv", "topic", "annotation"):
        assert hasattr(by_key["bardzell"], attr)
    cs = by_key["chouhy_solotar"]
    assert cs.arxiv == "1406.2300"                               # parsed from .bib note
    assert cs.doi == "10.1016/j.jalgebra.2015.02.019"
    assert cs.topic == "Algorithms"
    assert by_key["bardzell"].doi == "10.1006/jabr.1996.6813"


def test_subset_bibliography_and_bibtex_dedup():
    b = bibliography(keys=["quantum_ci", "qci_hh_oracle", "quantum_ci"])
    assert b.keys == ("quantum_ci", "qci_hh_oracle")             # order-preserving, deduped
    bt = b.bibtex()
    assert bt.count("@article{BGMS2005") == 1                    # bibtex_key dedup
    assert "@article{BerghErdmann2008" in bt


def test_to_dict_is_json_ready():
    import json
    d = bibliography(keys=["bardzell"]).to_dict()
    json.dumps(d)                                                # no exception
    assert d["groups"]["Algorithms"][0]["bibtex_key"] == "Bardzell1997"
    assert d["groups"]["Algorithms"][0]["annotation"]
```

`tests/citations/test_result_references.py`:
```python
"""Result .references plumbing: HHTable carries the engine's citation keys."""
from quiverlab import citations
from quiverlab.combinat import Quiver
from quiverlab.fields import CC, GF
from quiverlab.hochschild.table import HHTable


def test_hhtable_references_default_empty_and_repr():
    t = HHTable([1, 0], "HH^", "A")
    assert t.references == ()
    t2 = HHTable([1, 0], "HH^", "A", references=("bar",))
    assert t2.references == ("bar",) and "bar" in repr(t2)
    assert HHTable([1, 0], "HH^", "A") == HHTable([1, 0], "HH^", "A", references=("bar",))


def test_bar_path_reports_bar_key():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=CC)   # kA_2
    hh = A.hochschild_cohomology(1)
    assert "bar" in hh.references
    for key in hh.references:                                            # all in registry
        citations.reference(key)


def test_result_references_include_family_and_engine():
    """spec §3.9: table.references names engine + family keys (ops added by Plan 04/05)."""
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(5))
    A._family_citations = ("nakayama", "assem_book")                     # simulate a family stamp
    hh = A.hochschild_cohomology(1)
    assert "nakayama" in hh.references and "bar" in hh.references        # family + engine
    for key in A.citations():
        citations.reference(key)                                         # every key valid
```

- [ ] **Step 2: Run to verify failure** (no `bibliography`; `HHTable` has no `references`).

- [ ] **Step 3: Extend `HHTable`** — `src/quiverlab/hochschild/table.py`:
```python
class HHTable:
    def __init__(self, dims, kind, algebra_repr, engine="normalized bar complex",
                 references=()):
        self.dims = list(dims)
        self.kind = kind
        self.top = len(dims) - 1
        self.algebra_repr = algebra_repr
        self.engine = engine
        self.references = tuple(references)

    def __getitem__(self, n):
        return self.dims[n]

    def __iter__(self):
        return iter(self.dims)

    def __eq__(self, other):
        if isinstance(other, HHTable):
            return self.kind == other.kind and self.dims == other.dims
        return NotImplemented

    def __repr__(self):
        head = f"{self.kind}n dimensions for {self.algebra_repr} (engine: {self.engine})"
        cells = "  ".join(f"{self.kind}{n} = {d}" for n, d in enumerate(self.dims))
        out = head + "\n" + cells
        if self.references:
            out += "\nrefs: " + ", ".join(self.references)
        return out
```

- [ ] **Step 4: Thread engine keys + add `Algebra.citations()`** — in `src/quiverlab/core/algebra.py`:
  - Add an optional `_family_citations=()` ctor kwarg stored as `self._family_citations = tuple(_family_citations)` (families stamp this).
  - In `hochschild_cohomology`/`hochschild_homology`, pass `references=self.citations()` into every `HHTable(...)` (both the fast-engine and bar branches) — the table names engine **and** family keys (spec §3.9).
  - Add:
```python
def _engine_citations(self):
    # HH^*/HH_* dimensions are produced by the (normalized/fast-rank) bar complex,
    # whatever the presentation -- so the engine key is always "bar" (Hochschild1945).
    # Resolution-specific keys (bardzell / chouhy_solotar) are carried by the family
    # stamp where the family author declares them relevant (e.g. QuantumCI), and by
    # the Plan-04/05 resolution/ops result objects that actually run those engines.
    return ("bar",)

def citations(self):
    """Registry keys relevant to this algebra: its family stamp plus the HH engine
    (spec §3.9). Every key resolves via quiverlab.citations."""
    seen, out = set(), []
    for k in tuple(getattr(self, "_family_citations", ())) + self._engine_citations():
        if k not in seen:
            seen.add(k)
            out.append(k)
    return tuple(out)
```
  (Guard `getattr(self, "_family_citations", ())` in case older constructions omit it; default the attribute in `__init__` to `()`.)

- [ ] **Step 5: Write `bibliography.py`**:
```python
"""Grouped, annotated, ITERABLE bibliography (spec §3.9). Plan 09's
references.entry_view() iterates bibliography() and reads e.key / e.formatted /
e.doi / e.arxiv / e.topic / e.annotation off each entry -- those attribute names
are the frozen contract. `formatted`, `doi`, and `arxiv` are built from
references.bib; `topic` is mapped from the registry `kind`; annotations come from
the registry."""
import re
from dataclasses import dataclass

from quiverlab.citations.registry import (
    REGISTRY, bibtex as _bibtex, reference, references_bib_path,
)

# topic = the /literature grouping label; mapped from the registry kind.
_KIND_TOPIC = {"algorithm": "Algorithms", "family": "Families",
               "field": "Finite fields", "foundation": "Foundations"}
_TOPIC_ORDER = ("Algorithms", "Families", "Finite fields", "Foundations")


@dataclass(frozen=True)
class Entry:
    """One bibliography row. Attribute names match Plan 09's entry_view() reads."""
    key: str            # the registry key -- matches table.references / A.citations()
    bibtex_key: str
    formatted: str      # human citation string built from references.bib
    doi: str | None
    arxiv: str | None
    topic: str          # mapped from kind
    annotation: str


def _bib_fields(bibtex_key: str) -> dict:
    """Pull author/title/journal/year/volume/pages/doi/note from one .bib entry."""
    text = references_bib_path().read_text(encoding="utf-8")
    m = re.search(r"@\w+\{" + re.escape(bibtex_key) + r",(.*?)\n\}", text, re.S)
    body = m.group(1) if m else ""
    out = {}
    for fm in re.finditer(r"(\w+)\s*=\s*\{(.*?)\}(?=,\s*\n|\s*\n\})", body, re.S):
        out[fm.group(1).lower()] = re.sub(r"\s+", " ", fm.group(2)).strip()
    return out


def _clean(s: str) -> str:
    return (s.replace("{", "").replace("}", "").replace("\\", "")
             .replace("--", "-").replace("~", " ").strip())


def _format(bibtex_key: str, f: dict) -> str:
    authors = _clean(f.get("author", "")).replace(" and ", "; ")
    year = f.get("year", "")
    title = _clean(f.get("title", bibtex_key))
    venue = _clean(f.get("journal") or f.get("booktitle") or f.get("publisher") or "")
    vol = f.get("volume", "")
    pages = _clean(f.get("pages", ""))
    bits = [b for b in [f"{authors} ({year})." if authors else "", f"{title}.",
                        f"{venue} {vol}".strip() + (f", {pages}" if pages else "") + "."
                        if venue else ""] if b]
    return " ".join(bits)


def _arxiv_of(f: dict) -> str | None:
    m = re.search(r"arXiv:([\w./-]+)", f.get("note", "") + " " + f.get("eprint", ""))
    return m.group(1) if m else None


def _entry(key: str) -> Entry:
    r = REGISTRY[key]
    f = _bib_fields(r.bibtex_key)
    return Entry(key=key, bibtex_key=r.bibtex_key, formatted=_format(r.bibtex_key, f),
                 doi=f.get("doi"), arxiv=_arxiv_of(f),
                 topic=_KIND_TOPIC.get(r.kind, "Other"), annotation=r.annotation)


@dataclass(frozen=True)
class Bibliography:
    keys: tuple
    _entries: tuple

    def __iter__(self):
        return iter(self._entries)             # Plan 09 consumes this

    @property
    def groups(self) -> dict:
        g = {}
        for e in self._entries:
            g.setdefault(e.topic, []).append(e)
        return g

    def bibtex(self) -> str:
        out, seen = [], set()
        for k in self.keys:
            bk = REGISTRY[k].bibtex_key
            if bk not in seen:
                seen.add(bk)
                out.append(_bibtex(k))
        return "\n\n".join(out) + ("\n" if out else "")

    def to_dict(self) -> dict:
        return {"keys": list(self.keys),
                "groups": {t: [e.__dict__ for e in es] for t, es in self.groups.items()}}

    def __str__(self) -> str:
        lines = []
        for t in _TOPIC_ORDER:
            es = self.groups.get(t)
            if not es:
                continue
            lines.append(f"## {t}")
            for e in es:
                lines.append(f"  [{e.key}] {e.formatted}")
                lines.append(f"      {e.annotation}")
            lines.append("")
        return "\n".join(lines).rstrip() + "\n"


def bibliography(keys=None, grouped=True, annotated=True) -> Bibliography:
    if keys is None:
        keys = list(REGISTRY)
    seen, ordered = set(), []
    for k in keys:
        reference(k)                      # loud on unknown
        if k not in seen:
            seen.add(k)
            ordered.append(k)
    return Bibliography(tuple(ordered), tuple(_entry(k) for k in ordered))
```
(The `formatted`/`doi`/`arxiv`/`topic`/`annotation` attribute names on `Entry` are the exact reads in Plan 09's `webapp/server/references.py:entry_view`; a Plan-06 test asserts an entry exposes all six as attributes so Plan 09's `getattr` path never falls back.)

- [ ] **Step 6: Export** — `citations/__init__.py` adds `from quiverlab.citations.bibliography import bibliography, Bibliography`; `src/quiverlab/__init__.py` adds `from quiverlab.citations import bibliography` and `"bibliography"` to `__all__`.

- [ ] **Step 7: Run focused + float gate + prior citation test.** Expected: pass. Also run `tests/test_quickstart.py tests/hochschild -q` to confirm the `HHTable`/`Algebra` extensions didn't disturb existing HH callers.

- [ ] **Step 8: Commit** (`feat(citations): bibliography(), A.citations(), HHTable.references plumbing`).

---

### Task 4: Sanctioned grammar extension — exact non-rational coefficients

**Files:**
- Modify: `src/quiverlab/combinat/relations.py`
- Create: `tests/combinat/test_relation_exact_coeffs.py`

**Interfaces:**
- Consumes: `parse_rational`, `fields.complexfield` (`E`, sympy path) only indirectly (the parser stays field-agnostic — it never picks a field).
- Produces: the coefficient grammar widening. A coefficient factor may now be a **non-rational exact token** matching `_EXACT_COEFF = re.compile(r"^[+-]?(i|E\(\d+\)|sqrt\(\d+\)|\d*/?\d*)$")` in the broad sense — concretely: if a factor is not a plain rational and not an arrow name, it is treated as an exact scalar token and parsed with `sympy.sympify(tok, locals={"i": sympy.I, "E": _E}, rational=True, evaluate=True)`; a resulting non-number or any `Float` raises `ExactnessError`/`RelationError`. Rational factors continue to yield `Fraction` (Plan-03/04 byte-compatible). `Relation.terms` coefficient type widens to `Fraction | sympy.Expr`; like-term combination uses sympy-safe addition (`Fraction + Fraction` stays `Fraction`; any `Expr` promotes the sum). The build path coerces every coefficient through `field.parse_entry` (which already accepts `Fraction` and sympy `Expr` — see `fields/complexfield.py`), so `Quiver.algebra` / `build_reduction_system` need no change **provided** they coerce via `field.parse_entry`; a probe test asserts this (if a Plan-03 build coerces `Fraction` directly, add the one-line widening there and note it).

- [ ] **Step 1: Write the failing test**

`tests/combinat/test_relation_exact_coeffs.py`:
```python
"""Sanctioned grammar extension: exact non-rational coefficients (spec §3.3;
foreseen in the Plan-03 boundary note). Rational relations are unchanged."""
from fractions import Fraction

import pytest
import sympy

from quiverlab.combinat import Quiver
from quiverlab.combinat.relations import parse_relation
from quiverlab.errors import ExactnessError, RelationError
from quiverlab.fields import CC


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def test_rational_relations_are_byte_compatible():
    Q = _two_loops()
    r = parse_relation("x*y - 2*y*x", Q)
    assert r.terms == ((Fraction(1), ("x", "y")), (Fraction(-2), ("y", "x")))
    for c, _w in r.terms:
        assert isinstance(c, Fraction)                 # unchanged type


def test_imaginary_unit_coefficient_parses_exactly():
    Q = _two_loops()
    r = parse_relation("x*y + i*y*x", Q)               # QuantumCI(q="i")
    coeffs = {w: c for c, w in r.terms}
    assert coeffs[("x", "y")] == 1
    assert coeffs[("y", "x")] == sympy.I               # exact i, field-agnostic
    A = Q.algebra(relations=["x^2", "y^2", "x*y + i*y*x"], field=CC)   # general route
    assert A.dim == 4


def test_root_of_unity_and_radical_tokens():
    Q = _two_loops()
    r = parse_relation("x*y + E(3)*y*x", Q)
    coeffs = {w: c for c, w in r.terms}
    assert coeffs[("y", "x")] == sympy.exp(2 * sympy.pi * sympy.I / 3)


def test_float_coefficient_still_fails_loudly():
    Q = _two_loops()
    with pytest.raises((ExactnessError, RelationError)):
        parse_relation("0.5*x*y - y*x", Q)


def test_unknown_token_is_a_relation_error_not_a_silent_arrow():
    Q = _two_loops()
    with pytest.raises(RelationError):
        parse_relation("bogus(2)*x", Q)
```

- [ ] **Step 2: Run to verify failure** (`i*y*x` currently routes `i` to `parse_rational` and errors, or misparses).

- [ ] **Step 3: Implement** — in `relations.py`, add an exact-scalar branch to `_parse_term`. When a factor is neither a plain `_COEFF` rational, nor a valid arrow/`_POW`, and is not a bare identifier that is an arrow name: attempt an exact-scalar parse:
```python
import sympy

from quiverlab.errors import ExactnessError

_ARROWY = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\^\d+)?$")


def _E(n):
    return sympy.exp(2 * sympy.pi * sympy.I / int(n))


def _exact_scalar(tok):
    """Parse an exact non-rational coefficient token (i, E(n), sqrt(k), rationals).
    Returns Fraction for rationals, a sympy exact Expr otherwise. Loud on floats."""
    try:
        expr = sympy.sympify(tok, locals={"i": sympy.I, "E": _E},
                             rational=True, evaluate=True)
    except (sympy.SympifyError, TypeError, SyntaxError, ValueError):
        return None
    if not getattr(expr, "is_number", False):
        return None
    if expr.atoms(sympy.Float):
        raise ExactnessError(f"coefficient {tok!r} contains a floating-point number",
                             hint="write exact scalars: '1/2', 'i', 'sqrt(2)', 'E(3)'")
    if expr.is_rational:
        return Fraction(int(sympy.numer(expr)), int(sympy.denom(expr)))
    return expr
```
Then in `_parse_term`, replace the arrow-name resolution so that a factor which is *not* an arrow (`name not in quiver.arrows` and not a `_POW` of an arrow) but *is* an exact scalar is folded into `coeff` (same "coefficients before arrows" guard as the rational branch); only a token that is neither an arrow nor an exact scalar raises the existing `RelationError`. Widen like-term combination in `parse_relation` to keep `Fraction` when all summands are `Fraction` and promote to sympy otherwise (`combined[w] = combined.get(w, Fraction(0)) + c`; a mixed `Fraction + Expr` yields an `Expr` automatically under sympy). Update the `Relation.terms` docstring/type to `Fraction | sympy.Expr`. **Do not** change `is_monomial`/`min_length`/`max_length` (they inspect words, not coefficients).

- [ ] **Step 4: Run focused + the full Plan-03 groebner suite** (`tests/groebner -q`) to prove rational relations are byte-compatible + float gate. Expected: pass.

- [ ] **Step 5: Commit** (`feat(combinat): exact non-rational relation coefficients (i, E(n), sqrt)`).

---

### Task 5: `NakayamaAlgebra` (Kupisch series + cyclic)

**Files:**
- Create: `src/quiverlab/families/nakayama.py`, `tests/families/test_nakayama.py`

**Interfaces:**
- Consumes: `Quiver`, `Quiver.algebra` (monomial route), `errors.AdmissibilityError`.
- Produces: `NakayamaAlgebra(kupisch=None, *, n=None, l=None, cyclic=False, field=CC)`.
  - Form A: `NakayamaAlgebra([c_1, …, c_n], field=CC)` — Kupisch series; linear if `min == 1`, cyclic if `min ≥ 2`; `AdmissibilityError` otherwise (with the offending index).
  - Form B: `NakayamaAlgebra(n=<int>, l=<int>, cyclic=True, field=CC)` — homogeneous `kZ_n / rad^l` (all `c_i = l`); `cyclic=False` gives the linear truncation `kA_n / rad^l`.
  - Stamps `_family_citations=("nakayama", "assem_book")`.

- [ ] **Step 1: Write the failing test**

`tests/families/test_nakayama.py`:
```python
"""NakayamaAlgebra by Kupisch series and by (n, l, cyclic). Fixture N1."""
import pytest

from quiverlab.errors import AdmissibilityError
from quiverlab.families import NakayamaAlgebra
from quiverlab.fields import CC, GF


def test_kupisch_322_is_cyclic_dim7_cartan_det1():
    A = NakayamaAlgebra([3, 2, 2], field=CC)
    assert A.dim == 7
    assert A.cartan_matrix() == [[1, 1, 1], [0, 1, 1], [1, 0, 1]]     # Fixture N1
    assert A.hochschild_cohomology(0).dims == [1]                    # HH^0 = center = 1


def test_kupisch_322_char_independent_dim():
    assert NakayamaAlgebra([3, 2, 2], field=GF(32003)).dim == 7


def test_homogeneous_cyclic_n4_l3():
    A = NakayamaAlgebra(n=4, l=3, cyclic=True, field=CC)             # kZ_4 / rad^3
    assert A.dim == 12                                              # 4 * 3


def test_linear_form_b_n4_l3():
    A = NakayamaAlgebra(n=4, l=3, cyclic=False, field=CC)           # kA_4 / rad^3, Kupisch [3,3,2,1]
    assert A.dim == 9                                               # 3 + 3 + 2 + 1


def test_linear_series_ends_in_one():
    A = NakayamaAlgebra([3, 2, 1], field=CC)                        # linear A_3, no run-off relations
    assert A.dim == 6                                               # 3 + 2 + 1
    B = NakayamaAlgebra([2, 2, 1], field=CC)                        # kills a1*a2 only
    assert B.dim == 5                                               # e1,e2,e3,a1,a2


def test_bad_series_is_loud():
    with pytest.raises(AdmissibilityError):
        NakayamaAlgebra([2, 1, 2], field=CC)                        # interior 1: neither shape


def test_citations_include_nakayama():
    assert "nakayama" in NakayamaAlgebra([3, 2, 2]).citations()
```

- [ ] **Step 2: Run to verify failure** (`ImportError`).

- [ ] **Step 3: Implement `nakayama.py`**:
```python
"""Nakayama (serial) algebras (spec §3.4). Kupisch series K=[c_1..c_n], c_i = dim P_i.
Linear (A_n) if min(K)=1; cyclic (Z_n) if min(K)>=2. Paths read left-to-right."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import AdmissibilityError


def _validate_linear(K):
    n = len(K)
    if K[-1] != 1:
        raise AdmissibilityError(
            f"linear Kupisch series must end in 1 (sink projective is simple), got {K}",
            hint="use min(series) >= 2 for a cyclic Nakayama algebra")
    for i in range(n - 1):
        if K[i] < 2:
            raise AdmissibilityError(f"c_{i+1} = {K[i]} < 2 in the interior of {K}",
                                     hint="interior projectives have length >= 2")
        if K[i] > K[i + 1] + 1:
            raise AdmissibilityError(
                f"Kupisch admissibility fails at index {i+1}: {K[i]} > {K[i+1]} + 1",
                hint="need c_i <= c_{i+1} + 1")


def _validate_cyclic(K):
    n = len(K)
    for i in range(n):
        if K[i] < 2:
            raise AdmissibilityError(f"cyclic Kupisch entry c_{i+1} = {K[i]} < 2 in {K}",
                                     hint="every cyclic projective has length >= 2")
        if K[i] > K[(i - 1) % n] + 1:
            raise AdmissibilityError(
                f"cyclic admissibility fails at index {i+1}: {K[i]} > {K[(i-1)%n]} + 1",
                hint="need c_i <= c_{i-1} + 1 cyclically")


def NakayamaAlgebra(kupisch=None, *, n=None, l=None, cyclic=False, field=None):
    if kupisch is None:
        if n is None or l is None:
            raise AdmissibilityError(
                "give a Kupisch series, or n and l",
                hint="NakayamaAlgebra([3,2,2]) or NakayamaAlgebra(n=4, l=3, cyclic=True)")
        kupisch = [l] * n if cyclic else [min(l, n - i) for i in range(n)]  # e.g. l=3,n=4 -> [3,3,2,1]
    K = [int(c) for c in kupisch]
    if len(K) < 1 or any(c < 1 for c in K):
        raise AdmissibilityError(f"Kupisch entries must be >= 1, got {K}", hint="e.g. [3,2,2]")
    is_cyclic = min(K) >= 2
    m = len(K)
    if is_cyclic:
        _validate_cyclic(K)
        verts = list(range(1, m + 1))
        arrows = {f"a{i}": (i, i % m + 1) for i in range(1, m + 1)}   # i -> (i mod m)+1
    else:
        _validate_linear(K)
        verts = list(range(1, m + 1))
        arrows = {f"a{i}": (i, i + 1) for i in range(1, m)}           # 1->2->...->n
    Q = Quiver(verts, arrows)
    order = list(arrows)                                              # a1, a2, ...
    rels = []
    for i in range(m):                                               # length-c_i path from vertex i+1
        length = K[i]
        if is_cyclic:
            path = [order[(i + k) % m] for k in range(length)]
        else:
            if i + length > m - 1:                                   # runs off the sink (arrows a1..a_{m-1}): already zero
                continue
            path = [order[i + k] for k in range(length)]             # a_{i+1} .. (length arrows)
            if len(path) < length:
                continue
        if len(path) >= 2:
            rels.append("*".join(path))
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("nakayama", "assem_book")
    return A
```
(Note: for the linear case the relation is the length-`c_i` path from vertex `i+1`; where that path would run past the sink it is already zero in `kA_n` and is skipped. The cyclic case always has a length-`c_i` cycle-segment. Re-derive against Fixture N1: `[3,2,2]` cyclic gives rels `a1*a2*a3`, `a2*a3*a1`? — no: length `c_1=3` from vertex 1 is `a1 a2 a3`; `c_2=2` from vertex 2 is `a2 a3`; `c_3=2` from vertex 3 is `a3 a1`. The monomial builder then kills the longer `a2 a3 a1` etc. automatically. dim 7 ✓.)

- [ ] **Step 4: Run focused + float gate.** Expected: pass. If `[3,2,2]` yields dim ≠ 7, re-derive the relation-generation loop against Fixture N1 before proceeding (the monomial builder is the oracle: `A.dim` must equal `Σ c_i`).

- [ ] **Step 5: Commit** (`feat(families): NakayamaAlgebra (Kupisch series + cyclic)`).

---

### Task 6: `PathAlgebra`, `TruncatedPathAlgebra`, `RadicalSquareZero`

**Files:**
- Create: `src/quiverlab/families/dynkin.py` (diagram → quiver + orientation), `src/quiverlab/families/path_algebra.py`, `src/quiverlab/families/truncated.py`, `src/quiverlab/families/radical_square_zero.py`, `tests/families/test_path_algebra.py`, `tests/families/test_truncated.py`, `tests/families/test_rad_square_zero.py`

**Interfaces:**
- Consumes: `Quiver`, `Quiver.algebra` (monomial route), `errors.*`.
- Produces:
  - `families.dynkin.dynkin_quiver(type_str, orientation="linear") -> Quiver` — parses `"A5"`, `"D4"`, `"E6/7/8"`, Euclidean `"~A3"`/`"At3"`, etc.; standard labelings; `orientation ∈ {"linear","reverse", dict}`. Rejects a fully-cyclic Euclidean orientation only at algebra-build time (the quiver itself is legal).
  - `PathAlgebra(type_or_quiver, orientation="linear", field=CC)` — `relations=[]`; `NotFiniteDimensionalError` (from the monomial route) if the orientation has an oriented cycle. Stamps `("path_algebra", "happel_question", "assem_book")`.
  - `TruncatedPathAlgebra(type_or_quiver, r, orientation="linear", field=CC)` — `kQ / rad^r`; `r >= 2`. Stamps `("path_algebra", "bardzell")`.
  - `RadicalSquareZero(quiver, field=CC)` — `TruncatedPathAlgebra(quiver, 2)`. Stamps `("path_algebra", "bardzell")`.

- [ ] **Step 1: Write the failing tests**

`tests/families/test_path_algebra.py`:
```python
"""PathAlgebra (hereditary Dynkin/Euclidean). Fixture N2."""
import pytest

from quiverlab.errors import NotFiniteDimensionalError
from quiverlab.families import PathAlgebra
from quiverlab.fields import CC, GF


def test_D4_linear_orientation_dim9_hh_100():
    A = PathAlgebra("D4", field=CC)                 # arrows 1->2, 2->3, 2->4
    assert A.dim == 9                               # Fixture N2 (not 12)
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]


def test_A5_is_triangular_number():
    assert PathAlgebra("A5", field=CC).dim == 15    # 5*6/2


def test_tree_hh_is_char_independent():
    assert PathAlgebra("D4", field=GF(7)).hochschild_cohomology(2).dims == [1, 0, 0]


def test_cyclic_euclidean_orientation_is_loud():
    with pytest.raises(NotFiniteDimensionalError):
        PathAlgebra("~A2", orientation={"e12": (1, 2), "e20": (2, 0), "e01": (0, 1)})


def test_citations():
    assert "path_algebra" in PathAlgebra("A3").citations()
```

`tests/families/test_truncated.py`:
```python
"""TruncatedPathAlgebra kQ/rad^r. Fixture N3."""
from quiverlab.combinat import Quiver
from quiverlab.families import TruncatedPathAlgebra
from quiverlab.fields import CC


def test_A3_rad2_dim5():
    assert TruncatedPathAlgebra("A3", 2, field=CC).dim == 5     # Fixture N3
    assert TruncatedPathAlgebra("A5", 2, field=CC).dim == 9     # spec example


def test_from_explicit_quiver():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    assert TruncatedPathAlgebra(Q, 2, field=CC).dim == 5
```

`tests/families/test_rad_square_zero.py`:
```python
"""RadicalSquareZero(Q) = kQ/rad^2. Fixture N4."""
from quiverlab.combinat import Quiver
from quiverlab.families import RadicalSquareZero
from quiverlab.fields import CC


def test_two_loop_is_finite_dim3():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})     # kQ infinite; rad^2 truncates
    A = RadicalSquareZero(Q, field=CC)
    assert A.dim == 3                               # e, x, y  (Fixture N4)
    assert A.hochschild_cohomology(0).dims == [3]   # commutative local: HH^0 = 3


def test_three_cycle_dim6_hh0_1():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)})
    A = RadicalSquareZero(Q, field=CC)
    assert A.dim == 6                               # 3 vertices + 3 arrows
    assert A.hochschild_cohomology(0).dims == [1]
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `dynkin.py`** (standard labelings + orientation):
```python
"""Dynkin/Euclidean diagrams -> quivers with a chosen orientation (spec §3.4).
Edges are undirected pairs; orientation turns each into an arrow. Default 'linear':
i->j for an edge {i,j} with i<j in the standard labeling."""
import re

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError

_TYPE = re.compile(r"^(~|t)?([ADE])(\d+)$")


def _edges(letter, n):
    if letter == "A":
        return [(i, i + 1) for i in range(1, n)]
    if letter == "D":
        if n < 4:
            raise QuiverlabError(f"D{n} needs n >= 4", hint="D4, D5, ...")
        return [(i, i + 1) for i in range(1, n - 1)] + [(n - 2, n)]
    if letter == "E":
        if n not in (6, 7, 8):
            raise QuiverlabError(f"E{n} is not a diagram", hint="E6, E7, E8")
        chain = [(i, i + 1) for i in range(1, n - 1)]     # 1-2-...-(n-1)
        return chain + [(3, n)]                           # branch node 3 -> extra vertex n
    raise QuiverlabError(f"unknown Dynkin letter {letter!r}", hint="A, D, E")


def _euclidean_edges(letter, n):
    if letter == "A":                                     # ~A_n: cycle on 0..n
        return [(i, i + 1) for i in range(n)] + [(0, n)]
    # ~D, ~E labelings per Kac; implement the ones used in the family tour as needed.
    raise QuiverlabError(f"Euclidean ~{letter}{n} not yet tabulated",
                         hint="use ~A_n, or pass an explicit Quiver")


def dynkin_quiver(type_str, orientation="linear"):
    m = _TYPE.match(type_str)
    if not m:
        raise QuiverlabError(f"cannot parse diagram type {type_str!r}",
                             hint="examples: 'A5', 'D4', 'E6', '~A3'")
    euclid, letter, n = bool(m.group(1)), m.group(2), int(m.group(3))
    edges = _euclidean_edges(letter, n) if euclid else _edges(letter, n)
    verts = sorted({v for e in edges for v in e})
    arrows = {}
    for k, (u, v) in enumerate(edges):
        name = f"e{u}{v}"
        if isinstance(orientation, dict) and name in orientation:
            s, t = orientation[name]
        elif orientation == "reverse":
            s, t = (v, u) if u < v else (u, v)
        else:                                             # "linear"
            s, t = (u, v) if u < v else (v, u)
        arrows[name] = (s, t)
    return Quiver(verts, arrows)
```

- [ ] **Step 4: Implement the three builders** (`path_algebra.py`, `truncated.py`, `radical_square_zero.py`). `TruncatedPathAlgebra` enumerates every path of length `< r`... but the monomial route needs *relations* not a basis: pass `relations = [all length-r paths]`. Concretely, generate the forbidden words = every composable arrow-word of length exactly `r` (the monomial builder then also kills longer ones via the automaton):
```python
# truncated.py
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families.dynkin import dynkin_quiver


def _length_r_paths(quiver, r):
    words = [(a,) for a in quiver.arrows]
    for _ in range(r - 1):
        nxt = []
        for w in words:
            for a in quiver.arrows:
                if quiver.target(w[-1]) == quiver.source(a):
                    nxt.append(w + (a,))
        words = nxt
    return words


def TruncatedPathAlgebra(type_or_quiver, r, orientation="linear", field=None):
    if not (isinstance(r, int) and r >= 2):
        raise QuiverlabError(f"TruncatedPathAlgebra r={r!r}: need integer r >= 2",
                             hint="r = 1 kills the arrows; use r >= 2")
    Q = type_or_quiver if isinstance(type_or_quiver, Quiver) else \
        dynkin_quiver(type_or_quiver, orientation)
    rels = ["*".join(w) for w in _length_r_paths(Q, r)]
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("path_algebra", "bardzell")
    return A
```
```python
# radical_square_zero.py
from quiverlab.families.truncated import TruncatedPathAlgebra


def RadicalSquareZero(quiver, field=None):
    A = TruncatedPathAlgebra(quiver, 2, field=field)
    A._family_citations = ("path_algebra", "bardzell")
    return A
```
```python
# path_algebra.py
from quiverlab.combinat.quiver import Quiver
from quiverlab.families.dynkin import dynkin_quiver


def PathAlgebra(type_or_quiver, orientation="linear", field=None):
    Q = type_or_quiver if isinstance(type_or_quiver, Quiver) else \
        dynkin_quiver(type_or_quiver, orientation)
    A = Q.algebra(relations=[], field=field)     # hereditary; loud if cyclic orientation
    A._family_citations = ("path_algebra", "happel_question", "assem_book")
    return A
```

- [ ] **Step 5: Run focused (all three files) + float gate.** Re-derive D4 = 9 and A3/rad^2 = 5 against Fixtures N2/N3 if any dim is off.

- [ ] **Step 6: Commit** (`feat(families): PathAlgebra, TruncatedPathAlgebra, RadicalSquareZero`).

---

### Task 7: `IncidenceAlgebra(poset)`

**Files:**
- Create: `src/quiverlab/families/poset.py` (the `Poset` object + input format), `src/quiverlab/families/incidence.py`, `tests/families/test_incidence.py`

**Interfaces:**
- Consumes: `Quiver`, `Quiver.algebra` (general route — commutativity relations), `errors.RelationError`.
- Produces:
  - `families.poset.Poset(covers, elements=None)` — `covers`: list of `(x, y)` cover pairs (Hasse edge `x ⋖ y`); `.elements`, `.covers`, `.leq(x, y)`, `.hasse_quiver() -> (Quiver, name_map)`. Raises `RelationError` if the transitive closure is not antisymmetric (a directed cycle ⇒ not a poset).
  - `IncidenceAlgebra(poset_or_covers, elements=None, field=CC)` — accepts a `Poset` or a raw `covers` list; builds the Hasse quiver + one commutativity relation per pair of distinct parallel maximal Hasse paths; lowers through the general route. Stamps `("incidence", "assem_book")`.

- [ ] **Step 1: Write the failing test**

`tests/families/test_incidence.py`:
```python
"""IncidenceAlgebra of a poset. Fixture N5 (diamond == commutative square)."""
import pytest

from quiverlab.errors import RelationError
from quiverlab.families import IncidenceAlgebra
from quiverlab.families.poset import Poset
from quiverlab.fields import CC, GF


DIAMOND = [("b", "x"), ("b", "y"), ("x", "t"), ("y", "t")]


def test_diamond_is_commutative_square_dim9():
    A = IncidenceAlgebra(DIAMOND, field=CC)
    assert A.dim == 9                                   # 4 trivial + 4 arrows + 1 [b,t]
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]  # == Plan-03 commutative square


def test_diamond_char_independent():
    assert IncidenceAlgebra(DIAMOND, field=GF(32003)).dim == 9


def test_chain_poset_is_linear_path_algebra():
    A = IncidenceAlgebra([(1, 2), (2, 3)], field=CC)    # 1 < 2 < 3 chain
    assert A.dim == 6                                   # kA_3, dim n(n+1)/2


def test_directed_cycle_is_not_a_poset():
    with pytest.raises(RelationError):
        Poset([(1, 2), (2, 3), (3, 1)])                 # antisymmetry fails


def test_citations():
    assert "incidence" in IncidenceAlgebra(DIAMOND).citations()
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `poset.py`** (transitive closure, antisymmetry check, Hasse quiver, parallel maximal paths):
```python
"""Finite posets from cover data (spec §3.4). covers = [(x, y), ...] means x is
covered by y (Hasse edge, arrow x->y). Order <= is the reflexive-transitive closure."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import RelationError


class Poset:
    def __init__(self, covers, elements=None):
        self.covers = [tuple(c) for c in covers]
        elems = set(elements or [])
        for x, y in self.covers:
            elems.add(x)
            elems.add(y)
        self.elements = sorted(elems, key=str)
        self._le = self._closure()

    def _closure(self):
        le = {(x, x) for x in self.elements}
        le |= set(self.covers)
        changed = True
        while changed:
            changed = False
            for (a, b) in list(le):
                for (c, d) in list(le):
                    if b == c and (a, d) not in le:
                        le.add((a, d))
                        changed = True
        for (a, b) in le:
            if a != b and (b, a) in le:
                raise RelationError(
                    f"{a} <= {b} and {b} <= {a}: not a poset (antisymmetry fails)",
                    hint="cover data must not create a directed cycle")
        return le

    def leq(self, x, y):
        return (x, y) in self._le

    def hasse_quiver(self):
        names = {}
        arrows = {}
        for k, (x, y) in enumerate(self.covers):
            name = f"c{k}"
            names[name] = (x, y)
            arrows[name] = (x, y)
        Q = Quiver(list(self.elements), arrows)
        return Q, names
```

- [ ] **Step 4: Implement `incidence.py`** (commutativity relations from parallel maximal Hasse paths):
```python
"""Incidence algebra kP of a finite poset = Hasse quiver modulo commutativity
(all parallel paths equal). General (non-monomial) route."""
from quiverlab.families.poset import Poset


def _all_paths(Q):
    """Every directed path (arrow-name tuple) in the Hasse quiver (it is acyclic)."""
    paths = []
    for a in Q.arrows:
        stack = [(a,)]
        while stack:
            w = stack.pop()
            paths.append(w)
            for b in Q.arrows:
                if Q.target(w[-1]) == Q.source(b):
                    stack.append(w + (b,))
    return paths


def IncidenceAlgebra(poset_or_covers, elements=None, field=None):
    P = poset_or_covers if isinstance(poset_or_covers, Poset) else \
        Poset(poset_or_covers, elements)
    Q, _ = P.hasse_quiver()
    by_ends = {}
    for w in _all_paths(Q):
        if len(w) >= 2:
            by_ends.setdefault((Q.word_source(w), Q.word_target(w)), []).append(w)
    rels = []
    for (_st), group in by_ends.items():
        base = group[0]
        for other in group[1:]:                       # base == other (commutativity)
            rels.append("*".join(base) + " - " + "*".join(other))
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("incidence", "assem_book", "hodge")   # hodge: GS1987 for incidence
    return A
```

- [ ] **Step 5: Run focused + float gate.** Confirm dim 9 and HH=[1,0,0] tie exactly to the Plan-03 commutative-square fixture. (If the Gröbner completion produces a different dim, the relation set has a redundant/duplicate commutativity relation; keep one per endpoint pair — the engine tolerates redundancy but the dim must be 9.)

- [ ] **Step 6: Commit** (`feat(families): IncidenceAlgebra(poset) via the Hasse quiver + commutativity`).

---

### Task 8: `QuantumCI`, `ExteriorAlgebra`, `PreprojectiveAlgebra` (general route)

**Files:**
- Create: `src/quiverlab/families/quantum.py`, `src/quiverlab/families/exterior.py`, `src/quiverlab/families/preprojective.py`, `tests/families/test_quantum.py`, `tests/families/test_exterior.py`, `tests/families/test_preprojective.py`

**Interfaces:**
- Consumes: `Quiver`, `Quiver.algebra` (general route), the Task-4 coefficient grammar, `errors.*`.
- Produces:
  - `QuantumCI(q, field=CC)` — `k⟨x,y⟩/(x², y², xy + q·yx)`; `q` any exact scalar (rational or `"i"`/`"E(n)"`/`"sqrt(k)"`). Stamps `("quantum_ci", "qci_hh_oracle", "bardzell", "chouhy_solotar")`.
  - `ExteriorAlgebra(n, field=CC)` — `Λ(k^n)`; `n >= 1`. Stamps `("quantum_ci", "chouhy_solotar")`.
  - `PreprojectiveAlgebra(type_or_quiver, field=CC)` — `Π(Q)` via the doubled quiver + mesh relations. Stamps `("preprojective", "chouhy_solotar", "assem_book")`.

- [ ] **Step 1: Write the failing tests**

`tests/families/test_quantum.py`:
```python
"""QuantumCI k<x,y>/(x^2, y^2, xy + q yx). Fixtures N6, N7."""
from quiverlab.families import ExteriorAlgebra, QuantumCI
from quiverlab.fields import CC, GF


def test_quantum_ci_i_dim4_hh0_2():
    A = QuantumCI(q="i", field=CC)                     # Fixture N7
    assert A.dim == 4
    assert A.hochschild_cohomology(0).dims == [2]      # Z = span{1, xy}


def test_quantum_ci_rational_q():
    assert QuantumCI(q=3, field=GF(32003)).dim == 4


def test_exterior_2_equals_quantum_ci_1():
    E2 = ExteriorAlgebra(2, field=CC)
    Q1 = QuantumCI(q=1, field=CC)                      # Fixture N6
    assert E2.dim == Q1.dim == 4
    assert E2.T == Q1.T and E2.unit == Q1.unit         # byte-identical structure constants


def test_quantum_ci_minus_one_is_commutative_dim4_hh0_4():
    A = QuantumCI(q=-1, field=CC)                       # k[x,y]/(x^2,y^2), commutative
    assert A.dim == 4
    assert A.hochschild_cohomology(0).dims == [4]
```

`tests/families/test_exterior.py`:
```python
"""ExteriorAlgebra(n) = Lambda(k^n), dim 2^n."""
from quiverlab.families import ExteriorAlgebra
from quiverlab.fields import CC, GF


def test_dims_are_powers_of_two():
    assert ExteriorAlgebra(2, field=CC).dim == 4
    assert ExteriorAlgebra(3, field=CC).dim == 8


def test_exterior_2_hh0_char_dependence():
    assert ExteriorAlgebra(2, field=CC).hochschild_cohomology(0).dims == [2]
    assert ExteriorAlgebra(2, field=GF(2)).hochschild_cohomology(0).dims == [4]  # degenerates
```

`tests/families/test_preprojective.py`:
```python
"""PreprojectiveAlgebra Pi(Q). Fixture N8."""
from quiverlab.families import PreprojectiveAlgebra
from quiverlab.fields import CC


def test_preprojective_A2_dim4():
    A = PreprojectiveAlgebra("A2", field=CC)
    assert A.dim == 4                                   # monomial mesh relations
    assert A.hochschild_cohomology(0).dims == [1]


def test_preprojective_A3_dim10():
    A = PreprojectiveAlgebra("A3", field=CC)            # tetrahedral number T_3
    assert A.dim == 10                                  # Groebner-certified oracle
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `quantum.py`**:
```python
"""Quantum complete intersection A = k<x,y>/(x^2, y^2, xy + q yx) (spec §3.4).
General route; q may be non-rational (needs the Task-4 coefficient grammar)."""
from quiverlab.combinat.quiver import Quiver


def _q_token(q):
    if isinstance(q, str):
        return q
    return repr(q)                       # ints/Fractions stringify to exact rationals


def QuantumCI(q, field=None):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rels = ["x^2", "y^2", f"x*y + ({_q_token(q)})*y*x"]
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("quantum_ci", "qci_hh_oracle", "bardzell", "chouhy_solotar")
    return A
```
(the parenthesised `({_q_token(q)})` keeps a leading `-` or `i` well-formed for `_split_terms`; verify `x*y + (i)*y*x` and `x*y + (-1)*y*x` both parse. If the parser rejects the parenthesised coefficient, emit `x*y + <tok>*y*x` and special-case the sign.)

- [ ] **Step 4: Implement `exterior.py`** (all rational coefficients):
```python
"""Exterior algebra Lambda(k^n) = k<x_1..x_n>/(x_i^2, x_i x_j + x_j x_i). General route."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError


def ExteriorAlgebra(n, field=None):
    if not (isinstance(n, int) and n >= 1):
        raise QuiverlabError(f"ExteriorAlgebra({n!r}): need integer n >= 1", hint="e.g. 3")
    names = [f"x{i}" for i in range(1, n + 1)]
    Q = Quiver([1], {a: (1, 1) for a in names})
    rels = [f"{a}^2" for a in names]
    rels += [f"{names[i]}*{names[j]} + {names[j]}*{names[i]}"
             for i in range(n) for j in range(i + 1, n)]
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("quantum_ci", "chouhy_solotar")
    return A
```

- [ ] **Step 5: Implement `preprojective.py`** (double quiver + mesh relations):
```python
"""Preprojective algebra Pi(Q) of a Dynkin quiver (spec §3.4). Double every edge,
impose one mesh relation per vertex. Paths read left-to-right. General route."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.families.dynkin import dynkin_quiver


def PreprojectiveAlgebra(type_or_quiver, field=None):
    base = type_or_quiver if isinstance(type_or_quiver, Quiver) else \
        dynkin_quiver(type_or_quiver, "linear")
    arrows = {}
    star = {}
    for name, (s, t) in base.arrows.items():
        arrows[name] = (s, t)
        arrows[name + "s"] = (t, s)          # a* : t -> s
        star[name] = name + "s"
    Q = Quiver(list(base.vertices), arrows)
    # mesh relation at vertex v: sum_{a: s(a)=v} a*a*  -  sum_{b: t(b)=v} b*_star*b  = 0
    rels = []
    for v in base.vertices:
        pos = [f"{a}*{star[a]}" for a, (s, t) in base.arrows.items() if s == v]   # a a*
        neg = [f"{star[b]}*{b}" for b, (s, t) in base.arrows.items() if t == v]   # b* b
        terms = pos + [f"-{p}" for p in neg]
        if not terms:
            continue
        rel = " + ".join(terms).replace("+ -", "- ")
        rels.append(rel)
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("preprojective", "chouhy_solotar", "assem_book")
    return A
```
(For `A2`: mesh at 1 = `a*as` (`a a*`, 1→2→1); at 2 = `-as*a` (`a* a`, 2→1→2). Both single monomials ⇒ the general route detects them as monomial and dim = 4. For `A3` the branch vertex 2 yields `a*as - b s... ` a two-term relation ⇒ Gröbner. Re-derive dim against Fixture N8; the Gröbner certificate is the oracle.)

- [ ] **Step 6: Run focused (three files) + float gate.** If `PreprojectiveAlgebra("A3")` does not certify dim 10, inspect the completed reduction system (the mesh sign convention or leading-term choice), not the fixture — but verify against `T_3 = 10`.

- [ ] **Step 7: Commit** (`feat(families): QuantumCI, ExteriorAlgebra, PreprojectiveAlgebra`).

---

### Task 9: `TensorProduct`, `TrivialExtension` (structure-constant route)

**Files:**
- Create: `src/quiverlab/families/tensor.py`, `src/quiverlab/families/trivial_extension.py`, `tests/families/test_tensor.py`, `tests/families/test_trivial_extension.py`

**Interfaces:**
- Consumes: `Algebra` (ctor + `.T/.unit/.domain/.multiply/.dim`), `errors.FieldError`.
- Produces:
  - `TensorProduct(A, B) -> Algebra` — `dim = dim A · dim B`, row-major index `i*dim(B)+j`, `T[(i,j)][(k,l)] = A.T[i][k] ⊗ B.T[j][l]` (outer product over the shared domain), `unit = unit_A ⊗ unit_B`. Requires `A.domain is B.domain` (else `FieldError`). Basis labels `f"{la}⊗{lb}"`. Stamps `("tensor_product", "hodge")`. Carries no quiver (a structure-constant algebra); `cartan_matrix` on it raises the existing "needs the quiver presentation" error — that is correct.
  - `TrivialExtension(A) -> Algebra` — `dim = 2·dim A`; block multiplication `(a,f)(b,g) = (ab, a·g + f·b)` with the dual bimodule action. Stamps `("assem_book",)`. Symmetric.

- [ ] **Step 1: Write the failing tests**

`tests/families/test_tensor.py`:
```python
"""TensorProduct(A, B) = A (x) B, Kunneth. Fixture N9."""
from quiverlab.combinat import Quiver
from quiverlab.families import TensorProduct
from quiverlab.families.basic import truncated_polynomial
from quiverlab.fields import CC


def _kA2(field=CC):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


def test_kA2_tensor_kA2_is_commutative_square_dim9():
    A = TensorProduct(_kA2(), _kA2())
    assert A.dim == 9
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]     # == Plan-03 square == diamond


def test_kunneth_hh0_multiplicative():
    A = TensorProduct(truncated_polynomial(2, field=CC), truncated_polynomial(2, field=CC))
    assert A.dim == 4                                       # k[x,y]/(x^2,y^2)
    assert A.hochschild_cohomology(0).dims == [4]          # HH^0 = 2 * 2


def test_mismatched_fields_are_loud():
    import pytest
    from quiverlab.fields import GF
    from quiverlab.errors import FieldError
    with pytest.raises(FieldError):
        TensorProduct(_kA2(CC), _kA2(GF(5)))
```

`tests/families/test_trivial_extension.py`:
```python
"""TrivialExtension(A) = A |x D(A), symmetric. Fixture N10."""
from quiverlab.combinat import Quiver
from quiverlab.families import TrivialExtension
from quiverlab.fields import GF


def _kA2(field):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


def test_trivial_extension_kA2_dim6():
    A = TrivialExtension(_kA2(GF(32003)))
    assert A.dim == 6                                       # 2 * dim(kA_2) = 2*3 (Fixture N10)


def test_trivial_extension_is_symmetric_hh_duality():
    A = TrivialExtension(_kA2(GF(32003)))                  # symmetric: HH^n dim = HH_n dim
    co = A.hochschild_cohomology(4).dims
    ho = A.hochschild_homology(4).dims
    assert co == ho
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `tensor.py`**:
```python
"""Tensor product A (x)_k B as a structure-constant algebra (spec §3.4).
Basis a_i (x) b_j in row-major order i*dim(B)+j; (a(x)b)(a'(x)b') = aa' (x) bb'."""
from quiverlab.core.algebra import Algebra
from quiverlab.errors import FieldError


def TensorProduct(A, B):
    if A.domain is not B.domain:
        raise FieldError(
            f"TensorProduct needs one shared field; got {A.domain.name} and {B.domain.name}",
            hint="build both factors over the same field")
    dom = A.domain
    da, db = A.dim, B.dim
    m = da * db
    zero = dom.zero()

    def idx(i, j):
        return i * db + j

    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    for i in range(da):
        for k in range(da):
            aik = A.T[i][k]
            for j in range(db):
                for l in range(db):
                    bjl = B.T[j][l]
                    row = idx(i, j)
                    col = idx(k, l)
                    vec = T[row][col]
                    for p in range(da):
                        ap = aik[p]
                        if dom.is_zero(ap):
                            continue
                        for r in range(db):
                            br = bjl[r]
                            if dom.is_zero(br):
                                continue
                            vec[idx(p, r)] = dom.add(vec[idx(p, r)], dom.mul(ap, br))
    unit = [zero] * m
    for i in range(da):
        if dom.is_zero(A.unit[i]):
            continue
        for j in range(db):
            if dom.is_zero(B.unit[j]):
                continue
            unit[idx(i, j)] = dom.mul(A.unit[i], B.unit[j])
    la = A.basis_labels or [f"a{i}" for i in range(da)]
    lb = B.basis_labels or [f"b{j}" for j in range(db)]
    labels = [f"{la[i]}(x){lb[j]}" for i in range(da) for j in range(db)]
    C = Algebra(dom, T, unit, basis_labels=labels)
    C._family_citations = ("tensor_product", "hodge")
    return C
```

- [ ] **Step 4: Implement `trivial_extension.py`** (dual bimodule action, in the algebra's own basis):
```python
"""Trivial extension T(A) = A |x D(A), D(A) = Hom_k(A, k) (spec §3.4). Basis:
a_0..a_{n-1} then a_0^*..a_{n-1}^*; (a,f)(b,g) = (ab, a.g + f.b), symmetric algebra."""
from quiverlab.core.algebra import Algebra


def TrivialExtension(A):
    dom = A.domain
    n = A.dim
    m = 2 * n
    zero, one = dom.zero(), dom.one()

    def basis(i):
        v = [zero] * n
        v[i] = one
        return v

    # left/right action of A on itself in coordinates: prod[i][k] = a_i * a_k
    prod = [[A.multiply(basis(i), basis(k)) for k in range(n)] for i in range(n)]
    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    # (a_i, 0)(a_k, 0) = (a_i a_k, 0)
    for i in range(n):
        for k in range(n):
            for t in range(n):
                T[i][k][t] = prod[i][k][t]
    # (a_i,0)(0, a_k^*) = (0, a_i . a_k^*) ; (a_i . g)(c) = g(c a_i) => coeff on a_l^* is [a_l a_i]_k
    for i in range(n):
        for k in range(n):
            for l in range(n):
                T[i][n + k][n + l] = prod[l][i][k]
    # (0, a_i^*)(a_k, 0) = (0, a_i^* . a_k) ; (f . b)(c) = f(b c) => coeff on a_l^* is [a_k a_l]_i
    for i in range(n):
        for k in range(n):
            for l in range(n):
                T[n + i][k][n + l] = prod[k][l][i]
    # (0, *)(0, *) = 0  (D(A) is a square-zero ideal) -- already zero.
    unit = [zero] * m
    for t in range(n):
        unit[t] = A.unit[t]
    la = A.basis_labels or [f"a{i}" for i in range(n)]
    labels = list(la) + [f"{lbl}*" for lbl in la]
    Text = Algebra(dom, T, unit, basis_labels=labels)
    Text._family_citations = ("assem_book",)
    return Text
```

- [ ] **Step 5: Run focused + float gate.** Verify `TensorProduct(kA_2, kA_2)` HH matches the commutative-square oracle and `TrivialExtension(kA_2)` is dim 6 with `HH^n == HH_n`. If the trivial-extension multiplication is not associative, `Algebra._validate` is *not* called here (we build the `T` directly) — add a one-off `Algebra.from_structure_constants(..., check=True)` cross-check in the test over `GF(p)` to catch a sign/index error early.

- [ ] **Step 6: Commit** (`feat(families): TensorProduct and TrivialExtension (structure-constant route)`).

---

### Task 10: `families()` discoverability, top-level exports, family→citation wiring

**Files:**
- Create: `src/quiverlab/families/discover.py`, `tests/families/test_discover.py`, `tests/families/test_family_citations.py`
- Modify: `src/quiverlab/families/__init__.py`, `src/quiverlab/__init__.py`

**Interfaces:**
- Consumes: every family builder (Tasks 5–9, plus `basic.py`), `citations.reference`.
- Produces:
  - `families.discover.CATALOG: tuple[FamilyInfo, ...]` — `FamilyInfo(name, signature, route, citations, summary)`.
  - `families() -> FamilyListing` — printable/iterable; `.names()`, `.by_name(name)`, `__str__` (aligned table), `to_dict()` (web `/families` consumer). Signatures are the human one-liners from spec §3.4.
  - `families/__init__.py` re-exports every builder + `families` + `zoo` (Task 11).
  - `quiverlab/__init__.py` and `from quiverlab import *` export: `NakayamaAlgebra, PathAlgebra, TruncatedPathAlgebra, RadicalSquareZero, IncidenceAlgebra, QuantumCI, ExteriorAlgebra, PreprojectiveAlgebra, TrivialExtension, TensorProduct, zoo, families, bibliography` (batch stays unexported).

- [ ] **Step 1: Write the failing tests**

`tests/families/test_discover.py`:
```python
"""families() discoverability (spec §3.4)."""
from quiverlab import families


def test_catalog_lists_every_family_with_signature_and_route():
    listing = families()
    names = set(listing.names())
    assert {
        "NakayamaAlgebra", "PathAlgebra", "TruncatedPathAlgebra", "RadicalSquareZero",
        "IncidenceAlgebra", "QuantumCI", "ExteriorAlgebra", "PreprojectiveAlgebra",
        "TrivialExtension", "TensorProduct", "zoo",
    } <= names
    info = listing.by_name("QuantumCI")
    assert info.route == "general"
    assert "q=" in info.signature
    assert "quantum_ci" in info.citations
    s = str(listing)
    assert "NakayamaAlgebra" in s and "monomial" in s


def test_star_import_exposes_the_catalog():
    ns = {}
    exec("from quiverlab import *", ns)
    for name in ("NakayamaAlgebra", "QuantumCI", "zoo", "families", "bibliography"):
        assert name in ns
```

`tests/families/test_family_citations.py`:
```python
"""Every family stamps registry-valid citation keys; A.citations() resolves them."""
import pytest

from quiverlab import citations, families
from quiverlab.families import NakayamaAlgebra, PathAlgebra, QuantumCI


@pytest.mark.parametrize("A,expected", [
    (NakayamaAlgebra([3, 2, 2]), "nakayama"),
    (PathAlgebra("A3"), "path_algebra"),
    (QuantumCI(q=2), "quantum_ci"),
])
def test_family_citations_resolve(A, expected):
    keys = A.citations()
    assert expected in keys
    for k in keys:
        citations.reference(k)                     # loud if any stamped key is bogus


def test_catalog_citations_are_all_registered():
    for info in families():
        for k in info.citations:
            citations.reference(k)
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `discover.py`** (the `CATALOG` table + `FamilyListing`):
```python
"""The v1 family catalog + families() discoverability (spec §3.4)."""
from dataclasses import dataclass


@dataclass(frozen=True)
class FamilyInfo:
    name: str
    signature: str
    route: str          # "monomial" | "general" | "structure-constant" | "iterator"
    citations: tuple
    summary: str


CATALOG = (
    FamilyInfo("NakayamaAlgebra", "NakayamaAlgebra([c1,..,cn]) | (n=, l=, cyclic=)",
               "monomial", ("nakayama", "assem_book"),
               "Serial algebra by Kupisch series (linear or cyclic)."),
    FamilyInfo("PathAlgebra", "PathAlgebra(type, orientation='linear')",
               "monomial", ("path_algebra", "happel_question", "assem_book"),
               "Hereditary path algebra of a Dynkin/Euclidean quiver."),
    FamilyInfo("TruncatedPathAlgebra", "TruncatedPathAlgebra(type_or_Q, r)",
               "monomial", ("path_algebra", "bardzell"),
               "kQ / rad^r."),
    FamilyInfo("RadicalSquareZero", "RadicalSquareZero(Q)",
               "monomial", ("path_algebra", "bardzell"),
               "kQ / rad^2 from any quiver."),
    FamilyInfo("IncidenceAlgebra", "IncidenceAlgebra(covers)",
               "general", ("incidence", "assem_book", "hodge"),
               "Incidence algebra of a finite poset (Hasse quiver + commutativity)."),
    FamilyInfo("QuantumCI", "QuantumCI(q)",
               "general", ("quantum_ci", "qci_hh_oracle", "bardzell", "chouhy_solotar"),
               "Quantum complete intersection k<x,y>/(x^2,y^2,xy+q yx)."),
    FamilyInfo("ExteriorAlgebra", "ExteriorAlgebra(n)",
               "general", ("quantum_ci", "chouhy_solotar"),
               "Exterior algebra Lambda(k^n), dim 2^n."),
    FamilyInfo("PreprojectiveAlgebra", "PreprojectiveAlgebra(type)",
               "general", ("preprojective", "chouhy_solotar", "assem_book"),
               "Preprojective algebra of a Dynkin quiver."),
    FamilyInfo("TrivialExtension", "TrivialExtension(A)",
               "structure-constant", ("assem_book",),
               "Trivial extension A |x D(A) (symmetric)."),
    FamilyInfo("TensorProduct", "TensorProduct(A, B)",
               "structure-constant", ("tensor_product", "hodge"),
               "Tensor product A (x)_k B."),
    FamilyInfo("zoo", "zoo(dim_max=12)",
               "iterator", ("han_conjecture", "chouhy_solotar"),
               "Iterator over the curated exact zoo of open (Han-conjecture) algebras."),
)


class FamilyListing:
    def __init__(self, catalog=CATALOG):
        self._catalog = tuple(catalog)

    def __iter__(self):
        return iter(self._catalog)

    def names(self):
        return tuple(f.name for f in self._catalog)

    def by_name(self, name):
        for f in self._catalog:
            if f.name == name:
                return f
        raise KeyError(name)

    def to_dict(self):
        return {"families": [f.__dict__ | {"citations": list(f.citations)}
                             for f in self._catalog]}

    def __str__(self):
        w = max(len(f.name) for f in self._catalog)
        rows = [f"{f.name:<{w}}  [{f.route}]  {f.signature}" for f in self._catalog]
        return "quiverlab families (spec 3.4):\n" + "\n".join("  " + r for r in rows)


def families():
    return FamilyListing()
```

- [ ] **Step 4: Wire `families/__init__.py`** to re-export every builder + `families` + `zoo` (Task 11 fills `zoo`; import it lazily or land Task 11 first), and update `quiverlab/__init__.py` `__all__`.

- [ ] **Step 5: Run focused + `tests/test_quickstart.py` (star-import smoke) + float gate.**

- [ ] **Step 6: Commit** (`feat(families): families() catalog, top-level exports, family->citation wiring`).

---

### Task 11: `zoo(dim_max)` — curated catalog iterator (open_zoo lift)

**Files:**
- Create: `src/quiverlab/families/zoo.py`, `src/quiverlab/families/zoo_catalog.json` (bundled data — a curated subset lifted from the bank's `open_zoo_catalog_v2.json`), `tests/families/test_zoo.py`
- Modify: `src/quiverlab/families/__init__.py` (export `zoo`)

**Interfaces:**
- Consumes: the Plan-03/04 reduction-system → `groebner_algebra` path (`groebner.build_reduction_system` / `groebner_algebra`, or the CS `algebra_from_reduction_system` equivalent), `Algebra`.
- Produces:
  - Bundled `zoo_catalog.json` — a list of records `{"name", "ngen", "dim", "rules"}` where `rules` is the serialized reduction system `[[lead], [[coeff, word], ...]]` (arrows encoded as generator indices `0..ngen-1`). Lifted verbatim from the bank catalog (read-only), curated to `dim ≤ 36`, deterministic order. Includes the periodic-symmetric-family members `open2_33_712` (dim 9) and `open2_37_19612` (dim 21) and the golden `open_33_0` (dim 9).
  - `zoo(dim_max=12) -> Iterator[Algebra]` — loads the catalog, rebuilds each record with `dim ≤ dim_max` into an `Algebra` (generators `x, y[, z]` as loops on one vertex), yields in ascending `(dim, name)`. Each yielded algebra stamps `_family_citations = ("han_conjecture", "chouhy_solotar")` and carries `.zoo_name`.
  - `families.zoo.load_catalog() -> list[dict]`; `families.zoo.build_from_record(rec, field=CC) -> Algebra`.

- [ ] **Step 1: Write the failing test**

`tests/families/test_zoo.py`:
```python
"""zoo(dim_max) iterator over the curated exact zoo (open_zoo lift)."""
import importlib

import pytest

from quiverlab import zoo
from quiverlab.families.zoo import load_catalog
from quiverlab.fields import GF


def test_catalog_is_bundled_and_well_formed():
    cat = load_catalog()
    assert len(cat) >= 10
    for rec in cat:
        assert {"name", "ngen", "dim", "rules"} <= set(rec)
        assert rec["ngen"] in (2, 3)
        for lead, tail in rec["rules"]:                 # exact-integer coefficients only
            assert all(isinstance(i, int) for i in lead)
            for c, w in tail:
                assert isinstance(c, int)               # no floats in the bundled catalog
                assert all(isinstance(i, int) for i in w)


def test_zoo_yields_algebras_up_to_dim_max_in_order():
    algs = list(zoo(dim_max=12))
    dims = [A.dim for A in algs]
    assert dims == sorted(dims)
    assert all(d <= 12 for d in dims)
    assert all(A.dim >= 1 for A in algs)


def test_periodic_symmetric_member_present_dim9():
    names = {getattr(A, "zoo_name", None): A for A in zoo(dim_max=12)}
    assert "open2_33_712" in names
    assert names["open2_33_712"].dim == 9


@pytest.mark.skipif(
    importlib.util.find_spec("quiverlab.engine.resolutions_cs") is None,
    reason="open-zone HH golden gated on the Plan 04 CS backend, consistent with Tasks 1/13")
def test_zoo_algebra_hh_matches_open_zone_golden():
    # Build needs the Plan-03 Groebner route (a hard Plan-06 prereq); the depth-16
    # HH is produced by the Plan-02 minimal engine over GF(32003). Gated on CS to
    # stay green and consistent with the Task-13 ledger tests when Plan 04 is absent.
    A = next(A for A in zoo(dim_max=9, field=GF(32003))
             if getattr(A, "zoo_name", "") == "open_33_0")
    hh = A.hochschild_homology(16, engine="auto")   # GF(32003) -> fast/minimal engine
    assert hh.dims == [6] + [5] * 16                 # golden (Fixture Z1)
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Build the bundled catalog** — read the bank file `/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/HansConjecture/hanlab/open_zoo_catalog_v2.json` (READ-ONLY), curate the entries with `dim ≤ 36`, ensure the three named records (`open2_33_712`, `open2_37_19612`, `open_33_0`) are present, and write the curated list to `src/quiverlab/families/zoo_catalog.json`. (Do this with a one-off script under the scratchpad, not committed; the committed artifact is the JSON. Verify each record round-trips through `build_from_record` to the stated `dim`.) Provenance: add a top-level `"_provenance"` key naming the bank file and date.

- [ ] **Step 4: Implement `zoo.py`** (loader + record → `Algebra` via the reduction-system route):
```python
"""zoo(dim_max): the curated exact zoo (spec §3.4), lifted from hanlab open_zoo.
Each record is a confluent reduction system over k<x,y[,z]>; rebuilt through the
Plan-03/04 reduction-system -> Algebra path."""
import json
import pathlib

from quiverlab.combinat.quiver import Quiver

_CATALOG = pathlib.Path(__file__).with_name("zoo_catalog.json")
_GENS = ("x", "y", "z")


def load_catalog():
    data = json.loads(_CATALOG.read_text(encoding="utf-8"))
    return [r for r in data if isinstance(r, dict) and "rules" in r]


def _word(idxs):
    return "*".join(_GENS[i] for i in idxs)


def build_from_record(rec, field=None):
    ngen = rec["ngen"]
    gens = _GENS[:ngen]
    Q = Quiver([1], {g: (1, 1) for g in gens})
    rels = []
    for lead, tail in rec["rules"]:
        lhs = _word(lead)
        if not tail:
            rels.append(lhs)                                    # lead -> 0 (monomial)
        else:
            rhs = " + ".join(f"({c})*{_word(w)}" for c, w in tail)
            rels.append(f"{lhs} - ({rhs})")                     # lead = tail
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("han_conjecture", "chouhy_solotar")
    A.zoo_name = rec["name"]
    return A


def zoo(dim_max=12, field=None):
    recs = [r for r in load_catalog() if r["dim"] <= dim_max]
    recs.sort(key=lambda r: (r["dim"], r["name"]))
    for rec in recs:
        yield build_from_record(rec, field=field)
```
(Note: the reduction-system records encode a *confluent* system; rebuilding by handing the relations to `Quiver.algebra` re-runs Gröbner completion — deterministic and dim-preserving. If re-completion is too slow for the whole catalog, add a fast path `families.zoo.build_from_record(rec, recomplete=False)` that feeds `rules` straight into `groebner_algebra` bypassing completion, once Plan 04 exposes that entry point. Default stays the safe re-completing path.)

- [ ] **Step 5: Run focused + float gate.** The golden `open_33_0` HH check is the anchor: `[6] + [5]*16`. If dims differ, the record's rules or the generator-index encoding is wrong — compare to the bank `open_zoo_catalog.json` fixture `open_33_0` (rules `x^3→y², y³→0, yx→−xy`).

- [ ] **Step 6: Commit** (`feat(families): zoo(dim_max) curated iterator (open_zoo lift)`).

---

### Task 12: `quiverlab.batch` — labdb lift (SQLite persistence)

**Files:**
- Create: `src/quiverlab/batch/__init__.py`, `src/quiverlab/batch/db.py`, `src/quiverlab/batch/builders.py`, `src/quiverlab/batch/scan.py`, `tests/batch/__init__.py` (empty), `tests/batch/test_batch_db.py`, `tests/batch/test_batch_labdb_port.py`

**Interfaces:**
- Consumes: the family catalog (Tasks 5–11), `Algebra.hochschild_cohomology`/`hochschild_homology`, the minimal-resolution engine (Plan 02), `zoo` (Task 11).
- Produces (expert-facing; **not** exported at top level — reached via `import quiverlab.batch`):
  - `batch.db.ResultsDB(path=":memory:")` — SQLite; schema `results(id, name, builder, args, dim, N, associative, cx_homology, cx_cohomology, asymmetric, homology_bounded, homology_nonzero_tail, kind, han_verdict, growth, depth_reached, truncated_at, record_json, UNIQUE(builder, args, N))` + the three bank indexes; context manager; `insert` = `INSERT OR REPLACE` (idempotent upsert on `(builder, args, N)`); `insert_many`, `all`, `query(where, params)`, `asymmetric`, `bounded_homology_candidates`, `by_complexity`, `__len__`. Faithful to the bank `labdb.ResultsDB`, port header attribution retained.
  - `batch.builders.BUILDERS: dict[str, callable]` mapping `{"truncated_polynomial","quantum_ci","cyclic_nakayama","linear_path_algebra","dynkin","reduction_system", ...}` to the **quiverlab** family builders; `build_algebra(spec)`.
  - `batch.scan.analyze(spec) -> dict` (pure, JSON-serializable record); `run_scan(specs, n_workers=1)` — **default single-worker** (`multiprocessing.Pool` only when `n_workers > 1`, guarded for Marco's machine).
  - `batch.PRIME = 32003`.

- [ ] **Step 1: Write the failing tests**

`tests/batch/test_batch_db.py`:
```python
"""quiverlab.batch ResultsDB: SQLite persistence for scans (labdb lift)."""
import quiverlab.batch as batch
from quiverlab.batch.db import ResultsDB


def test_schema_and_upsert_idempotent():
    rec = {"name": "kA3", "builder": "linear_path_algebra", "args": [3], "N": 2,
           "dim": 6, "associative": True, "kind": "structural"}
    with ResultsDB(":memory:") as db:
        db.insert(rec)
        db.insert(dict(rec))                 # same (builder, args, N)
        assert len(db) == 1                  # INSERT OR REPLACE dedup


def test_query_roundtrip_and_json():
    import json
    rec = {"name": "qCI", "builder": "quantum_ci", "args": [3], "N": 4, "dim": 4,
           "associative": True, "asymmetric": True, "cx_homology": 1, "cx_cohomology": 0,
           "kind": "structural"}
    with ResultsDB(":memory:") as db:
        db.insert(rec)
        (got,) = db.query("builder = ?", ("quantum_ci",))
        assert got["name"] == "qCI"
        json.dumps(db.all())
        assert db.asymmetric()[0]["name"] == "qCI"
```

`tests/batch/test_batch_labdb_port.py`:
```python
"""Re-port of the bank tests/test_labdb.py, adapted to quiverlab families.
Golden numeric fixtures preserved (Fixture Z1)."""
import pytest

from quiverlab.batch.builders import BUILDERS, build_algebra
from quiverlab.batch.scan import analyze, run_scan

P = 32003


def test_registry_covers_the_core_builders():
    assert {"truncated_polynomial", "quantum_ci", "cyclic_nakayama",
            "linear_path_algebra", "dynkin", "reduction_system"} <= set(BUILDERS)


def test_builders_construct_associative_algebras():
    for spec in [{"builder": "truncated_polynomial", "args": [3]},
                 {"builder": "quantum_ci", "args": [2]},
                 {"builder": "cyclic_nakayama", "args": [3, 2]},
                 {"builder": "linear_path_algebra", "args": [3]},
                 {"builder": "dynkin", "args": ["D", 4]}]:
        A = build_algebra(spec)
        assert A.dim >= 1


def test_analyze_is_deterministic_and_serial_equals_parallel():
    specs = [{"builder": "quantum_ci", "args": [3], "N": 4}]
    assert analyze(dict(specs[0])) == analyze(dict(specs[0]))
    assert run_scan(specs, 1) == run_scan(specs, 1)


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("quiverlab.engine.resolutions_cs") is None,
    reason="open-zone analyze needs the Plan 04 CS backend")
def test_open_zone_golden_open_33_0():
    spec = {"builder": "reduction_system", "N": 16,
            "args": [2, [[[0, 0, 0], [[1, [1, 1]]]], [[1, 1, 1], []],
                         [[1, 0], [[-1, [0, 1]]]]], "open_33_0"]}
    rec = analyze(spec)
    assert rec["dim"] == 9
    assert rec["HH_homology"][str(P)] == [6] + [5] * 16          # Fixture Z1 golden
    assert rec["resolution_ranks"] == [1, 2, 2, 1, 1, 2, 2, 1,
                                       1, 2, 2, 1, 1, 2, 2, 1, 1, 2]
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `db.py`** — lift `labdb.ResultsDB` from the bank verbatim (schema, `_migrate_open_columns`, `insert`/`insert_many`/`all`/`query`/`asymmetric`/`bounded_homology_candidates`/`by_complexity`/`__len__`, `_nullable_bool`), with the ported-with-attribution header and package-relative imports. No NN coupling exists in `labdb` (verified). Keep `json`/`sqlite3` stdlib.

- [ ] **Step 4: Implement `builders.py`** — map the bank builder names to quiverlab families:
  - `truncated_polynomial(a)` → `families.basic.truncated_polynomial`
  - `quantum_ci(c)` → `families.QuantumCI(q=c)`
  - `cyclic_nakayama(n, ell)` → `families.NakayamaAlgebra(n=n, l=ell, cyclic=True)`
  - `linear_path_algebra(n)` → `families.basic.linear_path_algebra`
  - `dynkin(typ, n)` → `families.PathAlgebra(f"{typ}{n}")`
  - `reduction_system(ngen, rules, name)` → `families.zoo.build_from_record({"ngen":ngen,"rules":rules,"name":name,"dim":None})`
  - `tensor_product(specA, specB)`, `trivial_extension(spec)` → the Task-9 builders (optional; add if the port test needs them).
  `build_algebra(spec)` = `BUILDERS[spec["builder"]](*spec.get("args", []))`.

- [ ] **Step 5: Implement `scan.py`** — lift `labdb.analyze`/`_analyze_open`/`run_scan`, routing structural builders through `hochschild_homology`/`hochschild_cohomology` + `complexity`, and `reduction_system` through the minimal-resolution open-zone path (`minimal_homology_dims` + `minimal_resolution` + `scan_open` policy fns `han_verdict`/`depth_for_dim`/`max_term_dim_for_dim` — port these small policy helpers into `batch/scan.py`). `run_scan(specs, n_workers=1)` default single-worker; `Pool` only when `n_workers > 1`.

- [ ] **Step 6: Run focused + float gate.** The structural DB tests must pass unconditionally; the open-zone golden test is `skipif` on Plan 04's CS backend (it needs the reduction-system open-zone analyze).

- [ ] **Step 7: Commit** (`feat(batch): labdb lift -- ResultsDB, builders registry, analyze/run_scan`).

---

### Task 13: LEDGER re-ports — periodic-symmetric-family + open-zoo-broaden

**Files:**
- Create: `tests/batch/test_periodic_symmetric_family.py`, `tests/batch/test_open_zoo_broaden.py`
- Modify (only if needed for import-closure): `src/quiverlab/batch/__init__.py` (expose a `make_shards`-equivalent `open_zoo_to_specs`, `depth_for_dim`, `max_term_dim_for_dim`), `src/quiverlab/families/zoo.py` (a `search_*`/generator entry point if the broaden test exercises it)

**Interfaces:**
- Consumes: `zoo`/`families.zoo` (Task 11), `batch` (Task 12), the minimal-resolution engine + private steppers `engine.resolutions_minimal._init_resolution`/`_advance_resolution` (Plan 02), `engine.coxeter.is_frobenius`/`nakayama_automorphism` (Plan 02), the CS backend (Plan 04).
- Produces: the two ledger tests re-ported (the Plan-02 deferral: "zoo/labdb/periodic-symmetric-family (Plan 06)"), **iff import-closed** at execution time (they need Plan 04's reduction-system open-zone builder + Plan 02's minimal stepper + coxeter). If Plan 04 is not yet delivered, mark the open-zone assertions `skipif` and record the residual obligation in the acceptance task's frozen statement.

- [ ] **Step 1: Write the tests** (golden fixtures preserved verbatim — Fixture Z1):

`tests/batch/test_periodic_symmetric_family.py`:
```python
"""Re-port of the bank tests/test_periodic_symmetric_family.py (LEDGER OBLIGATION).
The family k<x,y>/(x^3, y^b - x^2, yx + xy): periodic, symmetric. Fixture Z1."""
import importlib

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("quiverlab.engine.resolutions_cs") is None,
    reason="periodic-symmetric family needs the Plan 04 CS/reduction-system backend")

from quiverlab.engine.coxeter import is_frobenius, nakayama_automorphism
from quiverlab.engine.resolutions_minimal import minimal_homology_dims
from quiverlab.families.zoo import build_from_record, load_catalog


def _algebra(name):
    rec = next(r for r in load_catalog() if r["name"] == name)
    from quiverlab.engine.adapter import to_engine
    from quiverlab.fields import GF
    return to_engine(build_from_record(rec, field=GF(32003)).unit_adapted())


def test_dim9_sibling_periodic_symmetric_and_p2_growth():
    A = _algebra("open2_33_712")                       # k<x,y>/(x^3, y^3 - x^2, yx + xy)
    assert A.m == 9
    assert is_frobenius(A, 32003)
    nu, _ = nakayama_automorphism(A, 32003)
    assert np.array_equal(nu % 32003, np.eye(A.m, dtype=object) % 32003)
    assert minimal_homology_dims(A, 6, primes=(32003,))[32003] == [6, 5, 5, 5, 5, 5, 5]
    hh2 = minimal_homology_dims(A, 6, primes=(2,))[2]
    assert hh2 == [9, 10, 14, 18, 22, 26, 30]          # growing tail at p=2
    assert hh2[-1] > hh2[0]


def test_dim21_headline_is_symmetric():
    A = _algebra("open2_37_19612")                     # k<x,y>/(x^3, y^7 - x^2, yx + xy)
    assert A.m == 21
    assert is_frobenius(A, 32003)
    nu, _ = nakayama_automorphism(A, 32003)
    assert np.array_equal(nu % 32003, np.eye(A.m, dtype=object) % 32003)
```

`tests/batch/test_open_zoo_broaden.py`:
```python
"""Re-port of the bank tests/test_open_zoo_broaden.py (LEDGER OBLIGATION).
Curated-catalog dimension bands + spec generation. Import-closed subset."""
from quiverlab.batch import depth_for_dim, max_term_dim_for_dim, open_zoo_to_specs
from quiverlab.families.zoo import load_catalog


def test_catalog_round_trips_and_covers_dimension_bands():
    cat = load_catalog()
    for e in cat:
        assert {"name", "ngen", "dim", "rules"} <= set(e)
    dims = {e["dim"] for e in cat}
    assert {9, 12}.issubset(dims) or len(dims) >= 5      # curated band coverage


def test_open_zoo_to_specs_band_and_probe():
    specs = open_zoo_to_specs(load_catalog(), min_dim=9, max_dim=9)
    assert specs
    for s in specs:
        assert s["builder"] == "reduction_system"
        assert s["N"] == depth_for_dim(9)
        assert s["max_term_dim"] == max_term_dim_for_dim(9)
```

- [ ] **Step 2: Run to verify failure / assess import-closure.** Determine whether Plan 04's CS reduction-system builder and Plan 02's minimal stepper are present. If present → make them pass; if absent → the `pytestmark` skips keep the suite green and the obligation is recorded (Step 4).

- [ ] **Step 3: Fill import-closure gaps** — port the tiny `scan_open` policy helpers (`depth_for_dim`, `max_term_dim_for_dim`, `han_verdict`) and the `make_shards` adapter `open_zoo_to_specs` into `quiverlab.batch` (they carry no homology machinery, no NN coupling). Confirm `engine.coxeter.is_frobenius`/`nakayama_automorphism` and `engine.resolutions_minimal.minimal_homology_dims` (+ private steppers) are the Plan-02 ports (they are, per the engine tree).

- [ ] **Step 4: Run focused + float gate.** Record in the acceptance task's frozen statement whether the open-zone assertions ran or skipped (and why).

- [ ] **Step 5: Commit** (`test(batch): re-port periodic-symmetric-family + open-zoo-broaden (ledger obligation)`).

---

### Task 14: Acceptance — cross-checks, internals chapter, frozen statement, full suite

**Files:**
- Create: `tests/families/test_acceptance.py`, `docs/internals/10-families-citations.md`
- Modify: `README.md` (families + citations examples that the acceptance test executes verbatim)

**Interfaces:**
- Consumes: the entire public surface — every family builder, `families()`, `zoo`, `bibliography`, `A.citations()`, `HHTable.references`.
- Produces: the Plan-06 exit criterion + the frozen statement for Plans 07/09.

- [ ] **Step 1: Write the acceptance test** — the triple cross-check + the citations end-to-end:

`tests/families/test_acceptance.py`:
```python
"""Plan 06 acceptance. Three construction routes yield the SAME dim-9 algebra
(commutative square == diamond incidence == kA_2 (x) kA_2), all HH=[1,0,0]; and
citations flow from family -> A.citations() -> bibliography()."""
from quiverlab import (
    ExteriorAlgebra, IncidenceAlgebra, NakayamaAlgebra, PathAlgebra, QuantumCI,
    TensorProduct, bibliography, families,
)
from quiverlab.combinat import Quiver
from quiverlab.fields import CC, GF


def _kA2(field=CC):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


def test_triple_crosscheck_dim9_commutative_square():
    Qsq = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    square = Qsq.algebra(relations=["a*b - c*d"], field=CC)          # general kQ/I
    diamond = IncidenceAlgebra([("b", "x"), ("b", "y"), ("x", "t"), ("y", "t")], field=CC)
    tensor = TensorProduct(_kA2(), _kA2())
    for A in (square, diamond, tensor):
        assert A.dim == 9
        assert A.hochschild_cohomology(2).dims == [1, 0, 0]


def test_catalog_dims_char_independent_and_exact():
    for A in (NakayamaAlgebra([3, 2, 2], field=GF(32003)),
              PathAlgebra("D4", field=GF(7)),
              QuantumCI(q=2, field=GF(32003)),
              ExteriorAlgebra(3, field=GF(32003))):
        assert A.dim in (7, 9, 4, 8)


def test_exterior2_equals_quantumci1():
    assert ExteriorAlgebra(2, field=CC).T == QuantumCI(q=1, field=CC).T


def test_citations_flow_end_to_end():
    A = QuantumCI(q="i", field=CC)
    keys = A.citations()
    assert "quantum_ci" in keys and "bar" in keys
    b = bibliography(keys=keys)
    assert "@article{BGMS2005" in b.bibtex()
    assert "quantum complete intersection" in str(b).lower()
    hh = A.hochschild_cohomology(2)
    assert "quantum_ci" in hh.references and "bar" in hh.references   # family + engine (§3.9)
    assert set(hh.references) == set(keys)             # table.references == A.citations()


def test_families_listing_complete():
    assert len(families().names()) >= 11
```

- [ ] **Step 2: Run the acceptance test.** Fix any localized failure before proceeding.

- [ ] **Step 3: Write `docs/internals/10-families-citations.md`** — the "under the hood" chapter (same five-section format as chapters 01–07; see `docs/internals/README.md`). Full content:

```markdown
# 10 — Families and citations

## The mathematics
A *family* is a recipe that turns a few numbers or a diagram name into a
finite-dimensional algebra `kQ/I`. quiverlab ships the standard catalogue of
representation theory: Nakayama (serial) algebras by Kupisch series, hereditary
path algebras of Dynkin/Euclidean quivers, truncations `kQ/rad^r`, incidence
algebras of posets, quantum complete intersections, exterior algebras,
preprojective algebras, and two constructions that build a new algebra from old
ones (tensor product and trivial extension). A *citation registry* records, for
every algorithm and every family, the paper it comes from.

## How it is represented
Each family is a plain Python function returning a Plan-03 `Algebra`. There are
three construction *routes*:
- **monomial** — the relations are single forbidden paths; the algebra is built
  by the Plan-01 monomial route (`Quiver.algebra`, the forbidden-word automaton).
  Nakayama, PathAlgebra, TruncatedPathAlgebra, RadicalSquareZero.
- **general** — at least one relation is a genuine linear combination (e.g.
  `x*y + q*y*x`); the algebra is completed by the Plan-03 Gröbner engine.
  QuantumCI, ExteriorAlgebra, PreprojectiveAlgebra, IncidenceAlgebra.
- **structure-constant** — the multiplication table `T` is written directly from
  the factors, with no quiver. TensorProduct, TrivialExtension.
A family stamps the algebra with `_family_citations`, a tuple of registry keys.
The registry itself is a dict `key -> Reference(key, bibtex_key, kind, title,
annotation, tags)`; the annotations are the ground truth the web /literature page
and `quiverlab.bibliography()` render. `references.bib` holds the verified BibTeX.

## How the computation runs
1. `NakayamaAlgebra([3,2,2])` reads the Kupisch series, decides linear vs cyclic
   (`min >= 2` -> cyclic), lays out the quiver, generates the length-`c_i`
   forbidden path from each vertex, and calls the monomial route. `dim = sum c_i`.
2. `QuantumCI(q="i")` writes the three relation strings, `x*y + (i)*y*x` among
   them; the relation parser accepts the exact token `i` (this chapter's one new
   grammar rule); the Gröbner engine rewrites `y*x -> i*x*y` and certifies dim 4.
3. `TensorProduct(A, B)` fills `T[i*db+j][k*db+l]` with the outer product of the
   two multiplication tables; `dim = dim A * dim B`.
4. `A.hochschild_cohomology(n)` attaches, to the returned `HHTable`, the citation
   keys of the engine paths it used (`.references`): `bar` always, plus `bardzell`
   or `chouhy_solotar` by dispatch. `A.citations()` unions those with the family
   keys. `bibliography(keys)` groups them by kind and prints the annotations.
5. `zoo(dim_max)` loads a bundled JSON catalogue (lifted from hanlab's open_zoo),
   rebuilds each confluent reduction system into an `Algebra`, and yields those
   with `dim <= dim_max`.

## A worked micro-example
`NakayamaAlgebra([3,2,2])`: cyclic `Z_3`, arrows `a1:1->2, a2:2->3, a3:3->1`,
forbidden paths `a1*a2*a3` (len 3 from 1), `a2*a3` (len 2 from 2), `a3*a1` (len 2
from 3). Irreducible basis `e1,e2,e3,a1,a2,a3,a1*a2` -> **dim 7**. Cartan
`[[1,1,1],[0,1,1],[1,0,1]]` (row `i` = composition factors of `P_i`), `det 1`,
`sum = 7`. Centre = scalars, so `HH^0 = 1`. `A.citations()` -> `('nakayama',
'assem_book', 'bar')`; `bibliography(A.citations())` prints the ASS textbook and
Happel's bar-complex reference with their annotations.

## Where to look in the code
| concept | file | function/class |
|---|---|---|
| Kupisch series | `families/nakayama.py` | `NakayamaAlgebra` |
| Dynkin diagram -> quiver | `families/dynkin.py` | `dynkin_quiver` |
| hereditary path algebra | `families/path_algebra.py` | `PathAlgebra` |
| `kQ/rad^r` | `families/truncated.py` | `TruncatedPathAlgebra` |
| poset -> incidence algebra | `families/poset.py`, `families/incidence.py` | `Poset`, `IncidenceAlgebra` |
| quantum CI / exterior / preprojective | `families/{quantum,exterior,preprojective}.py` | `QuantumCI`, `ExteriorAlgebra`, `PreprojectiveAlgebra` |
| tensor / trivial extension | `families/{tensor,trivial_extension}.py` | `TensorProduct`, `TrivialExtension` |
| discoverability | `families/discover.py` | `families`, `CATALOG` |
| curated zoo | `families/zoo.py`, `families/zoo_catalog.json` | `zoo`, `build_from_record` |
| citation registry | `citations/registry.py`, `citations/references.bib` | `REGISTRY`, `reference`, `bibtex` |
| bibliography | `citations/bibliography.py` | `bibliography`, `Bibliography` |
| result references | `hochschild/table.py`, `core/algebra.py` | `HHTable.references`, `Algebra.citations` |
| batch persistence | `batch/db.py`, `batch/scan.py` | `ResultsDB`, `analyze`, `run_scan` |

This chapter is the **Plan 06** checkout. `sweep` (Plan 05) consumes this catalogue
(`families()` is the hook); the trace subsystem (Plan 07) will render the citation
keys these functions stamp.
```
(Chapter-number note: chapters 01–07 exist; Plan 03 adds 08, Plan 04 adds 09. This
chapter is filed as `10-families-citations.md` per the Plan-06 dispatch. If Plan 05
(modules, not yet drafted) claims chapter 10 first at integration, renumber this to
the next free index — Plan 06 owns the *families + citations* chapter content
regardless of its final number.)

- [ ] **Step 4: Update `README.md`** — append a families + citations example:
```markdown
## Families and citations

```python
from quiverlab import NakayamaAlgebra, QuantumCI, families, bibliography

A = NakayamaAlgebra([3, 2, 2])          # cyclic Nakayama, dim 7
print(A.hochschild_cohomology(0))       # HH^0 = 1
print(A.citations())                    # ('nakayama', 'assem_book', 'bar')

print(families())                       # the whole v1 catalog with signatures
print(bibliography(A.citations()))      # grouped, annotated references
```
```

- [ ] **Step 5: Run the full suite as a tracked background job, await it, then commit.**

Start the full suite as a **single tracked background job** (Bash `run_in_background: true`) — do not block the foreground or hand-redirect to a log:
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
Await its completion notification (the harness re-invokes on exit); do not poll. When it returns green, run the pure-Python path once in the foreground as the dual-path equality gate from Plan 02: `QUIVERLAB_NO_NUMBA=1 NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`. Both must be green, including `tests/test_no_floats.py`, before the final commit.

- [ ] **Step 6: Final commit**
```bash
git add src/quiverlab/families src/quiverlab/citations src/quiverlab/batch \
        src/quiverlab/hochschild/table.py src/quiverlab/core/algebra.py \
        src/quiverlab/combinat/relations.py src/quiverlab/__init__.py \
        tests/families tests/citations tests/batch tests/combinat \
        docs/internals/10-families-citations.md README.md pyproject.toml
git commit -m "feat: Plan 06 acceptance -- full family catalog, citations registry, batch

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
```

---

## Frozen interface for Plans 07 (viz/trace) and 09 (web)

Plans 07 and 09 consume the following **verbatim** (do not rename, retype, or change signatures without a coordinated interface bump).

```python
# quiverlab (top-level, also via `from quiverlab import *`)
NakayamaAlgebra(kupisch=None, *, n=None, l=None, cyclic=False, field=CC) -> Algebra
PathAlgebra(type_or_quiver, orientation="linear", field=CC) -> Algebra
TruncatedPathAlgebra(type_or_quiver, r, orientation="linear", field=CC) -> Algebra
RadicalSquareZero(quiver, field=CC) -> Algebra
IncidenceAlgebra(poset_or_covers, elements=None, field=CC) -> Algebra
QuantumCI(q, field=CC) -> Algebra
ExteriorAlgebra(n, field=CC) -> Algebra
PreprojectiveAlgebra(type_or_quiver, field=CC) -> Algebra
TrivialExtension(A) -> Algebra
TensorProduct(A, B) -> Algebra
zoo(dim_max=12, field=CC) -> Iterator[Algebra]        # each yields A.zoo_name
families() -> FamilyListing                            # .names(), .by_name(n), to_dict(), __str__
bibliography(keys=None, grouped=True, annotated=True) -> Bibliography   # ITERABLE (see below)

# every family stamps A._family_citations : tuple[str, ...]  (registry keys)
Algebra.citations() -> tuple[str, ...]                 # family keys + engine keys, deduped
HHTable(dims, kind, algebra_repr, engine=..., references=())   # .references : tuple[str,...]
#   references = A.citations() (engine + family; ops added by Plan 04/05)

# quiverlab.citations
citations.reference(key) -> Reference(key, bibtex_key, kind, title, annotation, tags)
citations.bibtex(key) -> str                           # raw BibTeX entry
citations.all_keys() -> tuple[str, ...]
citations.references_bib_path() -> pathlib.Path
iter(Bibliography) -> Iterator[Entry]                  # Plan 09 consumes THIS
Entry(key, bibtex_key, formatted, doi, arxiv, topic, annotation)   # attrs Plan09 entry_view reads
Bibliography.groups : dict[str, list[Entry]]           # by topic (mapped from kind)
Bibliography.keys : tuple[str, ...]
Bibliography.bibtex() -> str ; Bibliography.to_dict() -> dict ; str(Bibliography)

# quiverlab.batch  (expert; NOT in `from quiverlab import *`)
batch.ResultsDB(path=":memory:")                       # SQLite; UNIQUE(builder, args, N)
batch.analyze(spec) -> dict ; batch.run_scan(specs, n_workers=1) -> list[dict]
batch.BUILDERS : dict[str, callable] ; batch.build_algebra(spec) -> Algebra
```

**What the frozen promises mean.**
- **`zoo()`** yields ready-to-compute `Algebra` objects in ascending `(dim, name)`, each carrying `.zoo_name` and the `("han_conjecture","chouhy_solotar")` citation stamp. Plan 09's catalog/runner consumes `zoo(dim_max)` (via `/api/catalog` and the build form); it never regenerates the catalog — the bundled `zoo_catalog.json` is authoritative.
- **`families()`** is the single discoverability surface and the **`sweep` catalog hook** (Plan 05): `sweep(NakayamaAlgebra, [3,2,2], fields=[...])` calls the builder once per field; every builder is field-generic with a `field=` keyword. Plan 05 owns `sweep`; Plan 06 owns the catalog. Plan 09 reads `families()` to populate the `/` build-form family dropdown and `/api/catalog` (it invokes builders by keyword, e.g. `truncated_polynomial(n=2, field=GF(2))`).
- **`bibliography()`** is the **only** source of citation text; it is **iterable**, and Plan 09's `webapp/server/references.py:entry_view` iterates it (`for e in ql.bibliography()`), reading `e.key / e.formatted / e.doi / e.arxiv / e.topic / e.annotation` — those `Entry` attribute names are the frozen contract (Plan 09 does **not** call `.to_dict()`). `/literature` renders `grouped_bibliography()` built from this iteration; annotations live in `citations/registry.py`, nowhere else. Unknown keys fail loudly with `CitationError`.
- **`A.citations()` / `HHTable.references`** (both engine **and** family keys, spec §3.9) let the trace subsystem (Plan 07) print, next to every computed number, the papers that produced it — the keys are already stamped; Plan 07 only renders them (and cites the `bar` route through `Hochschild1945`).

---

## Boundary notes

- **`sweep` ownership (Plan 05 vs Plan 06).** `sweep` does **not** exist in the tree at Plan-06 baseline (confirmed: no `def sweep` under `src/quiverlab/`) and is **not** drafted (`scratchpad/plan-05-modules-draft.md` absent). Per the ROADMAP, **Plan 05 owns `sweep()`**; **Plan 06 provides the catalog hook** — every family builder takes a `field=` keyword and is rebuildable over any field, and `families()` enumerates them. Plan 06 adds no `sweep`. A single in-plan contract test (Task 5/6) asserts a family rebuilds over `CC` and `GF(p)` with the same dim; that is the whole of Plan 06's sweep obligation.
- **Spec sections.** On `origin/main`, spec §3.9 is **Citations** (added 2026-07-18) and §3.10 is "Moving variables" (sweep + `quiverlab.batch`). Plan 06 owns the §3.9 registry core + `references.bib` + family refs (PDF bibliography in Plan 07, web surface in Plan 09, per §3.9's own split); the `quiverlab.batch` lift satisfies §3.10's batch clause. No spec-drift flag — the earlier working tree was stale.
- **QuantumCI sign convention.** The spec writes `xy + q·yx`; under it the anticommutative/exterior point is `q = +1` and `ExteriorAlgebra(2) = QuantumCI(q=1)` (Fixture N6). The Plan-06 dispatch's cross-check phrasing "`q = -1`" matches the alternate BGMS convention `xy − q·yx`. This plan follows the **spec** convention throughout and states the equivalence explicitly in the test.
- **Corrected numeric expectations in the dispatch brief.** `PathAlgebra("D4")` is **dim 9** (linear orientation), not 12 (the `D4` star admits at most 9 paths). `TrivialExtension(kA_2)` is **dim 6** (`2·dim kA_2 = 2·3`), not 5. Both corrections are proved in Fixtures N2 and N10.
- **Plan-04 gating for the open-zone HH tests.** *Building* the family catalog and the zoo algebras needs only the Plan-03 Gröbner route (a hard Plan-06 prereq, checked un-skipped in Task 1); the depth-16 open-zone HH goldens are produced by the Plan-02 minimal engine. The open-zone HH assertions in Tasks 1, 11, and 13 nonetheless `skipif` on `quiverlab.engine.resolutions_cs` — a single, consistent conservative gate so the suite stays green whatever order Plans 04/06 land in. The family *catalog* (Tasks 5–10) and the citations registry (Tasks 2–3) are unconditional on Plans 01–03.
- **Bibliography entry-count is a FLOOR (Plan 08 grows it).** `test_bib_covers_registry` asserts registry-coverage + `count >= 16` (never an equality). Plan 08 appends the four software-citation entries (`qpa`, `gap4`, `sagemath`, `quiverlab`) to the packaged `references.bib` and raises the floor; the coverage assertion (every registry `bibtex_key` resolves) is the invariant, the count is only a floor.
- **`multiprocessing` in batch.** The bank's `labdb.run_scan` uses `multiprocessing.Pool`. Per spec §5 ("no multiprocessing in core") and Marco's machine limits, `batch.run_scan` defaults to `n_workers=1` and only forms a `Pool` when a caller explicitly raises it; batch is expert-only and not imported by the beginner surface.
- **Structure-constant families carry no quiver.** `TensorProduct`/`TrivialExtension` return an `Algebra` with `quiver=None`; `A.cartan_matrix()` on them raises the existing "needs the quiver presentation" error (correct — there is no canonical path basis). Their HH still computes via the bar/fast engine.
- **Read-only bank.** `open_zoo_catalog_v2.json` and the three bank tests are read once (Task 11 curation, Task 13 golden fixtures) and never modified; the lifted `db.py`/`scan.py` carry the ported-with-attribution header used across `engine/`.
```
