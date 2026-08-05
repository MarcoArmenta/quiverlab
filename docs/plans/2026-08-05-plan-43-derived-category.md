# Plan 43: Derived-Category Surface — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The derived-category surface on top of P39's complex layer: reified
hyper-Hom classes (`Hom_{D^b}(X, Y[n])` as actual chain maps), the derived AR
translate `τ_{D^b} = ν∘[−1]` on perfect complexes (**loud refusal at infinite
global dimension**, per Happel), a **tilting-complex verifier** (rigidity via
hyper-Hom + K₀ unimodular generation) with `End(T)` recovered as a
structure-constant algebra (Rickard — the derived-equivalent algebra), and a
**derived fingerprint** panel comparing two algebras on Coxeter polynomial,
Cartan (det + Smith), HH/HC and center — in **necessary-condition language only**.
This is the metaplan's P43 card and the "derived categories" research axis.

**Architecture:** One new package `src/quiverlab/derived/` — `homs.py`
(hyper-Hom basis reification), `tau.py` (`τ_{D^b}`/`τ⁻_{D^b}`), `tilting.py`
(tilting/silting verifier + `End(T)`), `fingerprint.py` (the invariant tuple +
comparison), and the shared `block.py::derived_fingerprint_block` the GUI/webapp/
HPC runners call. Every routine is a thin exact-linear-algebra layer over
primitives that already exist: P39's `ChainComplex`/`ChainMap`/`hyper_hom_dims`
and its **private** `_hom_total_blocks`/`_delta_total` (reified here, never
recomputed — the basis-mismatch discipline), P37's `end_algebra`/`ModuleHom.then`
structure-constant machinery, P38's `cartan_matrix`/`coxeter_matrix`/
`coxeter_polynomial`/`euler_form_matrix`, P41's `nakayama_functor` + the
corner-transpose helper `duality.transpose_module` uses, and `modules/ext.py::
global_dimension` for the Happel gate. **No new math engine.** Every constructed
object self-certifies (chain-map squares commute, `τ_{D^b}` output pinned against
the trusted module `τ`, `End(T)` corner-Cartan pinned against the known Cartan)
and is cross-checked against QPA (`TauOfComplex`) and closed-form Dynkin ground
truth — correctness never rests on trusting a construction.

**Tech Stack:** `modules/linalg_mod` (`matmul`, `mat_rank`, `kernel_columns`,
`cols_to_matrix`, `solve_columns`, `independent_modulo`, `matvec`, `transpose`,
`col`, `identity`, `zeros`), `fields.linalg.reduce_mod_nullspace` (byte-stable
canonicalisation — the CS/Plan-17 precedent), `ModuleHom`/`end_algebra`/
`hom_basis`/`direct_sum` (**Plan 37**), `ChainComplex`/`ChainMap`/`hyper_hom_dims`/
`projective_model`/`from_projective_resolution` (**Plan 39**), `nakayama_functor`/
the corner-transpose (**Plan 41**), `invariants/cartan`+`invariants/forms`
(**Plan 38**), `sympy` (`smith_normal_form`/`invariant_factors`, `Poly`,
`charpoly`). No floats in `src/`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **Prerequisites: P39 AND P41 merged to `dev`** (P39 for the complex/hyper-Hom
  substrate; P41 for `nakayama_functor`/`nakayama_functor_minus` module `τ`
  cross-checks and the factored corner-transpose helper). P37 and P38 are already
  merged (this is W3; P37/P38 landed in W1). Branch `plan-43-derived-category`
  off `dev` **after P39 and P41 have merged**. If P41's `_corner_transpose_matrix`
  helper is present, **reuse it**; if P41 factored it differently or not at all,
  factor it out of `duality.transpose_module` (`duality.py:79-99`, the `d1star`
  construction) into one shared home and import it from both places — the task
  below detects and adapts.
- `tests/modules/` auto-assigns to the **deep** bucket, `tests/qpa/` to the
  **qpa** bucket, `tests/webapp/`+`tests/hpc/` are extras-gated (unmarked)
  (`tests/conftest.py`). New library tests live in
  `tests/modules/test_derived_*.py` — **do NOT create a new top-level tests
  dir**. Run new tests by path during development; finish each task with a
  `-m deep` (or `-m qpa`) spot-run of the touched files.
- **`decompose` refuses over GF(p) when char ≤ dim M** (the trace-form radical is
  unreliable there — `modules/decompose.py`). Any battery that decomposes a
  perfect complex's terms into indecomposable projectives, or certifies module
  indecomposability, or builds `End(T)` runs over **QQ or GF(32003)** so both the
  split search and the locality certificate decide. State this per battery.
- House conventions: homological complexes (`d_n: C_n → C_{n-1}`, rows=target —
  P39's module docstring); a module-map matrix is dense `tgt.dim × src.dim` over
  the shared Domain; composition is left-to-right, `f.then(g)` (P37 `ModuleHom`,
  extended to `ChainMap` here). Hom/comparisons refuse loudly across sides/
  algebras. All refusals are `QuiverlabError`; `check=False` fast paths only for
  internally-constructed data whose certificate is asserted separately.
- **Honest scope (metaplan §6, stated on the verification page as each piece
  lands):** *deciding* derived equivalence is not algorithmic — P43 ships
  **verifiers** and **necessary-condition** invariants, never a decider; the
  derived fingerprint's language is "distinguished / not distinguished by these
  invariants", **never** "(in)equivalent". Serre functor / AR triangles / `τ_{D^b}`
  refuse **loudly** at infinite global dimension (`k[x]/(x²)` = the pinned
  negative test). Classifying `D^b` indecomposables for wild algebras is out of
  scope.
- **The K₀-action Coxeter matrix is the CONJUGATE of P38's `coxeter_matrix`**
  (see Task 2): P38's `coxeter_matrix(A) = −C⁻ᵀC` has the correct **Coxeter
  polynomial** (char poly) but does **not** act on dimension vectors as `[τ]`; the
  dim-vector action is `c = −C·C⁻ᵀ` (metaplan §2's `−CᵀC⁻¹` in the ASS Cartan
  convention). Both are conjugate, same char poly; the fingerprint uses the
  polynomial (P38, unchanged), the `τ_{D^b}` K₀ oracle uses `c`, arbitrated by the
  concrete kA₂/kA₃ `τ`.
- Plan-32 markers: chain-map-square / certificate / round-trip / `d²=0` tests =
  `oracle_selfcert`; two-independent-route dim agreement (hyper-Hom ≡ module Ext,
  complex-`τ` ≡ module-`τ`) = `oracle_crossengine`; Dynkin closed-form values
  (mesh dim vectors, tilting Cartan, D₄/A₄ Coxeter distinction) =
  `oracle_literature`; QPA comparisons live in `tests/qpa/` (bucket = the class,
  never double-marked). Every plan merge updates `docs/verification.md` (new
  oracle rows + recounted class counts, `tests/release/test_oracle_classes.py`
  green) and adds citations to `citations/references.bib` + `registry.py`.
- **Merge-train count note:** the metaplan runs Plans 36–50 in parallel waves;
  the suite/oracle-class counts on the verification page drift as siblings merge.
  Recount from a live collection at *this* plan's merge (do not copy a stale
  number); note "counts as of the P43 merge; sibling plans in flight may shift
  them" beside the table, matching the P39/P41 recount discipline.
- Conventional commits; green tests at every commit.

---

### Task 1: `ChainMap.then` (P39 extension) + `hyper_hom_basis` reification

The net-new primitive P43 is built on: P39 computes hyper-Hom **dimensions**
(`hyper_hom_dims`) but exposes no way to get the classes as objects. Task 1
reifies a basis of `H^n(Hom^•(X, Y))` as actual `ChainMap`s `X → Y[n]` (a cocycle
in `Hom^n` **is** a degree-0 chain map into the shift, self-certified by
`ChainMap(..., check=True)`), and adds `ChainMap.then` so those maps compose
(needed for `End(T)` in Task 3).

**Files:**
- Modify: `src/quiverlab/modules/complexes.py` (add `ChainMap.then`; stash
  per-term projective-summand vertices in `from_projective_resolution` for Task 2)
- Create: `src/quiverlab/derived/__init__.py`, `src/quiverlab/derived/homs.py`
- Test: `tests/modules/test_derived_homs.py`

**Interfaces:**
- Consumes (P39, read first): `ChainComplex.{degrees,term,shift,is_perfect}`,
  `ChainMap.{component,cone}`, `hyper_hom_dims(X, Y, lo, hi)` (signature is
  `(X, Y, lo, hi)` — a range, not a single degree), and the **private**
  `_hom_total_blocks(X, Y, n, dom) -> (blocks, dim)` /
  `_delta_total(X, Y, n, dom) -> (matrix, src_dim, tgt_dim)` /
  `_combine_homs(H, coeffs, sdim, tdim, dom)`. `_hom_total_blocks`' block dicts are
  `{"p", "q"=p-n, "homs", "offset", "count"}` with `homs` the `hom_space(X_p, Y_q)`
  basis (`Y_q.dim × X_p.dim` matrices).
- Produces:
  ```python
  # modules/complexes.py -- ChainMap method (mirrors ModuleHom.then, morphism.py:67)
  def then(self, g: "ChainMap") -> "ChainMap":
      # self then g (left-to-right): (self.then(g)).component(n) = g_n @ self_n.
      # Requires g.src is self.tgt (same middle complex). check=False (each square
      # is the matmul of two commuting squares; the callers re-validate on demand).

  # derived/homs.py
  def hyper_hom_basis(X: ChainComplex, Y: ChainComplex, n: int) -> list[ChainMap]:
      # A basis of H^n(Hom^.(X, Y)) reified as chain maps X -> Y.shift(n). For X
      # PERFECT this is Hom_{D^b}(X, Y[n]) (P39 header). Each returned map is built
      # from a canonical coset representative of ker(delta^n)/im(delta^{n-1}) and is
      # a genuine chain map (ChainMap(..., check=True) is the reification self-cert).
      # SELF-CERT: len(result) == hyper_hom_dims(X, Y, n, n)[n] (asserted here).
      # QuiverlabError if X is not certified perfect (use projective_model first).
  ```

**Why a cocycle is a chain map `X → Y[n]` (state verbatim in the docstring).**
`Hom^n(X, Y) = ⊕_p Hom(X_p, Y_{p-n})`; a component `f_p: X_p → Y_{p-n}` lands in
`Y.shift(n).term(p) = Y.term(p-n)`. P39's cocycle condition (convention `(*)`,
Weibel 2.7.4) is `d^Y f_p = (-1)^n f_{p-1} d^X`, i.e.
`(-1)^n d^Y f_p = f_{p-1} d^X`; the shifted complex `Y[n]` carries differential
`(-1)^n d^Y`, so this is exactly the chain-map square for `X → Y[n]`. Hence
`ChainMap(X, Y.shift(n), comps, check=True)` **passes** on any cocycle and
**fails** on a non-cocycle — the construction is self-certifying.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_derived_homs.py
"""hyper_hom_basis reifies H^n(Hom^.(X, Y)) as chain maps X -> Y[n]. Self-cert:
every reified class is a valid chain map (check=True), the count equals
hyper_hom_dims, and degree-0 composition matches End on stalks. Cross-engine:
degree-n classes count Ext^n on a projective-resolution source."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import ChainComplex, identity_chain_map
from quiverlab.modules.complexes import hyper_hom_dims
from quiverlab.derived.homs import hyper_hom_basis

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


@xeng
def test_basis_count_equals_ext_on_resolution():
    A = _a3()
    for v in (1, 2, 3):
        M, N = A.simple(v), A.simple(1)
        X = ChainComplex.from_projective_resolution(M, length=5)
        Y = ChainComplex.stalk(N, 0)
        for n in range(0, 5):
            basis = hyper_hom_basis(X, Y, n)
            assert len(basis) == A.ext(M, N, n)                 # count = Ext^n
            assert len(basis) == hyper_hom_dims(X, Y, n, n)[n]  # count = dim H^n


@selfcert
def test_every_reified_class_is_a_chain_map():
    # ChainMap(check=True) inside hyper_hom_basis is the self-cert; re-assert the
    # target is exactly Y.shift(n) and each square commutes (rebuild with check).
    A = _a3()
    M, N = A.simple(2), A.simple(1)
    X = ChainComplex.from_projective_resolution(M, length=4)
    Y = ChainComplex.stalk(N, 0)
    for n in range(0, 4):
        for f in hyper_hom_basis(X, Y, n):
            assert f.tgt.degrees() == Y.shift(n).degrees()
            from quiverlab.modules.complexes import ChainMap
            ChainMap(f.src, f.tgt, {k: f.component(k) for k in f.src.degrees()},
                     check=True)                                # squares commute


@selfcert
def test_degree0_endo_composition_and_identity():
    # degree-0 hyper-Hom of a projective stalk with itself = End_A(P); ChainMap.then
    # composes them; the identity chain map is one of the classes (up to homotopy).
    A = _a3()
    P1 = A.projective(1)
    X = ChainComplex.stalk(P1, 0)
    basis = hyper_hom_basis(X, X, 0)
    assert len(basis) == A.hom(P1, P1)
    idX = identity_chain_map(X)
    # then() type-checks and stays a chain map:
    comp = basis[0].then(idX)
    assert comp.component(0) == basis[0].component(0)           # f . id = f


@selfcert
def test_nonperfect_source_refused():
    A = _a3()
    X = ChainComplex.stalk(A.simple(2), 0)      # simple: not projective
    with pytest.raises(QuiverlabError, match="perfect"):
        hyper_hom_basis(X, X, 0)
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_derived_homs.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.derived` / `AttributeError: then`.

- [ ] **Step 3: Implement**

`ChainMap.then` in `modules/complexes.py` (place it beside `cone`):

```python
def then(self, g):
    """``self`` then ``g`` (left-to-right, mirroring ``ModuleHom.then``):
    ``(self.then(g)).component(n) = g.component(n) @ self.component(n)``. Refuses
    if the middle complexes differ (``g.src is self.tgt``)."""
    if g.src is not self.tgt:
        raise QuiverlabError("ChainMap.then: middle complexes differ")
    dom = self.domain
    comps = {}
    degs = set(self.src.degrees()) | set(g.tgt.degrees())
    for n in degs:
        fn, gn = self.component(n), g.component(n)
        if self.tgt.term(n).dim == 0:                 # matmul is shapeless through 0
            comps[n] = lm.zeros(g.tgt.term(n).dim, self.src.term(n).dim, dom)
        else:
            comps[n] = lm.matmul(gn, fn, dom)
    return ChainMap(self.src, g.tgt, comps, check=False)
```

Stash the projective-summand vertices in `from_projective_resolution` (a one-line
provenance addition Task 2 reads — the resolution already knows them):

```python
# in ChainComplex.from_projective_resolution, after building `terms`/`dmats`:
X._proj_vertices = {n: list(terms_list[n].vertices)     # degree -> [summand vertices]
                    for n in terms if terms_list[n].module is not None}
```

`derived/homs.py`:

```python
"""Reified hyper-Hom: a basis of H^n(Hom^.(X, Y)) as chain maps X -> Y[n]
(Plan 43 / derived category). Thin accessor over P39's private Hom-total-complex
internals (_hom_total_blocks / _delta_total / _combine_homs) -- never recompute the
Hom complex a second way (the P41 basis-mismatch discipline). Canonical coset
representatives via fields.linalg.reduce_mod_nullspace (byte-reproducible, the CS
precedent)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.complexes import (ChainComplex, ChainMap, hyper_hom_dims,
                                         _hom_total_blocks, _delta_total,
                                         _combine_homs)


def _coset_reps(cocycle_cols, coboundary_cols, dom):
    """Canonical representatives of ker/im: choose a maximal subset of `cocycle_cols`
    independent modulo `coboundary_cols` (dim = dim H^n), then reduce each modulo the
    coboundary span for byte-stability (reduce_mod_nullspace against the coboundary
    matrix's row-solve -- the Plan-17 free-variables-zero coset canonicalisation)."""
    idx = lm.independent_modulo(cocycle_cols, coboundary_cols, dom)
    reps = [list(cocycle_cols[i]) for i in idx]
    if coboundary_cols:
        B = lm.cols_to_matrix(coboundary_cols)
        reps = [linalg.reduce_mod_nullspace(r, lm.transpose(B), dom) for r in reps]
    return reps


def hyper_hom_basis(X, Y, n):
    if not X.is_perfect():
        raise QuiverlabError(
            "hyper_hom_basis: X must be a perfect complex; resolve it with "
            "projective_model(X, ...) first (P39).")
    dom = X.domain
    blocks, cdim = _hom_total_blocks(X, Y, n, dom)
    dn, _s, _t = _delta_total(X, Y, n, dom)             # delta^n : Hom^n -> Hom^{n+1}
    dn1, _s1, _t1 = _delta_total(X, Y, n - 1, dom)      # delta^{n-1}: Hom^{n-1}->Hom^n
    if dn and dn[0]:
        cocycles = lm.kernel_columns(dn, dom)
    else:                                               # delta^n = 0 => every cochain
        ident = lm.identity(cdim, dom)
        cocycles = [lm.col(ident, j) for j in range(cdim)]
    cobounds = ([lm.col(dn1, j) for j in range(len(dn1[0]))]
                if (dn1 and dn1[0]) else [])
    reps = _coset_reps(cocycles, cobounds, dom)
    Yn = Y.shift(n)
    maps = []
    for rep in reps:
        comps = {}
        for b in blocks:
            coeffs = rep[b["offset"]: b["offset"] + b["count"]]
            comps[b["p"]] = _combine_homs(b["homs"], coeffs,
                                          X.term(b["p"]).dim, Y.term(b["q"]).dim, dom)
        maps.append(ChainMap(X, Yn, comps, check=True))  # cocycle => valid chain map
    # self-cert: the homotopy-quotient dimension equals P39's rank-formula dimension.
    if len(maps) != hyper_hom_dims(X, Y, n, n)[n]:
        raise QuiverlabError(
            "hyper_hom_basis: reified class count != hyper_hom_dims (coset/rank "
            "mismatch -- basis reification bug)")
    return maps
```

**Adjust to reality:** confirm `_combine_homs`' `(sdim, tdim)` order against
`complexes.py:646` (it returns a `tdim × sdim` = target×source matrix; here
`sdim = X_p.dim`, `tdim = Y_q.dim`). If `reduce_mod_nullspace` wants the *nullspace
generators*, feed `transpose(B)` so its nullspace is the coboundary span, or use
whatever coset-canonicaliser P39/CS already exposes — the byte-stability is what
matters; the **count** self-cert is the real guarantee. `hyper_hom_dims(X,Y,n,n)`
re-runs the two `_delta_total`s; if that double-compute shows in profiling, factor
a private `_delta_at(X,Y,m,dom,cache)`  — but never a second *independent* route.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/complexes.py src/quiverlab/derived/__init__.py \
        src/quiverlab/derived/homs.py tests/modules/test_derived_homs.py
git commit -m "feat(derived): ChainMap.then + hyper_hom_basis -- reified Hom_{D^b}(X, Y[n]) classes"
```

---

### Task 2: `tau.py` — the derived AR translate `τ_{D^b}` / `τ⁻_{D^b}`

**Files:**
- Create: `src/quiverlab/derived/tau.py`
- Modify: `src/quiverlab/core/algebra.py` (optional `Algebra.tau_complex(M, ...)`
  convenience — build the resolution, apply `τ_{D^b}`; read the delegation idiom)
- Test: `tests/modules/test_derived_tau.py`

**Math pinned (Happel).** For `A` of **finite global dimension**, `D^b(mod A)` has
a Serre functor `S = ν` (the derived Nakayama functor) and AR triangles, with
`τ_{D^b} = ν ∘ [−1]`; at infinite global dimension **there is no Serre functor** —
refuse loudly (`k[x]/(x²)` = the pinned negative case). On a **perfect** complex
`X` (bounded, projective terms) `ν` is applied **termwise**: each projective term
`P_v ↦ I_v = ν P_v` (`builders.injective`), each differential (a map between
projective terms) `↦` its `ν`-image (a map between injectives) via the
corner-transpose. `τ_{D^b}(X) = (ν X)[−1]`.

**Interfaces:**
- Consumes: `modules/ext.py::global_dimension` (`.exact` ⇔ finite here — `exact`
  is True only when every simple resolves within the bound); `builders.injective`
  (`I_v`, `builders.py:74`); `builders.projective`; **P41's corner-transpose
  helper** (reuse `modules/ar.py::_corner_transpose_matrix` if P41 factored it,
  else factor it out of `duality.transpose_module`'s `d1star` construction,
  `duality.py:79-99`, into a shared home); `duality.tau`/P41
  `nakayama_functor` (the trusted-`τ` arbiter); `ChainComplex.{shift,homology,
  homology_dims}`; `decompose`+`identify_standard` (the provenance fallback);
  `invariants/cartan::cartan_matrix` (the K₀ oracle's `c`).
- Produces:
  ```python
  def tau_Db(X: ChainComplex) -> ChainComplex        # (nu X)[-1]; X certified perfect
  def tau_Db_minus(X: ChainComplex) -> ChainComplex  # (nu^-1 X)[+1] via the A^op dual
  ```

**Per-term summand vertices.** `ν` needs each perfect term's projective-summand
vertex multiset. Fast path: `X._proj_vertices` (Task-1 provenance, set by
`from_projective_resolution` and carried through `shift`/`cone`). Fallback: recover
per term via `decompose(term)` + `identify_standard(summand) == ("projective", v)`
(loud if any summand is non-projective — `X` was not perfect; **char scope**
QQ/GF(32003) for `decompose`). The differential's block layout must match the
vertex order — assert it (a `ν(d)` that lands off the injective grading raises).

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_derived_tau.py
"""Derived AR translate tau_{D^b} = nu[-1] on perfect complexes. Self-cert:
d.d=0 in the output (ChainComplex check), and the round-trip tau^-_Db . tau_Db is
a quasi-iso to X. Cross-engine: on a projective resolution of a non-projective
indecomposable M over kA_n, homology of tau_Db is concentrated in degree 0 and
isomorphic to the trusted module tau(M). Literature: the K0 identity
chi(tau_Db X) = c . chi(X) with c the K0-action Coxeter matrix. Negative:
k[x]/(x^2) (infinite gl.dim) refuses."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.tau import tau_Db, tau_Db_minus

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


def _a3():
    return linear_path_algebra(3, field=QQ)     # kA3, 1->2->3, hereditary


@xeng
@pytest.mark.parametrize("n, v", [(2, 1), (2, 2), (3, 2)])
def test_tau_db_of_nonprojective_is_module_tau(n, v):
    A = linear_path_algebra(n, field=QQ)
    M = A.simple(v)
    if M.tau().dim == 0:                          # projective: skip (no module tau)
        pytest.skip("projective simple has no module tau")
    X = ChainComplex.from_projective_resolution(M, length=6)
    T = tau_Db(X)
    hd = T.homology_dims()
    # concentrated in degree 0 (nu M = 0 for non-projective interval modules over kA_n)
    assert all(d == 0 for k, d in hd.items() if k != 0)
    assert hd.get(0, 0) == M.tau().dim
    assert T.homology(0).dimension_vector() == M.tau().dimension_vector()


@lit
def test_k0_coxeter_bookkeeping():
    # chi(tau_Db X) = c . chi(X), c = -C * C^-T the K0-ACTION Coxeter matrix (NOT
    # P38's coxeter_matrix, which is the conjugate -C^-T C -- same char poly).
    import sympy as sp
    A = _a3()
    C = sp.Matrix(A.cartan_matrix())
    c = -C * C.inv().T
    verts = list(A.quiver.vertices)
    for v in (1, 2):                              # non-projective simples over kA3
        M = A.simple(v)
        X = ChainComplex.from_projective_resolution(M, length=6)
        T = tau_Db(X)
        chiX = _chi_vec(X, verts)
        chiT = _chi_vec(T, verts)
        assert list(c * sp.Matrix(chiX)) == [sp.Integer(x) for x in chiT]


@selfcert
def test_round_trip_is_quasi_iso():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(2), length=6)
    back = tau_Db_minus(tau_Db(X))
    # homology dimension vectors agree degreewise (a quasi-iso to X in D^b)
    assert back.homology_dims() == X.homology_dims()


@selfcert
def test_infinite_gldim_refused():
    from quiverlab.families import truncated_polynomial
    A = truncated_polynomial(2, field=GF(7))     # k[x]/(x^2): gl.dim = infinity
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    with pytest.raises(QuiverlabError, match="global dimension|Serre|gl.dim"):
        tau_Db(X)


def _chi_vec(Z, verts):
    out = [0] * len(verts)
    for k in Z.degrees():
        dv = Z.term(k).dimension_vector()
        for i, w in enumerate(verts):
            out[i] += (-1) ** k * dv.get(w, 0)
    return out
```

**Sharpening note (kA_n concentration).** `ν M = 0` for a non-projective interval
module over `kA_n` (verified by hand — kA₂: `ν S_1 = coker([1,2,3]→S_1) = 0`;
kA₃: `ν[1,2] = ν S_2 = 0`), so `tau_Db` homology is concentrated in degree 0.
State this scope in the test; on an algebra where `ν M ≠ 0` the homology is *not*
concentrated (that is correct, not a bug) — the concentration oracle is a **kA_n**
fact, the K₀ identity is general.

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
"""tau_{D^b} = nu[-1] and tau^-_{D^b} on perfect complexes (Plan 43 / Happel). The
Serre functor of D^b(mod A) exists iff gl.dim A < infinity; refuse loudly otherwise
(k[x]/(x^2) is the negative test). nu is applied termwise on a perfect (projective)
complex -- P_v -> I_v and the corner-transpose on the maps -- then shifted. The
output self-certifies (d.d=0 via ChainComplex check=True) and is pinned against the
trusted module tau; correctness never rests on the corner-transpose bookkeeping."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.builders import injective, projective
from quiverlab.modules.complexes import ChainComplex


def _require_finite_gldim(A):
    from quiverlab.modules.ext import global_dimension
    g = global_dimension(A)
    if not g.exact:                              # exact == resolved == finite here
        raise QuiverlabError(
            f"tau_Db: D^b(mod A) has no Serre functor -- global dimension is not "
            f"finite ({g!r}); the derived AR translate is undefined (Happel).",
            hint="k[x]/(x^2) is the pinned negative case; needs finite gl.dim")
    return g.value


def _term_vertices(X, n):
    """[summand vertex list] of the perfect term X_n: provenance fast path, else
    decompose + identify_standard (loud on a non-projective summand / char scope)."""
    prov = getattr(X, "_proj_vertices", None)
    if prov is not None and n in prov:
        return list(prov[n])
    from quiverlab.modules.decompose import decompose
    from quiverlab.modules.hom import identify_standard
    verts = []
    for s, _mult in decompose(X.term(n)):
        std = identify_standard(s)
        if not std or std[0] != "projective":
            raise QuiverlabError(
                f"tau_Db: degree-{n} term is not projective -- X is not perfect "
                f"(got summand {std})")
        verts.append(std[1])
    return verts


def _nu_map(d, src_verts, tgt_verts, A):
    """nu of a map d: (+)P_{src} -> (+)P_{tgt} between projective terms, as a matrix
    (+)I_{src} -> (+)I_{tgt} between the injective bases. d^* = corner-transpose of d
    (reuse P41's _corner_transpose_matrix if present, else factor it from
    duality.transpose_module's d1star); nu(d) = D(d^*) = its k-dual (transpose in the
    injective basis, exactly builders.injective's `transpose(Lb)` convention)."""
    try:
        from quiverlab.modules.ar import _corner_transpose_matrix       # P41 (reuse)
    except ImportError:                                                  # factor here
        from quiverlab.derived._corner import _corner_transpose_matrix
    dstar = _corner_transpose_matrix(d, tgt_verts, src_verts, A)        # (+)Ae_tgt -> (+)Ae_src
    return lm.transpose(dstar)                                          # D: right action = transpose


def _apply_nu_perfect(X):
    A, dom = X.algebra, X.domain
    terms, dmats = {}, {}
    vlists = {}
    for n in X.degrees():
        if X.term(n).dim == 0:
            continue
        vlists[n] = _term_vertices(X, n)
        blocks = [injective(A, v) for v in vlists[n]]
        terms[n] = _direct_sum_injectives(blocks, dom, name=f"nu_C_{n}")
    for n in X.degrees():
        d = X._dmats.get(n) if hasattr(X, "_dmats") else None
        if d and d[0] and n in terms and (n - 1) in terms:
            dmats[n] = _nu_map(d, vlists[n], vlists[n - 1], A)
    return ChainComplex(terms, dmats, check=True)          # d.d=0 self-cert


def tau_Db(X):
    _require_finite_gldim(X.algebra)
    return _apply_nu_perfect(X).shift(-1)


def tau_Db_minus(X):
    # tau^-_Db over A = D . tau_Db over A^op . D  (the P41 nu^- contravariant-dual
    # trick lifted to complexes); certified by the round-trip quasi-iso to X.
    _require_finite_gldim(X.algebra)
    return _op_dual(tau_Db(_op_dual(X))).shift(2)          # see note on the net shift
```

**Adjust to reality — the load-bearing bookkeeping (arbitrated by the oracle,
not asserted a priori):**
- `_corner_transpose_matrix(d, from_verts, to_verts, A)` returns the corner-transpose
  `d^*` of a map between projective terms. **The exact argument order, the row/col
  layout, and the final `transpose` are pinned by
  `test_tau_db_of_nonprojective_is_module_tau`** — if `homology(0)` disagrees with
  `M.tau()`, the corner-transpose direction or the `D`-transpose is flipped: fix it
  **once** against the kA₂/kA₃ hand computation and pin the passing identity with a
  derivation comment (the P39/P41 convention-flip discipline). Do **not** special-case
  per degree.
- `_direct_sum_injectives` mirrors `resolution._direct_sum` (the block-diagonal
  module builder duality already uses); the injective basis order must match
  `builders.injective` (paths **ending** at `v`).
- `tau_Db_minus`'s net shift: `[-1]` from `tau_Db` on `A^op`, plus `D`'s
  cohomological degree flip (`D` sends homological `n ↦ −n`) — the composite lands a
  `shift` off `X`; **derive the exact integer from the round-trip test**
  (`tau_Db_minus(tau_Db(X))` must be quasi-iso to `X`) and pin it. If a cleaner
  self-contained route (an injective co-model, dual to P39's `projective_model`) is
  wanted later, it is a strict superset — but the `A^op` dual reuses `tau_Db` and is
  enough for the round-trip + `ν⁻(I_v) ≅ P_v` self-cert. Keep the two **numeric**
  oracles (concentration ≅ module-`τ`, K₀ identity) on `tau_Db`.
- If `X._proj_vertices` is absent and `decompose` refuses (char ≤ dim), raise the
  `decompose` error unchanged (honest) — the batteries run over QQ/GF(32003).

- [ ] **Step 4: Run tests** — Expected: PASS (corner-transpose flip resolved)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/derived/tau.py src/quiverlab/core/algebra.py \
        tests/modules/test_derived_tau.py
git commit -m "feat(derived): tau_{D^b} = nu[-1] on perfect complexes -- Happel gate, module-tau + K0 oracles"
```

---

### Task 3: `tilting.py` — tilting/silting verifier + `End(T)` as an algebra

**Files:**
- Create: `src/quiverlab/derived/tilting.py`
- Test: `tests/modules/test_derived_tilting.py`

**Interfaces:**
- Consumes: `hyper_hom_dims` (P39 — rigidity), `hyper_hom_basis`+`ChainMap.then`
  (Task 1 — `End(T)`), `end_algebra`/`_structure_constants`/`regular_corner_dims`
  (P37 `modules/endomorphism.py` — the corner oracle template),
  `Algebra.from_structure_constants`, `fields.linalg.solve`+`reduce_mod_nullspace`,
  `cartan_matrix` (P38), `minimal_resolution` (the 2-term presentation).
- Produces:
  ```python
  @dataclass
  class TiltingReport:
      is_tilting: bool
      rigid: bool                 # hyper_hom_dims(T, T)[n] == 0 for all n != 0 in window
      generates: bool             # K0 g-matrix square & unimodular (det +-1)
      window: tuple[int, int]     # (n_min, n_max): the EXACT rigidity-check range
      g_matrix: list[list[int]]   # rows = summand K0 classes chi (Euler char per vertex)
      det: int

  def is_tilting_complex(summands: list[ChainComplex]) -> TiltingReport
      # Every summand must be certified perfect (loud otherwise). rigidity is DECIDED
      # (perfect => bounded => the window is finite and exact); generation is DECIDED
      # (integer determinant). No semi-decision here.

  def end_algebra_of_complex(summands: list[ChainComplex]) -> Algebra
      # End_{D^b}(T) = (+)_{i,j} Hom_{D^b}(T_i, T_j) with composition = ChainMap.then
      # reduced to canonical homotopy reps; structure constants -> from_structure_constants.

  def two_term_silting_from_presentation(M) -> tuple[ChainComplex, TiltingReport]
      # the 2-term complex [P_1 --d_1--> P_0] (degrees 1,0) of M's minimal projective
      # presentation + the is_tilting-style rigidity report (the AIR bridge P45 consumes).
  ```

**The rigidity window (decided, not semi-decided).** For perfect summands each
spanning degrees `[lo_i, hi_i]`, `Hom^n(T_i, T_j)` is the **zero** cochain group
unless some degree `p` has both `T_i` at `p` and `T_j` at `p−n` nonzero, i.e.
`n ∈ [lo_i − hi_j, hi_i − lo_j]`. Over all pairs the honest window is
`n_min = min_i lo_i − max_j hi_j`, `n_max = max_i hi_i − min_j lo_j`; **outside it
hyper-Hom is provably 0**, so rigidity is fully decided by scanning `n ∈
[n_min, n_max] \ {0}`. Report the window (metaplan honesty: name what was checked).

**K₀ generation.** Each summand `T_i`'s class in `K₀ = ℤ^{#vertices}` is its Euler
characteristic `χ(T_i) = Σ_n (−1)^n dim-vec(T_i,n)`. Stack these as the rows of the
`#summands × #vertices` **g-matrix**. `T` generates `K^b(proj)` (Rickard) iff this
matrix is **square** (`#summands = #simples`) and **unimodular** (`det = ±1`) — the
classes are then a `ℤ`-basis of `K₀`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_derived_tilting.py
"""Tilting-complex verifier + End(T). Self-cert: T=A (projective stalks) is tilting
with End(T) ~ A (corner-Cartan). Literature: the kA2 APR tilt T = P1 (+) S1 is
tilting, End(T) has the reoriented-A2 Cartan. Negative: a missing summand fails
generation; X (+) X[1] fails rigidity."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.tilting import (is_tilting_complex, end_algebra_of_complex,
                                       two_term_silting_from_presentation)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a2():
    return linear_path_algebra(2, field=QQ)      # kA2, 1->2, hereditary


@selfcert
def test_regular_module_is_tilting_and_end_is_A():
    # T = A = (+)_v stalk(P_v): the trivial tilting complex; End(T) ~ A.
    A = _a2()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    rep = is_tilting_complex(T)
    assert rep.is_tilting and rep.rigid and rep.generates
    assert rep.window == (0, 0)                  # width-0 summands: Ext^{!=0}(P,P)=0
    E = end_algebra_of_complex(T)
    from quiverlab.derived.tilting import corner_cartan_of_complex
    assert corner_cartan_of_complex(T) == A.cartan_matrix()   # End(A_A) ~ A oracle


@lit
def test_ka2_apr_tilt():
    # DERIVATION (kA2, arrow a:1->2; right modules). Indecomposables: S1=(1,0),
    # S2=P2=(0,1), P1=[1,2]=(1,1). Non-projective: S1 (pd 1: 0->P2->P1->S1->0).
    # APR tilt at the sink 2: T = P1 (+) tau^{-1}(P2) = P1 (+) S1 (the unique AR
    # sequence 0->S2->P1->S1->0 gives tau^{-1}(S2)=S1). Checks: pd P1=0, pd S1=1;
    # Ext^1(S1,P1)=0 (coker(Hom(P1,P1)->Hom(P2,P1)) = coker(k->k iso)=0),
    # Ext^1(S1,S1)=0 (Hom(P2,S1)=0); Ext^1(P1,-)=0. Summands=2=#simples. So T is
    # tilting. End(T): Hom(P1,P1)=Hom(S1,S1)=k, Hom(P1,S1)=k (P1->>S1), Hom(S1,P1)=0
    # => Cartan(End T) = [[1,1],[0,1]] (hereditary A2, reoriented) = A's Cartan.
    A = _a2()
    P1 = A.projective(1)
    S1 = A.simple(1)
    T = [ChainComplex.stalk(P1, 0),
         ChainComplex.from_projective_resolution(S1, length=2)]   # S1 as a perfect cx
    rep = is_tilting_complex(T)
    assert rep.is_tilting and rep.rigid and rep.generates
    from quiverlab.derived.tilting import corner_cartan_of_complex
    assert corner_cartan_of_complex(T) == [[1, 1], [0, 1]]        # = A.cartan_matrix()


@selfcert
def test_missing_summand_fails_generation():
    A = _a2()
    T = [ChainComplex.stalk(A.projective(1), 0)]     # one summand, two simples
    rep = is_tilting_complex(T)
    assert rep.generates is False and rep.is_tilting is False     # g-matrix not square


@selfcert
def test_shifted_copy_fails_rigidity():
    # T = X (+) X[1] with X = stalk(P1): Hom_{D^b}(X, X[1][-1]) = End(X) != 0, so
    # rigidity fails at n = -1 (and symmetrically at n = +1).
    A = _a2()
    X = ChainComplex.stalk(A.projective(1), 0)
    T = [X, X.shift(1)]
    rep = is_tilting_complex(T)
    assert rep.rigid is False and rep.is_tilting is False
    assert rep.window[0] <= -1 <= rep.window[1]


@selfcert
def test_two_term_silting_from_presentation():
    A = _a2()
    cx, rep = two_term_silting_from_presentation(A.simple(1))
    assert set(cx.degrees()) <= {0, 1} and cx.is_perfect()
    assert rep.rigid                                 # a 2-term silting object is rigid
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
def _span(cx):
    ds = cx.degrees()
    return (ds[0], ds[-1]) if ds else (0, 0)


def _chi(cx, verts):
    out = [0] * len(verts)
    for n in cx.degrees():
        dv = cx.term(n).dimension_vector()
        for i, v in enumerate(verts):
            out[i] += (-1) ** n * dv.get(v, 0)
    return out


def _direct_sum_complex(summands):
    """The perfect complex (+)_i T_i (block-diagonal terms, block-diagonal diffs),
    reusing complexes._block_diag_module degreewise. Perfect (each summand is)."""
    ...  # block-diagonal per degree; carry _perfect / _proj_vertices provenance


def is_tilting_complex(summands):
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.complexes import hyper_hom_dims
    import sympy as sp
    for T in summands:
        if not T.is_perfect():
            raise QuiverlabError("is_tilting_complex: every summand must be a "
                                 "certified perfect complex")
    A = summands[0].algebra
    verts = list(A.quiver.vertices)
    spans = [_span(T) for T in summands]
    n_min = min(lo for lo, _ in spans) - max(hi for _, hi in spans)
    n_max = max(hi for _, hi in spans) - min(lo for lo, _ in spans)
    Tsum = _direct_sum_complex(summands)
    hh = hyper_hom_dims(Tsum, Tsum, n_min, n_max)
    rigid = all(hh.get(n, 0) == 0 for n in range(n_min, n_max + 1) if n != 0)
    g = [_chi(T, verts) for T in summands]
    det = int(sp.Matrix(g).det()) if len(g) == len(verts) else 0
    generates = (len(g) == len(verts)) and det in (1, -1)
    return TiltingReport(is_tilting=rigid and generates, rigid=rigid,
                         generates=generates, window=(n_min, n_max),
                         g_matrix=g, det=det)


def end_algebra_of_complex(summands):
    """End_{D^b}(T) via degree-0 hyper-Hom classes between all summand pairs,
    composed with ChainMap.then and reduced to canonical homotopy reps, then handed
    to from_structure_constants (mirrors endomorphism._structure_constants)."""
    from quiverlab.derived.homs import hyper_hom_basis
    ...  # (a) collect basis = [class for i,j in pairs for class in
         #     hyper_hom_basis(T_i, T_j, 0)] with a block index (i,j);
         # (b) product basis[a] * basis[b] = basis[b].then(basis[a]) reduced to its
         #     canonical homotopy rep, expressed over `basis` by solve_columns (the
         #     coset canonicalisation of Task 1); composable only when the middle
         #     complex matches (else the structure constant is 0);
         # (c) T, unit -> Algebra.from_structure_constants(T, unit, field=dom, check=True).


def corner_cartan_of_complex(summands):
    """dim_k Hom_{D^b}(T_i, T_j) (degree 0) as an integer matrix -- the corner-Cartan
    of End(T). = cartan_matrix(A) when T=A (the regular_corner_dims analogue)."""
    from quiverlab.derived.homs import hyper_hom_basis
    return [[len(hyper_hom_basis(Ti, Tj, 0)) for Tj in summands] for Ti in summands]


def two_term_silting_from_presentation(M):
    from quiverlab.modules.complexes import ChainComplex
    from quiverlab.modules.resolution import minimal_resolution
    terms, dmats = minimal_resolution(M, 1)      # P_1 --d_1--> P_0 --> M
    cx = ChainComplex({0: terms[0].module, 1: terms[1].module}, {1: dmats[1]},
                      check=True)
    cx._perfect = True
    cx._proj_vertices = {0: list(terms[0].vertices), 1: list(terms[1].vertices)}
    rep = is_tilting_complex([cx])               # rigidity report (generation is n/a
                                                 # for a single object; report it honest)
    return cx, rep
```

**Adjust to reality:**
- `corner_cartan_of_complex` counts hyper-Hom **classes**; for `T = A` (projective
  stalks) `dim Hom_{D^b}(P_v, P_w) = dim e_v A e_w = C[v][w]` (higher Ext vanish), so
  it equals `cartan_matrix(A)` — the corner oracle. Whether the matrix comes out
  `C` or `Cᵀ` is fixed by the direction convention of `hyper_hom_basis`; the
  `test_ka2_apr_tilt` value `[[1,1],[0,1]]` and the `test_regular_module` equality
  pin it — if transposed, transpose the return **once** and comment.
- `end_algebra_of_complex`' composition must be reduced to canonical homotopy reps
  **before** reading structure constants (two chain maps equal in `D^b` differ by a
  coboundary); reuse Task 1's `_coset_reps`/`reduce_mod_nullspace` so `End(T)` is
  byte-reproducible. The `_direct_sum_complex` is only needed for the `is_tilting`
  rigidity scan; `end_algebra_of_complex` works pairwise (no direct sum required).
- `two_term_silting`'s `TiltingReport.generates` for a **single** object is `False`
  (one row, `#vertices` columns) — that is correct (a 2-term *silting* object need
  not be *tilting*); the AIR consumer (P45) reads `.rigid` + the g-matrix
  (a single g-vector), not `.is_tilting`. Document this in the return.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/derived/tilting.py tests/modules/test_derived_tilting.py
git commit -m "feat(derived): tilting-complex verifier (rigidity+K0 generation) + End(T) as an algebra"
```

---

### Task 4: `fingerprint.py` + `block.py` — the derived fingerprint

**Files:**
- Create: `src/quiverlab/derived/fingerprint.py`, `src/quiverlab/derived/block.py`
- Test: `tests/modules/test_derived_fingerprint.py`

**Interfaces:**
- Consumes (one call each): `A.coxeter_polynomial()` (P38 — exact sympy `Poly`),
  `A.cartan_matrix()` (P38), `sympy.matrices.normalforms.{smith_normal_form,
  invariant_factors}` (verified importable — sympy 1.14.0 in the venv; **unused in
  the repo, wired fresh here** for the Cartan integer-equivalence class),
  `A.hochschild_cohomology(top)` / `A.hochschild_homology(top)`,
  `A.cyclic_homology(top)`, `A.center()` (→ `(dim, basis)`), `A.global_dimension()`
  (honest `repr`). All exact over any Domain (Plans 19/35).
- Produces:
  ```python
  def derived_fingerprint(A, top=4) -> dict
      # {"coxeter_polynomial": str|None, "cartan_det": int|None,
      #  "cartan_smith": [int]|None,        # invariant factors of C (Smith normal form)
      #  "hh_cohomology_dims": [int], "hh_homology_dims": [int],
      #  "cyclic_dims": [int], "center_dim": int, "gl_dim": str}   # gl_dim = honest repr
      # Fields that RAISE on this input (singular Cartan -> coxeter/smith;
      # presentation-less -> cartan) are captured as {"error": <msg>} per field
      # (the recognizers per-flag-error precedent), NEVER a 500.
  def compare_fingerprints(fa, fb) -> dict
      # {"distinguished_by": [field names where fa != fb], "verdict":
      #  "distinguished" | "not distinguished by these invariants"}
      # NECESSARY-CONDITION language ONLY -- never "(in)equivalent". Only fields both
      # sides computed (no {"error"}) are compared.
  ```

**Honest-scope wording (binding, verification page).** The fingerprint is a tuple
of **derived-invariants-or-necessary-conditions**: the Coxeter polynomial and
`|det C|` are derived invariants; HH/HC/center are derived invariants (Rickard);
the Cartan Smith factors are an integer-equivalence invariant (a necessary
condition, coarser than `ℤ`-congruence). Equal fingerprints ⇏ derived equivalent —
`compare_fingerprints` says **"not distinguished by these invariants"**, the
cospectral-trees pin is the standing counterexample.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_derived_fingerprint.py
"""Derived fingerprint: a necessary-condition invariant tuple + honest comparison.
Literature: D4 vs A4 distinguished by the Coxeter polynomial; the 8-vertex
cospectral trees NOT distinguished (equal fingerprint) -- honest verdict."""
import pytest

from quiverlab import Quiver
from quiverlab.derived.fingerprint import derived_fingerprint, compare_fingerprints

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _d4():
    return Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra()


def _a4():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra()


def _spider():
    return Quiver(list(range(1, 9)),
                  {"a1": (1, 2), "a2": (2, 3), "a3": (3, 4),
                   "b": (1, 5), "c": (1, 6), "d": (1, 7), "e": (1, 8)}).algebra()


def _double_star():
    return Quiver(list(range(1, 9)),
                  {"m": (1, 2), "a": (1, 3), "b": (1, 4), "c": (1, 5),
                   "d": (2, 6), "e": (2, 7), "f": (2, 8)}).algebra()


@lit
def test_d4_a4_distinguished_by_coxeter():
    fa, fd = derived_fingerprint(_a4()), derived_fingerprint(_d4())
    res = compare_fingerprints(fa, fd)
    assert "coxeter_polynomial" in res["distinguished_by"]
    assert res["verdict"] == "distinguished"


@lit
def test_cospectral_trees_not_distinguished():
    fs = derived_fingerprint(_spider())
    fd = derived_fingerprint(_double_star())
    res = compare_fingerprints(fs, fd)
    assert res["distinguished_by"] == []
    assert res["verdict"] == "not distinguished by these invariants"
    # (these trees are NOT derived equivalent -- the honest-scope demonstration)


@selfcert
def test_fingerprint_never_claims_equivalent():
    res = compare_fingerprints(derived_fingerprint(_a4()), derived_fingerprint(_a4()))
    assert "equivalent" not in res["verdict"]        # necessary-condition language only


@selfcert
def test_singular_or_presentationless_fields_are_honest_errors():
    from quiverlab.families import truncated_polynomial
    from quiverlab.fields import GF
    fp = derived_fingerprint(truncated_polynomial(2, field=GF(7)))   # non-unimodular C
    # coxeter still computes (det C = 2 != 0 -> Phi integral, t+1); center/HH compute;
    # gl_dim is an honest infinite lower bound.
    assert isinstance(fp["cartan_det"], int) and fp["cartan_det"] == 2
    assert "gl.dim" in fp["gl_dim"] or ">=" in fp["gl_dim"]
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
"""Derived fingerprint (Plan 43): a tuple of derived-invariants / necessary
conditions and an honest comparison. NEVER a derived-equivalence decider -- the
metaplan honest-scope rule: 'distinguished / not distinguished by these
invariants', never '(in)equivalent'."""
from __future__ import annotations

import sympy as sp
from sympy.matrices.normalforms import invariant_factors

from quiverlab.errors import QuiverlabError


def _field(fn):
    try:
        return fn()
    except QuiverlabError as exc:
        return {"error": str(exc)}


def derived_fingerprint(A, top=4):
    fp = {}
    fp["coxeter_polynomial"] = _field(lambda: str(A.coxeter_polynomial().as_expr()))
    fp["cartan_det"] = _field(lambda: int(sp.Matrix(A.cartan_matrix()).det()))
    fp["cartan_smith"] = _field(
        lambda: [int(x) for x in invariant_factors(sp.Matrix(A.cartan_matrix()))])
    fp["hh_cohomology_dims"] = list(A.hochschild_cohomology(top).dims)
    fp["hh_homology_dims"] = list(A.hochschild_homology(top).dims)
    fp["cyclic_dims"] = list(A.cyclic_homology(top))
    fp["center_dim"] = A.center()[0]
    fp["gl_dim"] = repr(A.global_dimension())
    return fp


def compare_fingerprints(fa, fb):
    distinguished = []
    for key in fa:
        va, vb = fa.get(key), fb.get(key)
        if isinstance(va, dict) or isinstance(vb, dict):   # a field errored one side
            continue                                       # not comparable -> skip
        if va != vb:
            distinguished.append(key)
    verdict = ("distinguished" if distinguished
               else "not distinguished by these invariants")
    return {"distinguished_by": distinguished, "verdict": verdict}
```

`derived/block.py` (the shared GUI/webapp/HPC block builder — Task 6 calls it):

```python
"""Shared derived-fingerprint block for the three runners (Plan 43 / Task 6)."""
from quiverlab.derived.fingerprint import derived_fingerprint


def derived_fingerprint_block(A, top=4):
    fp = derived_fingerprint(A, top)
    keys = ["happel_triangulated", "rickard_derived", "lenzing_delapena_spectral"]
    return {"kind": "derived_fingerprint", "top": top, "fingerprint": fp,
            "latex": _fp_latex(fp), "references": keys}   # citations added by the caller
```

**Adjust to reality:** `A.hochschild_cohomology(top)` returns a table with `.dims`
(read `spec.py:1160`'s `A.cyclic_homology` handling for the exact accessor);
`A.cyclic_homology(top)` returns the dims list (or `(table, reps)` under
`with_reps=True` — use the plain form). For a tree path algebra HH is cheap
(`HH^0 = k`, `HH^{≥1} = 0`); keep `top` small (default 4). `cartan_smith` uses
`invariant_factors` (Smith normal form's diagonal) — verify the import path against
the installed sympy (`sympy.matrices.normalforms.invariant_factors`, confirmed
importable in the venv). Do **not** claim congruence: it is the two-sided
`GL_n(ℤ)`-equivalence class (a necessary condition), state so in the docstring.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/derived/fingerprint.py src/quiverlab/derived/block.py \
        tests/modules/test_derived_fingerprint.py
git commit -m "feat(derived): derived fingerprint + honest necessary-condition comparison (D4/A4, cospectral trees)"
```

---

### Task 5: QPA cross-oracle battery (`TauOfComplex`, probe-first)

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py`, `src/quiverlab/qpa/crosscheck.py`
- Test: `tests/qpa/test_derived_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.run`, `session.should_skip_qpa`,
  `session.libgap_handle`), `scripts.quiver_and_algebra_script`+`scripts.module_decl`
  (the module serialiser), the existing `crosscheck_tau` (`crosscheck.py:133` — the
  DTr/TrD template), QPA Ch. 10–11 `TauOfComplex`/`Complex`/`HomologyOfComplex` —
  **probe live for the exact names first** (the Plan-35/P39 `NamesGVars()`/
  `IsBoundGlobal` sweep precedent).
- Produces: `crosscheck_tau_complex(algebra, M)` — build `X` = projective
  resolution of `M` both sides, compare `tau_Db(X)` **homology dimension vectors**
  against QPA's `TauOfComplex`; wire `what == "tau_complex"` into the `crosscheck`
  dispatch ladder + hint list. **Honest documented fallback (P39 Ch.10 finding):**
  if QPA's `Complex`/`TauOfComplex` objects do not script cleanly through libgap
  (streams/category-object oddities), fall back to the **module-level** oracle —
  `tau_Db(resolution of M)` homology(0) vs QPA `DTr(M)` via the existing
  `crosscheck_tau` — and record exactly what was compared in the test docstring.
  Never a silent skip.

- [ ] **Step 1: Probe + write the battery**

```python
# tests/qpa/test_derived_qpa.py
"""QPA as the oracle for the derived AR translate (Plan 43). qpa-marked: skips
locally, mandatory under QUIVERLAB_REQUIRE_QPA=1. Probes TauOfComplex; documents
the module-level fallback (P39 Ch.10 complex-scripting hazard)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_derived_surface_probe():
    lg = session.libgap_handle()
    present = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
               for name in ("TauOfComplex", "HomologyOfComplex", "MappingCone")}
    # record what QPA exposes; the battery below uses TauOfComplex if present, else
    # the DTr module-level fallback (documented, not skipped).
    assert present["HomologyOfComplex"]              # Ch.10 is present in QPA 1.37


@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ),
                               linear_path_algebra(4, field=QQ)])
def test_tau_complex_homology_vs_qpa(A):
    for v in A.quiver.vertices:
        M = A.simple(v)
        if M.tau().dim == 0:                         # projective end: skip
            continue
        A.crosscheck("tau_complex", M).assert_agree()
```

- [ ] **Step 2: Implement the crosscheck** in `crosscheck.py` (mirror
  `crosscheck_tau`): build `M` as a QPA module; if `TauOfComplex` scripts, construct
  the projective-resolution complex in QPA, apply `TauOfComplex`, read
  `HomologyOfComplex` dimension vectors, compare to `tau_Db(...)` homology dim
  vectors (`_flat_dimvec_multiset`); otherwise compare `tau_Db(...).homology(0)` to
  QPA `DTr(M)` via `crosscheck_tau`'s DTr path. Append `what == "tau_complex"` to the
  dispatch ladder (`crosscheck.py:369+`) and its hint list.

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_derived_qpa.py -v` (the venv
  has the `[qpa]` extra). Expected: PASS live (or the documented DTr fallback).

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/qpa/ tests/qpa/test_derived_qpa.py
git commit -m "test(qpa): derived tau battery -- TauOfComplex homology vs tau_Db, documented DTr fallback"
```

---

### Task 6: GUI story — the `derived_fingerprint` scalar compute kind

One scalar kind (algebra-only, **schema v1** — no new request block), wired through
the standard **seven-touchpoint** pattern (P38 `cartan`/`coxeter_polynomial`/
`ext_algebra` precedent, `spec.py:1174-1214`). **The two-algebra compare panel is
DEFERRED to P50** (it needs a second-algebra request field — a schema change);
record the deferral on the verification page (the P39 GUI-deferral precedent). The
single-algebra `derived_fingerprint` block ships now.

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch` scalar branch + `_snip` recipe)
- Modify: `docs/gui/runner.py` (the byte-identical Pyodide twin handler + ETA)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox id
  `qlgui-derived_fingerprint`, `S.ids`, push list, block renderer)
- Modify: `webapp/static/app.js` (webapp block renderer)
- Modify: `webapp/server/i18n/en.json`, `es.json` (labels)
- Modify: `src/quiverlab/trace/results_html.py` (report renderer branch)
- Modify: `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` (ONE new
  fixture, existing entries byte-identical)
- Test: `tests/webapp/test_derived_fingerprint_p43.py`

**Interfaces:**
- `derived_fingerprint` is a **scalar kind on the algebra block** (schema v1). The
  `spec.py` branch mirrors `ext_algebra`/`recognizers`:
  ```python
  if kind == "derived_fingerprint":
      top = item.hi if item.hi is not None else 4
      from quiverlab.derived.block import derived_fingerprint_block
      block = derived_fingerprint_block(A, top)
      block["citations"] = _citation_pairs(block["references"])
      return block, None
  ```
  The Pyodide twin (`docs/gui/runner.py`) calls the SAME
  `derived.block.derived_fingerprint_block`, so the blocks are byte-identical.
- Block: `{"kind": "derived_fingerprint", "top": int, "fingerprint": {...all
  fields, errored ones as {"error": ...}...}, "latex": ..., "references":
  ["happel_triangulated", "rickard_derived", "lenzing_delapena_spectral"],
  "citations": ...}`. The renderer shows the Coxeter polynomial (MathJax), Cartan
  det + Smith factors, HH/HC/center dims, gl.dim (honest), and one line of
  necessary-condition scope text ("a derived-invariant fingerprint; equal values are
  a necessary condition for derived equivalence, not a proof").

- [ ] **Step 1: Write the failing cross-runner test** (unmarked, extras-gated dir —
  copy `tests/webapp/test_module_blocks_m0729.py`'s runner-pair fixture)

```python
# tests/webapp/test_derived_fingerprint_p43.py
"""derived_fingerprint scalar kind: schema-v1, served by hpc.spec, mirrored
byte-identically by the Pyodide twin."""
import json


def test_derived_fingerprint_block_shape(tmp_path):
    # request: kA4 (path quiver 1->2->3->4) via the standard algebra request dict,
    # compute ["derived_fingerprint"]. Assert:
    #   block["fingerprint"]["coxeter_polynomial"] == "t**4 + t**3 + t**2 + t + 1"
    #   block["fingerprint"]["cartan_det"] == 1
    #   "happel_triangulated" in [k for k, _ in block["citations"]]
    ...


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py; assert
    # json.dumps(sort_keys=True) equality on the derived_fingerprint block.
    ...
```

(Fill the `...` by copying the m0729 fixture; the assertions are the contract.)

- [ ] **Step 2: Implement** the `spec.py` branch + the `docs/gui/runner.py` twin
  (shape-identical), the two `gui.js` renderers + `app.js`, the ETA entry
  (`"derived_fingerprint": 1.0` — HH to `top=4` dominates), i18n keys
  (`inv.derived_fingerprint`, `block.derived_fingerprint.title`,
  `block.derived_fingerprint.scope`, field labels — EN and ES), the `_snip` recipe
  (`"derived_fingerprint": "derived_fingerprint(A)"`), and the `results_html.py`
  branch.

- [ ] **Step 3: Add ONE golden fixture** (`derived_fingerprint_a4`) to
  `_runner_goldens.json`; note it in the `test_runner_delegation.py` docstring
  change-log. Verify existing goldens stay byte-identical BEFORE adding.

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_derived_fingerprint_p43.py tests/webapp/test_runner_delegation.py tests/hpc tests/gui -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): derived_fingerprint scalar kind (schema v1); compare panel deferred to P50"
```

---

### Task 7: citations, verification page, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/registry.py` (add `happel_triangulated ->
  Happel1988`; **reuse** `rickard_derived -> Rickard1989`, already registered)
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-05-metaplan-v0.2.0.md` (tick the P43 card delivery note)
- Modify: `src/quiverlab/__init__.py` / `core/algebra.py` (public exports — see below)
- Test: existing release gates (`tests/release/test_oracle_classes.py`) + one
  citation-presence assertion.

**Citations (BibTeX-verified — both entries already in `references.bib`):**
- `Happel1988` — Happel, *Triangulated Categories in the Representation Theory of
  Finite Dimensional Algebras*, LMS Lecture Note Series **119**, CUP, 1988
  (`references.bib:82`, DOI `10.1017/CBO9780511629228`) — the Serre-functor /
  `τ_{D^b}` / infinite-gl.dim refusal ground truth. Add the registry key
  `happel_triangulated -> Happel1988` (mirroring the `_r(...)` idiom,
  `registry.py:24`); **NOTE** the existing `happel_trivial_extension` also points at
  `Happel1988` but is semantically the trivial-extension use — add the clean derived
  alias, do not overload it.
- `Rickard1989` — Rickard, *Morita theory for derived categories*, J. London Math.
  Soc. **39** (1989) (`references.bib:250`, already keyed `rickard_derived`,
  `registry.py:165`) — `End(T)` = the derived-equivalent algebra. Reuse the key.

**Public surface:** export `from quiverlab.derived import (hyper_hom_basis, tau_Db,
tau_Db_minus, is_tilting_complex, end_algebra_of_complex,
two_term_silting_from_presentation, derived_fingerprint, compare_fingerprints)` via
`derived/__init__.py`; read `src/quiverlab/__init__.py` for the top-level export
idiom and add `derived_fingerprint`/`compare_fingerprints` + the `tau_Db`/tilting
names as the public derived surface (keep the engine-internal discipline — no
reaching into `quiverlab.engine.*`).

- [ ] **Step 1: Add the `happel_triangulated` registry key** (BibTeX already
  verified present):

```python
_r("happel_triangulated", "Happel1988", "foundation",
   "Triangulated Categories in the Representation Theory of Finite Dimensional Algebras",
   "The derived-category reference: the Serre functor / AR triangles of D^b(mod A) "
   "exist iff gl.dim < infinity, tau_{D^b} = nu[-1] -- the ground truth for Plan 43.",
   "book"),
```

- [ ] **Step 2: Verification page** — add the Plan-43 subsystem row
  (`derived/` | oracles: self-cert chain-map/`d²=0`/round-trip + cross-engine
  hyper-Hom≡Ext & complex-`τ`≡module-`τ` + literature kA_n tilting/mesh & D₄/A₄
  Coxeter distinction + live QPA `TauOfComplex`), the **honest-scope entries**
  (deciding derived equivalence is not algorithmic — verifiers + necessary-condition
  fingerprint only; `τ_{D^b}` refuses loudly at infinite gl.dim; the compare panel is
  deferred to P50), and recount the class table (`tests/release/
  test_oracle_classes.py` drives the numbers — run collection, paste, re-run to
  green; note "counts as of the P43 merge, sibling plans in flight"). Add the
  Happel1988/Rickard1989 pins to the Class-1 literature list. README: one line in the
  features list ("derived category: hyper-Hom classes, the derived AR translate
  `τ_{D^b}`, a tilting-complex verifier with `End(T)`, and a derived fingerprint").

- [ ] **Step 3: Full gate:**
  `... -m pytest tests/modules -q` (deep, the touched files),
  `... -m pytest -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release -q`,
  and a citation-presence check (`happel_triangulated`/`rickard_derived` resolve; the
  `derived_fingerprint` block carries them) — all green.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-43 derived-category oracle rows + Happel1988 citation + recounted classes"
```

---

## Acceptance (Plan-43 definition of done)

1. `hyper_hom_basis`, `ChainMap.then`, `tau_Db`/`tau_Db_minus`,
   `is_tilting_complex`/`TiltingReport`/`end_algebra_of_complex`/
   `corner_cartan_of_complex`/`two_term_silting_from_presentation`,
   `derived_fingerprint`/`compare_fingerprints` all public in
   `src/quiverlab/derived/` (with the `Algebra`/`__init__` exports named in the
   tasks), loudly validated, self-certified.
2. `hyper_hom_basis` reifies genuine chain maps `X → Y[n]` (`ChainMap(check=True)`),
   its count equals `hyper_hom_dims` **and** module `Ext^n` on a resolution source
   (cross-engine); `ChainMap.then` composes them.
3. `τ_{D^b} = ν∘[−1]` returns a `ChainComplex` (`d²=0` self-cert); on a
   non-projective indecomposable over `kA_n` its homology is concentrated in
   degree 0 and `≅` the trusted module `τ` (cross-engine), the K₀ identity
   `χ(τ_{D^b} X) = c·χ(X)` (`c = −C·C⁻ᵀ`, the K₀-action Coxeter matrix) holds
   (literature), and **infinite global dimension refuses loudly** (`k[x]/(x²)`
   pinned). `τ⁻_{D^b}` certified by the round-trip quasi-iso.
4. `is_tilting_complex` decides rigidity on the **exact reported window** and
   generation by `ℤ`-unimodularity; `T = A` is tilting with `End(T) ≅ A`
   (corner-Cartan = `cartan_matrix`), the kA₂ APR tilt `P₁ ⊕ S₁` is tilting with the
   reoriented-A₂ Cartan (derivation in the test comment), a missing summand fails
   generation and `X ⊕ X[1]` fails rigidity (both concretely derived).
5. `derived_fingerprint` distinguishes D₄ from A₄ (Coxeter polynomial) and
   **does not** distinguish the 8-vertex cospectral trees — `compare_fingerprints`
   in necessary-condition language only (never "(in)equivalent"); singular-Cartan /
   presentation-less fields captured as honest per-field errors.
6. Live QPA battery green (`-m qpa`): `TauOfComplex` homology vs `tau_Db`, or the
   documented module-level `DTr` fallback (P39 Ch.10 hazard).
7. The `derived_fingerprint` compute kind clickable end-to-end (GUI canvas → block →
   report) in EN+ES, schema-v1 scalar kind, one golden added with a documented
   reason, both runners byte-identical; the two-algebra compare panel deferral to
   P50 recorded on the verification page.
8. `docs/verification.md` updated (new oracle rows, recounted classes with the
   merge-train note, honest-scope entries); Happel1988 registry key added and
   BibTeX-verified, `rickard_derived` reused; `tests/release/` green; deep + qpa +
   fast buckets green on the touched surface.
