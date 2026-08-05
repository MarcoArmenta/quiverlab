# Plan 41: C3 Auslander–Reiten Completion — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finish the Auslander–Reiten surface. Almost-split sequences with the
middle term built and certified (`0 → τM → E → M → 0`); irreducible maps and
`rad(M,N)/rad²`; AR-quiver **knitting** with an honest semi-decision contract
(complete iff rep-finite, loud budget cap otherwise); the Nakayama functor ν /
ν⁻ named as a first-class functor; and stable Hom. The τ translate itself is
already delivered (`Module.tau()/tau_minus()`, Plan 23/24) — this plan adds the
sequences, maps, and quiver that τ anchors.

**Architecture:** One new module, `src/quiverlab/modules/ar.py`, a thin
exact-linear-algebra layer over primitives that already exist. Almost-split
sequences reuse `yoneda.baer_extension` (the pushout middle term is already
built and self-certifies exactness) fed a socle cocycle of `Ext¹(M, τM)`; the
socle is picked by the **left End(M)-action** on `Ext¹(M, τM)` annihilated by
`rad End(M)` (ARS IV.1–IV.3). The End(M)-action needs a general chain-map lift
(`ext_algebra.py::_lift` is simple-target-only — this plan generalises its
per-generator solve to any module endomorphism). ν is the cokernel of the
induced injective map from the minimal presentation, tied to the trusted τ by a
`ker ≅ τM` certificate. Stable Hom is the cokernel of "compose with the cover".
Knitting is a BFS from the projectives using `almost_split_sequence` +
`decompose` + `identify_standard`. No new math engines; every result
self-certifies (exact, non-split, ends indecomposable) and is cross-checked
against QPA and closed-form Dynkin/Nakayama ground truth.

**Tech Stack:** pure exact linear algebra over `Domain` (`modules/linalg_mod`,
`fields/linalg`), `ModuleHom`/`ShortExactSequence`/`direct_sum`/`hom_basis`/
`is_direct_summand`/`end_algebra` from **Plan 37** (this plan REQUIRES P37
merged), `yoneda.baer_extension`, `complex_reps.ext_reps`, `decompose`,
`duality`. No floats in `src/`.

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- `tests/modules/` auto-assigns to the **deep** bucket; `tests/qpa/` to the
  **qpa** bucket (`tests/conftest.py`). Run new tests by path during
  development; finish each task with a `-m deep` (or `-m qpa`) spot-run of the
  touched files.
- **P37 is the prerequisite** — consume `ModuleHom`, `ShortExactSequence`,
  `direct_sum`, `is_direct_summand`, `hom_basis`, `Module.identity_hom`,
  `Module.end_algebra` exactly as their signatures read in
  `docs/plans/2026-08-05-plan-37-categorical-glue.md`. Branch
  `plan-41-ar-completion` off `dev` **after P37 has merged**.
- **`decompose` refuses over GF(p) when char ≤ dim M** (the trace-form radical
  is unreliable there — see `modules/decompose.py`). Every AR battery that
  decomposes middle terms or certifies indecomposability runs over **QQ or
  GF(32003)** so both the split search and the locality certificate decide.
  `rad End(M)` (Task 5) reuses the same trace-form radical, same char scope.
- House conventions: a module-map matrix is dense `tgt.dim × src.dim` over the
  shared Domain (`hom_space`'s layout). Composition is left-to-right, `f.then(g)`
  (P37). Hom/comparisons refuse loudly across sides/algebras
  (`hom._assert_comparable`). All refusals are `QuiverlabError`; `check=False`
  fast paths only for internally-constructed data whose certificate is asserted
  separately.
- **Semi-decision honesty (metaplan §6):** AR knitting is complete **iff**
  rep-finite; a wild/large algebra hits the budget and refuses **loudly** with a
  `"budget"` status — never a silently truncated "AR quiver". State this on the
  verification page as each piece lands.
- Plan-32 markers: exactness/certificate/identity tests = `oracle_selfcert`;
  the AR formula (`Ext¹ ≅ stable Hom`) and any two-independent-route dim
  agreement = `oracle_crossengine`; Dynkin/Nakayama closed-form vertex counts and
  mesh structure = `oracle_literature`; QPA comparisons live in `tests/qpa/`
  (bucket = the class, never double-marked).
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted
  class counts, `tests/release/test_oracle_classes.py` green) and adds its
  citations to `citations/references.bib` + `registry.py`.
- Conventional commits; green tests at every commit.

---

### Task 1: the general chain-map lift `lift_endomorphism_along_resolution`

The net-new primitive P41 needs and `ext_algebra.py::_lift` cannot supply
(it lifts a summand-projection cocycle into the resolution of a **simple**;
here we lift an arbitrary module **endomorphism** `φ: M → M` over `M`'s own
minimal resolution).

**Files:**
- Create: `src/quiverlab/modules/ar.py`
- Test: `tests/modules/test_ar_lift.py`

**Interfaces:**
- Consumes: `modules/resolution.py::minimal_resolution(M, length) -> (terms, dmats)`
  (`terms[n].module`, `terms[n].vertices`, `terms[n].dim`; `dmats[n]: P_n → P_{n-1}`
  shape `dim P_{n-1} × dim P_n`; **`dmats[0]` is the augmentation `P_0 → M`**);
  `ModuleHom` (P37); `linalg_mod` (`matmul`, `matvec`, `zeros`, `col`);
  `fields.linalg.solve` + `reduce_mod_nullspace` (byte-reproducible solve, the
  Plan-17 canonicalisation `ext_algebra._solve_in_component` already uses).
- Produces:
  ```python
  def lift_endomorphism_along_resolution(phi, res=None, degrees=1):
      # phi: a ModuleHom M -> M (or its raw M.dim x M.dim matrix + a module M).
      # res: (terms, dmats) from minimal_resolution(M, degrees+1); built if None.
      # Returns [phi_0, phi_1, ..., phi_degrees], phi_n a dim P_n x dim P_n matrix
      #   with  eps . phi_0 == phi . eps           (phi_0 covers phi; eps = dmats[0])
      #   and   d_n . phi_n == phi_{n-1} . d_n     (a chain map over the resolution).
      # Canonical (reduce_mod_nullspace) so the lift is byte-reproducible.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_lift.py
"""General chain-map lift of a module endomorphism over its minimal resolution.
Self-certifying: every square d_n . phi_n == phi_{n-1} . d_n is asserted exactly,
and identity lifts to identities."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.ar import lift_endomorphism_along_resolution
from quiverlab.modules.morphism import hom_basis
from quiverlab.modules.resolution import minimal_resolution

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(32003))


def _squares_commute(phis, dmats, dom):
    # eps . phi_0 == phi . eps is the caller's cover condition; here check the
    # interior chain-map squares d_n . phi_n == phi_{n-1} . d_n.
    for n in range(1, len(phis)):
        dn = dmats[n]
        lhs = lm.matmul(dn, phis[n], dom)
        rhs = lm.matmul(phis[n - 1], dn, dom)
        if lhs != rhs:
            return False
    return True


def test_identity_lifts_to_identities():
    A = _a3()
    for v in (1, 2, 3):
        M = A.projective(v) if v == 1 else A.simple(v)
        idM = M.identity_hom()
        phis = lift_endomorphism_along_resolution(idM, degrees=3)
        for n, ph in enumerate(phis):
            assert ph == lm.identity(len(ph), M.domain)


def test_endomorphism_squares_commute():
    A = _a3()
    S2 = A.simple(2)
    terms, dmats = minimal_resolution(S2, 4)
    for phi in hom_basis(S2, S2):                 # every endo of S2
        phis = lift_endomorphism_along_resolution(phi, (terms, dmats), degrees=3)
        # cover condition eps . phi_0 == phi . eps
        eps = dmats[0]
        assert lm.matmul(eps, phis[0], A.domain) == lm.matmul(phi.matrix, eps,
                                                              A.domain)
        assert _squares_commute(phis, dmats, A.domain)


def test_lift_is_deterministic():
    A = _a3()
    P1 = A.projective(1)
    phi = hom_basis(P1, P1)[0]
    a = lift_endomorphism_along_resolution(phi, degrees=3)
    b = lift_endomorphism_along_resolution(phi, degrees=3)
    assert a == b                                 # reduce_mod_nullspace => byte-stable
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_ar_lift.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.modules.ar`

- [ ] **Step 3: Implement** — the cover step solves against the augmentation
  `eps = dmats[0]`; the interior steps mirror
  `ext_algebra._build_module_map`/`_solve_in_component` (read them first — the
  per-generator, per-vertex-component solve is verbatim reusable, only the "seed"
  differs: `_lift` seeds a generator-to-generator map, we seed the cover).

```python
"""Auslander-Reiten completion (Plan 41 / C3): the general chain-map lift, the
Nakayama functor, stable Hom, the End(M)-action on Ext^1, almost-split
sequences, irreducible maps, and AR-quiver knitting. Right modules; exact over
any Domain. Every constructed object self-certifies -- correctness never rests
on trusting a construction."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.resolution import minimal_resolution


def _summand_layout(term, algebra):
    """[(vertex, offset, dim, basis_labels)] for the direct-sum blocks of a
    resolution term P_n = (+)_s P_{v_s} -- rebuilt from projective(A, v), whose
    ordered path basis matches the block layout the resolution used (the exact
    idiom of ext_algebra._SimpleResolution.info / complex_reps._terms_info)."""
    from quiverlab.modules.builders import projective
    out, off = [], 0
    for v in term.vertices:
        Pv = projective(algebra, v)
        out.append((v, off, Pv.dim, Pv._pv_basis_labels))
        off += Pv.dim
    return out


def _solve_generator_image(dtarget, target_mod, vertex, rhs, dom):
    """Solve d^target . y = rhs for a generator image y in the vertex-`vertex`
    component of the target term (the exact ext_algebra._solve_in_component
    routine: restrict to the idempotent-`e_vertex` columns, solve, canonicalise
    with reduce_mod_nullspace, scatter back). Empty component => rhs must be 0."""
    ev = target_mod.action[f"e_{vertex}"]
    dimT = target_mod.dim
    cols_v = [c for c in range(dimT) if not dom.is_zero(ev[c][c])]
    if not cols_v:
        return [dom.zero()] * dimT
    Dred = [[dtarget[r][c] for c in cols_v] for r in range(len(dtarget))]
    z = linalg.solve(Dred, rhs, dom)
    if z is None:
        raise QuiverlabError(
            "ar: chain-map lift solve is inconsistent",
            hint="on a genuine minimal resolution the obstruction is a boundary; "
                 "report the presentation")
    z = linalg.reduce_mod_nullspace(z, Dred, dom)
    y = [dom.zero()] * dimT
    for c, zc in zip(cols_v, z):
        y[c] = zc
    return y


def _map_from_generator_images(layout, target_mod, images, dom):
    """The full matrix P_src -> target sending each summand generator (the block's
    first basis vector) to images[s], extended A-linearly on the path basis:
    (generator_s . path) |-> images[s] . path = action[path] @ images[s]."""
    G = lm.zeros(target_mod.dim, sum(dv for _, _, dv, _ in layout), dom)
    for (v, off, dv, labels), y in zip(layout, images):
        for k, plabel in enumerate(labels):
            col = lm.matvec(target_mod.action[plabel], y, dom)
            for r in range(target_mod.dim):
                G[r][off + k] = col[r]
    return G


def lift_endomorphism_along_resolution(phi, res=None, degrees=1):
    from quiverlab.modules.morphism import ModuleHom
    if isinstance(phi, ModuleHom):
        M, phimat = phi.src, phi.matrix
        if phi.tgt is not phi.src:
            raise QuiverlabError("ar: lift needs an endomorphism M -> M")
    else:                                         # (matrix, M) -- read from res owner
        raise QuiverlabError("ar: pass a ModuleHom endomorphism")
    dom = M.domain
    if res is None:
        res = minimal_resolution(M, degrees + 1)
    terms, dmats = res
    layout0 = _summand_layout(terms[0], M.algebra)
    eps = dmats[0]                                # P_0 ->> M augmentation
    # phi_0 covers phi: eps . phi_0 = phi . eps, solved per P_0 generator against eps.
    imgs0 = []
    for (v, off, dv, labels) in layout0:
        gen_col = [eps[r][off] for r in range(len(eps))]        # eps(generator_s) in M
        rhs = lm.matvec(phimat, gen_col, dom)                   # phi(eps(gen_s))
        imgs0.append(_solve_generator_image(eps, M, v, rhs, dom))
    phis = [_map_from_generator_images(layout0, terms[0].module, imgs0, dom)]
    # phi_{n}: d_n . phi_n = phi_{n-1} . d_n, per generator via the component solve.
    for n in range(1, degrees + 1):
        layout = _summand_layout(terms[n], M.algebra)
        dn = dmats[n]
        imgs = []
        for (v, off, dv, labels) in layout:
            gencol = [dn[r][off] for r in range(len(dn))]       # d_n(generator_s)
            rhs = lm.matvec(phis[n - 1], gencol, dom)           # phi_{n-1}(d_n gen)
            imgs.append(_solve_generator_image(dmats[n], terms[n - 1].module, v,
                                               rhs, dom))
        phis.append(_map_from_generator_images(layout, terms[n].module, imgs, dom))
    return phis
```

**Adjust to reality:** `_solve_generator_image` assumes the target term's
`e_v` action is a diagonal idempotent on the `P_w` block basis — verify against
`builders.projective` (it is, by construction). If a term's `module` is `None`
(resolution exhausted) before `degrees`, raise loudly naming the depth.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar.py tests/modules/test_ar_lift.py
git commit -m "feat(ar): general chain-map lift of a module endomorphism over its minimal resolution"
```

---

### Task 2: the Nakayama functor ν / ν⁻ named

**Files:**
- Modify: `src/quiverlab/modules/ar.py`
- Modify: `src/quiverlab/modules/module.py` (add `Module.nakayama()`, `Module.nakayama_minus()`)
- Test: `tests/modules/test_ar_nakayama.py`

**Interfaces:**
- Consumes: `builders.injective(A, v)` (`ν(P_v) = I_v = D(A e_v)`, `builders.py:74`),
  `builders.projective`, `duality.transpose_module` (the corner-transpose `d_1^*`
  matrix, `duality.py:43-99` — reuse its corner-slicing, do NOT rederive),
  `duality.dualize`, `duality.tau` (the TRUSTED translate — the `ker ≅ τM`
  arbiter), `ModuleHom.cokernel`/`.kernel` (P37), `is_isomorphic`,
  `identify_standard`, `minimal_resolution(M, 1)`.
- Produces:
  ```python
  def nakayama_functor(M) -> Module          # ν M = D Hom_A(M, A) = coker(νP_1 -> νP_0)
  def nakayama_functor_minus(M) -> Module    # ν^- M = Hom_A(DA, M); ν^-(I_v) = P_v
  ```
  `Module.nakayama()` / `.nakayama_minus()` delegate.

**Math pinned (ARS):** for `M` with minimal presentation `P_1 --d_1--> P_0 -> M -> 0`,
applying `ν = D Hom_A(-, A)` gives the induced map `g: νP_1 -> νP_0` between
injectives with **`ker g = τM`** and **`coker g = νM`** (the two ends of
`0 -> τM -> νP_1 -> νP_0 -> νM -> 0`). `g` is `D` of the corner-transpose `d_1^*`
that `transpose_module` builds. `ν^-` is implemented as the opposite-algebra
dual `ν^-_A M = D( ν_{A^op}( D M ) )` (contravariant D each side of the
A^op Nakayama functor) — no separate injective-copresentation code, and
`ν^-(I_v) = P_v` self-certifies it.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_nakayama.py
"""Nakayama functor. Self-cert: nu(P_v) ~ I_v, nu^-(I_v) ~ P_v, and the ker of
the induced injective map ~ tau M (ties nu to the trusted translate). Additive."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import nakayama_functor, nakayama_functor_minus
from quiverlab.modules.hom import is_isomorphic

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def test_nu_of_projective_is_injective():
    A = _a3()
    for v in (1, 2, 3):
        assert is_isomorphic(nakayama_functor(A.projective(v)), A.injective(v))


def test_nu_minus_of_injective_is_projective():
    A = _a3()
    for v in (1, 2, 3):
        assert is_isomorphic(nakayama_functor_minus(A.injective(v)), A.projective(v))


def test_nu_ties_to_tau_on_nonprojectives():
    # ker(nu P_1 -> nu P_0) ~ tau M is the internal certificate; here we check the
    # PUBLIC consequence: nu preserves the dimension vector of the projective cover
    # data and tau M is recovered by the module's own tau() (arbiter).
    A = _a3()
    S1 = A.simple(1)
    # S1 is non-projective over kA3 with a*b=0? S1 = simple at the source: it IS
    # not projective here; tau S1 is nonzero. The nu certificate is asserted inside
    # nakayama_functor via an internal assert; this test pins the public tie:
    assert S1.tau().dim == nakayama_functor(S1)  # SHARPEN in Step 3 -- see note
    # (the real assertion is dim/iso of ker g vs tau S1; write it against the
    #  internal handle the implementation exposes for the test, e.g. a debug tuple)


def test_nu_additive():
    A = _a3()
    from quiverlab.modules.morphism import direct_sum
    D, _, _ = direct_sum(A.simple(1), A.simple(2))
    nD = nakayama_functor(D)
    n1, n2 = nakayama_functor(A.simple(1)), nakayama_functor(A.simple(2))
    assert nD.dimension_vector() == {v: n1.dimension_vector().get(v, 0)
                                     + n2.dimension_vector().get(v, 0)
                                     for v in A.quiver.vertices}
```

**Sharpening rule for `test_nu_ties_to_tau_on_nonprojectives`:** the placeholder
line is a stand-in. During Step 3 expose the induced map `g` (or its kernel) so
the test can assert `is_isomorphic(ker_g, S1.tau())` exactly — the internal
`assert` in `nakayama_functor` is the real certificate; the test pins the same
identity through a small debug handle (e.g. `nakayama_functor(M, _return_ker=True)`
or a module-private `_nu_kernel`). A test that cannot be sharpened means ν is
wrong — stop and debug, do not weaken.

- [ ] **Step 2: Run to verify failure** — `ImportError: nakayama_functor`

- [ ] **Step 3: Implement** — build `νP_1`, `νP_0` as direct sums of injectives
  over the presentation's summand vertices; build the induced injective map `g`
  by reusing `transpose_module`'s corner-slicing (`duality.py:79-99`) and
  dualising; `νM = ModuleHom(νP_1, νP_0, g).cokernel()`. **Certify per instance:**
  `is_isomorphic(ModuleHom(νP_1, νP_0, g).kernel()[0], duality.tau(M))` (raise
  `"nakayama: ker g not isomorphic to tau M"` otherwise) and, for each `P_v`,
  `identify_standard(νP_v) == ("injective", v)`. `ν^-` via
  `dualize(nakayama_functor(dualize(M).with_side(...)))` on `A^op` — read
  `duality.dualize` to get the side bookkeeping right (D flips the side), and
  self-cert `ν^-(I_v) ≅ P_v`.

  **Adjust to reality:** `transpose_module` returns the *cokernel module* (Tr M),
  not the raw induced map — factor its corner-slicing (lines 79-99, up to and
  including `d1star`) into a shared helper `_corner_transpose_matrix(M)` in
  `ar.py` (or lift it into `duality.py` and import), so both `transpose_module`
  and `nakayama_functor` build `d1_star` from one place. Then `g = D(d1_star)`
  is a matrix between the injective bases; confirm the injective-basis ordering
  matches `builders.injective`'s (paths ending at `v`) — if not, permute once and
  pin it with the `ker ≅ τM` arbiter.

- [ ] **Step 4: Run tests** (+ the trusted-τ neighbours):

Run: `... -m pytest tests/modules/test_ar_nakayama.py tests/modules/test_duality*.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar.py src/quiverlab/modules/module.py tests/modules/test_ar_nakayama.py
git commit -m "feat(ar): Nakayama functor nu/nu^- -- certified nu(P_v)=I_v, ker(induced)=tau M"
```

---

### Task 3: stable Hom `stable_hom_dim` + `hom_factors_through_projective`

**Files:**
- Modify: `src/quiverlab/modules/ar.py`
- Modify: `src/quiverlab/modules/morphism.py` (add `ModuleHom.factors_through_projective()`)
- Test: `tests/modules/test_ar_stable_hom.py`

**Interfaces:**
- Consumes: `hom_basis` (P37), `Module.projective_cover_hom()` (P37 Task 6 —
  `P(N) ->> N`), `linalg_mod.solve_columns`/`cols_to_matrix`/`mat_rank`.
- Produces:
  ```python
  def stable_hom_dim(M, N) -> int
      # dim underline{Hom}(M, N) = dim Hom(M,N) - dim P(M,N),
      # P(M,N) = image of  Hom(M, P(N)) --(compose with cover pi_N)--> Hom(M, N).
      # (Any map factoring through a projective factors through the cover of N.)
  def hom_factors_through_projective(f) -> bool   # f a ModuleHom: f in P(M,N)?
  ```
  `ModuleHom.factors_through_projective()` delegates.

**Cross-engine oracle (AR formula, ARS):**
`dim Ext¹_A(M, N) = dim underline{Hom}(τ⁻N, M)` — pinned degreewise on Dynkin
simples where both sides are cheap and closed-form.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_stable_hom.py
"""Stable Hom (mod projectives). Self-cert: underline Hom(P, -) = 0 (source
projective), and a cover composite factors through a projective. Cross-engine:
the Auslander-Reiten formula dim Ext^1(M,N) = dim underline Hom(tau^- N, M)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import hom_factors_through_projective, stable_hom_dim
from quiverlab.modules.morphism import hom_basis

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _a4():
    return linear_path_algebra(4, field=QQ)


@selfcert
def test_stable_hom_from_projective_vanishes():
    A = _a4()
    P1 = A.projective(1)
    for v in (1, 2, 3, 4):
        assert stable_hom_dim(P1, A.simple(v)) == 0     # source projective


@selfcert
def test_cover_composite_factors_through_projective():
    A = _a4()
    S1 = A.simple(1)
    pi = S1.projective_cover_hom()                       # P(S1) ->> S1, factors trivially
    # any g: S1 -> P(S1) composed with pi factors through the projective P(S1)
    for g in hom_basis(S1, pi.src):
        assert hom_factors_through_projective(g.then(pi))


@xeng
def test_auslander_reiten_formula():
    A = _a4()
    simples = [A.simple(v) for v in (1, 2, 3, 4)]
    for M in simples:
        for N in simples:
            lhs = A.ext(M, N, 1)
            rhs = stable_hom_dim(N.tau_minus(), M)
            assert lhs == rhs, (M.name, N.name, lhs, rhs)
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
def _projective_maps_columns(M, N):
    """Columns (over the vec of the tgt.dim x src.dim matrix space) spanning
    P(M,N) = { (compose with cover) . g : g in Hom(M, P(N)) }."""
    dom = M.domain
    pi = N.projective_cover_hom()                 # P(N) ->> N  (P37 Task 6)
    cols = []
    for g in hom_basis(M, pi.src):                # g: M -> P(N)
        comp = g.then(pi).matrix                  # M -> N through the projective
        cols.append([comp[i][j] for j in range(M.dim) for i in range(N.dim)])
    return cols


def stable_hom_dim(M, N):
    from quiverlab.modules.hom import _assert_comparable
    _assert_comparable(M, N, "stable Hom")
    dom = M.domain
    homs = hom_basis(M, N)
    hom_cols = [[f.matrix[i][j] for j in range(M.dim) for i in range(N.dim)]
                for f in homs]
    proj_cols = _projective_maps_columns(M, N)
    total = lm.mat_rank(lm.cols_to_matrix(hom_cols), dom) if hom_cols else 0
    if not proj_cols:
        return total
    both = lm.mat_rank(lm.cols_to_matrix(hom_cols + proj_cols), dom)
    # P(M,N) subset Hom(M,N), so rank(hom_cols + proj_cols) == rank(hom_cols):
    #   dim P(M,N) = rank(proj_cols); underline dim = total - dim P(M,N).
    dim_proj = lm.mat_rank(lm.cols_to_matrix(proj_cols), dom)
    assert both == total, "stable_hom_dim: P(M,N) is not inside Hom(M,N) (bug)"
    return total - dim_proj


def hom_factors_through_projective(f):
    dom = f.domain
    proj_cols = _projective_maps_columns(f.src, f.tgt)
    fvec = [[f.matrix[i][j] for j in range(f.src.dim) for i in range(f.tgt.dim)]]
    if not proj_cols:
        return f.is_zero()
    sol = lm.solve_columns(lm.cols_to_matrix(proj_cols),
                           lm.cols_to_matrix(fvec), dom)
    return sol is not None
```

**Adjust to reality:** `Module.projective_cover_hom` is P37 Task 6 — confirm its
matrix orientation (`P(N) ->> N`, `tgt = N`) so `.then(pi)` type-checks. The
vec convention (`[... for j in range(dim_src) for i in range(dim_tgt)]`) must be
consistent between the two helpers and with any coordinate you compare — keep the
single order used above.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar.py src/quiverlab/modules/morphism.py tests/modules/test_ar_stable_hom.py
git commit -m "feat(ar): stable Hom (mod projectives) + factors-through-projective; AR-formula cross-engine pin"
```

---

### Task 4: the End(M)-action on `Ext¹(M, N)`

**Files:**
- Modify: `src/quiverlab/modules/ar.py`
- Test: `tests/modules/test_ar_end_action.py`

**Interfaces:**
- Consumes: `Task 1` (`lift_endomorphism_along_resolution`), `hom_basis(M, M)`
  (P37 — the End(M) basis; or `Module.end_algebra()` when the algebra structure is
  wanted), `complex_reps.ext_reps(A, M, N, 1)` +
  `complex_reps._reconstruct_cocycle(col, homs_n, n_dim, width, dom)` (the tested
  cocycle-reconstruction the Yoneda pipeline `_one_ext_sequence` mirrors),
  `hom_basis(P_0, N)` for the coboundaries, `solve_columns`.
- Produces:
  ```python
  def end_action_on_ext1(M, N):
      # Returns (basis, action) with:
      #   basis   = list of End(M) basis endomorphisms (ModuleHom, from hom_basis(M,M))
      #   action  = list of e x e matrices (e = dim Ext^1(M,N)); action[i] is the
      #             matrix of  [f] |-> [f . (phi_i)_1]  (left End(M)-action by
      #             precomposition of the lifted degree-1 chain map) in the
      #             ext_reps cohomology basis of Ext^1(M,N).
  ```

**Convention (arbitrated, not assumed):** the action is `[f] ↦ [f ∘ φ₁]`
(precompose the degree-1 lift). Whether this is the socle-carrying action for the
almost-split socle (Task 5) is **decided by Task 5's non-split certificate**, not
asserted here — if the socle came out wrong the AR sequence would split and the
test would fail (the house convention-flip pattern: End-composition order, cup
sign). Task 4 only certifies the representation axioms below.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_end_action.py
"""End(M) acts on Ext^1(M, N) by precomposition of the lifted chain map.
Self-cert: id acts as identity; the action is a representation of End(M)
(phi.psi acts as the product of the action matrices, left-to-right)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.ar import end_action_on_ext1

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def test_identity_acts_as_identity():
    A = _a3()
    M, N = A.simple(2), A.simple(2).tau()          # Ext^1(M, tau M) is the AR case
    basis, action = end_action_on_ext1(M, N)
    e = len(action[0]) if action else 0
    idx = next(i for i, f in enumerate(basis) if f.is_iso() and f.matrix
               == M.identity_hom().matrix)
    assert action[idx] == lm.identity(e, A.domain)


def test_action_is_a_representation():
    A = _a3()
    M, N = A.simple(2), A.simple(2).tau()
    basis, action = end_action_on_ext1(M, N)
    dom = A.domain
    # pick two basis endos; their composite (left-to-right, f.then(g)) must act as
    # the product of action matrices (in the matching order).
    from quiverlab.modules.morphism import ModuleHom
    for i in range(len(basis)):
        for j in range(len(basis)):
            comp = basis[i].then(basis[j])          # M -> M
            # express comp over the End basis, combine action matrices linearly,
            # compare to action of comp computed directly (helper in ar.py).
            direct = _action_of(comp, M, N, dom)    # test-local: recompute via ar
            combined = lm.matmul(action[j], action[i], dom)  # precompose order
            assert direct == combined
```

(`_action_of` is a one-line test helper that calls the same internal
`_ext1_action_matrix(phi, ...)` `end_action_on_ext1` uses — expose it privately;
the representation-order assertion (`action[j] @ action[i]` vs `f.then(g)`) is the
arbiter that fixes the left/right bookkeeping. If it fails, swap the product
order once and document, per house style.)

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
def _ext1_data(A, M, N):
    """(cocycle_mats, coboundary_mats, cocycle_cols, homs1) for Ext^1(M, N).
    cocycle_mats[j] = f_j: P_1 -> N (N.dim x P_1.dim) the reconstructed class reps;
    coboundary_mats = { psi . d_1 : psi in Hom(P_0, N) } (the delta^0 image)."""
    from quiverlab.modules.complex_reps import _reconstruct_cocycle, ext_reps
    from quiverlab.modules.morphism import hom_basis
    dom = A.domain
    terms, dmats = minimal_resolution(M, 2)
    dims, payload = ext_reps(A, M, N, 1)
    # rebuild the ambient Hom(P_1, N) basis matrices exactly as ext_reps did:
    from quiverlab.modules.complex_reps import _hom_adjunction_basis, _terms_info
    tinfo = _terms_info(M.algebra, terms)
    homs1 = _hom_adjunction_basis(N, tinfo[1], dom, {})[1]
    width = len(dmats[1][0]) if (dmats[1] and dmats[1][0]) else 0
    cocycle_mats = [_reconstruct_cocycle([int(coeff) and coeff for coeff in _col_of(cls)],
                                         homs1, N.dim, width, dom)
                    for cls in payload["basis_classes"]["1"]]
    d1 = dmats[1]
    coboundary_mats = [lm.matmul(psi.matrix, d1, dom)
                       for psi in hom_basis(terms[0].module, N)]
    return cocycle_mats, coboundary_mats, terms, dmats


def _ext1_action_matrix(phi, cocycle_mats, coboundary_mats, res, dom):
    """The e x e matrix of [f] |-> [f . phi_1] on the Ext^1 basis, phi an
    endomorphism of M, phi_1 its degree-1 lift (Task 1)."""
    phis = lift_endomorphism_along_resolution(phi, res, degrees=1)
    phi1 = phis[1]
    e = len(cocycle_mats)
    basis_cols = [_vec(f) for f in cocycle_mats]
    bnd_cols = [_vec(b) for b in coboundary_mats]
    B = lm.cols_to_matrix(basis_cols + bnd_cols) if (basis_cols or bnd_cols) else []
    out = lm.zeros(e, e, dom)
    for j, f in enumerate(cocycle_mats):
        composed = lm.matmul(f, phi1, dom)        # f . phi_1 : P_1 -> N
        sol = lm.solve_columns(B, lm.cols_to_matrix([_vec(composed)]), dom)
        if sol is None:
            raise QuiverlabError("ar: f.phi_1 is not a cocycle mod coboundaries "
                                 "(the lift or the class is inconsistent)")
        for i in range(e):
            out[i][j] = sol[0][i]                 # cocycle-basis part only
    return out
```

(`_vec(mat)` = column-stacked flatten; `_col_of(cls)` reconstructs the dense
coefficient column from a `basis_classes` entry's sparse `vector` field — read
`complex_reps.serialize_class` for the sparse shape and invert it. **Adjust to
reality:** the cleaner route is to skip the serialized reps entirely and pull
the raw cohomology columns from `ext_reps`'s internal `cols_by_deg` — if that is
not returned, either reconstruct from `basis_classes[...]["vector"]` or add a
thin internal accessor to `complex_reps.ext_reps`; do NOT recompute the Ext
complex a second way, that would risk a basis mismatch.)

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar.py tests/modules/test_ar_end_action.py
git commit -m "feat(ar): End(M)-action on Ext^1(M,N) via the general lift; representation-axiom certified"
```

---

### Task 5: almost-split sequences `almost_split_sequence`

**Files:**
- Modify: `src/quiverlab/modules/ar.py`
- Modify: `src/quiverlab/modules/module.py` (add `Module.almost_split_sequence()`)
- Test: `tests/modules/test_ar_almost_split.py`

**Interfaces:**
- Consumes: `is_indecomposable`/`decompose` (Plan 30, budget 512, char scope),
  `identify_standard`, `duality.tau`, `Task 4` (`end_action_on_ext1` with
  `N = τM`), `decompose._trace_form_rank`-style trace-form radical for
  `rad End(M)`, `complex_reps._reconstruct_cocycle`, `yoneda.baer_extension`,
  `ShortExactSequence` + `ModuleHom` (P37), `linalg_mod.kernel_columns`.
- Produces:
  ```python
  def almost_split_sequence(M) -> ShortExactSequence
      # For M certified-indecomposable and NON-projective: the almost-split (AR)
      # sequence  0 -> tau M -> E -> M -> 0.  Self-certified: exact + non-split +
      # both ends indecomposable.  Loud refusal for projective M, decomposable M,
      # or an undecidable (budget/char) input.
  ```
  `Module.almost_split_sequence()` delegates.

**Algorithm (ARS IV.1–IV.3):** for `M` indecomposable non-projective, `End(M)` is
local (Fitting), so `rad End(M)` = the trace-form radical (char 0 or char > dim M).
`End(M)` acts on `Ext¹(M, τM)` (Task 4); the elements annihilated by `rad End(M)`
form the **socle**, which is nonzero and (over `End/rad = k`) 1-dimensional; ANY
nonzero socle class `ξ` gives THE almost-split sequence. Reconstruct `ξ`'s cocycle
`f: P_1 → τM`, `baer_extension(M, τM, f)` builds the pushout `E` and the exact
`0 → τM → E → M → 0`; certify.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_almost_split.py
"""Almost-split sequences. Self-cert: exact + non-split + indecomposable ends.
Literature: the kA_n AR sequences have the mesh middle-term dimension vectors."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a3():
    return linear_path_algebra(3, field=QQ)


@selfcert
def test_almost_split_is_exact_nonsplit_indecomposable_ends():
    A = _a3()
    # the injective-but-not-projective / interior indecomposable: S1 over kA3
    M = A.simple(1)
    assert not M.is_indecomposable() is False       # (S1 is indecomposable)
    ses = M.almost_split_sequence()                 # 0 -> tau M -> E -> M -> 0
    assert ses.N.is_isomorphic(M) or ses.N.dim == M.dim   # right end is M
    assert ses.M.dim == ses.L.dim + M.dim           # exact (P37 SES certifies)
    assert ses.is_split() is False                  # almost-split => non-split
    from quiverlab.modules.decompose import is_indecomposable
    assert is_indecomposable(ses.L)                 # tau M indecomposable
    assert is_indecomposable(M)


@selfcert
def test_projective_input_refused():
    A = _a3()
    with pytest.raises(QuiverlabError, match="projective"):
        A.projective(1).almost_split_sequence()     # no AR sequence ends at a projective


@selfcert
def test_decomposable_input_refused():
    A = _a3()
    from quiverlab.modules.morphism import direct_sum
    D, _, _ = direct_sum(A.simple(1), A.simple(2))
    with pytest.raises(QuiverlabError, match="indecomposable"):
        # almost_split_sequence needs an indecomposable end
        from quiverlab.modules.ar import almost_split_sequence
        almost_split_sequence(D)


@lit
def test_ka3_middle_terms_match_the_mesh():
    # kA3 (1->2->3, linear) AR quiver is the ZA3 slice; the AR sequence ending at
    # the interval module [1..3] has the two intervals [1..2],[2..3] as middle,
    # etc. Pin the middle-term dimension vectors against the closed-form mesh.
    A = _a3()
    # interval modules realised as indecomposable projectives/simples over kA3:
    M = A.simple(2)                                 # interval [2]
    ses = M.almost_split_sequence()
    got = sorted(ses.M.dimension_vector().values())
    # SHARPEN: fill the exact mesh dimension vector from the ZA3 picture in Step 3.
    assert sum(got) == ses.M.dim
```

(`test_ka3_middle_terms_match_the_mesh` is written against the closed-form ZAₙ
mesh — replace the placeholder `sum(...)` with the exact middle-term dimension
vector during Step 3, deriving it from the kA₃ AR quiver picture and pinning it.)

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
def _rad_end_basis(M):
    """A basis of rad End(M) as End-basis endomorphisms. M indecomposable =>
    End local => rad = trace-form radical (Dickson/CIW: char 0 or char > dim M).
    Reuses decompose's trace-form; refuses loudly outside that char scope."""
    from quiverlab.modules.hom import hom_space
    dom = M.domain
    char = dom.characteristic
    if not (char == 0 or char > M.dim):
        raise QuiverlabError(
            f"almost_split: rad End(M) unreliable in char {char} <= dim M = {M.dim}",
            hint="run over QQ or a characteristic > dim M (e.g. GF(32003))")
    H = hom_space(M, M)                              # End basis (matrices)
    # Gram matrix T[i][j] = tr(H_i H_j); rad End = kernel of T (as End-coordinates).
    r = len(H)
    T = lm.zeros(r, r, dom)
    for i in range(r):
        for j in range(r):
            prod = lm.matmul(H[i], H[j], dom)
            T[i][j] = _trace(prod, dom)
    return H, lm.kernel_columns(T, dom)              # columns = rad basis in End coords


def almost_split_sequence(M):
    from quiverlab.modules.complex_reps import _reconstruct_cocycle
    from quiverlab.modules.decompose import is_indecomposable
    from quiverlab.modules.morphism import ModuleHom
    from quiverlab.modules.ses import ShortExactSequence
    from quiverlab.modules.duality import tau as tau_of
    from quiverlab.modules.builders import _require_provenance
    dom = M.domain
    if not is_indecomposable(M):                     # raises loudly if undecidable
        raise QuiverlabError("almost_split: M must be indecomposable")
    if identify_standard(M) and identify_standard(M)[0] == "projective":
        raise QuiverlabError("almost_split: no AR sequence ends at a projective M")
    tauM = tau_of(M)
    if tauM.dim == 0:
        raise QuiverlabError("almost_split: tau M = 0 -- M is projective")
    # End(M) acts on Ext^1(M, tau M); the socle = annihilator of rad End(M).
    A = M.algebra
    basis, action = end_action_on_ext1(M, tauM)
    Hbasis, rad_coords = _rad_end_basis(M)           # rad basis in End coords
    # action of a rad element = the same End-coord combination of the action matrices
    stacked = []
    e = len(action[0]) if action else 0
    for rc in rad_coords:
        Arad = lm.zeros(e, e, dom)
        for i, c in enumerate(rc):
            if dom.is_zero(c):
                continue
            for a in range(e):
                for b in range(e):
                    Arad[a][b] = dom.add(Arad[a][b], dom.mul(c, action[i][a][b]))
        stacked.append(Arad)
    socle_cols = _joint_kernel(stacked, e, dom)      # cap_r ker(A_rad)
    if not socle_cols:
        raise QuiverlabError("almost_split: empty socle of Ext^1(M, tau M) -- "
                             "the End(M)-action or the class basis is inconsistent")
    xi = socle_cols[0]                               # any nonzero socle class
    cocycle_mats, _cob, terms, dmats = _ext1_data(A, M, tauM)
    # f = xi-combination of the Ext^1 basis cocycles f_j : P_1 -> tau M
    f = _combine_matrices(cocycle_mats, xi, tauM.dim, len(dmats[1][0]), dom)
    seq = baer_extension(M, tauM, f, terms, dmats)   # 0 -> tau M -> E -> M -> 0
    seq.assert_exact()
    ses = ShortExactSequence(ModuleHom(seq.modules[0], seq.modules[1], seq.maps[0]),
                             ModuleHom(seq.modules[1], seq.modules[2], seq.maps[1]))
    if ses.is_split():
        raise QuiverlabError("almost_split: constructed sequence splits -- the socle "
                             "class was not the almost-split class (bug/convention)")
    if not (is_indecomposable(seq.modules[0]) and is_indecomposable(M)):
        raise QuiverlabError("almost_split: an end is decomposable (bug)")
    return ses
```

**Certification basis (state in the docstring):** by ARS IV.1–IV.3 a sequence
`0 → τM → E → M → 0` whose class lies in `soc Ext¹(M, τM)` (annihilated by
`rad End(M)`), with both ends indecomposable, **is** almost split. The returned
object additionally carries the computational self-certs — exact
(`YonedaSequence.assert_exact` + P37 `ShortExactSequence`), non-split (P37
`is_split() is False`), ends indecomposable (Plan 30) — so nothing rests on
trusting the socle pick: a wrong pick splits and is refused loudly.

**Adjust to reality:** `_baer` cocycle shape is `τM.dim × P_1.dim` (`yoneda`
checks it). `_joint_kernel(mats, e, dom)` stacks the `A_rad` matrices vertically
and returns `kernel_columns`. `_combine_matrices` linearly combines the cocycle
matrices by `xi`. If `end_action_on_ext1` returned the socle-carrying action
under the *opposite* order, the `is_split()` guard fires — flip Task 4's product
convention once and re-run (the sequence is the arbiter).

- [ ] **Step 4: Run tests** — Expected: PASS (with the mesh pin sharpened)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar.py src/quiverlab/modules/module.py tests/modules/test_ar_almost_split.py
git commit -m "feat(ar): almost-split sequences 0->tauM->E->M->0 -- socle pick, certified exact/nonsplit"
```

---

### Task 6: irreducible maps `irreducible_maps`

**Files:**
- Modify: `src/quiverlab/modules/ar.py`
- Test: `tests/modules/test_ar_irreducible.py`

**Interfaces:**
- Consumes: `hom_basis`, `is_isomorphic` (M ≅ N branch), `_rad_end_basis`
  (Task 5, for `rad End(M)` when `M ≅ N`), `linalg_mod` rank/solve.
- Produces:
  ```python
  def irreducible_maps(M, N, within) -> int
      # dim rad(M,N)/rad^2(M,N) -- the AR-quiver arrow multiplicity M -> N, with
      # rad^2 computed by composing through the finite indecomposable set `within`
      # (the knitted component's modules).  M, N indecomposable.
      # rad(M,N) = Hom(M,N) when M !~ N; = rad End(M) when M ~ N.
      # rad^2(M,N) = sum_{L in within} span{ h.then(g) : h in rad(M,L), g in rad(L,N) }.
      # Exact once `within` contains all indecomposables (rep-finite).
  ```

**Spec-ambiguity resolution (recorded):** the metaplan card says "returns the
arrow-multiplicity matrix". `irreducible_maps(M, N, within)` returns the single
**integer** multiplicity for the pair `(M, N)`; the full AR-quiver arrow **matrix**
is assembled by `knit_ar_quiver` (Task 7) by calling this over `within × within`.
This keeps the signature honest (three positional args = a pair query) and the
matrix a derived artefact.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_irreducible.py
"""Irreducible maps = dim rad(M,N)/rad^2(M,N). Literature: kA_n has exactly the
mesh arrows (each interior irreducible-map multiplicity 1)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import irreducible_maps

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def _indecs(A):
    # kA3 indecomposables = the 6 interval modules; realise the small ones here.
    return [A.simple(v) for v in (1, 2, 3)] + [A.projective(v) for v in (1, 2)]


@selfcert
def test_no_irreducible_self_map_for_bricks():
    A = _a3()
    S2 = A.simple(2)
    within = _indecs(A)
    assert irreducible_maps(S2, S2, within) == 0     # rad End(S2) = 0 (a brick)


@lit
def test_ka3_mesh_arrow_multiplicities():
    A = _a3()
    within = _indecs(A)
    P1, P2, S2 = A.projective(1), A.projective(2), A.simple(2)
    # in the ZA3 mesh the irreducible maps have multiplicity 1 along the mesh edges
    # and 0 off it -- pin a couple of edges (SHARPEN the exact pairs in Step 3).
    assert irreducible_maps(S2, P1, within) in (0, 1)
    assert all(irreducible_maps(X, X, within) == 0 for X in within)  # no loops
```

(The `in (0, 1)` placeholders are sharpened in Step 3 to the exact kA₃ mesh
multiplicities read off the ZA₃ picture — the literature ground truth.)

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```python
def _rad_basis(M, N):
    """A spanning column set (vec) of rad(M,N) in the tgt.dim x src.dim vec space:
    all of Hom(M,N) when M !~ N; rad End(M) when M ~ N."""
    from quiverlab.modules.hom import is_isomorphic
    if not is_isomorphic(M, N):                      # non-iso indecs: every map radical
        return [_vec(f.matrix) for f in hom_basis(M, N)]
    H, rad_coords = _rad_end_basis(M)                # rad End(M) in End coords
    return [_vec(_combine_end(H, rc, M.domain)) for rc in rad_coords]


def irreducible_maps(M, N, within):
    dom = M.domain
    rad_MN = _rad_basis(M, N)
    if not rad_MN:
        return 0
    dim_rad = lm.mat_rank(lm.cols_to_matrix(rad_MN), dom)
    # rad^2: compose radical maps M -> L -> N through every L in `within`.
    rad2 = []
    for L in within:
        left = _rad_basis(M, L)                       # vecs of rad(M,L)
        right = _rad_basis(L, N)                       # vecs of rad(L,N)
        for hvec in left:
            h = _unvec(hvec, L.dim, M.dim)             # M -> L
            for gvec in right:
                g = _unvec(gvec, N.dim, L.dim)         # L -> N
                rad2.append(_vec(lm.matmul(g, h, dom)))  # g . h : M -> N
    if not rad2:
        return dim_rad
    both = lm.mat_rank(lm.cols_to_matrix(rad_MN + rad2), dom)
    assert both == dim_rad, "irreducible_maps: rad^2 not inside rad (bug)"
    dim_rad2 = lm.mat_rank(lm.cols_to_matrix(rad2), dom)
    return dim_rad - dim_rad2
```

**Adjust to reality:** `_vec`/`_unvec` are the column-stack pair (shared with
Tasks 3–5 — one definition at the top of `ar.py`). Guard the `M ≅ N` branch: if
`is_isomorphic` refuses (large GF(p^n)), let it raise — `within` is built over
QQ/GF(32003) by the knitting caller, so this is not hit on the oracles.

- [ ] **Step 4: Run tests** — Expected: PASS (mesh pins sharpened)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar.py tests/modules/test_ar_irreducible.py
git commit -m "feat(ar): irreducible-map multiplicity dim rad(M,N)/rad^2 relative to a finite indecomposable set"
```

---

### Task 7: AR-quiver knitting `knit_ar_quiver` + `ARQuiver`

**Files:**
- Modify: `src/quiverlab/modules/ar.py`
- Modify: `src/quiverlab/core/algebra.py` (add `Algebra.ar_quiver(...)` delegating)
- Test: `tests/modules/test_ar_knit.py`

**Interfaces:**
- Consumes: `almost_split_sequence` (Task 5), `duality.tau_minus`, `decompose`,
  `identify_standard`, `is_isomorphic`, `irreducible_maps` (Task 6),
  `builders.projective`/`.radical()`.
- Produces:
  ```python
  class ARQuiver:
      vertices    # list of dicts: {"name": "S_2"|"P_1"|None, "dimvec": {...}, "module": Module}
      arrows      # dict {(i, j): multiplicity} over vertex indices (irreducible maps i->j)
      tau_orbits  # list of vertex-index lists (tau-linked classes)
      is_complete # bool: True iff closure reached (rep-finite), False iff budget-capped
      status      # "complete" | "budget"
  def knit_ar_quiver(A, budget_modules=256, budget_dim=4096) -> ARQuiver
      # BFS from the indecomposable projectives via almost-split sequences;
      # STOPS LOUDLY at the budget with status "budget" (never a silent partial quiver).
  ```
  `Algebra.ar_quiver(budget_modules=256, budget_dim=4096)` delegates.

**Algorithm (classical knitting, honest semi-decision):** seed the discovered set
with `decompose(projective(A, v))`'s indecomposable summands for every `v` (each
`P_v` is indecomposable; also add the arrows into `P_v` from
`decompose(rad P_v)`). Then BFS: for each discovered indecomposable `X` that is
**non-injective** (`identify_standard(X)` not an injective and
`tau_minus(X).dim > 0`), let `Y = tau_minus(X)` (always non-projective
indecomposable), compute `almost_split_sequence(Y)` = `0 → X → E → Y → 0`
(`τY = X`), decompose `E`, register `Y` and every summand `E_i`, and add mesh
arrows `X → E_i` and `E_i → Y` with multiplicities from `irreducible_maps`. New
modules are matched to discovered ones by `is_isomorphic` (dim-vector prefilter
first). Continue until the frontier empties (**complete** — the connected
AR-component closed, so the algebra is rep-finite on this component) or the
budget trips (**budget** — refuse loudly). Every module carries `identify_standard`
naming (unnamed shown by dim vector). `tau_orbits` are read off the `X ↔ tau_minus X`
links.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_ar_knit.py
"""AR-quiver knitting. Literature ground truth: kA_n has n(n+1)/2 indecomposables
(A2->3, A3->6, A4->10), D4 has 12; Nakayama is serial. Contract: a wild algebra
trips the budget LOUDLY (status 'budget', not a silent partial quiver)."""
import pytest

from quiverlab import Quiver, dynkin_quiver, linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import knit_ar_quiver

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_linear_a_n_indecomposable_count(n, count):
    A = linear_path_algebra(n, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    assert len(ar.vertices) == count


@lit
def test_d4_has_twelve_indecomposables():
    A = dynkin_quiver("D4").algebra(relations=[], field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete and len(ar.vertices) == 12


@lit
def test_nakayama_serial_ar_quiver():
    # a self-injective / serial Nakayama algebra: the AR quiver count is closed-form.
    A = NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    # SHARPEN: the exact indecomposable count from the Kupisch series in Step 3.
    assert len(ar.vertices) >= 3


@selfcert
def test_mesh_relations_hold_on_a3():
    # every non-projective vertex X: sum of multiplicities into X (from E summands)
    # equals sum out of tau X -- the mesh/additivity self-check.
    A = linear_path_algebra(3, field=QQ)
    ar = knit_ar_quiver(A)
    for j, vj in enumerate(ar.vertices):
        if vj["name"] and vj["name"].startswith("P_"):
            continue
        into = sum(m for (i, k), m in ar.arrows.items() if k == j)
        assert into >= 1                       # every non-projective is a mesh sink


@selfcert
def test_wild_algebra_trips_budget_loudly():
    # 3-Kronecker (wild): knitting cannot close -- must stop with status 'budget'.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2), "c": (1, 2)}).algebra(
        relations=[], field=QQ)
    ar = knit_ar_quiver(A, budget_modules=40)
    assert ar.is_complete is False and ar.status == "budget"
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** the BFS per the algorithm above. Key points:
  - `_find(discovered, X)`: dim-vector prefilter then `is_isomorphic`; returns the
    index or `None`. New modules get appended with their `identify_standard` name.
  - Budget: trip when `len(discovered) > budget_modules` or any module's
    `dim > budget_dim`; set `status="budget"`, `is_complete=False`, and RETURN
    (do not raise inside the BFS — the ARQuiver object records the loud status;
    `knit_ar_quiver` returns it, callers read `.is_complete`). A hard-refusal
    variant (`Algebra.ar_quiver(..., strict=True)` raising) is optional; the test
    reads the status field.
  - `arrows` accumulate multiplicities via `irreducible_maps(E_i, Y, within)` where
    `within` = the current discovered list (exact once complete; the multiplicity
    is re-read at the end over the final `within` so it is not undercounted by a
    mid-BFS partial set — recompute the arrow matrix in a final pass).
  - `tau_orbits`: union-find over the `X — tau_minus(X)` links.

  **Adjust to reality:** guard `almost_split_sequence` refusals (a module the AR
  machinery cannot certify) — surface them as a loud `status="error"` with the
  offending dim vector, never a silent skip that would under-report the quiver.
  For a self-injective Nakayama algebra the component is periodic (a tube), not a
  postprojective slice — knitting still closes (finitely many indecomposables) but
  seed from projectives AND their `tau_minus` orbit; verify the count oracle.

- [ ] **Step 4: Run tests** — Expected: PASS (counts sharpened)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/modules/ar.py src/quiverlab/core/algebra.py tests/modules/test_ar_knit.py
git commit -m "feat(ar): AR-quiver knitting -- BFS from projectives, honest budget cap, ARQuiver object"
```

---

### Task 8: QPA cross-oracle battery

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py` (AR-sequence / predecessor script builders)
- Modify: `src/quiverlab/qpa/crosscheck.py` (`crosscheck_almost_split`, `crosscheck_predecessors`)
- Test: `tests/qpa/test_ar_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.run`, `session.libgap_handle`),
  `scripts.module_decl`/`quiver_and_algebra_script` (the existing module
  serialiser — read `crosscheck_tau`, `crosscheck.py:133`, it is the template),
  QPA GAP functions `AlmostSplitSequence(M)` and `PredecessorsOfModule(M, n)` —
  **probe live for the exact names first** (`NamesGVars()` grep, the Plan-35
  precedent in `tests/qpa/test_products_qpa.py`): the QPA manual documents both
  but they are UNSCRIPTED here.
- Produces: `crosscheck_almost_split(algebra, M)` — compare the middle-term
  **dimension vector** and its **summand-dimension-vector multiset** (via
  `_flat_dimvec_multiset`, the Plan-30 order-independent invariant already in
  `crosscheck.py`) against QPA's `AlmostSplitSequence(M)` middle; and
  `crosscheck_predecessors(algebra, M)` spot-checking a knitted neighbourhood
  against `PredecessorsOfModule`. `tests/qpa/test_ar_qpa.py` runs them on kA₃ and
  a Nakayama algebra; standard skipif header
  `pytest.mark.skipif(session.should_skip_qpa(), reason=...)`.

- [ ] **Step 1: Probe live QPA** for the AR surface names (mirror
  `tests/qpa/test_products_qpa.py`'s `NamesGVars()` sweep — but here the names are
  expected PRESENT, so the probe confirms them and the battery proceeds; if a name
  is ABSENT the test skips honestly naming what QPA lacks).

```python
# tests/qpa/test_ar_qpa.py
"""QPA as the oracle for almost-split sequences + AR predecessors (Plan 41).
qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1."""
import pytest

from quiverlab import linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ
from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def test_qpa_exposes_almost_split_surface():
    lg = session.libgap_handle()
    for name in ("AlmostSplitSequence", "PredecessorsOfModule"):
        assert bool(lg.eval(f'IsBoundGlobal("{name}")')), name


@pytest.mark.parametrize("A", [linear_path_algebra(3, field=QQ),
                               NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)])
def test_almost_split_middle_vs_qpa(A):
    for v in A.quiver.vertices:
        M = A.simple(v)
        if M.tau().dim == 0:                     # projective end: no AR sequence
            continue
        A.crosscheck("almost_split", M).assert_agree()
```

- [ ] **Step 2: Implement the crosschecks** in `crosscheck.py` (mirror
  `crosscheck_tau`: build `M` as a QPA module via `module_decl`, call
  `AlmostSplitSequence(M)`, read the middle object's `DimensionVector` and
  `DecomposeModule` dimension vectors; compare to
  `M.almost_split_sequence().M` and `decompose` of it via `_flat_dimvec_multiset`),
  and wire the `what="almost_split"`/`"predecessors"` cases into the `crosscheck`
  dispatch (append to the `if what == ...` ladder + the hint list).
  **If `AlmostSplitSequence` returns a QPA complex/streams object that does not
  script cleanly through libgap** (the Plan-39 QPA-Ch.10 hazard), fall back to
  comparing the middle-term dimension vector via
  `ObjectOfComplex(...,1)`/`RangeOfMorphism` and record exactly what was compared
  in the test docstring — honest coverage, never a silent skip.

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_ar_qpa.py -v` (the venv
  has the `[qpa]` extra per project memory). Expected: PASS live.

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/qpa/ tests/qpa/test_ar_qpa.py
git commit -m "test(qpa): AR battery -- AlmostSplitSequence middle-term + PredecessorsOfModule (live)"
```

---

### Task 9: GUI story — the `almost_split` module compute kind

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (`MODULE_KINDS`, `_MOD_REFS`, `_dispatch_module`, `_MODULE_TRACE_KINDS`?)
- Modify: `docs/gui/runner.py` (the byte-identical twin: `_MODULE_KINDS`, `_MOD_REFS`, `_module_block`)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (pick-list entry + block renderer)
- Modify: `webapp/static/app.js` (webapp block renderer)
- Modify: `webapp/server/i18n/en.json`, `es.json` (labels)
- Modify: `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py` (ONE new fixture, existing byte-identical)
- Test: `tests/webapp/test_almost_split_p41.py`

**Interfaces (the tau kind is the closest template — follow its schema-2 gating
exactly, `spec.py:1557`):**
- `almost_split` is a MODULE kind: add it to `MODULE_KINDS`, so a request naming
  it **requires schema 2** with a `module` block (the same guard as `tau`,
  `spec.py:626-631`). `_MOD_REFS["almost_split"] = ["assem_book", "ars_book"]`.
- `_dispatch_module` branch: `M = the built module`; `ses = M.almost_split_sequence()`;
  build the block
  ```python
  {"kind": "almost_split", "exists": true, "side": M.side,
   "tau": _mod_repr(ses.L),                       # tau M as a full representation
   "middle": {"summands": [summand_blocks(s, m) for (s, m) in decompose(ses.M)]},
   "M": _mod_view(M), "indecomposable": True,
   "latex": r"0 \to \tau M \to E \to M \to 0"}
  ```
  A **projective or decomposable or undecidable** input is caught into an honest
  refusal block `{"kind": "almost_split", "exists": false, "reason": <the loud
  message>}` (the `_target_translates` per-entry-error precedent — never a 500,
  the request stays a clean typed 4xx-worthy result).
- Both runners byte-identical on this block (`summand_blocks` is the shared
  serialiser — a summand isomorphic to a standard indecomposable is NAMED, others
  ship full matrices).
- GUI: an `almost_split` entry in the module-panel pick-list (alongside
  `tau`/`tau_minus`); the renderer prints the sequence `0 → τM → E → M → 0`, `E`'s
  summand decomposition, the certified-indecomposability note, and the honest
  refusal text for projective/decomposable input; MathJax + citations like the
  existing module blocks.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked, extras-gated dir
  — copy `tests/webapp/test_module_blocks_m0729.py`'s runner-pair fixture):

```python
# tests/webapp/test_almost_split_p41.py
"""The almost_split module kind: schema-2 gated, served by hpc.spec, mirrored by
the Pyodide twin, honest refusal for projective input."""
import json


def test_almost_split_block_shape(tmp_path):
    # schema-2 request: kA3 over QQ (or GF(32003)), module = builtin simple S_2,
    # compute ["almost_split"]. Assert:
    #   block["exists"] is True
    #   block["tau"] is a full {dims, maps} representation
    #   block["middle"]["summands"] is a non-empty list
    #   "assem_book" in [k for k, _ in block["citations"]]
    ...


def test_almost_split_projective_refused(tmp_path):
    # module = builtin projective P_1: block["exists"] is False, "reason" present,
    # NO exception escapes (clean typed result).
    ...


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py; assert
    # json.dumps(sort_keys=True) equality on the almost_split block.
    ...
```

(Fill the `...` by copying the m0729 fixture; the assertions listed are the
contract. `almost_split` needs `schema: 2` and a `module` block — follow the tau
kind's request shape verbatim.)

- [ ] **Step 2: Implement** the `spec.py` branch + the `docs/gui/runner.py` twin
  (keep them shape-identical), the two `gui.js` renderers + `app.js`, the ETA
  entry (`"almost_split": 0.3` beside `"tau": 0.1`), i18n keys
  (`inv.almost_split`, `block.almost_split.title`, `block.almost_split.sequence`,
  `block.almost_split.middle`, `block.almost_split.refused` — EN and ES),
  `_snip` recipe (`"almost_split": "M.almost_split_sequence()"`), and the
  `results_html.py` branch if the report renders module blocks.

- [ ] **Step 3: Add ONE golden fixture** (`almost_split_a3_s2`) to
  `_runner_goldens.json`; note it in the `test_runner_delegation.py` docstring
  change-log. Verify existing goldens stay byte-identical BEFORE adding.

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_almost_split_p41.py tests/webapp/test_runner_delegation.py tests/hpc tests/gui -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): almost_split module kind -- 0->tauM->E->M->0 clickable, honest refusal"
```

---

### Task 10: verification page, citations, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` (add `@book{ARS1995,...}`)
- Modify: `src/quiverlab/citations/registry.py` (register `ars_book -> ARS1995`)
- Modify: `docs/verification.md`
- Modify: `README.md`
- Modify: `docs/plans/2026-08-05-metaplan-v0.2.0.md` (tick the P41 card delivery note)
- Test: existing release gates (`tests/release/test_oracle_classes.py`), plus one
  citation-presence assertion.

- [ ] **Step 1: Add the ARS citation** (BibTeX-verified — Auslander, Reiten,
  Smalø, *Representation Theory of Artin Algebras*, Cambridge Studies in Advanced
  Mathematics **36**, Cambridge University Press, 1995):

```bibtex
@book{ARS1995,
  author    = {Auslander, Maurice and Reiten, Idun and Smal{\o}, Sverre O.},
  title     = {Representation Theory of {A}rtin {A}lgebras},
  series    = {Cambridge Studies in Advanced Mathematics},
  volume    = {36},
  publisher = {Cambridge University Press},
  year      = {1995},
}
```

and in `registry.py` (mirroring `_r("assem_book", "ASS2006", "foundation", ...)`):

```python
_r("ars_book", "ARS1995", "foundation",
   "Representation Theory of Artin Algebras",
   "The Auslander-Reiten theory reference: almost-split sequences, irreducible "
   "maps, the AR quiver, and the Nakayama functor -- the ground truth for Plan 41.",
   "book"),
```

  **Spec-ambiguity resolution (recorded):** the metaplan says "add ARS1995 if
  missing". The BibTeX **entry key** is `ARS1995`; the **registry/citation key**
  the module blocks reference is snake-case `ars_book` (the house convention —
  `assem_book` → `ASS2006`). `_MOD_REFS["almost_split"]` uses `ars_book`.

- [ ] **Step 2: Verification page** — add the Plan-41 subsystem row
  (`modules/ar.py` | oracles: self-cert lift/nakayama/stable/end-action/almost-
  split certificates + cross-engine AR formula + literature Dynkin/Nakayama counts
  + live QPA `AlmostSplitSequence`/`PredecessorsOfModule`), the honest-scope entry
  (**AR knitting is complete iff rep-finite; loud budget cap otherwise** — the
  metaplan §6 semi-decision ledger), and recount the class table
  (`tests/release/test_oracle_classes.py` drives the numbers — run collection,
  paste, re-run to green). Add the ARS1995 pins to the Class-1 literature list.
  README: one line in the features list ("almost-split sequences, irreducible
  maps, AR-quiver knitting, the Nakayama functor, stable Hom").

- [ ] **Step 3: Full gate:**
  `... -m pytest tests/modules -q` (deep, the touched files),
  `... -m pytest -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release -q`,
  and a citation-presence check (`ars_book` resolves; the `almost_split` block
  carries it) — all green.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-41 AR-completion oracle rows + ARS1995 citation + recounted classes"
```

---

## Acceptance (Plan-41 definition of done)

1. `lift_endomorphism_along_resolution`, `nakayama_functor`/`_minus`,
   `stable_hom_dim` + `hom_factors_through_projective`, `end_action_on_ext1`,
   `almost_split_sequence`, `irreducible_maps`, `knit_ar_quiver`/`ARQuiver` all
   public in `modules/ar.py` (with the `Module`/`Algebra` delegations named in the
   tasks), all loudly-validated, all self-certified.
2. Almost-split sequences returned only when certified: exact (Yoneda + P37 SES),
   non-split (P37 `is_split` False), ends indecomposable (Plan 30); the socle pick
   is arbitrated by the non-split guard, not assumed.
3. ν tied to the trusted τ by `ker(induced) ≅ τM` and `ν(P_v) ≅ I_v`; the AR
   formula `dim Ext¹(M,N) = dim underline Hom(τ⁻N, M)` pinned cross-engine.
4. AR knitting reproduces the closed-form indecomposable counts (kA₂/kA₃/kA₄ =
   3/6/10, D₄ = 12, Nakayama serial) and refuses **loudly** with `status="budget"`
   on a wild algebra — the honest semi-decision contract.
5. Live QPA battery green (`-m qpa`): `AlmostSplitSequence` middle-term multiset +
   `PredecessorsOfModule` spot-checks, or the documented honest-fallback comparison.
6. The `almost_split` compute kind clickable end-to-end (GUI canvas → block →
   report) in EN+ES, schema-2 gated like `tau`, one golden re-frozen with a
   documented reason, honest refusal rendered for projective/decomposable input.
7. `docs/verification.md` updated (new oracle rows, recounted classes, semi-
   decision honest-scope entry); ARS1995 citation added and BibTeX-verified;
   `tests/release/` green; deep + qpa + fast buckets green on the touched surface.
