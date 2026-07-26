# Plan 32 — Oracle-class test markers (v0.1.0 release gate)

**Date:** 2026-07-26  **Branch:** `plan-32-oracle-markers`  **Backlog item:** Tier 1a
"Oracle-class test markers (v0.1.0 release gate)".

## Goal

Ship the four-way oracle taxonomy as **orthogonal pytest markers**, so a JOSS
reviewer can run each oracle class as a one-liner
(`pytest -m oracle_literature`, …). This is a *classification*, not a file
reorganization: **zero file moves**, and the runtime buckets
(`fast`/`deep`/`qpa`) stay **byte-identical** (804 / 1180 / 112). Tests overlap
classes on purpose — a battery pins a literature value *and* asserts cross-engine
agreement in the same test, so it carries both markers.

The four classes:

- **`qpa`** (the existing marker — no new marker): our value ≡ live GAP/QPA.
- **`oracle_literature`** (new): the pass criterion is a value/identity that
  exists **outside the library** — a paper-pinned dimension, a closed form, a
  theorem identity (Coxeter tables, Happel trace, Cartan identities,
  `dim Tor = dim Ext`, `kZ_n/J^L` symmetry classifications, the read-only bank's
  closed forms, CRS/BGMS/BE/LMR pins), frozen as a constant the engine must
  reproduce.
- **`oracle_crossengine`** (new): two **independent implementations** are run and
  asserted equal **live** — CS ≡ bar ≡ Bardzell ≡ minimal degreewise, numba ≡ pure
  parity, presented ≡ ⋉ iso-invariance, native ≡ transported cup/cap, engine ≡
  generic-Domain path.
- **`oracle_selfcert`** (new): an **internal mathematical certificate** is the
  assertion — `d∘d = 0`, the CS order condition, dimension/iso certificates,
  byte-reproducibility / canonicality pins, adversarial-solver determinism,
  arbitration identities (unit / Leibniz / module identities used as sign
  arbiters).
- **UNMARKED** = contract & infrastructure: refusal honesty / error surfaces, API
  contracts, the no-floats AST gate, freshness/interface gates, GUI/webapp/HPC
  plumbing, release/docs/packaging tests, collection guards, foundational datatype
  correctness.

## Crisp class boundary (the rule that keeps the sweep consistent)

The literature/crossengine line is drawn by **frozen vs live**:

- A **frozen** value/identity from an external source (published number, closed
  form, theorem identity, a hand-derived bank differential frozen as a golden) →
  **`oracle_literature`**.
- Two computations run **in the same test** and asserted equal (two live code
  paths) → **`oracle_crossengine`**.
- The assertion **is** an internal axiom the object must satisfy by construction
  (no external number, no second engine) → **`oracle_selfcert`**.

A test can be several at once (bank byte-level oracle: the bank closed form is
`literature`, and it runs the bank code live against the CS engine, so also
`crossengine`).

## Edge-case rulings (recorded so the sweep is reproducible)

1. **QPA tests stay `qpa`-only.** The `qpa` marker *is* the QPA oracle class. A
   `tests/qpa/` test is **not** additionally marked `oracle_crossengine` merely
   for comparing against QPA, nor `oracle_literature` for referencing a value —
   the QPA agreement is the operative oracle. `tests/qpa/` is untouched by the
   sweep. (The two explicitly-`fast` files there — `test_scripts.py`,
   `test_session_guard.py` — are pure string-builder/guard infrastructure and stay
   UNMARKED.)
2. **numba ≡ pure and sparse ≡ dense are `oracle_crossengine`.** They are two
   genuinely different implementations of the same kernel
   (`test_kernels.py`, `test_linalg_fast.py`, `test_cohomology_fast_rank.py`).
3. **CS ≡ bar / CS ≡ Bardzell / minimal ≡ bar / minimal-coh ≡ CS-coh →
   `oracle_crossengine`.** If the same file *also* freezes a closed form, it also
   carries `oracle_literature`.
4. **The read-only bank byte-level oracle** (`test_battery_bank_oracle.py`) is
   `oracle_literature` **+** `oracle_crossengine`: the bank differentials are
   closed forms (literature-class in this repo's taxonomy, per
   `docs/verification.md`) *and* the bank is a wholly separate implementation run
   live against CS. Bank **frozen** golden fixtures (the batch labdb re-ports) are
   `oracle_literature` only.
5. **Derived-invariance / theorem-identity oracles are `oracle_literature`, not
   `oracle_crossengine`.** When HH^• agrees across two quiver orientations because
   a *theorem* (derived invariance, Happel's trace identity, incidence ≅ nerve,
   symmetric ⇒ `HH^n = HH_n`) says it must, the theorem is the oracle. The
   "two computations agree" is the theorem's content, not two independent engines.
6. **`d∘d=0`, order condition, canonicalization, adversarial-solver,
   dimension/iso certificates, unit/Leibniz/module arbitration identities →
   `oracle_selfcert`.**
7. **"Same engine behind an indirection" is NOT `oracle_crossengine` — it is
   plumbing (UNMARKED).** CLI ≡ public-API parity, webapp runner byte-stable
   delegation, GUI runner parity: both sides call the *same* library, so this is a
   delegation/plumbing contract, not two independent mathematical engines. The
   task's own list files GUI/webapp/HPC as infrastructure.
8. **Foundational datatype layers are UNMARKED (contract).** `tests/fields/`
   (field arithmetic axioms, base-change smoke), `tests/core/` +
   `tests/combinat/` (structure-constant identities, path law, relation parsing)
   verify the datatypes' API contract, not a homological/representation-theoretic
   computation. Left unmarked.
9. **Gröbner admissibility / finiteness certificates are UNMARKED (structural
   gate).** `docs/verification.md` files admissibility under "Other structural
   gates," separate from the two oracle classes; `oracle_selfcert` is reserved for
   certificates *internal to a (co)homology computation*. `tests/groebner/` stays
   unmarked (engine-mechanics correctness + loud refusal).
10. **Worked-steps trace tests are UNMARKED (self-consistency plumbing).** They
    assert printed claims equal what the *same* engine computed (golden-file
    binding discipline) — infrastructure, not an external/second oracle.
11. **`tests/hochschild/` (the bar complex) is `oracle_literature`.** The bar
    tests pin *known closed forms* (dual numbers `HH_0=2, HH_n=1`, hereditary
    vanishing, symmetric `HH^n = HH_n`) — the bar engine reproducing a literature
    value. (The bar complex is the base oracle *for the deeper engines*; where a
    deeper engine is measured against it live, that deeper-engine test is
    `oracle_crossengine`.)
12. **Acceptance/exit-criteria files** take the marker(s) of the oracle they
    actually assert (usually `oracle_crossengine` — they run engines and compare —
    sometimes `+ oracle_selfcert` for structural certification). Pure
    docs/quickstart acceptance (`test_quickstart.py`) is UNMARKED.

## Audit mechanism decision

The audited table on `docs/verification.md` is gated by a new release test,
`tests/release/test_oracle_classes.py`, following the existing **badge==page**
doctrine (cf. `tests/release/test_markers.py`, `test_readme.py`). Mechanism:
**out-of-process `pytest --collect-only -m <expr>`** (subprocess, exactly as
`test_markers.py` already does for the buckets) for each of `oracle_literature`,
`oracle_crossengine`, `oracle_selfcert`, and `qpa`, plus the union
`oracle_literature or oracle_crossengine or oracle_selfcert or qpa`; the counts
are parsed out of the `docs/verification.md` table and asserted equal to the live
collection. Subprocess (not in-process nested `pytest.main`) is chosen because
`test_markers.py` already proves it robust in this repo, it inherits the installed
extras environment, and the whole test stays well under 30 s. Marked `deep` (like
`test_markers.py`) so it shells out once on the deep leg, not on all fast matrix
cells. The bucket-invariance numbers (804/1180/112) remain gated by the existing
`test_markers.py` partition test.

## Sweep mechanics

- Module-level `pytestmark = [pytest.mark.oracle_…, …]` where a whole file shares
  one oracle profile (most battery files). Per-test / per-class decorators only
  where a file genuinely splits classes across distinct tests.
- Existing `pytestmark` (bucket or skipif) is **preserved** — merged into a list,
  never replaced. **No `fast`/`deep`/`qpa`/`slow` marker is ever added by this
  plan** (buckets stay auto-assigned by `tests/conftest.py`), which is what keeps
  the buckets byte-identical.
- `tests/qpa/` and `tests/conftest.py` are not edited.

## Forward-going convention

New battery files declare their oracle class at module level
(`pytestmark = [pytest.mark.oracle_literature, …]`) as part of the plan that adds
them, and update the audited counts on `docs/verification.md` (the release gate
fails otherwise) — the same standing rule that already requires every plan to add
its oracles to the verification page.

## Result (delivered 2026-07-26)

Audited global counts (`pytest --collect-only -m <expr>`, badge==page):

| Marker expression | Tests |
|---|---:|
| `-m oracle_literature` | 670 |
| `-m oracle_crossengine` | 396 |
| `-m oracle_selfcert` | 604 |
| `-m qpa` | 112 |
| `-m "oracle_literature or oracle_crossengine or oracle_selfcert or qpa"` | 1227 |

Per-directory contributions reconcile exactly (e.g. literature = engine 233 +
resolutions_cs/hochschild 54 + modules 181 + invariants/families/batch 202 = 670).
`tests/modules/` carries no `oracle_crossengine` (its cross-checks are QPA, in the
`qpa` bucket, plus internal certificates) — expected.

**Bucket reconciliation (the one honest delta).** The marker sweep is
byte-identical: with the new gate file excluded, `fast`/`deep`/`qpa` collect at
**804 / 1180 / 112** exactly (verified `--ignore`-ing
`tests/release/test_oracle_classes.py`). The gate itself is deliverable #5 — 3 new
`deep` tests (like `test_markers.py`) — so the *post-addition* totals are
**804 / 1183 / 112**, suite **2099**. `docs/verification.md` (bucket table, intro
count, subsystem release row) and the README tests badge are bumped 2096 → 2099 /
1180 → 1183 accordingly (the standing-rule count bump; `test_readme.py` keeps
badge == page green). `test_markers.py` pins the *partition*, not absolute counts,
so it is unaffected.

## Acceptance (all met)

- `pyproject.toml` registers `oracle_literature` / `oracle_crossengine` /
  `oracle_selfcert`; `qpa` documented as the QPA oracle class.
- The sweep marks `tests/` completely: every test is either oracle-marked or
  deliberately contract/infra. **89 test files** carry an oracle marker
  (engine 33, resolutions_cs 17, hochschild 2, modules 13, invariants 9,
  families 13, batch 2); `tests/qpa/` and every contract/infra directory
  (fields, core, combinat, groebner, citations, viz, trace, gui, hpc, webapp,
  release, top-level) are untouched.
- `docs/verification.md` gains the "Oracle classes as runnable markers" section +
  the audited class × marker-expression × count table + the intro sentence.
- `tests/release/test_oracle_classes.py` gates page == live collection (3 tests,
  ~21 s, subprocess collection — well under 30 s; `deep`).
- `-m fast` green (801 passed / 3 skipped); the marker sweep leaves
  `fast`/`deep`/`qpa` byte-identical (804/1180/112); whole-repo `--collect-only`
  clean (2099).
