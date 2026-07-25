# Plan 24 — left modules alongside right, right as default (Tier 1b)

Date: 2026-07-25. Branch: `plan-24-left-modules` (from `plan-23-module-surface`,
NOT from `main` — this plan builds directly on the Plan-23 op+D engine). Backlog:
Tier 1b (`docs/plans/DEEPER-ENGINES-BACKLOG.md`, Marco 2026-07-25).

**Directive (Marco, 2026-07-25):** "We should allow both left and right modules,
right as default."

## The idea in one line

A left `A`-module **is** a right `A^op`-module. Plan 23 made `A^op` a first-class
`Algebra` (`modules/opposite.py`) precisely so this becomes a *presentation-only*
wrapper: the left side reuses **every** existing right-module algorithm, run over
the representation algebra `A^op`, with no duplicated mathematics.

## Representation

`Module` keeps storing its action as a **right module over a representation
algebra `self.algebra`**. A new attribute `self.side ∈ {"right","left"}` records
how the user reads that representation:

| user's view | `self.algebra` (representation) | `self.side` | `self.base_algebra` |
|---|---|---|---|
| right `A`-module | `A` | `"right"` (default) | `A` |
| left `A`-module | `A^op` | `"left"` | `A` |

`base_algebra` is a derived property: `self.algebra` when right, `self.algebra.opposite()`
when left (the algebra the user thinks of the module as living over). Because
`opposite()` is a cached involution, `base_algebra` is exact and free.

**Right stays byte-unchanged.** `side` defaults to `"right"`; for a right module
`base_algebra is self.algebra`, so `__repr__` (`f"{name}: {side} {base_algebra}
module, dim {dim}, dimvec {dv}"`) reproduces the old string character-for-character.
`side` is invisible to `dimension_vector`, `check_module`, Hom/Ext, resolutions —
all of which read only `self.algebra` and `self.action`. The entire committed
right-module + qpa suite passes untouched (that IS the right-as-default gate).

## Why the mathematics is not duplicated

Every algorithm consumes `(self.algebra, self.action)` and is blind to `side`.
A left `A`-module carries `self.algebra = A^op`, so:

- `radical/top/soc`, `Hom/End/Ext`, `is_isomorphic`, projective covers &
  resolutions, injective resolutions & injective dimension, `τ/τ⁻` — all compute
  the right-`A^op`-module invariant, which is by definition the left-`A`-module
  invariant. Not one line of homological code is forked.

`side` is threaded only so that *derived* modules of a left module are themselves
presented as left modules (an honest repr, never a right `A^op` surprise):
`submodule`/`quotient` inherit the side of their ambient module; `_direct_sum`
and `projective_cover` tag their output with the covered module's side; resolution
terms inherit through the syzygy chain. These are label-only edits; the numbers
are identical.

## The duality functor `D` and the API decision for `dualize`

Classical `D = Hom_k(−,k)` is contravariant and **exchanges the two sides over the
same algebra**: `D(right A-mod) = left A-mod` and `D(left A-mod) = right A-mod`.

At the representation level D is unchanged from Plan 23 (transpose every action
matrix, reverse every label, representation algebra `R → R^op`). The Plan-24 change
is a **one-line side flip**: `dualize(M).side = flip(M.side)`. This makes
`base_algebra` invariant under `D` (right `A` ↔ left `A`, same `A`) — the
mathematically honest classical statement.

**API decision.** `M.dualize()` now returns a module carrying the **OTHER side flag
over the SAME base algebra**. For a right `A`-module it is a **left `A`-module**
(previously it was surfaced as a right `A^op`-module). The stored representation
algebra (`A^op`) and every action matrix are **byte-identical to Plan 23**; only
the `side` tag and the repr change. Consequently every downstream user of D
(`injective.py`, `τ`, `τ⁻`) is numerically unchanged: those consumers read the
representation, and D flips a tag they do not consult. `τ = D∘Tr` and `τ⁻ = Tr∘D`
each flip the side twice, so `τ`/`τ⁻` **preserve** side (τ of a right module is a
right module; τ of a left module is a left module) — as they must.

`Tr` (`transpose_module`) is likewise a right-`A` → right-`A^op` functor at the
representation level and gets the same side flip: `Tr M` is presented as the left
`A`-module it classically is (`Hom_A(P_•,A)` lands in `A`-mod = left modules).

### Plan-23 backward compatibility (two tests change semantics, loudly)

Two Plan-23 tests expressed the identity `I_v = D(A e_v)` through the idiom
`A.opposite().projective(v).dualize()`, relying on the *old* presentation of D as
"right `A^op`-mod → right `A`-mod" so the result compared equal to the right
`A`-module `injective(A,v)`:

- `tests/modules/test_duality_tau.py::test_D_of_opposite_projective_is_injective`
- `tests/modules/test_module_iso.py::test_iso_certificate_over_gfp_square`

Under the honest new semantics `dualize` flips the side, so that idiom yields a
*left*-tagged module, and `is_isomorphic` now **refuses to compare across sides**
(a category error). Both tests are updated to the honest form
`A.projective(v, side="left").dualize()` — literally "`I_v = D` of the **left**
projective `A e_v`" — which produces a right `A`-module equal to `injective(A,v)`.
This is a strictly more faithful statement of the theorem, not a weakening. The
representation compared is identical; only the construction is spelled in the
side-aware surface. (Recorded here and in the commit message per the standing
"do not silently change committed expectations" rule.)

## The side-aware surface

Constructors gain `side="right"` (default), routing left construction through
`A^op` and re-tagging:

- `A.simple(v, side=…)`, `A.projective(v, side=…)`, `A.injective(v, side=…)`:
  `side="left"` returns `opposite_algebra(A).<builder>(v).with_side("left")`.
- `A.module(dimension_vector, arrow_action, side=…, name=…)` (**new**): a thin
  convenience over `Module.from_arrow_action`. Right builds over `A`; left builds
  over `A^op` (the user supplies the `A^op`-representation, i.e. the arrow matrices
  of the opposite quiver — the honest data of a left module = right `A^op`-module)
  and re-tags left. `Module.from_arrow_action` itself is unchanged (right only).
- `M.side` is exposed; `M.with_side(s)` returns a re-tagged twin (same
  representation, other categorical label) — the explicit "side translation"
  between a right `A^op`-module and a left `A`-module.

`is_isomorphic`, `A.hom`, `A.ext` **refuse loudly** (`QuiverlabError`, a
category-error message) when the two arguments are not the same side over the same
base algebra — comparing a left to a right module is a category error, not `False`.
The guard runs *before* the dim/dim-vector fast-paths, so a left `S_v` and a right
`S_v` sharing a dimension vector still refuse rather than falsely report iso.

## Edge cases (explicit tests)

- **Multi-vertex asymmetry (kA₂, `a:1→2`).** Right `P(1)=e_1A={e_1,a}` has dimvec
  `{1:1,2:1}`; left `P(1)=Ae_1={e_1}` has dimvec `{1:1,2:0}` — the honest smallest
  witness that sides are not conflated. Both sides' S/P/I dimension vectors pinned
  on kA₂ and on the commutative square (`tests/modules/test_left_modules.py`).
- **τ commutes with the side translation.** `A.simple(v,side="left").tau()` and
  `A.opposite().simple(v).tau()` have *identical* representation (same `A^op`,
  same action bytes) and differ only in the side tag — pinned degreewise.
- **`D` exchanges the translates.** `τ(DM) = D(τ⁻M)` (both equal `D∘Tr∘D M`) and
  `D(τM) = τ⁻(DM)` (both equal `Tr M`) — tested at dim-vector + `is_isomorphic`
  level on kA₂/kA₃ (AR duality of translates, ASS Ch. IV/VIII).
- **`injective_dimension` of a left module = `pd` of its `D`-dual right module.**
  Definitional identity, tested directly (the `D`-dual of a left `A`-module is a
  right `A`-module over the same `A`).
- **Local (single-vertex) algebras.** `k[x]/(x^n)` is self-opposite; left and
  right S/P/I dimension vectors and pd/inj.dim agree.
- **Iso across sides refuses loudly** (`pytest.raises(QuiverlabError)`).

## Oracles

- **Literature (no `[qpa]`):** left/right duality `D` and the AR translates are
  ASS (Assem–Simson–Skowroński, *Elements of the Representation Theory of
  Associative Algebras* Vol. 1, 2006 — bib key `ASS2006`, already in
  `references.bib`). Cited in test docstrings at the chapter granularity we can
  verify (Ch. II projectives/injectives `e_vA` vs `Ae_v`; Ch. III duality `D`;
  Ch. IV/VIII AR translates) — no invented theorem numbers.
- **QPA (`-m qpa`, `[qpa]` live):** QPA is right-module native; left-side
  quantities are crosschecked by **feeding QPA the opposite algebra**
  (`A.opposite()`), against which the left module's underlying right-`A^op`
  representation is the natural QPA input. `tests/qpa/test_left_modules_qpa.py`
  ties the left surface (`A.simple(v,side="left").tau()` etc.) to the
  QPA-validated `A.opposite()` right-module computation.
- **Right-suite-unchanged** is itself the primary gate.

## Acceptance

- New src: `Module.side`/`base_algebra`/`with_side`, side-aware `__repr__`;
  `dualize`/`transpose_module` side flip; side inheritance in `radtopsoc`,
  `resolution` (`submodule`/`quotient`/`_direct_sum`/`projective_cover`); the
  category guard in `hom.py` (reused by `ext`, `is_isomorphic`); `Algebra.simple/
  projective/injective/module` gain `side=`. No floats in `src/` (AST gate).
  Domain-generic. No new `quiverlab.engine.*` imports.
- Tests (deep bucket): `tests/modules/test_left_modules.py`; qpa bucket:
  `tests/qpa/test_left_modules_qpa.py`. Two Plan-23 tests updated (above).
- Docs: this plan doc; `docs/internals/10-modules.md` gains a "Left modules"
  section; `CLAUDE.md` status line; ROADMAP row; backlog checkbox.
- **Verification page:** `docs/verification.md` does NOT exist on this branch (it
  lives on the concurrent plan-22 branch). Plan 24's oracles (this section — the
  left-module subsystem row maps to `test_left_modules.py` +
  `test_left_modules_qpa.py`) MUST be added to that page **at merge time**.

## Deferred / out of scope

- No-code left/right picker in GUI/webapp — Tier 1b item 3 (next slice); the
  `side` field enters that schema then.
- General `A.module(side="left")` with fully custom matrices is supported, but the
  worked-example battery focuses on the S/P/I builders (where the asymmetry is
  sharpest and QPA-checkable).
