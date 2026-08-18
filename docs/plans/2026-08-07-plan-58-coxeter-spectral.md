# Plan 58: Coxeter Spectral Analysis (record R20) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the existing exact Coxeter surface (`coxeter_polynomial`,
`spectral_radius`, `mahler_measure`, `is_cyclotomic_product`) into a **certified
spectral report** on `Algebra`: an exact `ℤ[x]` cyclotomic factorization with `Φ_n`
labels, a quasi-unipotence / finite-order verdict, an exact count of roots outside
the unit circle, and the **spectral radius `ρ` and Mahler measure `M` as certified
algebraic numbers** (minimal polynomial over `ℚ` + rational isolating interval +
root index) — never a float — plus the per-class Lehmer-dichotomy note. Clickable
end-to-end (GUI canvas → block → report), all four locales (EN/ES/FR/ZH), both
runners byte-identical.

**Architecture — CRITICAL scoping facts (2026-08-07 audit).** Almost every
primitive already ships; **P58 adds NO new root-finding or spectral machinery**, it
adds a *certification + labelling + assembly* layer and the no-code surface:

- `A.coxeter_polynomial()` → exact `sympy.Poly` in `t`, `charpoly(-C^{-T}C)`
  (`src/quiverlab/invariants/cartan.py`). Loud on singular Cartan. **Untouched.**
- `A.coxeter_matrix()` → exact integer (or rational) matrix `-C^{-T}C`
  (`cartan.py`). Used by P58 only to decide **finite order** (`Φ^m = I`). **Untouched.**
- `invariants/spectral.py::spectral_radius / mahler_measure` — **already exact**,
  already SOUND for complex off-circle roots (the `z + 1/z` Dickson/Sturm branch:
  `real_roots` when they provably suffice, `all_roots` otherwise), already
  short-circuits cyclotomic input to the exact integer `1`, already returns `None`
  on degree `< 1`. **Their tests are pinned (`tests/invariants/test_spectral.py`)
  — do NOT change their behaviour;** P58 only *wraps* their output and *reuses*
  their private off-circle helper `_off_circle_roots` / `_noncyclotomic_part` for
  the exact count.
- `engine/coxeter_spectrum.py::is_cyclotomic_product` — **already exact** (factor
  over `ℚ`, match each factor against `Φ_n` for `n ≤ 2·deg²`). Ported read-only
  from hanlab (MIT, Marco Armenta). P58 **imports** it; it does **not** touch that
  module's operator-level `σ`/paracyclic-monodromy internals (`twisted_bimodule`,
  `column_degeneration`, … — the arXiv:2606.15595 "operator level" program). See
  Design decision **§D7**.

New library code: extend `invariants/spectral.py` with the certified-number +
labelling helpers, and one new assembler module `invariants/coxeter_spectral.py`
holding `coxeter_spectral(A)` + the shared block builder `coxeter_spectral_block(A)`.
One new scalar compute kind `coxeter_spectral` wired through the standard
seven-touchpoint pattern (spec dispatch, Pyodide runner twin, `gui.js` ×2, i18n
**all four locales EN/ES/FR/ZH**, ETA, `_snip`, `results_html`, golden). The
`derived_fingerprint` block is
**explicitly out of scope** (see §D5) so its goldens and cross-runner byte-identity
stay untouched.

**Tech Stack:** sympy exact only — `factor_list`, `cyclotomic_poly`, `totient`,
`minimal_polynomial`, `Poly.intervals()` (public root isolation, **rational**
endpoints), `real_roots`, matrix powers over `ℚ`. **No floats in `src/`** — the AST
gate `tests/test_no_floats.py` scans `src/` only; `sympy` `Integer`/`Rational`/
`CRootOf` and rational interval endpoints all pass (confirmed against the gate's
rules — no `float(`, no decimal literal, no `.evalf()` in `src/`; the bank-float
values survive only as **test** oracles, which are exempt).

## Global Constraints

- Python is always `.venv/bin/python`; tests via
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`.
- Bucketing (`tests/conftest.py` auto-assign by directory): `tests/invariants/` →
  **fast**; `tests/hpc/` → **fast**; `tests/qpa/` → **qpa**; `tests/webapp/` →
  unmarked (extras-gated dir, Plan-32 rule).
- Plan-32 oracle-class markers on every non-webapp test: `oracle_literature`
  (frozen theory/paper values), `oracle_crossengine` (two independent
  implementations agree), `oracle_selfcert` (internal certificates), `qpa`.
- The new compute kind acts on the **existing algebra block** ⇒ **schema stays v1**
  (Plan-26/38 precedent: a version bump is only for a new request *block*).
- **i18n is FOUR locales**: `webapp/server/i18n/{en,es,fr,zh}.json` all carry the
  same key set, gated by `tests/webapp/test_i18n.py::test_key_parity_with_english[fr]`
  and `[zh]`. Every new key MUST be added to all four (supply FR and ZH strings, not
  just EN/ES) or the parity gate fails.
- Every new citation key is BibTeX-verified before use (Plan-29 rule).
- All refusals loud (`QuiverlabError`); presentation-less algebras (no path basis)
  already refuse at `cartan_matrix` — that refusal propagates and is captured as a
  per-field `{"error": …}` in the block, never a crash, never a silent default.
- Conventional commits; green at every commit; branch `plan-58-coxeter-spectral`
  off `dev`. Do NOT commit unless asked.

---

## The record (verbatim — R20)

> **R20 — Coxeter spectral analysis.** [C-scout P8; keep-with-corrections] Object:
> cyclotomic test (exact ℤ[x] factorization + Φ_n recognition); Mahler measure and
> spectral radius as CERTIFIED ALGEBRAIC NUMBERS (minimal polynomial + rational
> isolating interval — never a float; sympy CRootOf/Sturm). The Lehmer dichotomy
> (M = 1 or ≥ μ₀ ≈ 1.17628) is proven ONLY for restricted classes (accessible
> algebras) — an open problem in general; report per-class only. Refs: de la Peña
> arXiv:1310.1910, arXiv:1310.1557 (both sole-author, verified); de la Peña–Takane
> Arch. Math. 55 (1990) 120–134. Oracles: kA₂ χ = Φ₃; 3-Kronecker χ = x²−7x+1,
> ρ = (7+3√5)/2; T_{2,3,7} = E₁₀ realizes Lehmer's degree-10 polynomial (verified).
> Size S–M.

Metaplan card: `docs/plans/2026-08-07-metaplan-v1.0.0.md` §5 P58.

---

## Reference re-verification (mandatory — findings)

Fetched/verified 2026-08-07. Record the outcome; `# PIN` anything untranscribable.

- **arXiv:1310.1910** — *"On the Mahler measure of the Coxeter polynomials of
  algebras"*, **José-Antonio de la Peña (sole author)**, Adv. Math. 2014. **VERIFIED.**
  Existing bib key **`dlPena2014mahler`** (`references.bib` — no new key needed).
  Central theorem transcribed verbatim from the abstract: *"for any accessible
  algebra `A` either `M(χ_A) = 1` or `M(χ_B) ≥ μ₀` for some convex subcategory `B`
  of `A`."* Lehmer (1933): the polynomial `T^10 + T^9 − T^7 − T^6 − T^5 − T^4 − T^3
  + T + 1` has Mahler measure `μ₀ = 1.176280…`. **Honest-scope note (drives §D4):**
  the dichotomy is *class-conditional* (accessible algebras) **and** phrased over
  convex subcategories, NOT a statement about `M(χ_A)` of an arbitrary `A`. The
  general Lehmer problem is open. quiverlab must NOT emit a verdict on it.
- **arXiv:1310.1557** — *"Algebras whose Coxeter polynomials are products of
  cyclotomic polynomials"*, **de la Peña (sole author)**, 2013. **VERIFIED.** Defines
  *algebras of cyclotomic type* (χ_A a product of cyclotomics). Results used by
  §D1: (i) fractional-Calabi–Yau ⇒ **periodic** Coxeter transformation ⇒ cyclotomic
  type; (ii) non-negative homological form ⇒ cyclotomic type. This is the source
  that separates **periodic (finite order)** from merely **cyclotomic type
  (quasi-unipotent)** — the exact distinction P58 reports. **New bib key required:
  `dlPena2013cyclotomic`** (BibTeX-verify at Task 1 Step 3 — arXiv preprint form; if
  a journal version is not verifiable, cite the eprint).
- **de la Peña–Takane, Arch. Math. 55 (1990) 120–134** — title corroborated by web
  search as *"Spectral properties of Coxeter transformations and applications"*
  (volume/pages/year match the record). Full text not fetched. Used by §D3 for the
  *reality of the dominant eigenvalue on the wild-hereditary locus* (`ρ` a real
  eigenvalue). **New bib key `dlPenaTakane1990spectral`** — `# PIN`: transcribe the
  title/pages from the record + this search; **BibTeX-verify at Task 1 Step 3**, and
  if the exact title cannot be confirmed, cite by author/journal/volume/pages only.
- **Lehmer / E₁₀ realization** — cross-checked against the codebase's own pin: the
  wild hereditary star `[2,3,7]` = `T_{2,3,7}` = **E₁₀** (10 vertices; the wild
  extension of `E₈ = T_{2,3,5}`, past the affine `Ẽ₈ = T_{2,3,6}`) has Coxeter
  polynomial exactly Lehmer's polynomial `t^10 + t^9 − t^7 − t^6 − t^5 − t^4 − t^3 +
  t + 1` (`tests/invariants/test_coxeter_literature.py::test_lehmer_polynomial`,
  already green). ρ = M = μ₀ (the smallest known Salem number). **VERIFIED.**
- **QPA surface** — live `NamesGVars()` probe (2026-08-07): QPA **has**
  `CoxeterMatrix` / `CoxeterPolynomial` as attributes on quiver algebras, and
  `CartanMatrix`. Probe on `kA₂`: `CoxeterPolynomial(A) = x_1^2 + x_1 + 1` — matches
  our convention **exactly** (variable rename `x_1 → t`; the matrix representative
  differs by basis convention but the polynomial agrees). QPA has **no** Mahler /
  spectral-radius / cyclotomic-recognition surface (`Mahler → []`, `Spectral → []`;
  `Cyclotomic*` are GAP's polynomial *constructors*, not algebra predicates). Drives
  Task 4: crosscheck the *polynomial* against QPA; honest scope for the spectral
  analysis. **VERIFIED.**

---

## Design decisions settled

### D1 — Cyclotomic test, Φ_n recognition, and periodicity link

- **Verdict:** reuse `is_cyclotomic_product(χ)` (single source of truth) for the
  boolean `cyclotomic`. All factors cyclotomic ⟺ every eigenvalue of `Φ` is a root
  of unity ⟺ `ρ = 1` ⟺ `M = 1`.
- **Φ_n labelling (new):** `cyclotomic_factorization(χ)` runs `sp.factor_list(χ)`
  over `ℤ`; for each monic irreducible factor `f` of degree `d`, search
  `n ∈ {1, …, 2d²}` for `totient(n) == d` **and** `sympy.Poly(f) == cyclotomic_poly(n)`
  (normalise sign; a factor equal to `t` — root 0 — is non-cyclotomic). The bound
  `n ≤ 2d²` is exact: `φ(n) ≥ √(n/2)` for all `n ≥ 1`, so `φ(n) = d ⇒ n ≤ 2d²`
  (consistent with, and slightly tighter than, the range `is_cyclotomic_product`
  already searches — `range(1, 2·d·d + 3)`, i.e. `n ≤ 2d² + 2`; both are correct).
  Each returned entry is `{factor, latex, multiplicity, cyclotomic_index: n|None}`.
  The search starts at `n = 1`, so `Φ_1 = t − 1` and `Φ_2 = t + 1` are labelled like
  any other index; the bare factor `t` (root 0) is not a root of unity and is
  labelled `None` (the `is_cyclotomic_product` `t`-guard). **Repeated cyclotomic
  factors keep their multiplicity** — `factor_list` returns `(factor, mult)` pairs
  and `multiplicity` carries `mult` straight through, so a `D_4` Coxeter polynomial
  `(t+1)²(t²−t+1) = Φ_2²·Φ_6` reports `{Φ_2, mult 2}, {Φ_6, mult 1}` and the affine
  `(t−1)²·Φ_3` reports `{Φ_1, mult 2}, {Φ_3, mult 1}`.
- **Periodicity link — report BOTH, precisely:**
  - `quasi_unipotent` = `cyclotomic` (all eigenvalues are roots of unity; equivalently
    some power of `Φ` is unipotent). This is the eigenvalue-level statement and is
    n&s from the factorization — no diagonalizability needed.
  - `coxeter_order` = the finite multiplicative order `m` with `Φ^m = I`, **or `None`**.
    Finite order additionally requires `Φ` diagonalizable. Computed **exactly and
    cheaply** without `is_diagonalizable`: if not `quasi_unipotent`, `order = None`
    (reason "not quasi-unipotent ⇒ ρ > 1 ⇒ infinite order"); else set
    `m = lcm(recognised n's)` (each `Φ_n` factor contributes a primitive `n`-th root
    whose order is exactly `n`, so the lcm is the *only candidate* minimal period)
    and **verify** `A.coxeter_matrix()^m == I` by exact matrix power. If it holds,
    `order = m` (minimal by construction); if not, `order = None` with reason
    "quasi-unipotent but the Coxeter transformation has a nontrivial Jordan block ⇒
    infinite order (typical of tame/affine type)".
  - **Worked check of the design (both regimes):** self-injective/symmetric `Φ = −I`
    ⇒ χ = `(t+1)^v = Φ₂^v`, `m = lcm{2} = 2`, `Φ² = I` ⇒ order 2. `kA₂` ⇒ χ = Φ₃,
    order 3. **2-Kronecker (affine `Ã₁`)** ⇒ χ = `(t−1)² = Φ₁²`, all eigenvalues 1 so
    `quasi_unipotent = True`, but `Φ = [[−1,−2],[2,3]]`, `Φ − I` has rank 1 ⇒ Jordan
    block ⇒ `Φ^1 ≠ I` ⇒ `order = None` (infinite). This is exactly de la Peña
    1310.1557's periodic ⊊ cyclotomic-type separation.

### D2 — ρ and M as certified algebraic numbers (never a float)

The record's representation: minimal polynomial over `ℚ` (as primitive `ℤ[x]`),
rational isolating interval, Sturm-certified unique root, plus a "which root" index.

**`certify_real_algebraic(alpha) -> dict`** (new, in `spectral.py`). `alpha` is an
exact non-negative real sympy value (what `spectral_radius`/`mahler_measure`
return: `Integer`, `Rational`, a radical `Add`/`Pow`, a `CRootOf`, or an `Abs` of a
complex `CRootOf`). Validated end-to-end (see the sanity runs below):

1. **rational** (`alpha.is_Rational`, includes cyclotomic ⇒ 1): `minpoly = [q, −p]`
   for `p/q`, `degree = 1`, `interval = (alpha, alpha)`, `root_index = 0`,
   `is_rational = True`.
2. **algebraic** (everything else): `m = sp.minimal_polynomial(alpha, x, polys=True)`
   (primitive `ℤ[x]`); `ivs = sp.Poly(m, x).intervals()` (ascending, **disjoint,
   isolating**, **rational** endpoints); locate the *unique* `k` with `alpha` in the
   `k`-th interval by the exact signs `(alpha − a_k).is_nonnegative` and
   `(b_k − alpha).is_nonnegative`. `interval = (a_k, b_k)`, `root_index = k`,
   `minpoly = m.all_coeffs()`, `degree = m.degree()`.
   **Refuse loudly** (`QuiverlabError`) if `minimal_polynomial` raises or the
   location is not unique (0 or >1 hits) — never guess a root.

**Caller contract — `certify_real_algebraic` is ONLY ever called on a
real-dominant `alpha`.** `minimal_polynomial` is intractable on the
complex-modulus form `sqrt(CRootOf·CRootOf)` (measured: **120.78 s** on the §D2
table's dim-5 example, `alpha.is_real == None` — it *succeeds slowly*, so an
exception-based guard would NEVER fire). Therefore the assembler (§D3) gates on the
shipped `_real_roots_suffice` predicate **before** any `minimal_polynomial` call and
never hands a complex-modulus value to this function. `certify_real_algebraic`
itself assumes its input is a real radical/`CRootOf`/rational; its only loud refusal
is the non-unique-location guard above.

Serialisation in the block: `interval` as `[str(a), str(b)]` (exact rational
strings like `"6"`, `"7/2"`), `minpoly` as `list[int]` (high→low degree), `latex`
via `sp.latex(alpha)` (radical for low degree, `CRootOf` notation otherwise),
`value` as `str(alpha)`. **No decimal ever** — the rational interval *is* the
numeric localisation.

**Validated (2026-08-07):**

| input χ | ρ = M | minpoly (ℤ[x]) | interval | root_index |
|---|---|---|---|---|
| `t²−7t+1` (3-Kronecker) | `(7+3√5)/2` | `[1,−7,1]` | `(6,7)` | 1 |
| Lehmer (E₁₀) | Salem `μ₀` | Lehmer `[1,1,0,−1,−1,−1,−1,−1,0,1,1]` | `(1,2)` | 1 |

Both real-dominant rows certify in **sub-second** time. The **complex-dominant**
case `(t+1)(t⁴−7t³+16t²−7t+1)` (a non-hereditary dim-5 `rad²=0` algebra) is
**NOT** in this table by design: there `ρ = |complex root|` returns as
`sqrt(CRootOf·CRootOf)` with `is_real == None`, and `minimal_polynomial` on it takes
**120.78 s** (measured 2026-08-07) — it does not raise, it hangs. §D3 refuses ρ/M
for that regime *deterministically* (never calling `minimal_polynomial`), while
still reporting the outside-circle count.

### D3 — Spectral radius: a deterministic real-dominant gate (never a hang)

- `ρ` and `M` are, by definition, **non-negative reals** and always **exactly
  computed** by the shipped `spectral_radius`/`mahler_measure` (sound for complex
  off-circle roots). P58 never re-derives them; it certifies them via §D2.
- **The gate (the fix for the measured 121 s hang).** An exception-based
  "try `minimal_polynomial`, catch on failure" guard is WRONG: on a complex-modulus
  `alpha` the call does not raise, it runs for ~121 s and then succeeds — so nothing
  catches it, and the offline desktop tier (no wall-clock kill) hangs indefinitely
  on a *small* quiver. Sizing cannot save us either: the trigger is
  complex-dominance, not dimension (the example is dim 5, which routes to the instant
  tier). Instead, `coxeter_spectral_block` **decides real-dominance up front,
  deterministically and cheaply, with a shipped predicate**, and only then calls
  `certify_real_algebraic`:

  > **Complex-dominant gate.** Let `q = _noncyclotomic_part(χ)`. If `q is None`
  > (all-cyclotomic short-circuit, §D4) then `ρ = M = 1` (rational, certified
  > trivially). Otherwise evaluate the shipped predicate `_real_roots_suffice(q)`
  > (`spectral.py`; measured **0.0004 s**, and it does **no** complex-root
  > isolation — it counts distinct real roots of the `z + 1/z` reduction against its
  > squarefree degree). If `True` (**real-dominant**: every off-circle root is real,
  > so `ρ` is a real radical/`CRootOf` and `minimal_polynomial` is sub-second — all
  > literature oracles: 3-Kronecker, Lehmer, m-Kronecker, Salem stars), certify `ρ`
  > and `M` via §D2. If `False` (**complex-dominant**: `ρ = |complex root|`), do
  > **NOT** call `minimal_polynomial` — set both `spectral_radius` and
  > `mahler_measure` to `{"error": "complex-dominant spectrum: ρ = |complex root|; "
  > "certified-algebraic-number extraction is out of scope; outside-circle count "
  > "reported"}` and continue. The `off_circle_root_count` is reported in **either**
  > branch.

  This is deterministic, size-independent, reproducible (no wall-clock, no
  threads/timeouts), and reuses machinery already in the module — nothing new to
  implement for the decision.
- **`ρ` and `M` share one cost path.** `M = |lc|·∏|off-circle roots|`, so the
  off-circle roots being all-real (or not) governs both. They are certified together
  when the gate is `True` and refused together (behind the *same* gate) when it is
  `False` — never one without the other.
- Real-dominant reality on the wild-hereditary / tree-bipartite locus (`ρ` a simple
  real eigenvalue `> 1`) is the classical A'Campo / de la Peña–Takane 1990 fact; the
  gate is exactly the shipped, sound decision of that reality on the actual `χ`.
- We do **not** invent a Möbius/Routh root counter: the shipped `_off_circle_roots`
  (the `z + 1/z` Dickson + Sturm branch) already counts outside-circle roots
  soundly for complex roots — §D4 reuses it. Respecting the existing engine beats
  reimplementing it.

### D4 — Outside-unit-circle root count (always reported, exact)

`off_circle_root_count(χ)` (new, in `spectral.py`) = `0` on cyclotomic input (all
roots on the circle), else `len(_off_circle_roots(_noncyclotomic_part(χ)))` — the
shipped, sound, exact helper (counts complex off-circle roots too).
**`_noncyclotomic_part(χ)` contract:** it returns the product (with multiplicity) of
the non-cyclotomic irreducible factors of `χ`, **or `None`** when `χ` is a product
of cyclotomics (the all-cyclotomic short-circuit). `off_circle_root_count` and the
§D3 gate both special-case that `None` first (count `0`; `ρ = M = 1`). Coxeter
polynomials are (anti)reciprocal, so outside-count = inside-count and
on-circle = deg − 2·outside; the block reports the outside count and notes the
reciprocal symmetry as prose (never relied on for correctness — the helper is
sound on non-reciprocal input via its `all_roots` fallback). Validated: 3-Kronecker
→ 1, Lehmer → 1 (Salem), complex-dominant example → 2 (a conjugate pair).

### D5 — API / payload (the `coxeter_spectral` block; schema v1)

New scalar compute kind `coxeter_spectral` on the algebra block. Shared builder
`coxeter_spectral_block(A)` drives **both** runners byte-identically. Shape:

```python
{
  "kind": "coxeter_spectral",
  "coxeter_polynomial": "<str expr in t>",     # str(χ)
  "coxeter_latex": "<sympy.latex(χ)>",
  "cyclotomic": true|false,                     # is_cyclotomic_product(χ)
  "factorization": [                            # per irreducible ℤ[x] factor
    {"factor": "t^2 + t + 1", "latex": "...", "multiplicity": 1, "cyclotomic_index": 3},
    {"factor": "t^2 - 7*t + 1", "latex": "...", "multiplicity": 1, "cyclotomic_index": null}
  ],
  "quasi_unipotent": true|false,                # == cyclotomic
  "coxeter_order": 3|null,                      # finite order Φ^m=I, or null
  "coxeter_order_reason": "<str>",              # "Φ^3 = I" / Jordan-block / ρ>1
  "outside_unit_circle_count": 1,
  "spectral_radius": <certified-dict> | {"error": "<loud reason>"},
  "mahler_measure":  <certified-dict> | {"error": "<loud reason>"},
  "lehmer_class_note": "<fixed documentation string; i18n>",
  "latex": "<\\begin{aligned}…\\end{aligned} summary>",
  "scope": "<honest-scope string>",
  "references": ["dlPena2014mahler", "dlPena2013cyclotomic",
                 "dlPenaTakane1990spectral", "lenzing_delapena_spectral", "assem_book"],
  # "citations": added by the dispatch caller from "references" (the standard pattern)
}
```

`<certified-dict>` = `{"minpoly": [ints], "degree": int, "interval": ["a","b"],
"root_index": int, "is_rational": bool, "latex": "...", "value": "..."}`. Any field
that raises on this input (singular/absent Cartan → `coxeter_polynomial`;
uncertifiable ρ/M → §D3) is captured as `{"error": <the loud message>}` per field
(Plan-30 τ-block / P43 fingerprint precedent), never a crash.

`lehmer_class_note` is a **fixed documentation string**, not a computed verdict:
*"de la Peña (arXiv:1310.1910): for an accessible algebra, M(χ) is either 1 or ≥ μ₀,
where μ₀ = 1.17628… is the Mahler measure of Lehmer's polynomial (minimal polynomial
t^10 + t^9 − t^7 − t^6 − t^5 − t^4 − t^3 + t + 1). quiverlab does not decide
accessibility and does not rule on the open general Lehmer problem; this documents
the class-conditional dichotomy only."* — μ₀ referenced by its minimal polynomial,
never a float.

`scope`: *"exact spectral analysis of the Coxeter polynomial χ = charpoly(−C^{-T}C);
ρ and M are certified algebraic numbers (minimal polynomial + rational isolating
interval). The finite/tame/wild reading is a theorem on the hereditary locus (see
form_type)."*

**Deviation from the metaplan card (stated honestly).** The metaplan card §5 P58
ends: *"GUI: extend `derived_fingerprint`/Coxeter block."* **This plan deviates**:
it adds a **new** `coxeter_spectral` compute kind and leaves the existing
`derived_fingerprint` block **byte-identical** rather than mutating it. The reasons
are concrete, not cosmetic:

- **Golden byte-stability.** `derived_fingerprint`/`derived_compare` ship frozen
  runner goldens (`tests/webapp/_runner_goldens.json`) and a cross-runner
  byte-identity contract (`test_runner_delegation.py`). Adding spectral fields would
  re-freeze those goldens and perturb every existing fingerprint/compare consumer
  for no new capability the standalone kind doesn't already deliver.
- **Block size / separation of concerns.** The spectral report (factorization with
  `Φ_n` labels, certified ρ/M, order, Lehmer note) is substantial; folding it into
  the derived-equivalence fingerprint conflates two distinct user intents (a
  necessary-condition fingerprint vs. a spectral analysis) and bloats a block whose
  point is a compact invariant tuple.
- **Existing consumers.** The fingerprint block already carries `coxeter_polynomial`
  as one of its derived invariants; the *spectral analysis of* that polynomial is a
  natural separate clickable, not a mutation of the fingerprint's contract.

The orchestrator has accepted this deviation (a new kind is cleaner than mutating
the byte-frozen fingerprint block and keeps existing goldens untouched). The record
R20 itself does not prescribe the GUI shape; the GUI clause lives only in the
metaplan card, quoted above.

### D6 — GUI + report

`gui.js` ×2 (docs + webapp, byte-identical twins): a `qlgui-coxeter_spectral`
checkbox in the invariants row; `S.ids` entry; push-list scalar entry; a block
renderer that shows the factorization with `Φ_n` labels, the cyclotomic /
quasi-unipotent / order line, the outside-count, and ρ/M as **minpoly + isolating
interval** (MathJax) with the `latex` exact form; the Lehmer-class note as an honest
footnote. `results_html.py` gets `_coxeter_spectral_html(b)` mirroring
`_recognizers_html`/`_derived_fingerprint_html`. i18n keys in **all four locales**
(`webapp/server/i18n/{en,es,fr,zh}.json`). ETA `scalars` entry `"coxeter_spectral":
0.5`: the §D3 gate makes cost *bimodal and deterministic* — real-dominant certifies
sub-second (the `0.5` estimate is comfortable), and complex-dominant is **refused
without computing** (never the 121 s path), so the estimate is honest in both
branches and needs no size-based routing hack. `_snip` recipe. One golden fixture.

### D7 — Relationship to `engine/coxeter_spectrum.py`

`engine/coxeter_spectrum.py` is the **operator-level** hanlab port (σ / paracyclic
monodromy / Serre-twist on `HH_•` — Marco's arXiv:2606.15595 "Coxeter spectra at
the operator level" program). P58 is the **K0-level** (classical) spectral analysis
of `χ_A = charpoly(−C^{-T}C)`, and is **independent** of that operator-level
machinery. The only shared symbol is `is_cyclotomic_product`, which P58 **imports**
(already public, already reused by `spectral.py`). P58 modifies **none** of that
module's internals (read-only bank attribution preserved).

---

### Task 1: certified-algebraic + labelling + count layer in `spectral.py`

**Files:**
- Modify: `src/quiverlab/invariants/spectral.py` (append helpers; existing
  `spectral_radius`/`mahler_measure`/private helpers UNCHANGED)
- Modify: `src/quiverlab/citations/references.bib` (add `dlPena2013cyclotomic`,
  `dlPenaTakane1990spectral` — BibTeX-verified)
- Test: `tests/invariants/test_coxeter_certified.py`

**Interfaces (new public functions in `spectral.py`):**
```python
def certify_real_algebraic(alpha) -> dict
    # exact non-negative real sympy value -> {minpoly, degree, interval,
    # root_index, is_rational, latex, value}; QuiverlabError on uncertifiable /
    # non-unique location (§D2). Never a float.
def cyclotomic_factorization(poly) -> list[dict]
    # [{factor, latex, multiplicity, cyclotomic_index: int|None}, ...] over ℤ (§D1).
def off_circle_root_count(poly) -> int
    # exact count of roots with |z| > 1 (0 on cyclotomic), via the shipped
    # sound _off_circle_roots/_noncyclotomic_part helpers (§D4).
```

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_coxeter_certified.py
"""Certified algebraic ρ/M, Φ_n labelling, off-circle count — self-certifying
(minpoly annihilates, interval brackets, count sound) + the record's exact
literature values."""
import sympy as sp
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.spectral import (
    certify_real_algebraic, cyclotomic_factorization, off_circle_root_count,
    spectral_radius, mahler_measure)

t, x = sp.Symbol("t"), sp.Symbol("x")
LEHMER = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1


@pytest.mark.oracle_selfcert
def test_certificate_is_self_consistent():
    for alpha in (spectral_radius(t**2 - 7*t + 1), spectral_radius(LEHMER)):
        c = certify_real_algebraic(alpha)
        m = sp.Poly(c["minpoly"], x)
        assert m.eval(sp.nsimplify(alpha)) == 0                    # minpoly annihilates
        a, b = sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])
        assert (alpha - a).is_nonnegative and (b - alpha).is_nonnegative  # brackets
        ivs = m.intervals()
        assert (a, b) == ivs[c["root_index"]][0]                   # index consistent
        assert not any(isinstance(z, float) for z in c["minpoly"])  # never a float


@pytest.mark.oracle_literature
def test_three_kronecker_certified_value():
    c = certify_real_algebraic(spectral_radius(t**2 - 7*t + 1))
    assert c["minpoly"] == [1, -7, 1] and c["degree"] == 2
    assert (sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])) == (6, 7)
    assert sp.simplify(sp.sympify(c["value"]) - (7 + 3*sp.sqrt(5))/2) == 0


@pytest.mark.oracle_literature
def test_lehmer_minpoly_and_interval():
    c = certify_real_algebraic(spectral_radius(LEHMER))
    assert sp.Poly(c["minpoly"], x) == sp.Poly(LEHMER.subs(t, x), x)
    assert (sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])) == (1, 2)


@pytest.mark.oracle_literature
def test_cyclotomic_labelling():
    facs = cyclotomic_factorization(t**2 + t + 1)                  # Φ_3
    assert facs == [{"factor": "t**2 + t + 1", "latex": facs[0]["latex"],
                     "multiplicity": 1, "cyclotomic_index": 3}]
    d4 = cyclotomic_factorization((t + 1)**2 * (t**2 - t + 1))     # Φ_2^2 · Φ_6
    idx = sorted((f["cyclotomic_index"], f["multiplicity"]) for f in d4)
    assert idx == [(2, 2), (6, 1)]
    mixed = cyclotomic_factorization(LEHMER)                       # irreducible, none
    assert [f["cyclotomic_index"] for f in mixed] == [None]


@pytest.mark.oracle_selfcert
def test_rational_and_count():
    c = certify_real_algebraic(sp.Integer(1))
    assert c["is_rational"] and c["degree"] == 1 and c["minpoly"] == [1, -1]
    assert off_circle_root_count(t**3 - 1) == 0                    # cyclotomic
    assert off_circle_root_count(t**2 - 7*t + 1) == 1              # one Salem/Pisot root
    assert off_circle_root_count(LEHMER) == 1
```

- [ ] **Step 2: Run to verify failure** — `ImportError` on the three new names.

- [ ] **Step 3: Implement the three helpers** at the end of `spectral.py`,
  exactly per §D1/D2/D4. Reuse `is_cyclotomic_product` (import already present) and
  the module's own `_off_circle_roots`/`_noncyclotomic_part`. `_locate` asserts a
  unique interval hit and raises `QuiverlabError` otherwise. **Confirm the no-floats
  gate stays green** (`tests/test_no_floats.py`) — `minimal_polynomial`,
  `intervals()`, `cyclotomic_poly`, `totient` are all exact, no `float(`, no
  `.evalf()`, no decimal literal.

- [ ] **Step 4: Add the two citations** to `references.bib`
  (`dlPena2013cyclotomic` = arXiv:1310.1557, de la Peña, 2013; and
  `dlPenaTakane1990spectral` = de la Peña–Takane, Arch. Math. 55 (1990) 120–134).
  BibTeX-verify; run `... -m pytest tests/trace/test_references.py -q` to confirm the
  bib still parses.

- [ ] **Step 5: Run tests** — Expected: PASS. Commit.

```bash
git add src/quiverlab/invariants/spectral.py src/quiverlab/citations/references.bib \
        tests/invariants/test_coxeter_certified.py
git commit -m "feat(invariants): certified algebraic ρ/M + Φ_n labelling + off-circle count"
```

---

### Task 2: the assembler `coxeter_spectral.py` + `Algebra.coxeter_spectral()`

**Files:**
- Create: `src/quiverlab/invariants/coxeter_spectral.py`
- Modify: `src/quiverlab/invariants/__init__.py` (export `coxeter_spectral`,
  `coxeter_spectral_block`), `src/quiverlab/core/algebra.py`
  (`Algebra.coxeter_spectral()`, beside `coxeter_polynomial` ~line 577)
- Test: `tests/invariants/test_coxeter_spectral.py`

**Interfaces:**
```python
def coxeter_spectral(A) -> dict          # the certified spectral report dict (no
                                         # "kind"/"references"/"citations"); the raw
                                         # analysis. quasi_unipotent + coxeter_order
                                         # use A.coxeter_matrix() (§D1).
def coxeter_spectral_block(A) -> dict    # the runner block: coxeter_spectral(A) plus
                                         # {"kind","references","scope",
                                         #  "lehmer_class_note","latex"}; per-field
                                         # {"error": …} on any loud refusal.
```
`Algebra.coxeter_spectral()` lazy-imports and returns `coxeter_spectral(self)`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_coxeter_spectral.py
"""The assembled Coxeter spectral report: verdicts, order, cross-checked against
the shipped spectral_radius/mahler_measure."""
import sympy as sp
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.invariants.spectral import spectral_radius, mahler_measure


def _kron(m, field=QQ):
    return Quiver([1, 2], {f"a{i}": (1, 2) for i in range(m)}).algebra(field=field)


@pytest.mark.oracle_literature
def test_ka2_is_phi3_order_3():
    r = linear_path_algebra(2, field=QQ).coxeter_spectral()
    assert r["cyclotomic"] and r["quasi_unipotent"]
    assert [f["cyclotomic_index"] for f in r["factorization"]] == [3]
    assert r["coxeter_order"] == 3
    assert r["outside_unit_circle_count"] == 0
    assert r["spectral_radius"]["value"] == "1" and r["mahler_measure"]["value"] == "1"


@pytest.mark.oracle_literature
def test_two_kronecker_cyclotomic_type_but_infinite_order():
    r = _kron(2).coxeter_spectral()
    assert r["cyclotomic"] and r["quasi_unipotent"]               # (t-1)^2 = Φ_1^2
    assert r["coxeter_order"] is None                            # Jordan block (affine)
    assert "Jordan" in r["coxeter_order_reason"] or "infinite" in r["coxeter_order_reason"]


@pytest.mark.oracle_literature
def test_three_kronecker_wild_certified():
    r = _kron(3).coxeter_spectral()
    assert not r["cyclotomic"] and r["coxeter_order"] is None
    assert r["outside_unit_circle_count"] == 1
    assert r["spectral_radius"]["minpoly"] == [1, -7, 1]         # x^2 - (m^2-2)x + 1
    assert sp.simplify(sp.sympify(r["spectral_radius"]["value"]) - (7 + 3*sp.sqrt(5))/2) == 0


@pytest.mark.oracle_crossengine
def test_report_matches_shipped_primitives():
    for A in (linear_path_algebra(4, field=QQ), _kron(3),
              _kron(5, field=GF(32003))):
        chi = A.coxeter_polynomial().as_expr()
        r = A.coxeter_spectral()
        if "value" in r["spectral_radius"]:
            assert sp.simplify(sp.sympify(r["spectral_radius"]["value"])
                               - spectral_radius(chi)) == 0
            assert sp.simplify(sp.sympify(r["mahler_measure"]["value"])
                               - mahler_measure(chi)) == 0


@pytest.mark.oracle_selfcert
def test_singular_cartan_refuses_per_field():
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(5))  # k[x]/(x^2)
    r = A.coxeter_spectral()
    assert "error" in r["coxeter_polynomial"]                    # singular Cartan, captured
```

- [ ] **Step 2: Run to verify failure** — `AttributeError: coxeter_spectral`.

- [ ] **Step 3: Implement `coxeter_spectral.py`.** Compute `χ = A.coxeter_polynomial()`
  inside a `_field(...)`-style guard (P43 fingerprint pattern) so a singular/absent
  Cartan is captured per field, never raised. From `χ`: `cyclotomic`,
  `factorization`, `off_circle_root_count`, `certify_real_algebraic(spectral_radius(χ))`
  and `…(mahler_measure(χ))` (each guarded → `{"error": …}` on §D2 refusal).
  `quasi_unipotent = cyclotomic`; `coxeter_order` via §D1 (`lcm(n's)` + exact
  `A.coxeter_matrix()` power verify, itself guarded on singular Cartan). Add the
  `Algebra.coxeter_spectral()` method and the `invariants/__init__.py` exports.

- [ ] **Step 4: Run tests** — Expected: PASS. Commit.

```bash
git add src/quiverlab/invariants/coxeter_spectral.py src/quiverlab/invariants/__init__.py \
        src/quiverlab/core/algebra.py tests/invariants/test_coxeter_spectral.py
git commit -m "feat(invariants): coxeter_spectral report — cyclotomic verdict, order, certified ρ/M"
```

---

### Task 3: literature-oracle battery

**Files:**
- Test: `tests/invariants/test_coxeter_spectral_literature.py`

**Interfaces:** consumes `Algebra.coxeter_spectral()`, the `families/dynkin.py`
generators, `linear_path_algebra`, `TruncatedPathAlgebra`, and the `_star`/`_canonical`
builders (copy from `test_coxeter_literature.py`). All `oracle_literature`.

- [ ] **Step 1: Write the tests** — the record's oracle set plus the sweeps:
  - **kA₂ = Φ₃** (order 3), **kA₄ = Φ₅** (`v₅`, prime ⇒ single cyclotomic factor,
    order 5), **kA₅ = Φ₂Φ₃Φ₆** (`v₆`, three cyclotomic factors — labelling test).
  - **3-Kronecker**: χ = `x²−7x+1`, ρ = M = `(7+3√5)/2` (minpoly `[1,−7,1]`, interval
    `(6,7)`), outside count 1, non-cyclotomic, order `None`.
  - **m-Kronecker ladder** `m ∈ {1,2,3,4,5}`: m=1 → A₂ = Φ₃ (cyclotomic, order 3);
    m=2 → `(t−1)²` (cyclotomic type, order `None`, Jordan/affine); m≥3 → wild
    `x²−(m²−2)x+1`, one root outside, ρ = M = `(m²−2+√((m²−2)²−4))/2`. Reconciles the
    P33 pin (m=3 ⇒ `x²−7x+1`).
  - **T_{2,3,7} = E₁₀ = Lehmer** (`_star([1,2,6])`, 10 vertices): χ = Lehmer
    polynomial (transcribe `t^10+t^9−t^7−t^6−t^5−t^4−t^3+t+1`), non-cyclotomic,
    outside count 1 (Salem), ρ = M with `minpoly == Lehmer`, interval `(1,2)`.
  - **Dynkin/Euclidean cyclotomic sweep**: A₄/D₅/E₆/E₇/E₈ all `cyclotomic == True`,
    finite `coxeter_order` equal to the Coxeter number (E₈ ⇒ 30, E₆ ⇒ 12, E₇ ⇒ 18;
    cross-check against `test_coxeter_literature`'s `_coxeter_number`), ρ = M = 1;
    affine `Ã_{1,3}`/`D̃₄`/`Ẽ₆` all `cyclotomic == True` but `coxeter_order is None`
    (defect/Jordan block), ρ = M = 1.
  - **Lehmer-class note present** and mentions μ₀ / Lehmer minimal polynomial; assert
    it never contains a float and never claims a general verdict.

- [ ] **Step 2: Run** — Expected: PASS. If a value mismatches, fix the **quiver /
  transcription**, never the pin (P38 rule). Commit.

```bash
git add tests/invariants/test_coxeter_spectral_literature.py
git commit -m "test(invariants): Coxeter spectral literature battery (Φ_n, Lehmer/E₁₀, m-Kronecker, ADE sweep)"
```

---

### Task 4: QPA crosscheck (Coxeter polynomial) + honest scope

**Files:**
- Test: `tests/qpa/test_coxeter_spectral_qpa.py`

**Interfaces:** `quiverlab.qpa.session` (`should_skip_qpa`, `run`, `libgap_handle`),
`quiverlab.qpa.scripts.quiver_and_algebra_script`. Probe confirmed QPA binds
`CoxeterPolynomial(A)` on quiver algebras and it agrees with ours after `x_1 → t`.

- [ ] **Step 1: Write the tests** (marker `qpa`, skip guard `should_skip_qpa()`):
  - **`test_coxeter_polynomial_matches_qpa`**: for kA₂, kA₄, 3-Kronecker, D₄, E₆
    (over `QQ`), parse QPA's `CoxeterPolynomial(A)` (a GAP polynomial in `x_1`),
    substitute `x_1 → t`, and assert exact equality with `A.coxeter_polynomial()`.
    (This anchors the *foundation* of the whole spectral surface cross-engine.)
  - **`test_qpa_has_no_spectral_surface`** (honest-scope guard, Plan-35 pattern): scan
    `NamesGVars()` **before** any probe; assert there is **no** name matching
    `Mahler`/`Spectral`/`Lehmer` and no algebra-level *cyclotomic-type recognizer*
    (`Cyclotomic*` are GAP polynomial constructors, not predicates on `A`). The test
    **FAILS** (not skips) if QPA ever grows such a surface — the standing signal to
    add a real crosscheck. Documents that ρ/M/cyclotomic-verdict/Lehmer-class are
    covered by literature + cross-engine (bank floats) + self-cert oracles.

- [ ] **Step 2: Run live** (`... -m pytest tests/qpa/test_coxeter_spectral_qpa.py -v -m qpa`;
  venv has `[qpa]`). Expected: PASS. Commit.

```bash
git add tests/qpa/test_coxeter_spectral_qpa.py
git commit -m "test(qpa): CoxeterPolynomial crosscheck + honest-scope guard for the spectral surface"
```

---

### Task 5: the no-code surface (`coxeter_spectral` compute kind, seven touchpoints)

**Files:**
- Modify: `src/quiverlab/hpc/spec.py` (`_dispatch` new kind ~line 1520; `_snip`
  scalar map ~line 2289; `ETA_MODEL["scalars"]` ~line 1255)
- Modify: `docs/gui/runner.py` (twin handler ~line 903; ETA `scalars` ~line 1255)
- Modify: `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox ~line 107; `S.ids`
  ~line 224; push-list ~line 858; a `renderBlock` branch)
- Modify: `webapp/static/app.js` (webapp block renderer, if it has its own branch set)
- Modify: `webapp/server/i18n/{en,es,fr,zh}.json` (**all four** — the parity gate
  `test_i18n.py::test_key_parity_with_english[fr]/[zh]` fails on a two-language add)
- Modify: `src/quiverlab/trace/results_html.py` (`_coxeter_spectral_html` + kind map)
- Modify: `tests/webapp/_runner_goldens.json` + `test_runner_delegation.py` docstring
- Test: `tests/webapp/test_coxeter_spectral_p58.py` (unmarked), `tests/hpc/test_p58_kinds.py`

**Interfaces:** consumes `coxeter_spectral_block(A)`; both runners share it, so the
blocks are byte-identical. Dispatch shape copies the `recognizers`/`derived_fingerprint`
handlers (`block = coxeter_spectral_block(A); block["citations"] =
_citation_pairs(block["references"]); return block, None`).

- [ ] **Step 1: Write the failing cross-runner test** (copy the runner-pair fixture
  from `tests/webapp/test_module_blocks_m0729.py`):

```python
# tests/webapp/test_coxeter_spectral_p58.py  (unmarked — extras-gated dir)
def test_coxeter_spectral_block_shape(tmp_path):
    # request: 3-Kronecker over QQ, compute ["coxeter_spectral"]; assert
    #   block["cyclotomic"] is False
    #   block["outside_unit_circle_count"] == 1
    #   block["spectral_radius"]["minpoly"] == [1, -7, 1]
    #   "dlPena2014mahler" in block["references"]
    ...
def test_ka2_block_is_phi3(tmp_path):
    # kA₂: block["factorization"][0]["cyclotomic_index"] == 3; coxeter_order == 3
    ...
def test_twin_parity(tmp_path):
    # SAME request through hpc.spec and docs/gui/runner.py; assert
    # json.dumps(sort_keys=True) equality on the coxeter_spectral block.
    ...
```

- [ ] **Step 2: Implement the dispatch handler** in `spec.py`, then the twin in
  `runner.py` (both call `coxeter_spectral_block`), then the GUI wiring: checkbox
  `qlgui-coxeter_spectral` ("Coxeter spectral analysis") in both `gui.js`, `S.ids`
  entry, push-list scalar entry (join the `["cartan","coxeter_polynomial",…]` loop),
  `renderBlock` branch (factorization with `Φ_n`, cyclotomic/order/outside line,
  ρ/M as minpoly + interval, Lehmer footnote), ETA `"coxeter_spectral": 0.5`,
  `_snip` recipe (`"coxeter_spectral": lambda it: "A.coxeter_spectral()"`), and
  i18n keys in **all four locales `{en,es,fr,zh}.json`** (`inv.coxeter_spectral`,
  `pick.kind.coxeter_spectral`, `block.coxeter_spectral.title`, `.cyclotomic_yes/_no`,
  `.quasi_unipotent`, `.order_finite/_infinite`, `.outside_count`, `.rho`, `.mahler`,
  `.lehmer_note`, `.complex_dominant` (the §D3 refusal message), `.scope`) — supply
  FR and ZH strings alongside EN/ES; run `test_i18n.py` to confirm key parity across
  all four. `results_html.py`: `_coxeter_spectral_html(b)` + the kind→renderer map
  entry + `_KIND_TITLES["coxeter_spectral"]`.

- [ ] **Step 3: Add the golden fixture** `coxeter_spectral_3kronecker_qq` to
  `_runner_goldens.json`; note it in the `test_runner_delegation.py` docstring
  change-log. Run the delegation test BEFORE adding to confirm existing goldens are
  byte-identical.

- [ ] **Step 4: Run the gates**

```
... -m pytest tests/webapp/test_coxeter_spectral_p58.py tests/webapp/test_runner_delegation.py tests/hpc -q
... -m pytest tests/gui -q      # the Pyodide twin's own suite, if present
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): coxeter_spectral compute kind — certified spectral analysis clickable end-to-end"
```

---

### Task 6: verification page, README, suite gate

**Files:**
- Modify: `docs/verification.md`, `README.md`
- Test: existing release gates

- [ ] **Step 1:** Verification page: add P58 rows to the subsystem→oracles→tests
  table (`invariants/` gains certified spectral analysis with its oracle classes;
  the `qpa/` row count grows; add the **honest-scope** entry: QPA has no
  Mahler/spectral/cyclotomic-type surface — the spectral layer is covered by
  literature pins (dlP 2014/2013, Lehmer), cross-engine (bank floats + the shipped
  primitives), and self-cert (minpoly annihilates + Sturm interval + ρ = M when one
  root outside)). Recount the class table (`tests/release/test_oracle_classes.py`
  drives the numbers — collect, paste, re-run green). Add the arXiv:1310.1910 /
  1310.1557 / dlP–Takane 1990 pins to the Class-1 literature list. README: one line
  in the features list ("certified Coxeter spectral analysis: cyclotomic Φ_n test,
  spectral radius & Mahler measure as certified algebraic numbers, Lehmer-class
  note").
- [ ] **Step 2:** Full gate:
  `... -m pytest -q -m fast` green; `... -m pytest tests/qpa -q -m qpa` green;
  `... -m pytest tests/release -q` green; `... -m pytest tests/test_no_floats.py -q`
  green.
- [ ] **Step 3: Commit**

```bash
git add docs/verification.md README.md
git commit -m "docs(verification): P58 Coxeter-spectral oracle rows + recounted classes"
```

---

## Acceptance (P58 definition of done)

1. `Algebra.coxeter_spectral()` public, exact, loud on refusal: cyclotomic verdict
   + Φ_n-labelled `ℤ[x]` factorization, quasi-unipotent flag, finite
   `coxeter_order` (or `None` with an honest reason), exact outside-unit-circle
   count, and ρ + M as **certified algebraic numbers** (minpoly + rational isolating
   interval + root index) or per-field loud refusal. **No floats in `src/`** (gate
   green).
2. Record oracles green: kA₂ = Φ₃ (order 3); 3-Kronecker χ = x²−7x+1,
   ρ = (7+3√5)/2 (minpoly [1,−7,1], interval (6,7)); T_{2,3,7} = E₁₀ = Lehmer
   (minpoly = the degree-10 polynomial, interval (1,2)); plus the m-Kronecker ladder
   and the ADE/affine cyclotomic sweep. Two new citation keys BibTeX-verified.
3. QPA `CoxeterPolynomial` crosscheck green live; the honest-scope guard FAILS if
   QPA ever grows a spectral surface.
4. `coxeter_spectral` clickable end-to-end (GUI canvas → block → report) in **all
   four locales EN/ES/FR/ZH** (i18n key parity gate green), both runners
   byte-identical, golden added with a documented change-log entry, schema still v1;
   `derived_fingerprint` byte-identical (untouched).
5. Verification page recounted with the honest-scope entry; fast + qpa + release +
   no-floats suites green.

## Explicit non-goals

- No change to `spectral_radius`/`mahler_measure`/`is_cyclotomic_product` behaviour
  (their pins stay green); P58 only wraps and reuses them.
- No new root-finding / Möbius / Routh counter (§D3/D4 reuse the shipped sound
  helpers).
- No extension of `derived_fingerprint`/`derived_compare` (§D5).
- No verdict on the general (open) Lehmer problem; the Lehmer note is
  class-conditional documentation only (§D4).
- No touching `engine/coxeter_spectrum.py` internals (§D7).
- No decimal/float display anywhere in `src/` — the rational isolating interval is
  the numeric localisation; bank floats remain test-only oracles.

---

## Change log

- **2026-08-07** — plan drafted (P58 / R20).
- **2026-08-07 adversarial review: 1 blocking + 2 majors + 2 minors applied**
  (complex-dominant deterministic `_real_roots_suffice` gate replacing the
  non-firing exception guard — 121 s hang measured; i18n ×4; the
  `derived_fingerprint` deviation argued honestly against the metaplan card).
