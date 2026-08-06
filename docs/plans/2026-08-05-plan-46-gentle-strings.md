# Plan 46: C5 Gentle / String Subsystem — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the special-biserial / string / gentle world first-class and no-code.
Butler–Ringel **string and band module classification** (fully algorithmic — the
indecomposables of a string algebra are exactly the string modules plus the band
families); string-module **τ** by the hook/cohook combinatorics, **cross-engine
arbitrated** against the trusted Plan-41 AR engine; the **Avella-Alaminos–Geiss
(AG) derived invariant** for gentle algebras — computed, **derived-invariant, and
labelled honestly as NOT complete**; and a **Brauer-graph-algebra constructor**
(presented `kQ/I`, per-instance dimension-certified, its Brauer-star special case
byte-checked against the existing Nakayama family). Every constructed module
self-certifies (`check_module` + indecomposability spot-checks); every construction
is dimension-certified or refuses loudly; the completeness edge (bands ⇒
rep-infinite ⇒ census is a budget-capped sample, never a silent "complete list")
is an honest semi-decision contract.

**Architecture:** One new package `src/quiverlab/strings/` and one new family
module, each a thin combinatorial / exact-linear-algebra layer over primitives
that already exist (no new math engines):

- `src/quiverlab/strings/walks.py` — the `Letter` alphabet (a direct arrow or a
  formal inverse), reduced-walk validity for string algebras (no forbidden
  length-2 subpath via the reduction system, no immediate backtrack, endpoint
  bookkeeping), the string-algebra sign functions `sigma`/`epsilon` (built from
  `_len2_in_ideal` + `_in_degree`/`_out_degree` — the AAG-gap data that lives
  nowhere else), `enumerate_strings(A, max_length, budget)` with the honest
  classification contract, and `find_bands(A, max_length)`.
- `src/quiverlab/strings/modules.py` — `string_module(A, walk)` and
  `band_module(A, walk, eigenvalue, mult=1)`, materialised through
  `Module.from_arrow_action` (0/1 arrow matrices from the walk; a band's closing
  letter carries an invertible companion block over the field — loud if the
  eigenvalue is not in the field).
- `src/quiverlab/strings/ar_strings.py` — `string_tau(A, walk)` /
  `string_tau_minus(A, walk)` by the Butler–Ringel hook/cohook rule, **arbitrated
  cross-engine**: on every test case `string_module(string_tau(w))` must be
  `is_isomorphic` to `string_module(w).tau()`; a convention disagreement (hook vs
  cohook side) flips ONCE against that arbiter and is documented.
- `src/quiverlab/strings/ag.py` — the permitted-thread / forbidden-thread data for
  gentle algebras and `ag_invariant(A) -> AGInvariant` (a multiset of `(n, m)`
  pairs, Avella-Alaminos–Geiss 2008). Honest: docstring + verification page state
  it is a derived **invariant** but **not complete**.
- `src/quiverlab/strings/block.py` — `strings_block(A)`, the algebra-only
  compute-kind builder (recognizer verdicts + a string census up to a default
  length + band presence ⇒ a rep-type verdict for string algebras + the AG
  invariant when gentle), mirroring `invariants/recognizers.py::recognizers_block`.
- `src/quiverlab/families/brauer.py` — `BrauerGraphAlgebra(graph, multiplicities,
  field)`: a presented `kQ/I` from a ribbon graph, per-instance dimension-certified
  by `dim = sum_v m_v * val(v)^2`, with the Brauer-star special case pinned
  `≅ NakayamaAlgebra(n, mn+1, cyclic=True)` (byte-equal Cartan matrices).

The whole thing composes on the P38 recognizers (`is_string`/`is_gentle`/
`is_special_biserial`, `reduction_system_of`, `_len2_in_ideal`,
`_in_degree`/`_out_degree`), the module materialiser (`Module.from_arrow_action` +
`check_module`), the P30 `decompose`/`is_indecomposable`, the P23/24
`Module.tau()`/`.tau_minus()`, the P41 AR engine (`knit_ar_quiver`,
`modules/ar.py` — the string-τ / census cross-engine oracle), and the presented
backbone (`Quiver.algebra`, the `families/trivial_extension.py` length-lex kernel
idiom, `families/nakayama.py`'s loud-hint validation template).

**Tech Stack:** pure exact combinatorics + exact linear algebra over `Domain`
(`modules/linalg_mod`, `fields/linalg`); no sympy beyond what `decompose` already
uses. **No floats in `src/`** (AST-gated by `tests/test_no_floats.py`, which scans
`src/` — the new `strings/` package and `families/brauer.py` are algebra and must
be exact: arrow matrices are `0`/`1`/int-mod-`p`/`Fraction`, band eigenvalues are
sanctioned Domain scalars, never Python floats).

## Global Constraints

- Python is always `.venv/bin/python`; tests run
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q ...`.
- **P38 (forms/recognizers) and P41 (AR completion) are MERGED to `dev` and are
  hard prerequisites.** This plan consumes `invariants/recognizers.py`
  (`is_string`/`is_gentle`/`is_special_biserial`/`is_nakayama`, `_len2_in_ideal`,
  `_in_degree`/`_out_degree`, `recognizers_block`), `resolutions_cs/build.py::
  reduction_system_of`, and `modules/ar.py` (`knit_ar_quiver`, `ARQuiver`) exactly
  as their signatures read on `dev`. Branch `plan-46-gentle-strings` off `dev`
  **after P38 and P41 have merged** (both are present on `dev` at authoring time —
  verify with `python -c "import quiverlab.modules.ar, quiverlab.invariants.recognizers"`
  before starting). **No dependency on P39/P40/P42–P45** — this branch merges
  independently.
- **Test buckets are auto-assigned by directory** (`tests/conftest.py`):
  `tests/modules/` and `tests/families/` → **deep**; `tests/invariants/`,
  `tests/webapp/`, `tests/gui/`, `tests/hpc/` → **fast**; `tests/qpa/` → **qpa**.
  Per the metaplan surfaces decision, the string / band / string-τ batteries live
  in **`tests/modules/test_strings_*.py`** (deep — they materialise and decompose
  modules), the Brauer batteries in **`tests/families/test_brauer*.py`** (deep),
  and the AG-invariant battery in **`tests/invariants/test_ag*.py`** (fast — pure
  combinatorics on the reduction system, no module materialisation). Run new tests
  by path during development; finish each task with a `-m deep` / `-m fast` /
  `-m qpa` spot-run of the touched files.
- **The `decompose` char-caveat is load-bearing.** `decompose`/`is_indecomposable`
  refuse over GF(p) when `char <= dim M` (the trace-form radical is unreliable —
  `modules/decompose.py`). Every battery that decomposes a materialised string /
  band module, or spot-checks indecomposability, runs over **QQ or GF(32003)** so
  both the split search and the locality certificate decide. The pure combinatorial
  layers (`walks.py`, `ag.py`) are field-agnostic and run over any sanctioned field.
- **Loud refusals, never silent wrong answers** (the `GlobalDimension` /
  `TrivialExtension` precedent): a `band_module` whose eigenvalue is not in the
  field, a `string_module` off a walk that fails `check_module`, a
  `BrauerGraphAlgebra` whose dimension certificate fails, a census asked of a
  non-string algebra, or an AG invariant asked of a non-gentle algebra — each
  raises `QuiverlabError` with a `hint=`, never returns a plausible-looking wrong
  object.
- **Honest semi-decision contract (metaplan §6):** for a **string** algebra the
  Butler–Ringel classification is complete **iff there are no bands** (rep-finite);
  when bands exist the algebra is rep-infinite and `enumerate_strings` returns a
  **length-capped sample** with `status="budget"`, never a `status="complete"`
  list. State this on the verification page.
- House conventions: path composition is **left-to-right** (`a*b` = first `a` then
  `b`, `target(a)=source(b)`); modules are **right** modules by default
  (`Module.from_arrow_action` uses the anti-homomorphism `action[x*y] =
  action[y] @ action[x]`); `is_isomorphic`/`decompose` refuse loudly across
  sides/algebras. All refusals are `QuiverlabError`.
- **Convention flips are ARBITRATED, not assumed** (the house cup-sign /
  composition-order precedent): the string-τ hook-vs-cohook side, the AG `(n, m)`
  vs `(m, n)` ordering, and the Brauer arrow direction are each pinned by an
  oracle test (engine-τ isomorphism; AAG paper values; the Nakayama Cartan
  byte-equality). Write both candidates, let the oracle decide, record the outcome
  in the docstring.
- Plan-32 markers: `check_module`/partition/dimension certificates and the
  string-module self-identities = `oracle_selfcert`; string-τ ≡ engine-τ, census
  count ≡ `knit_ar_quiver` vertex count, Brauer-star ≡ Nakayama = `oracle_crossengine`;
  Butler–Ringel / AAG worked values, the `n(n+1)/2` interval count, the Brauer
  dimension formula, `is_symmetric(BGA)` = `oracle_literature`; QPA comparisons live
  in `tests/qpa/` (bucket = the class, never double-marked).
- **Mid-merge-train counts:** v0.2.0 lands ~15 subplans in overlapping waves, so
  the absolute suite counts drift between this plan's authoring and its merge.
  **Task 8 recounts the oracle-class table at merge time by running
  `tests/release/test_oracle_classes.py`** (paste the live numbers, never a
  guessed-at-authoring count) and claims only the deltas this plan adds.
- Every plan merge updates `docs/verification.md` (new oracle rows + recounted
  class table green, `tests/release/test_oracle_classes.py` green) and adds its
  citations to `citations/references.bib` + `registry.py`. Conventional commits;
  green tests at every commit.

---

### Task 1: the `Letter` alphabet, walk validity, `sigma`/`epsilon`, `enumerate_strings`, `find_bands`

The combinatorial core. A **string** (walk) is a reduced sequence of letters — each
a direct arrow or a formal inverse — that avoids the relations; `enumerate_strings`
lists them up to a length bound with the honest completeness contract, and
`find_bands` decides rep-finiteness.

**Files:**
- Create: `src/quiverlab/strings/__init__.py`, `src/quiverlab/strings/walks.py`
- Test: `tests/modules/test_strings_walks.py`

**Interfaces:**
- Consumes: `resolutions_cs/build.py::reduction_system_of(A) -> ReductionSystem`
  (`.rules` — each `ReductionRule` has `.lead` a tuple of arrow-name strings,
  `.source`, `.target`; `.leading_words()`; `.normal_form(word) -> dict`, the dict
  **empty iff `word` is in `I`**; `.domain`); `invariants/recognizers.py::
  _len2_in_ideal(rs, a, b)` (`a`, `b` arrow **name strings** with
  `target(a)=source(b)`; returns `True` iff the length-2 path `a*b` lies in `I`,
  i.e. `rs.normal_form((a, b))` is empty), `_in_degree(Q)`/`_out_degree(Q)` (each a
  `Quiver -> {vertex: int}`), `is_string(A)`; `combinat/quiver.py::Quiver`
  (`.vertices`, `.arrows` a `{name: (src, tgt)}` dict, `.source(name)`,
  `.target(name)`, `.word_source(word)`, `.word_target(word)`, `.compose_ok(word)`).
- Produces:
  ```python
  # A letter is (arrow_name: str, direction: +1|-1). +1 = direct arrow, -1 = formal
  # inverse. A walk is a tuple of letters. Left-to-right: letter i ends where
  # letter i+1 starts, in the WALK sense (inverses reverse endpoints).
  Letter = tuple            # (str, int)  -- int in (+1, -1)
  def letter_source(Q, ell) -> vertex      # source(name) if +1 else target(name)
  def letter_target(Q, ell) -> vertex      # target(name) if +1 else source(name)
  def invert(ell) -> Letter                # (name, -dir)

  def string_signs(A) -> tuple[dict, dict]
      # (sigma, epsilon): Q_1 -> {+1, -1} satisfying the Butler-Ringel string
      # conditions (S1) arrows with equal SOURCE get distinct sigma; (S2) arrows
      # with equal TARGET get distinct epsilon; (S3) for beta*gamma NOT in I with
      # target(beta)=source(gamma), sigma(gamma) = -epsilon(beta). Built greedily
      # from _in_degree/_out_degree (each <= 2 for a string algebra) + _len2_in_ideal.
      # Loud if A is not a string algebra, or the constraints are inconsistent.

  def is_valid_walk(A, walk, rs=None) -> bool
      # composable (letter_target(i) == letter_source(i+1)); reduced
      # (walk[i+1] != invert(walk[i])); no forbidden direct/inverse length-2
      # subword (a direct pair a,b with a*b in I, or an inverse pair a^-1,b^-1 with
      # b*a in I); and the sigma/epsilon side-consistency at every internal vertex.

  def enumerate_strings(A, max_length=8, budget=4096) -> StringCensus
      # StringCensus(walks=[canonical walk reprs], status="complete"|"budget",
      #   max_length, count, has_bands). Canonicalised up to inversion (w ~ w^-1 give
      # isomorphic modules -- keep the lexicographically-smaller). HONEST CONTRACT:
      # status == "complete" ONLY when find_bands(A, max_length) is empty AND the
      # DFS closed within budget below max_length (rep-finite: the list is ALL
      # indecomposable non-projective... = all string modules). Bands present, or a
      # budget/length cut => status == "budget" (a sample, never claimed complete).

  def find_bands(A, max_length=8) -> list        # canonical band walks (cyclic strings)
      # a band: a cyclic reduced walk of length >= 1, not a proper power, containing
      # BOTH a direct and an inverse letter (a pure directed/inverse cycle is not a
      # band), every rotation a valid walk. Non-empty => A is rep-INFINITE.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_strings_walks.py
"""Reduced walks, string signs, string enumeration, band detection (Plan 46 / C5).
Self-cert: every enumerated walk is valid and canonicalised; the sign functions
satisfy the Butler-Ringel conditions. Literature: kA_n has n(n+1)/2 strings (the
interval modules) and NO bands (rep-finite); the 2-cycle kQ/(ab,ba) HAS a band."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.fields import QQ
from quiverlab.strings.walks import (enumerate_strings, find_bands, invert,
                                     is_valid_walk, string_signs)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _kA3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=QQ)


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def _two_cycle():
    # kQ/(ab, ba): a: 1->2, b: 2->1, self-injective gentle, dim 4. HAS a band.
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


@selfcert
def test_signs_satisfy_string_conditions():
    A = _gentle_a3()
    sigma, epsilon = string_signs(A)
    assert set(sigma) == set(A.quiver.arrows) == set(epsilon)
    assert all(v in (1, -1) for v in sigma.values())
    assert all(v in (1, -1) for v in epsilon.values())
    # (S3): for a*b not in I (a: 1->2, b: 2->3, a*b NOT in I over kA3-no-rel), the
    # only composable nonzero pair -- sigma(b) == -epsilon(a).
    assert sigma["b"] == -epsilon["a"]


@selfcert
def test_enumerated_walks_are_valid_and_canonical():
    A = _gentle_a3()
    cen = enumerate_strings(A, max_length=6)
    for w in cen.walks:
        assert is_valid_walk(A, w)
        assert tuple(w) <= tuple(invert(e) for e in reversed(w))  # canonical rep


@lit
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_linear_a_n_string_count_is_interval_count(n, count):
    # kA_n indecomposables = the n(n+1)/2 interval modules, ALL string modules; no
    # bands (rep-finite). Butler-Ringel; the interval classification.
    arrows = {chr(ord("a") + i): (i + 1, i + 2) for i in range(n - 1)}
    A = Quiver(list(range(1, n + 1)), arrows).algebra(relations=[], field=QQ)
    cen = enumerate_strings(A, max_length=2 * n)
    assert cen.status == "complete"
    assert cen.count == count          # n(n+1)/2
    assert find_bands(A, max_length=2 * n) == []


@lit
def test_two_cycle_has_a_band_and_is_not_complete():
    A = _two_cycle()
    bands = find_bands(A, max_length=6)
    assert bands != []                                   # rep-infinite
    cen = enumerate_strings(A, max_length=6)
    assert cen.status == "budget" and cen.has_bands      # never claim "complete"


@selfcert
def test_non_string_algebra_refused():
    from quiverlab.errors import QuiverlabError
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (1, 4)}).algebra(
        relations=[], field=QQ)                          # 3 arrows out of 1: not SB
    with pytest.raises(QuiverlabError):
        string_signs(A)
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/modules/test_strings_walks.py -v`
Expected: FAIL — `ModuleNotFoundError: quiverlab.strings.walks`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/strings/walks.py
"""Reduced walks / strings / bands over a string algebra (Plan 46 / C5).

A string algebra is kQ/I with I monomial admissible, at most 2 arrows in and 2 out
at each vertex, and the gentle-style branching condition (for each arrow at most one
nonzero and one zero length-2 continuation on each side). Its indecomposable modules
are the STRING modules (from reduced walks avoiding I) and the BAND modules (from
cyclic reduced walks); Butler-Ringel 1987. Everything here is exact combinatorics on
the reduction system -- float-free, field-agnostic."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.recognizers import (_in_degree, _len2_in_ideal,
                                              _out_degree, is_string)
from quiverlab.resolutions_cs.build import reduction_system_of


def invert(ell):
    name, d = ell
    return (name, -d)


def letter_source(Q, ell):
    name, d = ell
    return Q.source(name) if d > 0 else Q.target(name)


def letter_target(Q, ell):
    name, d = ell
    return Q.target(name) if d > 0 else Q.source(name)


def _require_string(A):
    if not is_string(A):
        raise QuiverlabError(
            "strings: this algebra is not a string algebra",
            hint="strings/bands are defined for special-biserial monomial kQ/I "
                 "with <=2 arrows in/out per vertex (is_string(A) is False)")


def string_signs(A):
    """(sigma, epsilon): Q_1 -> {+1,-1}. (S1) equal-source arrows get distinct
    sigma; (S2) equal-target arrows get distinct epsilon; (S3) beta*gamma not in I
    => sigma(gamma) = -epsilon(beta). Greedy assignment; a string algebra always
    admits one (<=2 branches). Loud if inconsistent or A is not a string algebra."""
    _require_string(A)
    Q = A.quiver
    if Q is None:
        raise QuiverlabError("strings: need a quiver-presented algebra",
                             hint="structure-constant algebras carry no walks")
    rs = reduction_system_of(A)
    arrows = list(Q.arrows)
    sigma, epsilon = {}, {}
    # (S1)/(S2): at each vertex, the (<=2) out-arrows split +1/-1 for sigma; the
    # (<=2) in-arrows split +1/-1 for epsilon.
    for v in Q.vertices:
        outs = [a for a in arrows if Q.source(a) == v]
        for i, a in enumerate(outs):
            sigma[a] = 1 if i == 0 else -1
        ins = [a for a in arrows if Q.target(a) == v]
        for i, a in enumerate(ins):
            epsilon[a] = 1 if i == 0 else -1
    # (S3): reconcile via the nonzero composable pairs. Propagate; if a forced value
    # contradicts an assigned one, raise loudly (the string-condition inconsistency).
    for b in arrows:
        for g in arrows:
            if Q.target(b) == Q.source(g) and not _len2_in_ideal(rs, b, g):
                want = -epsilon[b]
                if sigma[g] != want:
                    # flip the sigma-class at source(g) once to honour (S3), then
                    # re-check (S1) at that vertex; if still inconsistent, raise.
                    _flip_sigma_class(Q, sigma, Q.source(g))
                    if sigma[g] != want:
                        raise QuiverlabError(
                            "strings: sigma/epsilon constraints are inconsistent",
                            hint=f"arrows {b!r}, {g!r} force a contradictory sign; "
                                 "the algebra may violate the string branch condition")
    return sigma, epsilon


def _flip_sigma_class(Q, sigma, v):
    for a in Q.arrows:
        if Q.source(a) == v:
            sigma[a] = -sigma[a]


def is_valid_walk(A, walk, rs=None):
    """Composable, reduced, relation-avoiding, sign-consistent."""
    Q = A.quiver
    if rs is None:
        rs = reduction_system_of(A)
    if not walk:
        return True                                   # trivial walk e_v is valid
    for i in range(len(walk) - 1):
        if letter_target(Q, walk[i]) != letter_source(Q, walk[i + 1]):
            return False                              # not composable
        if walk[i + 1] == invert(walk[i]):
            return False                              # not reduced (backtrack)
        if not _pair_ok(rs, walk[i], walk[i + 1]):
            return False                              # hits a relation
    return _signs_consistent(A, walk, rs)


def _pair_ok(rs, ell1, ell2):
    (n1, d1), (n2, d2) = ell1, ell2
    if d1 > 0 and d2 > 0:
        return not _len2_in_ideal(rs, n1, n2)         # direct a*b must be nonzero
    if d1 < 0 and d2 < 0:
        return not _len2_in_ideal(rs, n2, n1)         # inverse pair: reverse b*a
    return True                                       # a direct-then-inverse turn is
                                                      # governed by sign-consistency


def _signs_consistent(A, walk, rs):
    """A direct->inverse or inverse->direct turn at a peak/valley must use two
    DISTINCT arrows on the same side (sigma/epsilon distinct); this is exactly the
    Butler-Ringel condition that the string module is well-defined."""
    sigma, epsilon = string_signs(A)
    Q = A.quiver
    for i in range(len(walk) - 1):
        (n1, d1), (n2, d2) = walk[i], walk[i + 1]
        if d1 > 0 and d2 < 0:                          # peak: ... a  b^-1 ...
            if n1 == n2:
                return False
            if epsilon[n1] == epsilon[n2]:             # equal-target side clash
                return False
        elif d1 < 0 and d2 > 0:                        # valley: ... a^-1 b ...
            if n1 == n2:
                return False
            if sigma[n1] == sigma[n2]:
                return False
    return True


def _canonical(walk):
    rev = tuple(invert(e) for e in reversed(walk))
    return min(tuple(walk), rev)


@dataclass(frozen=True)
class StringCensus:
    walks: tuple
    status: str          # "complete" | "budget"
    max_length: int
    count: int
    has_bands: bool


def enumerate_strings(A, max_length=8, budget=4096):
    _require_string(A)
    rs = reduction_system_of(A)
    Q = A.quiver
    alphabet = [(a, +1) for a in Q.arrows] + [(a, -1) for a in Q.arrows]
    seen = set()
    frontier = [(v,) if False else () for v in ()]     # see note: seed per vertex
    # seed: the trivial walk at each vertex is a string module (the simple), plus
    # single letters; then extend. (Canonical seeding written in the DFS below.)
    out = set()
    truncated = [False]

    def extend(walk):
        cw = _canonical(walk)
        if cw in seen:
            return
        seen.add(cw)
        out.add(cw)
        if len(seen) > budget:
            truncated[0] = True
            return
        if len(walk) >= max_length:
            truncated[0] = True
            return
        end = letter_target(Q, walk[-1]) if walk else None
        for ell in alphabet:
            nxt = walk + (ell,)
            if is_valid_walk(A, nxt, rs):
                extend(nxt)

    for v in Q.vertices:
        out.add(())                                    # trivial walk (simple S_v)
        # (represent trivial walks per-vertex; see "Adjust to reality" for the exact
        #  bookkeeping that keeps one census entry per vertex.)
    for ell in alphabet:
        extend((ell,))
    bands = find_bands(A, max_length)
    complete = (not bands) and (not truncated[0])
    walks = tuple(sorted(out))
    return StringCensus(walks, "complete" if complete else "budget",
                        max_length, len(walks), bool(bands))


def find_bands(A, max_length=8):
    _require_string(A)
    rs = reduction_system_of(A)
    Q = A.quiver
    alphabet = [(a, +1) for a in Q.arrows] + [(a, -1) for a in Q.arrows]
    bands, seen = [], set()

    def closes(walk):
        return (letter_target(Q, walk[-1]) == letter_source(Q, walk[0])
                and is_valid_walk(A, walk + (walk[0],), rs)            # cyclic-reduced
                and walk[0] != invert(walk[-1]))

    def is_band(walk):
        dirs = {d for _, d in walk}
        if dirs != {+1, -1}:                            # pure directed/inverse cycle
            return False
        if _is_proper_power(walk):
            return False
        return all(is_valid_walk(A, walk[k:] + walk[:k], rs)          # every rotation
                   for k in range(len(walk)))

    def grow(walk):
        if len(walk) > max_length:
            return
        if len(walk) >= 1 and closes(walk) and is_band(walk):
            cw = _canonical_cyclic(walk)
            if cw not in seen:
                seen.add(cw)
                bands.append(cw)
        end = letter_target(Q, walk[-1])
        for ell in alphabet:
            if letter_source(Q, ell) == end and is_valid_walk(A, walk + (ell,), rs):
                grow(walk + (ell,))

    for ell in alphabet:
        grow((ell,))
    return bands
```

**Adjust to reality (Task 1):**
- `is_string`/`_len2_in_ideal`/`_in_degree`/`_out_degree` are P38 — **read their exact
  bodies first** (`invariants/recognizers.py`; `_len2_in_ideal(rs, a, b)` takes
  arrow-name strings and assumes `target(a)=source(b)`, so guard the call with the
  endpoint check). `reduction_system_of` lives in `resolutions_cs/build.py`.
- **The trivial-walk bookkeeping** (one census entry per vertex = the simple module)
  needs a concrete representation for a length-0 walk *at a specific vertex*; a bare
  `()` cannot name its vertex. Represent a trivial walk as `((None, v),)` or carry a
  parallel `{vertex}` set — pick one and keep `string_module` (Task 2) consuming the
  same shape. The count oracle `n(n+1)/2` **includes** the `n` simples, so the trivial
  walks must be counted.
- `_signs_consistent` recomputes `string_signs(A)` per call — memoise it on `A`
  (a module-level `functools.lru_cache` keyed by `id(A)` or an `A`-attribute cache)
  so enumeration is not quadratic.
- `_is_proper_power(walk)` / `_canonical_cyclic(walk)`: a band is defined up to
  rotation AND inversion and must not be `u^k` for a shorter `u`; write both helpers
  (rotation-minimal canonical form; period check) at the top of the module.
- **The count arbiter is Task 3's `knit_ar_quiver` tie** — if `enumerate_strings`
  over/under-counts vs the AR-quiver vertex count on a rep-finite string algebra, the
  bug is here (a missed sign turn or a mis-canonicalised inversion). Do NOT weaken the
  count assertion; fix the walk validity.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/strings/__init__.py src/quiverlab/strings/walks.py tests/modules/test_strings_walks.py
git commit -m "feat(strings): Letter alphabet, sigma/epsilon signs, reduced-walk validity, enumerate_strings + find_bands (honest rep-finite-iff-no-bands contract)"
```

---

### Task 2: `string_module` / `band_module` via `from_arrow_action`

**Files:**
- Create: `src/quiverlab/strings/modules.py`
- Test: `tests/modules/test_strings_modules.py`

**Interfaces:**
- Consumes: `walks.py` (`Letter`, `letter_source`/`letter_target`,
  `is_valid_walk`, `find_bands`), `Module.from_arrow_action(algebra,
  dimension_vector, arrow_action, name)` (**classmethod**; vertex-ordered ambient
  basis; **RIGHT** modules; anti-homomorphism `action[x*y] = action[y] @ action[x]`;
  calls `check_module` and raises `QuiverlabError` on a relation violation),
  `modules/linalg_mod` (`zeros`, `identity`), `modules/decompose.py::
  is_indecomposable` (char-scoped), `fields/linalg` for the companion block.
- Produces:
  ```python
  def string_module(A, walk, name=None) -> Module
      # dim = (#vertices on the walk) = len(walk) + 1. Basis z_0..z_n, one per
      # visited vertex; each DIRECT letter alpha: z_{i-1} <- z_i (arrow acts as the
      # 0/1 partial map), each INVERSE letter the transpose. dimension_vector counts
      # z_i's by their vertex. Self-certifies via from_arrow_action's check_module.
  def band_module(A, walk, eigenvalue, mult=1, name=None) -> Module
      # a cyclic string; dim = mult * len(walk). 0/1 blocks for all but the closing
      # letter, which carries the mult x mult invertible Jordan block J_mult(lambda).
      # eigenvalue MUST be a nonzero field element (coerced via A.domain) -- loud
      # otherwise. mult >= 1. Self-certifies via check_module.
  ```

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_strings_modules.py
"""String and band modules (Plan 46 / C5). Self-cert: every materialised module
passes check_module (inside from_arrow_action) and is indecomposable; the simple /
projective string modules match the builtin ones. Loud: a band eigenvalue outside
the field. Over QQ / GF(32003) so is_indecomposable decides (char > dim)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.decompose import is_indecomposable
from quiverlab.modules.hom import is_isomorphic
from quiverlab.strings.modules import band_module, string_module

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def _two_cycle(field=QQ):
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=field)


@selfcert
def test_string_modules_are_indecomposable():
    A = _gentle_a3()
    for w in [((None, 1),), (("a", +1),), (("b", +1),),
              (("a", +1), ("b", -1))]:                 # peak: a then b^-1 at vertex 2
        M = string_module(A, w)
        assert is_indecomposable(M)


@xeng
def test_length_one_direct_string_is_the_projective_arm():
    # the walk (a) over gentle kA3/(ab): the module [1;2] = P1 (dim vector {1:1,2:1}).
    A = _gentle_a3()
    M = string_module(A, (("a", +1),))
    assert M.dimension_vector() == {1: 1, 2: 1, 3: 0}
    assert is_isomorphic(M, A.projective(1))


@xeng
def test_trivial_walk_is_the_simple():
    A = _gentle_a3()
    for v in (1, 2, 3):
        assert is_isomorphic(string_module(A, ((None, v),)), A.simple(v))


@selfcert
def test_band_module_dim_and_indecomposable():
    A = _two_cycle(field=QQ)
    band = _some_band(A)                               # find_bands(A)[0]
    M = band_module(A, band, eigenvalue=1, mult=1)
    assert M.dim == len(band)                          # mult=1
    assert is_indecomposable(M)


@selfcert
def test_band_eigenvalue_must_be_in_the_field():
    A = _two_cycle(field=GF(5))
    band = _some_band(A)
    with pytest.raises(QuiverlabError):
        band_module(A, band, eigenvalue="1/2", mult=1)   # 1/2 not a GF(5) literal? use
                                                         # an unrepresentable token here
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.strings.modules`

- [ ] **Step 3: Implement**

```python
# src/quiverlab/strings/modules.py
"""String and band module materialisation (Plan 46 / C5). Right modules built by
Module.from_arrow_action (check_module self-certifies). Float-free: matrices are
0/1 + exact field scalars."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module
from quiverlab.strings.walks import (invert, is_valid_walk, letter_source,
                                     letter_target)


def _walk_vertices(A, walk):
    """The n+1 vertices z_0..z_n visited by `walk` (z_i between letters). A trivial
    walk ((None, v),) yields the single vertex [v]."""
    Q = A.quiver
    if len(walk) == 1 and walk[0][0] is None:
        return [walk[0][1]]
    verts = [letter_source(Q, walk[0])]
    for ell in walk:
        verts.append(letter_target(Q, ell))
    return verts


def string_module(A, walk, name=None):
    if not (len(walk) == 1 and walk[0][0] is None) and not is_valid_walk(A, walk):
        raise QuiverlabError(f"strings: {walk!r} is not a valid string over A",
                             hint="use enumerate_strings(A) to get valid walks")
    Q, dom = A.quiver, A.domain
    verts = _walk_vertices(A, walk)
    n = len(verts)                                     # basis z_0..z_{n-1}
    dimvec = {v: 0 for v in Q.vertices}
    for v in verts:
        dimvec[v] += 1
    # per-arrow 0/1 action on the local basis (z_i indexed 0..n-1). A DIRECT letter
    # a at position i (a: verts[i] -> verts[i+1]) contributes, for a RIGHT module,
    # the map sending z_{i+1} |-> z_i (arrow acts "downward" on the string picture).
    action = {a: lm.zeros(n, n, dom) for a in Q.arrows}
    for i, ell in enumerate(walk):
        name_i, d = ell
        if name_i is None:
            continue
        if d > 0:                                       # direct: z_{i+1} <- z_i? fix by
            action[name_i][i][i + 1] = dom.one()        # the engine-tau arbiter (Task 3)
        else:                                           # inverse: transpose slot
            action[name_i][i + 1][i] = dom.one()
    return Module.from_arrow_action(A, dimvec, action, name=name or _walk_name(walk))


def band_module(A, walk, eigenvalue, mult=1, name=None):
    if mult < 1:
        raise QuiverlabError("strings: band multiplicity must be >= 1")
    dom = A.domain
    try:
        lam = dom.coerce(eigenvalue)
    except Exception as exc:
        raise QuiverlabError(
            f"strings: band eigenvalue {eigenvalue!r} is not in the field",
            hint="the eigenvalue must be a nonzero element of A.domain") from exc
    if dom.is_zero(lam):
        raise QuiverlabError("strings: band eigenvalue must be nonzero")
    Q = A.quiver
    L = len(walk)
    verts = [letter_source(Q, walk[k]) for k in range(L)]   # cyclic: L vertices
    block = mult
    dim = L * block
    dimvec = {v: 0 for v in Q.vertices}
    for v in verts:
        dimvec[v] += block
    # 0/1 blocks for letters 0..L-2; the closing letter L-1 carries J_mult(lambda).
    # (block placement + which letter closes is arbitrated by check_module +
    #  is_indecomposable; write the L-1 identity blocks and the one Jordan block.)
    action = {a: lm.zeros(dim, dim, dom) for a in Q.arrows}
    ...   # place identity/Jordan blocks per the walk directions; see note
    return Module.from_arrow_action(A, dimvec, action, name=name or "band")
```

**Adjust to reality (Task 2):**
- **The direct-letter matrix orientation (`z_{i+1} |-> z_i` vs `z_i |-> z_{i+1}`) is a
  CONVENTION** — write one, and let Task 3's `is_isomorphic(string_module(w),
  A.projective/simple)` and the string-τ ≡ engine-τ arbiter decide it. If
  `test_length_one_direct_string_is_the_projective_arm` fails (the walk `(a)` should be
  `P1 = [1;2]` over gentle kA3/(ab)), transpose the two slot assignments once and
  document. This is the house "oracle decides the orientation" pattern.
- `Module.from_arrow_action`'s `check_module` **is** the self-certificate: if the 0/1
  matrices do not satisfy the relations of `A`, it raises `QuiverlabError` — so a walk
  that (wrongly) crossed a relation is caught at materialisation, not silently
  accepted. This makes Task 1's `is_valid_walk` and Task 2 mutually checking.
- The band block placement (which single letter carries the Jordan block, closing the
  cycle) is the one structural choice; `is_indecomposable(band_module(...))` over QQ is
  the arbiter (a mis-placed block splits). For `mult=1`, `J_1(lambda) = [lambda]`.
- `dom.coerce`/`dom.one`/`dom.is_zero` are the `Domain` API (`fields/`); read one
  existing caller (`Module.from_arrow_action` itself) for the exact spelling. A
  GF(p) literal that is not a residue, or a `Fraction` over GF(p), must raise — the
  loud-eigenvalue test pins it.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/strings/modules.py tests/modules/test_strings_modules.py
git commit -m "feat(strings): string_module + band_module via from_arrow_action -- check_module self-certified, indecomposable spot-checks"
```

---

### Task 3: string-module τ by hooks/cohooks — cross-engine arbitrated

**Files:**
- Create: `src/quiverlab/strings/ar_strings.py`
- Test: `tests/modules/test_strings_ar.py`

**Interfaces:**
- Consumes: `walks.py` (`string_signs`, `letter_source`/`letter_target`,
  `invert`, `is_valid_walk`), `modules.py::string_module`, `Module.tau()`/
  `.tau_minus()` (the TRUSTED P23/24 translate — the arbiter),
  `modules/hom.py::is_isomorphic`, `modules/ar.py::knit_ar_quiver`/`ARQuiver`
  (P41 — the census / count cross-oracle).
- Produces:
  ```python
  def string_tau(A, walk) -> walk           # the walk w' with string_module(w') ~ tau(string_module(w))
  def string_tau_minus(A, walk) -> walk     # tau^-1, dual
      # Butler-Ringel: tau adds a cohook on one end and deletes a hook on the other
      # (adding/removing the maximal direct/inverse "arm" allowed by the string
      # conditions). A projective string returns the empty/None marker (tau = 0).
  ```

**Convention (arbitrated, not assumed):** the hook/cohook *side* (which end gets the
cohook added, which gets the hook removed) is fixed by the engine-τ isomorphism
below, exactly as the cup sign and composition order are arbitrated in P20/P37. Write
one side; if the isomorphism fails, swap ONCE and document.

- [ ] **Step 1: Write the failing tests**

```python
# tests/modules/test_strings_ar.py
"""String-module tau by hooks/cohooks (Plan 46 / C5), ARBITRATED against the trusted
Plan-41/23 engine tau: on every walk, string_module(string_tau(w)) ~ string_module(w).tau().
Cross-engine: the string census count equals knit_ar_quiver's vertex count on a
rep-finite string algebra."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.ar import knit_ar_quiver
from quiverlab.modules.hom import is_isomorphic
from quiverlab.strings.ar_strings import string_tau, string_tau_minus
from quiverlab.strings.modules import string_module
from quiverlab.strings.walks import enumerate_strings

xeng = pytest.mark.oracle_crossengine


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def _kA3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=[], field=QQ)


@xeng
@pytest.mark.parametrize("factory", [_kA3, _gentle_a3])
def test_string_tau_matches_engine_tau(factory):
    A = factory()
    cen = enumerate_strings(A, max_length=6)
    for w in cen.walks:
        M = string_module(A, w)
        if M.tau().dim == 0:                          # projective: no non-proj tau
            wt = string_tau(A, w)
            assert wt is None
            continue
        wt = string_tau(A, w)
        assert is_isomorphic(string_module(A, wt), M.tau())


@xeng
def test_tau_and_tau_minus_are_inverse_on_non_projective_injective():
    A = _gentle_a3()
    for w in enumerate_strings(A, max_length=6).walks:
        M = string_module(A, w)
        if M.tau().dim == 0 or M.tau_minus().dim == 0:
            continue
        assert is_isomorphic(string_module(A, string_tau_minus(A, string_tau(A, w))), M)


@xeng
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_census_count_equals_ar_quiver_vertex_count(n, count):
    arrows = {chr(ord("a") + i): (i + 1, i + 2) for i in range(n - 1)}
    A = Quiver(list(range(1, n + 1)), arrows).algebra(relations=[], field=QQ)
    cen = enumerate_strings(A, max_length=2 * n)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    assert cen.count == len(ar.vertices) == count
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** — the hook/cohook rule. A **hook** on a string `w` is the
  extension by the maximal *inverse* arm the string conditions allow at one end; a
  **cohook** is the maximal *direct* arm. Butler–Ringel: `tau(M(w))` is obtained by
  simultaneously *removing a hook* at one end of `w` and *adding a cohook* at the
  other (up to the projective/injective edge cases where an end has no hook to
  remove or no cohook to add — those are the ends that touch a projective/injective).

```python
# src/quiverlab/strings/ar_strings.py
"""String-module AR translate by hooks/cohooks (Butler-Ringel 1987). Combinatorial;
the trusted engine tau (Plan 23/24) is the arbiter of the side convention."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.strings.walks import (invert, is_valid_walk, letter_source,
                                     letter_target, string_signs)


def _add_cohook(A, walk, at_end):        # extend by the maximal DIRECT arm allowed
    ...                                  # append/prepend arrows while is_valid_walk
def _remove_hook(A, walk, at_end):       # delete the maximal INVERSE arm at an end
    ...


def string_tau(A, walk):
    """Remove a hook at the START, add a cohook at the END (the sign of "start"/"end"
    is arbitrated by the engine-tau isomorphism test). None when M(walk) is
    projective (no hook to remove there)."""
    if _is_projective_string(A, walk):
        return None
    w1 = _remove_hook(A, walk, at_end="start")
    if w1 is None:
        return None
    return _add_cohook(A, w1, at_end="end")


def string_tau_minus(A, walk):
    if _is_injective_string(A, walk):
        return None
    w1 = _remove_hook(A, walk, at_end="end")     # dual sides
    if w1 is None:
        return None
    return _add_cohook(A, w1, at_end="start")
```

**Adjust to reality (Task 3):**
- Write `_add_cohook`/`_remove_hook` as greedy maximal extensions/deletions gated by
  `is_valid_walk` and the `string_signs` side data. The **canonical reference is
  Butler–Ringel 1987 §3** (hooks/cohooks) — cite it; reproduce its rule, do not
  invent one.
- **The arbiter is `test_string_tau_matches_engine_tau`.** If the isomorphism fails
  on the very first walk, the hook/cohook SIDES are swapped: swap the `at_end`
  arguments in `string_tau`/`string_tau_minus` ONCE (start<->end) and re-run. If it
  fails only on projective/injective ends, the projective-string / injective-string
  predicate is off — fix `_is_projective_string` (a string module is projective iff it
  is `P_v` for some v, checkable by `identify_standard` on the materialised module, or
  by the walk being a maximal "direct then inverse" arm — use `identify_standard` as
  the ground truth). NEVER weaken the isomorphism assertion; the engine tau is trusted.
- `M.tau().dim == 0` is the projective test at the module level — cheap and exact;
  use it to branch, and cross-check against the combinatorial `_is_projective_string`
  (they must agree — assert it in an internal `oracle_selfcert` if convenient).
- The census-count tie uses `knit_ar_quiver` (P41). If the counts disagree on a
  rep-finite string algebra, the bug is Task 1's enumeration (a missed inversion
  canonicalisation) — the AR quiver is the trusted count.

- [ ] **Step 4: Run tests** — Expected: PASS (side convention flipped if needed)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/strings/ar_strings.py tests/modules/test_strings_ar.py
git commit -m "feat(strings): string_tau/string_tau_minus by hooks/cohooks -- arbitrated string_module(string_tau(w)) ~ M(w).tau(); census count = AR-quiver count"
```

---

### Task 4: the AG (Avella-Alaminos–Geiss) derived invariant

**Files:**
- Create: `src/quiverlab/strings/ag.py`
- Test: `tests/invariants/test_ag_invariant.py`

**Interfaces:**
- Consumes: `walks.py` (`string_signs`, `letter_source`/`letter_target`),
  `resolutions_cs/build.py::reduction_system_of`, `recognizers.py::is_gentle`,
  `_len2_in_ideal`, `combinat/quiver.py::Quiver`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class AGInvariant:
      pairs: tuple             # sorted tuple of (n, m), a MULTISET
      def as_multiset(self) -> dict     # {(n, m): multiplicity}
  def permitted_threads(A) -> list      # maximal NONZERO directed paths (partition Q_1)
  def forbidden_threads(A) -> list      # maximal paths of RELATIONS (partition Q_1; may be cyclic)
  def ag_invariant(A) -> AGInvariant
      # Avella-Alaminos-Geiss 2008: the alternating permitted/forbidden thread walk,
      # output the multiset of (n, m) = (# permitted arrows, # forbidden arrows) per
      # cycle, with the cyclic-thread special pairs (n, 0)/(0, m). GENTLE only (loud
      # otherwise). DERIVED-INVARIANT but NOT complete -- stated in the docstring.
  ```

**Mathematical grounding (derived in this plan; the paper is the binding oracle).**
For a gentle algebra `A = kQ/I` (I generated by length-2 monomials):
- A **permitted thread** is a maximal nonzero path — a directed path `c_1...c_k`
  with `c_i c_{i+1} not in I`, maximal on both ends. Because gentle, the permitted
  successor of each arrow (the unique `b` with `a*b not in I`, `source(b)=target(a)`)
  is unique when it exists, so **permitted threads partition `Q_1`**: every arrow lies
  in exactly one.
- A **forbidden thread** is a maximal path of relations — `c_1...c_l` with every
  `c_i c_{i+1} in I`, maximal; a single arrow in no relation is its own length-1
  forbidden thread. Gentle ⇒ the forbidden successor is unique, so **forbidden
  threads partition `Q_1`** too (a self-injective gentle algebra has a *cyclic*
  forbidden thread — the relations close up, e.g. `...abab...`).
- The AAG algorithm alternates permitted and forbidden threads around the quiver and
  records, per closed cycle, a pair `(n, m)` = (number of permitted arrows,
  number of forbidden arrows) traversed; a cyclic forbidden thread contributes an
  isolated `(0, m)`, a cyclic permitted thread an isolated `(n, 0)`.

**Two worked thread structures (derived here; the `(n, m)` counting convention is
arbitrated against AAG 2008 §3 — see the arbiter note):**

*Example A — `A = kA_3/(a*b)`, `Q: 1 -a-> 2 -b-> 3`, `I = (a*b)`, gentle, dim 5.*
Nonzero paths of length ≥1: `a`, `b` (a*b = 0). Permitted threads: `{a}` (cannot
extend — after `a` at vertex 2 the only out-arrow `b` has `a*b in I`; before `a` at
vertex 1 there is no in-arrow) and `{b}` (dually) — **two permitted threads, each
length 1** (`sum n = 2 = |Q_1|`). Forbidden threads: `a*b in I`, maximal (no in-arrow
at 1, no out-arrow at 3) — **one forbidden thread `a*b`, length 2** (`sum m = 2 =
|Q_1|`). The single non-self-injective cycle gives **one pair**; with `sum n = sum m
= 2` the candidate is `AG(kA_3/(a*b)) = {(2, 2)}` — *arbitrated*: if AAG's
trivial-thread bookkeeping splits it, the paper's value (e.g. `{(1,1),(1,1)}`)
wins; reproduce the paper's own tabulated `kA_3/rad^2` entry.

*Example B — `A = kQ/(a*b, b*a)`, the 2-cycle `1 -a-> 2 -b-> 1`, gentle self-injective,
dim 4.* Permitted threads: `{a}`, `{b}` (both length 1; `a*b, b*a in I` block
extension) — **two permitted threads**, `sum n = 2`. Forbidden threads: `a*b in I`
and `b*a in I` close up into a single **cyclic forbidden thread `(a b)` of period 2**
(`...abab...`) covering both arrows — `sum m = 2`. The cyclic forbidden thread
contributes `(0, 2)`; the permitted arms contribute the complementary `(2, 0)` (or
`(1,0)+(1,0)`) — candidate `AG(2-cycle) = {(2, 0), (0, 2)}` — *arbitrated* against
AAG 2008 (the self-injective / cyclic-thread case).

Because I can rigorously derive the **thread structure** but the exact `(n, m)`
counting + trivial/cyclic-thread bookkeeping is the paper's convention, this task's
**primary literature oracle is AAG 2008's own worked examples** (reproduce the
paper's tabulated invariants verbatim); the two examples above are pinned only after
the implementation reproduces the paper, flipping the pair order / split once if the
convention differs. The unconditional self-cert (below) does not depend on the
convention.

- [ ] **Step 1: Write the failing tests**

```python
# tests/invariants/test_ag_invariant.py
"""AG (Avella-Alaminos-Geiss 2008) derived invariant for gentle algebras.
Self-cert (convention-INDEPENDENT): permitted & forbidden threads each PARTITION the
arrows -- every arrow in exactly one of each; sum of permitted lengths = sum of
forbidden lengths = |Q_1|. Literature: the AAG 2008 worked values (reproduced
verbatim; the two hand-derived thread structures pinned after the convention is fixed).
Cross-engine: derived-invariance -- a gentle algebra and a KNOWN derived-equivalent
one share the invariant."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.strings.ag import (ag_invariant, forbidden_threads,
                                  permitted_threads)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def _two_cycle():
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


def _arrows_covered(threads, Q):
    seen = []
    for t in threads:
        seen += list(t)                              # t = tuple of arrow names
    return sorted(seen)


@selfcert
@pytest.mark.parametrize("factory", [_gentle_a3, _two_cycle])
def test_threads_partition_the_arrows(factory):
    A = factory()
    Q = A.quiver
    arrows = sorted(Q.arrows)
    # each arrow appears exactly once across permitted threads, once across forbidden
    assert _arrows_covered(permitted_threads(A), Q) == arrows
    # forbidden threads may be cyclic; count each arrow-slot once per period
    assert _arrows_covered(forbidden_threads(A), Q) == arrows


@selfcert
def test_invariant_sums_equal_arrow_count():
    A = _gentle_a3()
    inv = ag_invariant(A)
    assert sum(n for n, _ in inv.pairs) == len(A.quiver.arrows)
    assert sum(m for _, m in inv.pairs) == len(A.quiver.arrows)


@lit
def test_ag_of_gentle_a3():
    # SHARPEN in Step 3 to the AAG 2008 tabulated value for kA_3/rad^2 (derived
    # thread structure: permitted {a},{b}; forbidden {ab}). Pin after reproducing
    # the paper's own examples.
    A = _gentle_a3()
    inv = ag_invariant(A)
    assert inv.as_multiset() == {(2, 2): 1}          # arbitrated -- see plan Example A


@lit
def test_ag_of_two_cycle():
    A = _two_cycle()
    inv = ag_invariant(A)
    assert inv.as_multiset() == {(2, 0): 1, (0, 2): 1}   # arbitrated -- see Example B


@xeng
def test_derived_invariance_under_opposite():
    # A and A^op are derived equivalent (D is a duality); the AG invariant is a
    # DERIVED invariant, so ag_invariant(A) == ag_invariant(A^op) up to the pair
    # order convention (n,m)<->(m,n) that D induces -- assert the multiset match
    # under whichever orientation the implementation fixes (documented).
    A = _gentle_a3()
    invA = ag_invariant(A)
    invOp = ag_invariant(A.opposite())
    ms, msop = invA.as_multiset(), invOp.as_multiset()
    assert ms == msop or ms == {(m, n): c for (n, m), c in msop.items()}


@selfcert
def test_non_gentle_refused():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 3), "c": (3, 4)}).algebra(
        relations=["a*b", "a*b*c"], field=QQ)         # relation length 3: string not gentle
    with pytest.raises(QuiverlabError):
        ag_invariant(A)
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement** the thread enumeration (deterministic permitted/forbidden
  successors from `_len2_in_ideal`) and the AAG alternating walk. **Reproduce AAG
  2008 §3's algorithm and its worked examples first**; only then sharpen the two
  `@lit` pins. Refuse non-gentle inputs loudly (`is_gentle(A)` is False).

**Adjust to reality (Task 4):**
- Build the permitted/forbidden successor maps once from `reduction_system_of(A)` +
  `_len2_in_ideal`: `perm_succ[a]` = the unique arrow `b`, `source(b)=target(a)`,
  `a*b not in I` (else `a` ends a permitted thread); `forb_succ[a]` = the unique `b`
  with `a*b in I`. Threads = orbits of these partial permutations. The **cyclic**
  forbidden thread (self-injective case) is a full orbit of `forb_succ` with no start
  — detect it (no arrow in the orbit is a thread-start) and emit its `(0, m)` pair.
- **The partition self-cert (`test_threads_partition_the_arrows`) is convention-free
  and unconditional — it must pass regardless of the AAG turn rule.** If it fails, the
  successor maps are wrong (a missed gentle branch). This is the anchor that catches
  enumeration bugs independent of the `(n, m)` counting.
- **The `(n, m)` pairs are arbitrated against AAG 2008.** Sharpen the two `@lit`
  multiset assertions to the paper's convention: reproduce the paper's own tabulated
  algebras first (the primary pins), then confirm `kA_3/rad^2` and the 2-cycle match
  the derived thread structures. If the paper orders `(m, n)` or splits the single
  cycle differently, adopt the paper's convention and flip the two pins ONCE +
  document. **Never keep a `(n, m)` pin the implementation cannot reproduce from AAG.**
- Honesty: the module docstring AND the verification page state ag_invariant is a
  derived **invariant, not a complete** one (completeness needs the graded
  Opper–Plamondon–Schroll geometric data — out of scope). Cite AAG 2008.
- `A.opposite()` is the P23 opposite algebra; confirm it returns a gentle presented
  algebra (it does — reversing a gentle quiver is gentle). If `A.opposite()` yields a
  structure-constant (presentation-less) algebra for which `reduction_system_of`
  refuses, skip the opposite cross-check on that input and note it (the derived-
  invariance oracle then rests on the AAG worked pairs + a tilting-equivalent example).

- [ ] **Step 4: Run tests** — Expected: PASS (pins sharpened to AAG)
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/strings/ag.py tests/invariants/test_ag_invariant.py
git commit -m "feat(strings): AG (Avella-Alaminos-Geiss) derived invariant -- permitted/forbidden threads, honest not-complete; partition self-cert + AAG-2008 pins"
```

---

### Task 5: `BrauerGraphAlgebra(graph, multiplicities, field)`

**Files:**
- Create: `src/quiverlab/families/brauer.py`
- Modify: `src/quiverlab/families/__init__.py` (export `BrauerGraph`,
  `BrauerGraphAlgebra`), `src/quiverlab/families/discover.py` (CATALOG `FamilyInfo`
  + `_iter_families` skip-set for the non-scalar constructor)
- Test: `tests/families/test_brauer.py`

**Interfaces:**
- Consumes: `combinat/quiver.py::Quiver` + `Quiver.algebra(relations, field,
  degree_bound)`; the `families/nakayama.py` loud-hint validation template (read
  `_validate_*` there — mirror the `QuiverlabError(..., hint=...)` shape); the
  `families/trivial_extension.py` length-lex kernel idiom
  (`_relation_string`/`extract_relations`/`_solve_combo`) for emitting the
  Brauer-graph relations; `families/nakayama.py::NakayamaAlgebra(n, l, cyclic=True,
  field=...)` (the Brauer-star construction cross-oracle); `Algebra.cartan_matrix`,
  `Algebra.is_symmetric` (Plan-29/31 trace-form certifier), `Algebra.dim`.
- Produces:
  ```python
  @dataclass(frozen=True)
  class BrauerGraph:
      # a connected ribbon multigraph: edges (each an (u, v) unordered pair, loops
      # allowed), and a CYCLIC ordering of the edge-ends around each vertex.
      edges: tuple            # tuple of frozenset/tuple endpoint pairs, indexed 0..E-1
      cyclic_order: dict      # vertex -> tuple of (edge_index, end_tag) in ribbon order
      def valency(self, v) -> int          # # edge-ends at v (a loop counts twice)
      def validate(self) -> None           # connected, every end ordered, loud hints

  def BrauerGraphAlgebra(graph, multiplicities, field=None) -> Algebra
      # quiver: one VERTEX per EDGE of `graph`; arrows: for each graph-vertex v and
      # each consecutive pair (e, e') in the cyclic order around v, an arrow e -> e'
      # (a cycle of length valency(v)). relations: the standard Brauer presentation
      # (special-cycle powers agree at shared edges; C_v^{m_v} * (next) = 0; length-2
      # paths not around a common vertex are 0). Presented via Quiver.algebra.
      # CERTIFIED per instance: dim == sum_v multiplicities[v] * valency(v)^2 (loud
      # QuiverlabError otherwise). BGAs are symmetric (Plan-29 certifier -> True).
  ```

**Dimension formula (derived here).** For a Brauer graph algebra, the
projective-injective `P_e` of an edge `e = {u, v}` is the amalgam of two uniserial
"arms" (around `u`, around `v`) glued at a common top `S_e` and common socle `S_e`;
the `u`-arm has `m_u * val(u)` composition factors, the `v`-arm `m_v * val(v)`, so
`dim P_e = m_u*val(u) + m_v*val(v)` (a truncated leaf — valency 1, multiplicity 1 —
contributes `1`). Summing over edges and regrouping by vertex (each vertex `v` is an
end of `val(v)` edge-ends, each contributing `m_v*val(v)`):
`dim A = sum_e dim P_e = sum_{v in V} m_v * val(v)^2`. This is the per-instance
certificate. **Brauer-star check:** center `c` valency `n` multiplicity `m`, `n`
truncated leaves ⇒ `dim = m*n^2 + n*(1*1^2) = n*(m*n + 1)` = `n` projectives each
uniserial of length `mn+1` = the symmetric Nakayama `NakayamaAlgebra(n, l=m*n+1,
cyclic=True)`. For `m = 1`: `dim = n(n+1)`, `≅ kZ_n/J^{n+1} ≅ T(kA_n)` (Plan-31).

- [ ] **Step 1: Write the failing tests**

```python
# tests/families/test_brauer.py
"""Brauer graph algebras (Plan 46 / C5). Literature: dim = sum_v m_v*val(v)^2; BGAs
are symmetric. Cross-engine: the Brauer STAR (one central vertex, n truncated leaves)
is NakayamaAlgebra(n, mn+1, cyclic=True) -- byte-equal Cartan matrices; m=1 recovers
T(kA_n) = kZ_n/J^{n+1}."""
import pytest

from quiverlab import GF
from quiverlab.families.brauer import BrauerGraph, BrauerGraphAlgebra
from quiverlab.families.nakayama import NakayamaAlgebra
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _star(n):
    # center 0, leaves 1..n; edges e_i = {0, i}, i=1..n; cyclic order at 0 = e_1..e_n.
    edges = tuple((0, i) for i in range(1, n + 1))
    cyclic = {0: tuple((i, "at0") for i in range(n))}
    cyclic.update({i: ((i, "atleaf"),) for i in range(1, n + 1)})
    return BrauerGraph(edges=edges, cyclic_order=cyclic)


@lit
@pytest.mark.parametrize("n, m", [(3, 1), (4, 1), (3, 2)])
def test_star_dimension_formula(n, m):
    G = _star(n)
    mult = {0: m, **{i: 1 for i in range(1, n + 1)}}   # center mult m, leaves mult 1
    A = BrauerGraphAlgebra(G, mult, field=QQ)
    assert A.dim == n * (m * n + 1)                    # sum_v m_v*val(v)^2 = m*n^2 + n


@xeng
@pytest.mark.parametrize("n, m", [(3, 1), (4, 1), (5, 1), (3, 2)])
def test_brauer_star_is_symmetric_nakayama(n, m):
    G = _star(n)
    mult = {0: m, **{i: 1 for i in range(1, n + 1)}}
    A = BrauerGraphAlgebra(G, mult, field=QQ)
    N = NakayamaAlgebra(n=n, l=m * n + 1, cyclic=True, field=QQ)
    assert A.dim == N.dim
    assert A.cartan_matrix() == N.cartan_matrix()      # byte-equal Cartan (arbiter)


@lit
def test_brauer_graph_algebras_are_symmetric():
    G = _star(3)
    A = BrauerGraphAlgebra(G, {0: 2, 1: 1, 2: 1, 3: 1}, field=QQ)
    assert A.is_symmetric() is True                    # Plan-29 trace-form certifier


@lit
def test_dimension_certificate_refuses_on_bad_wiring():
    # a mis-specified cyclic order that breaks the presentation must raise, not
    # silently return a wrong-dim algebra.
    from quiverlab.errors import QuiverlabError
    bad = BrauerGraph(edges=((0, 1),), cyclic_order={0: (), 1: ()})   # no ends ordered
    with pytest.raises(QuiverlabError):
        BrauerGraphAlgebra(bad, {0: 1, 1: 1}, field=QQ)
```

- [ ] **Step 2: Run to verify failure** — `ModuleNotFoundError: quiverlab.families.brauer`

- [ ] **Step 3: Implement** — build the quiver (vertices = edges; arrows = consecutive
  pairs around each graph-vertex), emit the Brauer relations, present via
  `Quiver.algebra`, and certify `dim == sum_v m_v*val(v)^2`. Mirror
  `families/trivial_extension.py::_presented_trivial_extension` end-to-end (the
  quiver-build + `extract_relations` + widen-once dimension certificate) and
  `families/nakayama.py`'s validation shape.

```python
# src/quiverlab/families/brauer.py  (skeleton)
"""Brauer graph algebras (Plan 46 / C5). A ribbon multigraph with vertex
multiplicities presents a symmetric special-biserial kQ/I: vertices = edges, arrows =
consecutive edge pairs around graph-vertices, relations = the standard Brauer
presentation. Per-instance dimension certificate dim = sum_v m_v*val(v)^2; loud
otherwise. Schroll 2018 (survey); Wald-Waschbuesch 1985. Float-free / exact."""
from __future__ import annotations
from dataclasses import dataclass

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError


def BrauerGraphAlgebra(graph, multiplicities, field=None):
    graph.validate()                                   # loud on a broken ribbon graph
    for v in _vertices(graph):
        if v not in multiplicities or multiplicities[v] < 1:
            raise QuiverlabError(
                f"BrauerGraphAlgebra: missing/invalid multiplicity at vertex {v!r}",
                hint="every graph vertex needs an integer multiplicity >= 1")
    Q, rels = _present(graph, multiplicities)          # vertices=edges, Brauer arrows+rels
    A = Q.algebra(relations=rels, field=field)
    expected = sum(multiplicities[v] * graph.valency(v) ** 2 for v in _vertices(graph))
    if A.dim != expected:
        raise QuiverlabError(
            f"BrauerGraphAlgebra: dim {A.dim} != certificate {expected} "
            f"(= sum_v m_v*val(v)^2)",
            hint="the ribbon ordering or the multiplicities do not present a Brauer "
                 "graph algebra; check the cyclic edge order around each vertex")
    A._family_citations = ("schroll_brauer", "wald_waschbusch", "assem_book")
    return A
```

**Adjust to reality (Task 5):**
- **The special cycles + gluing relations are the mathematical crux.** Around each
  graph-vertex `v`, the arrows form a cycle `C_v` of length `val(v)`; the relations
  are (i) for each edge `e = {u, v}` the "special cycle" identity
  `C_u^{m_u}` (starting at `e`) `= C_v^{m_v}` (starting at `e`) as elements ending at
  the socle; (ii) `C_v^{m_v} * alpha = 0` for the arrow `alpha` leaving the last
  vertex; (iii) any length-2 path `e -> e' -> e''` that does not stay within one
  graph-vertex's cycle is `0`. Encode (ii)+(iii) as monomial relations directly; the
  gluing (i) is what makes `dim P_e = m_u*val(u) + m_v*val(v)` — the **dimension
  certificate is the arbiter** that the relation set is right (an over/under-relation
  changes `dim`). Read Schroll 2018 §2 for the exact presentation before writing
  `_present`; do not guess the relation set — the certificate will catch a wrong one,
  but start from the paper.
- **Truncated edges** (a valency-1, multiplicity-1 leaf) carry **no loop** — the leaf
  contributes no arrow; the edge's projective is uniserial. Guard `valency(v) == 1 and
  m_v == 1` specially (no `C_v` cycle). A valency-1 vertex with `m_v > 1` DOES carry a
  loop — handle both.
- **Arrow direction** (`e -> e'` vs `e' -> e'` for consecutive edges around `v`) is a
  CONVENTION pinned by `test_brauer_star_is_symmetric_nakayama` (Cartan byte-equality
  against the known Nakayama family). Write one direction; if the Cartan disagrees,
  reverse the cyclic-order reading once and document.
- `Quiver.algebra` may need a `degree_bound` for the higher-multiplicity /
  higher-valency cases so the Gröbner route terminates — pass a generous bound derived
  from the max projective length `max_e (m_u*val(u) + m_v*val(v)) + 1` (mirror
  `preprojective.py`'s auto-bound). If the presented backbone cannot certify
  finiteness it raises `NotFiniteDimensionalError` — but a genuine Brauer graph algebra
  IS finite-dimensional, so this signals a wiring bug, not an honest infinite case.
- `families/__init__.py` export + `discover.py` CATALOG entry (mirror an existing
  `FamilyInfo`); add `BrauerGraphAlgebra` to `_iter_families`'s skip-set (like the
  `"zoo"` precedent) since it takes a `BrauerGraph` + dict, not scalar args, so the
  webapp form-builder must not introspect it. Surface it as a drawable **preset**
  (a small star) via `scripts/gui_build_hook.py::_preset_algebras()` if that hook
  exists — otherwise defer the preset to Task 6.

- [ ] **Step 4: Run tests** — Expected: PASS
- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families/brauer.py src/quiverlab/families/__init__.py \
        src/quiverlab/families/discover.py tests/families/test_brauer.py
git commit -m "feat(families): BrauerGraphAlgebra(graph, multiplicities) -- presented kQ/I, dim=sum m_v*val(v)^2 certified, Brauer-star = symmetric Nakayama cross-oracle, symmetric"
```

---

### Task 6: `strings_block(A)` + the `strings` GUI scalar kind

**Files:**
- Create: `src/quiverlab/strings/block.py`
- Modify: `src/quiverlab/hpc/spec.py` (scalar-kind dispatch + `_snip`),
  `docs/gui/runner.py` (the Pyodide twin: matching dispatch + `_snip` + ETA),
  `docs/gui/gui.js` + `webapp/static/gui/gui.js` (checkbox, `S.ids`, push-list,
  `renderBlock`, `scheduleProbe`), `webapp/templates/index.html` (checkbox),
  `webapp/server/i18n/en.json` + `es.json` (`inv.strings`, `block.strings.*`),
  `src/quiverlab/trace/results_html.py` (`_HEADINGS` + a `_strings_html` branch),
  `tests/webapp/_runner_goldens.json` + `tests/webapp/test_runner_delegation.py`
  (ONE golden, existing byte-identical)
- Test: `tests/webapp/test_strings_block_p46.py`, `tests/gui/test_strings_runner_twin.py`

**Interfaces (the `recognizers` kind is the exact template — algebra-only, schema v1,
shared builder + per-runner `citations`):**
- `strings_block(A)` mirrors `invariants/recognizers.py::recognizers_block(A)`:
  returns a dict; each runner stamps `block["citations"] =
  _citation_pairs(block["references"])`. Shape:
  ```python
  {"recognizers": {"is_special_biserial": bool|{"error":...},
                   "is_string": ..., "is_gentle": ...},
   "strings": {"count": int, "status": "complete"|"budget", "max_length": int,
               "sample": [walk_repr, ...]},                     # census (string algebras)
   "bands": {"exist": bool, "sample": [band_repr, ...]},
   "rep_type": "finite"|"infinite"|"unknown",                   # string algebras only
   "ag_invariant": [[n, m], ...] | None,                        # only when gentle
   "references": ["butler_ringel", "avella_geiss", "assem_book"]}
  # A non-string algebra: strings/bands/rep_type honestly null-ish + a "note"; the
  # recognizer verdicts still populate. Guarded per-field like recognizers_block.
  ```
- `strings` is an **algebra-only scalar kind** (schema v1, NO module block), routed by
  `spec.py::_dispatch` (NOT `_dispatch_module`); `estimator.sizing_dim` is
  algebra-dim based, so **no estimator edit** (the recognizers-kind precedent).
- The **seven GUI touchpoints** follow the `recognizers` kind verbatim (agent-verified
  anchors): checkbox `qlgui-strings` (both `gui.js` + `index.html`), `S.ids`,
  the request push-list, the `renderBlock` branch (census table + AG multiset + a
  rep-type line), `scheduleProbe`, the i18n `inv.strings`/`block.strings.*` chain
  (`en.json` + `es.json`), and the `results_html.py` heading + branch. `app.js` (the
  server webapp renderer) is OPTIONAL — the recognizers kind renders through the
  shared GUI path; match whatever `recognizers` did.

- [ ] **Step 1: Write the failing cross-runner test** (unmarked — extras-gated dir; copy
  the recognizers-kind runner-pair fixture):

```python
# tests/webapp/test_strings_block_p46.py
"""The `strings` algebra-only scalar kind: served by hpc.spec, mirrored by the
Pyodide twin, byte-identical block."""


def test_strings_block_shape_gentle_a3(tmp_path):
    # schema-1 request: gentle kA3/(ab) over QQ, compute ["strings"]. Assert:
    #   block["recognizers"]["is_gentle"] is True
    #   block["strings"]["status"] == "complete" and count == 6
    #   block["bands"]["exist"] is False and block["rep_type"] == "finite"
    #   block["ag_invariant"] is not None                 # gentle => AG present
    #   "avella_geiss" in [k for k, _ in block["citations"]]
    ...


def test_strings_block_two_cycle_reports_infinite(tmp_path):
    # 2-cycle kQ/(ab,ba): bands exist, status "budget", rep_type "infinite".
    ...


def test_twin_parity(tmp_path):
    # run the same request through docs/gui/runner.py; json.dumps(sort_keys=True)
    # equality on the strings block (both runners byte-identical).
    ...
```

- [ ] **Step 2: Implement** `strings/block.py` + the spec.py branch + the
  `docs/gui/runner.py` twin (keep them shape-identical), the two `gui.js` touchpoints
  + `index.html` checkbox, the ETA entry (`"strings": <small>` beside the recognizers
  entry), the i18n keys (`inv.strings`, `block.strings.title`, `block.strings.census`,
  `block.strings.bands`, `block.strings.reptype`, `block.strings.ag` — EN and ES), the
  `_snip` recipe (`"strings": "strings_block(A)"` or the public
  `A.strings()`/`enumerate_strings(A)` spelling you expose), and the
  `results_html.py` `_HEADINGS["strings"] = "Strings & bands"` + `_block_html` branch.

- [ ] **Step 3: Add ONE golden fixture** (`strings_gentle_a3`) to `_runner_goldens.json`;
  note it in `test_runner_delegation.py`'s docstring change-log. Run the delegation
  test BEFORE adding to confirm existing goldens stay byte-identical.

- [ ] **Step 4: Run the gates**

Run: `... -m pytest tests/webapp/test_strings_block_p46.py tests/webapp/test_runner_delegation.py tests/gui/test_strings_runner_twin.py tests/hpc -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(gui,webapp,hpc): `strings` algebra-only scalar kind -- census table + AG multiset + honest rep-type verdict, both runners byte-identical"
```

---

### Task 7: QPA cross-oracle battery

**Files:**
- Modify: `src/quiverlab/qpa/scripts.py` / `src/quiverlab/qpa/crosscheck.py` (only if a
  new helper is needed; the recognizer + decompose crosschecks already exist)
- Test: `tests/qpa/test_strings_qpa.py`

**Interfaces:**
- Consumes: the QPA session (`session.should_skip_qpa()`, `session.run`,
  `session.libgap_handle`), `qpa/scripts.py::quiver_and_algebra_script(A)` (the
  recognizer crosscheck template — `IsGentleAlgebra`/`IsSpecialBiserialAlgebra` are
  **already** live-crosschecked in `tests/qpa/test_recognizers_qpa.py`),
  `qpa/scripts.py::module_decl(...)` + `qpa/crosscheck.py::_flat_dimvec_multiset` /
  `crosscheck_decompose` (the Plan-30 module-level decompose oracle — QPA's
  `DecomposeModuleWithMultiplicities`), `strings/modules.py::string_module`,
  `qpa_module.py::graded_form` (our module -> QPA row convention).
- Produces: `tests/qpa/test_strings_qpa.py`, standard skipif header. **Honest scope:
  QPA has NO string/band enumeration verb and NO AG invariant** — the crosschecks are
  therefore (a) the recognizer verdicts (already live) re-run on the string-family zoo,
  and (b) **module-level** validation: `decompose` of a direct sum of `string_module`s
  agrees with QPA's `DecomposeModuleWithMultiplicities` dimension-vector multiset.

- [ ] **Step 1: Probe live QPA + write the battery** (mirror
  `tests/qpa/test_recognizers_qpa.py` for the verdicts and the Plan-30 tor/decompose
  QPA style for the modules; probe `NamesGVars()` for any string surface and FAIL if
  one ever appears — the Plan-35 skip-that-fails-if-appears precedent):

```python
# tests/qpa/test_strings_qpa.py
"""QPA as the oracle for the string family (Plan 46). QPA has recognizers
(IsGentleAlgebra / IsSpecialBiserialAlgebra) but NO string/band enumerator and NO AG
invariant -- so the crosschecks are: recognizer PARITY on the string zoo, and
MODULE-level decompose of a sum of string modules. qpa-marked."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.morphism import direct_sum
from quiverlab.qpa import session
from quiverlab.strings.modules import string_module
from quiverlab.strings.walks import enumerate_strings

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _gentle_a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)


def test_qpa_has_no_string_enumeration_surface():
    lg = session.libgap_handle()
    for name in ("StringModules", "BandModules", "AGInvariant"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), \
            f"QPA now ships {name} -- add a real crosscheck (honest-scope changed)"


def test_recognizer_parity_on_string_zoo():
    A = _gentle_a3()
    A.crosscheck("is_gentle").assert_agree()             # IsGentleAlgebra (already live)
    A.crosscheck("is_special_biserial").assert_agree()


def test_string_module_sum_decompose_vs_qpa():
    # a direct sum of a few string modules must decompose (QPA) back to the same
    # dimension-vector multiset -- validates our materialisation at module level.
    A = _gentle_a3()
    walks = [w for w in enumerate_strings(A, max_length=6).walks][:4]
    mods = [string_module(A, w) for w in walks]
    D = mods[0]
    for M in mods[1:]:
        D, _, _ = direct_sum(D, M)
    A.crosscheck("decompose", D).assert_agree()          # DecomposeModuleWithMultiplicities
```

- [ ] **Step 2: Wire any missing crosscheck case.** `is_gentle`/`is_special_biserial`
  and `decompose` crosschecks already exist (P38 / P30) — if `A.crosscheck("decompose",
  M)` is not yet a dispatch case for an arbitrary materialised module, add it mirroring
  `crosscheck_decompose` (`crosscheck.py`); otherwise no `src/` change. **The
  string/band/AG crosschecks are honest skips** (no QPA surface) — the
  `test_qpa_has_no_string_enumeration_surface` probe is the standing guard.

- [ ] **Step 3: Run live** `... -m pytest tests/qpa/test_strings_qpa.py -v` (the venv has
  the `[qpa]` extra). Expected: PASS live; the enumeration probe confirms absence.

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/qpa tests/qpa/test_strings_qpa.py
git commit -m "test(qpa): string family -- recognizer parity (IsGentle/IsSpecialBiserial) + module-level decompose of string-module sums; honest no-enumeration-surface guard"
```

---

### Task 8: verification page, citations, README, suite gate

**Files:**
- Modify: `src/quiverlab/citations/references.bib` + `src/quiverlab/citations/registry.py`
- Modify: `docs/verification.md`, `README.md`
- Modify: `docs/plans/2026-08-05-metaplan-v0.2.0.md` (tick the P46 card delivery note)
- Test: existing release gates (`tests/release/test_oracle_classes.py`,
  `tests/citations/`)

- [ ] **Step 1: Citations** (BibTeX-VERIFIED only; `_r(key, bibtex_key, kind, title,
  annotation, *tags)` registry precedent). Add:

```bibtex
@article{ButlerRingel1987,
  author  = {Butler, M. C. R. and Ringel, Claus Michael},
  title   = {Auslander-{R}eiten sequences with few middle terms and applications
             to string algebras},
  journal = {Communications in Algebra},
  volume  = {15}, number = {1-2}, pages = {145--179}, year = {1987},
}
@article{AvellaAlaminosGeiss2008,
  author  = {Avella-Alaminos, Diana and Gei{\ss}, Christof},
  title   = {Combinatorial derived invariants for gentle algebras},
  journal = {Journal of Pure and Applied Algebra},
  volume  = {212}, number = {1}, pages = {228--243}, year = {2008},
}
@incollection{Schroll2018,
  author    = {Schroll, Sibylle},
  title     = {Brauer graph algebras},
  booktitle = {Homological Methods, Representation Theory, and Cluster Algebras},
  series    = {CRM Short Courses}, publisher = {Springer}, pages = {177--223},
  year      = {2018},
}
@article{WaldWaschbusch1985,
  author  = {Wald, Burkhard and Waschb{\"u}sch, Josef},
  title   = {Tame biserial algebras},
  journal = {Journal of Algebra},
  volume  = {95}, number = {2}, pages = {480--500}, year = {1985},
}
```

and in `registry.py` (mirror the `_r("assem_book", "ASS2006", "foundation", ...)`
shape):

```python
_r("butler_ringel", "ButlerRingel1987", "algorithm",
   "Auslander-Reiten sequences for string algebras",
   "Butler-Ringel: the string/band module classification and the hook/cohook "
   "description of the AR translate -- the ground truth for the string subsystem.",
   "modules"),
_r("avella_geiss", "AvellaAlaminosGeiss2008", "algorithm",
   "Combinatorial derived invariants for gentle algebras",
   "The AG-invariant: a multiset of (n,m) pairs from permitted/forbidden threads; "
   "a DERIVED invariant, provably NOT complete.", "invariants"),
_r("schroll_brauer", "Schroll2018", "family",
   "Brauer graph algebras (survey)",
   "The presentation of a Brauer graph algebra from a ribbon graph + multiplicities; "
   "the dimension and symmetric structure.", "families"),
_r("wald_waschbusch", "WaldWaschbusch1985", "foundation",
   "Tame biserial algebras",
   "Biserial / special-biserial structure underlying string and Brauer graph "
   "algebras.", "families"),
```

  **Spec-ambiguity resolution (recorded):** the metaplan card says "AAG 2008 worked
  examples, String-Applet / SBStrips spot values". SBStrips (a GAP package) and
  String-Applet are **not** installed in this venv and are not a scriptable oracle
  here — the honest oracles are AAG 2008 (literature pins), our own bar/CS/AR engines
  (cross-engine), and QPA recognizers. State this on the verification page; do NOT add
  a skipped SBStrips battery pretending to be an oracle.

- [ ] **Step 2: Verification page.** Add the Plan-46 subsystem rows:
  - `strings/walks.py` + `strings/modules.py` — `oracle_selfcert` (check_module +
    indecomposable materialisation, thread/walk partition); `oracle_literature`
    (Butler–Ringel `n(n+1)/2` interval count, 2-cycle band existence);
    `oracle_crossengine` (census count ≡ `knit_ar_quiver` vertex count).
  - `strings/ar_strings.py` — `oracle_crossengine` (string-τ ≡ engine τ; τ∘τ⁻ = id).
  - `strings/ag.py` — `oracle_selfcert` (arrow-partition + sum identities);
    `oracle_literature` (AAG 2008 worked values, `kA_3/rad^2`, 2-cycle);
    `oracle_crossengine` (derived-invariance under opposite).
  - `families/brauer.py` — `oracle_literature` (dim `= sum_v m_v*val(v)^2`,
    `is_symmetric`); `oracle_crossengine` (Brauer-star ≡ symmetric Nakayama Cartan).
  - `qpa` — recognizer parity + module-level decompose of string-module sums.
  Add the **honest-scope entries**: (a) for a **string** algebra the Butler–Ringel
  classification is complete **iff there are no bands** (rep-finite); bands ⇒
  rep-infinite ⇒ `enumerate_strings` returns a length-capped **sample**
  (`status="budget"`), never a "complete" list; (b) the **AG invariant is a derived
  invariant but NOT complete** (completeness needs the graded OPS geometric data — out
  of scope) — never claim completeness; (c) QPA has **no** string/band enumeration and
  **no** AG surface (crosschecks are recognizer-level + module-level decompose only);
  SBStrips / String-Applet are not installed and are not oracles here; (d) band modules
  need the eigenvalue in the field (loud otherwise); (e) `decompose`-based
  indecomposability certificates carry the char ≤ dim caveat (batteries over
  QQ / GF(32003)). **Recount the class table** (`tests/release/test_oracle_classes.py`
  drives the numbers — run collection, paste the LIVE counts, re-run to green; do NOT
  guess an at-authoring number given the mid-merge-train drift).

- [ ] **Step 3: README.** One features line: "string & band module classification,
  string-module τ (hooks/cohooks), the Avella-Alaminos–Geiss derived invariant for
  gentle algebras (honest: invariant, not complete), and a Brauer-graph-algebra
  constructor — the C5 gentle/string subsystem."

- [ ] **Step 4: Full gate:**
  `... -m pytest tests/modules tests/families -q` (deep, touched dirs),
  `... -m pytest tests/invariants tests/webapp tests/gui tests/hpc -q -m fast`,
  `... -m pytest tests/qpa -q -m qpa`,
  `... -m pytest tests/release tests/citations -q`,
  and a citation-presence check (`butler_ringel`/`avella_geiss`/`schroll_brauer`
  resolve; the `strings` block carries `avella_geiss`) — all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs(verification): Plan-46 gentle/string oracle rows + honest scope (bands=>rep-infinite sample, AG not complete, no QPA enumeration) + citations + recounted classes"
```

---

## Acceptance (Plan-46 definition of done)

1. `enumerate_strings`/`find_bands`/`string_signs`/`is_valid_walk` (walks),
   `string_module`/`band_module` (modules), `string_tau`/`string_tau_minus`
   (ar_strings), `permitted_threads`/`forbidden_threads`/`ag_invariant` (ag),
   `strings_block`, and `BrauerGraph`/`BrauerGraphAlgebra` all public in
   `quiverlab.strings` / `quiverlab.families`, all loudly-validated, all
   self-certified (materialised modules pass `check_module`; constructions carry a
   dimension certificate).
2. Every materialised string / band module self-certifies (`check_module` inside
   `from_arrow_action`) and is spot-checked indecomposable over QQ / GF(32003); the
   trivial walk = the simple, a length-1 direct string = the projective arm (the
   orientation convention arbitrated by these ties).
3. **String-τ is arbitrated cross-engine:** `string_module(string_tau(w))` is
   `is_isomorphic` to `string_module(w).tau()` on every enumerated walk (the
   hook/cohook side flipped once if needed, documented); `string_tau_minus∘string_tau
   = id` off the projective/injective edges; and the string census count equals the
   Plan-41 `knit_ar_quiver` vertex count (`kA_n = n(n+1)/2`).
4. **Honest classification contract:** `enumerate_strings` claims `status="complete"`
   ONLY when `find_bands` is empty and the DFS closed (rep-finite); the 2-cycle
   `kQ/(ab,ba)` has a band, is reported rep-infinite, and its census is a
   `status="budget"` sample — never a silent "complete" list.
5. **AG invariant is derived-invariant and honestly not complete:** permitted &
   forbidden threads partition the arrows (convention-free self-cert), the sums equal
   `|Q_1|`, the `(n, m)` multiset reproduces AAG 2008's worked values (the two
   hand-derived thread structures pinned to the paper's convention, flipped once if
   it orders `(m, n)` / splits cycles differently), and `ag_invariant(A) =
   ag_invariant(A^op)` up to the documented orientation; non-gentle input refuses
   loudly.
6. **Brauer graph algebras** satisfy `dim = sum_v m_v*val(v)^2` (per-instance
   certificate, loud on bad wiring), are certified `is_symmetric` (Plan-29), and the
   Brauer-star special case is byte-Cartan-equal to `NakayamaAlgebra(n, mn+1,
   cyclic=True)` (and `T(kA_n) = kZ_n/J^{n+1}` at `m=1`) — the construction
   cross-oracle; the arrow-direction convention arbitrated by that Cartan equality.
7. The `strings` algebra-only scalar kind is clickable end-to-end (GUI canvas →
   block → report) in EN+ES, schema v1, both runners byte-identical, ONE golden added
   with a documented change-log entry; the census table, AG multiset, and honest
   rep-type verdict render.
8. Live QPA battery green (`-m qpa`): recognizer parity (`IsGentleAlgebra`/
   `IsSpecialBiserialAlgebra`) on the string zoo + module-level decompose of
   string-module sums, with the standing `NamesGVars()` guard that FAILS if QPA ever
   ships a string/band/AG surface (honest scope: QPA has none today).
9. `docs/verification.md` recounted (live numbers, mid-merge-train honest) with the
   four honest-scope entries (bands ⇒ rep-infinite sample; AG not complete; QPA no
   enumeration + SBStrips not an oracle here; band-eigenvalue/decompose-char caveats);
   `butler_ringel`/`avella_geiss`/`schroll_brauer`/`wald_waschbusch` citations added
   and BibTeX-verified; README line added; deep (touched dirs) + fast + qpa + release
   + citations suites green. No dependency taken on P39/P40/P42–P45 (merged
   independently to `dev`).
