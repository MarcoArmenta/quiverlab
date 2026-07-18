# quiverlab Plan 03 — Gröbner Bases + General kQ/I Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A user hands `Quiver.algebra(relations=[...])` a *non-monomial* relation (`"a*b - c*d"`) and gets a certified finite-dimensional `kQ/I` as a Plan-01 `Algebra`, or a loud `AdmissibilityError`/`NotFiniteDimensionalError` — never a hang, never a guess. Under the hood: noncommutative Gröbner (Buchberger–Mora overlap) completion with a degree bound, exact reduction over any `Domain`, and a `ReductionSystem` object frozen for Chouhy–Solotar (Plan 04) to consume verbatim.

**Architecture:** A new `src/quiverlab/groebner/` package. A length-lexicographic admissible path order (refined by fixed arrow order) selects leading words; overlap/inclusion ambiguity detection (written in the **left-to-right** composition convention) drives Buchberger–Mora completion of the relation rules to a confluent reduction system, capped by a degree bound; the certificate reuses the Plan-01 irreducible-word automaton (`core.monomial.irreducible_paths`) to decide finiteness. `Quiver.algebra` keeps routing monomial inputs through the existing Plan-01 path and dispatches everything else through completion → irreducible-path basis → structure constants → `Algebra`.

**Tech Stack:** Python ≥ 3.10, no new hard dependencies (pure Python + the Plan-01 stack: sympy for exact ℂ). pytest. Reduction is combinatorial (no linear algebra); structure constants and HH reuse the frozen Plan-01 `Algebra`/bar machinery.

---

## Global Constraints

- **Repo root:** `/Users/marco/Desktop/HomologicalNetworks/quiverlab`. All paths below are relative to it.
- **Interpreter:** use the project venv **`/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`** (Python 3.12). The *system* `python`/`python3` is 3.8 and MUST NOT be used — it will fail on 3.10+ syntax (`X | None`, `list[int]`, etc.).
- **Thread throttle:** prefix **every** test command with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2` (Marco's machine has crashed under agent fleets; keep thread/memory pressure low). No new parallelism is introduced in this plan; core code is single-threaded.
- **Exact arithmetic only, through the `Domain` protocol.** All coefficient arithmetic in this plan goes through `dom.add/sub/neg/mul/inv/is_zero/eq/coerce/zero/one`. Never compare or combine coefficients with raw Python operators. Floats fail loudly (Plan-01 `ExactnessError`); relation coefficients are `fractions.Fraction` from the Plan-01 parser and are coerced into the chosen `Domain`.
- **Float ban (AST gate).** `tests/test_no_floats.py` scans all of `src/quiverlab/` for float/complex literals and `float()` calls and MUST stay green. Write **no** float or complex literals anywhere under `src/` (use ints; all "coefficients" are domain elements, never Python numbers). Run it as part of every task's suite.
- **Left-to-right path composition** (Assem–Simson–Skowroński): `a*b` = "first `a`, then `b`", requiring `target(a) == source(b)`. A word is a tuple of arrow names read left to right. **Every** overlap/ambiguity/S-polynomial/reduction definition in this plan is written in this convention (an overlap of leading words `w1, w2` is a word `b` that is a SUFFIX of `w1` and a PREFIX of `w2`). Document it in every relevant docstring.
- **Read-only banks.** `/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/` and every other bank are never read or written in this plan. Plan 03 depends only on Plan-01 source already in `src/quiverlab/`.
- **Full suite green at every commit.** Run `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q` from the repo root before each commit; it must pass.
- **Commits:** conventional prefixes (`feat:`, `test:`, `chore:`); every commit message ends with the trailer line
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- **Frozen Plan-01 interfaces consumed verbatim** (do not modify their signatures): `Domain` (`coerce/add/sub/neg/mul/inv/is_zero/eq/zero/one/characteristic/name`, hooks `parse_entry`/`make_domain`); `fields.CC`, `fields.GF`, `fields.QQ`; `combinat.Quiver` (`.vertices`, `.arrows`, `.source/.target`, `.word_source/.word_target`, `.compose_ok`, `.is_acyclic`); `combinat.Relation`/`parse_relations` (`.terms: tuple[tuple[Fraction, tuple[str,...]], ...]`, `.source`, `.target`, `.is_monomial`, `.max_length`, `.min_length`); `core.Algebra` (`.domain, .T, .unit, .dim, .basis_labels, .is_unit_adapted`, `.multiply`, `.unit_adapted`, ctor `Algebra(domain, T, unit, basis_labels=, _quiver=, _relations=)`); `core.monomial.build_monomial_algebra`, `core.monomial.irreducible_paths`; `errors.*`. **One compatible extension** is made: `Quiver.algebra` gains optional keyword args `degree_bound=None, trace=None` (existing calls `Q.algebra(relations=..., field=...)` are unaffected).

---

## Left-to-right Gröbner theory (read this before Task 1)

The implementer knows Python, not Gröbner theory. This section carries the mathematics; the code steps below implement exactly what is stated here.

**Path algebra.** For a quiver `Q`, `kQ` has as `k`-basis the paths (including a trivial path `e_v` per vertex). A path is a word `w = (a_1, …, a_n)` of arrows read left-to-right with `target(a_i) = source(a_{i+1})`. Product of paths `p·q` = concatenation `p⌢q` if `target(p) = source(q)`, else `0`. An element of `kQ` is a finite `k`-linear combination of paths.

**Admissible order (length-lex).** Fix the total order on arrows = their **insertion order in `Quiver.arrows`** (arrow at position 0 is smallest). For words `u, v`:

```
u < v  iff  len(u) < len(v),
         or len(u) == len(v) and at the first index i where they differ,
            rank(u[i]) < rank(v[i]).
```

This is a *well-founded, two-sided multiplicative* order (proof in the `PathOrder` docstring, Task 1): shorter-first makes it well-founded; the tie-break lex on ranks makes it total; and `u < v ⟹ w·u·x < w·v·x` for all composable `w, x` (length dominates, and on equal lengths the differing index survives concatenation with the same rank comparison). By Bergman's Diamond Lemma, a reduction system all of whose ambiguities resolve to `0` is **confluent**, and the irreducible words (those containing no leading word as a contiguous factor) form a `k`-basis of `kQ/I`.

**Reduction rules.** A relation `r = Σ_w c_w·w` (all `w` parallel: same source, same target) has a leading word `LM(r) = max_< {w : c_w ≠ 0}` with leading coefficient `LC(r)`. The rule is

```
lead := LM(r),   tail := -(1/LC(r))·(r - LC(r)·lead)  ==  { w : -(1/LC(r))·c_w  for w ≠ lead }
```

so `lead ≡ Σ tail` in `kQ/I`, and every word of `tail` is `< lead`. A monomial relation (single word) gives `lead → 0` (empty tail) — exactly the Plan-01 "forbidden word".

**Reducing.** To reduce an element, repeatedly find a word `W = u⌢lead⌢v` containing some rule's `lead` as a contiguous factor and replace that occurrence: `W ↦ Σ_{(c_t,t)∈tail} u⌢t⌢v` scaled by `W`'s coefficient. Because every `tail` word is `< lead`, each step strictly lowers the multiset of words in the well-order, so reduction **terminates**; if the system is confluent the result (the **normal form**) is independent of the order of application.

**Ambiguities (left-to-right).** Confluence is decided by overlap and inclusion ambiguities between leading words `w1, w2`:

- **Overlap.** A nonempty `b` that is a proper SUFFIX of `w1` and a proper PREFIX of `w2`: `w1 = a⌢b`, `w2 = b⌢c`, `b, a, c` possibly making the *ambiguity word* `W = a⌢b⌢c = w1⌢c = a⌢w2` (length `|w1|+|c| ≤ |w1|+|w2|-1`). (In right-to-left literature `b` would be a prefix of `w1` and suffix of `w2`; we use left-to-right throughout.) In code these are the `Ambiguity` fields `left` (`w1`, lead `a⌢b`) and `right` (`w2`, lead `b⌢c`).
- **Inclusion.** One leading word — the **outer** (`left` in code) — has another — the **inner** (`right`) — as a proper contiguous factor: `left.lead = a⌢(right.lead)⌢c` with `(a,c) ≠ ((),())`; ambiguity word `W = left.lead`. (Beware: this is the reverse of some texts that name the inner word first — here `left` is always the outer/containing rule and `right` the inner, matching the `Ambiguity` dataclass and the frozen Plan-04 interface.)

**S-polynomial.** Reduce `W` two ways and subtract (`left`/`right` are the `Ambiguity` fields):

- overlap: `S = (Σ tail(left))⌢c  −  a⌢(Σ tail(right))` (apply the `left` rule at the front, minus the `right` rule at the back);
- inclusion: `S = Σ tail(left)  −  a⌢(Σ tail(right))⌢c` (apply the outer `left` rule, minus the inner `right` rule in context).

`S` is a combination of words `< W`, parallel to `W`. **Completion (Buchberger–Mora):** reduce each `S` by the current rules; if the normal form `g ≠ 0`, add the rule from `g` (`lead = LM(g)`, tail from `g`), which introduces new ambiguities; repeat until a full pass adds nothing.

**Degree bound and certificate.** Only ambiguity words of length `≤ D` (the degree bound) are formed. Over a finite arrow set there are finitely many words of length `≤ D`, so completion under a fixed `D` **always terminates**. After completion let `L = max leading-word length`. Every overlap word has length `≤ 2L−1` and every inclusion word length `≤ L`; therefore if `2L−1 ≤ D` **all** ambiguities were formed and reduced to `0`, so the system is confluent and the certificate is valid. If `2L−1 > D` we refuse to certify → **`AdmissibilityError`** (raise `degree_bound`). If the system is confluent, enumerate irreducible words with the Plan-01 automaton: finite ⟹ certified finite-dimensional with irreducible-path basis; a cycle ⟹ **`NotFiniteDimensionalError`** naming the cycle.

---

## Hand-verified fixtures (expected values worked out by hand)

These five (plus one bonus) are the acceptance oracles; every worked number below is re-derived in the tasks that use it. Arrow order is always the insertion order shown.

### Fixture 1 — loop `x`, relation `x^3` (monomial; both routes must match)
`Q = Quiver([1], {"x": (1,1)})`, `relations=["x^3"]`. Monomial ⇒ Plan-01 route; also forced through Gröbner for the equivalence test.
Rule `xxx → 0`. Self-overlaps of `xxx` give ambiguity words `xxxxx`, `xxxx`; both reduce to `0` (front/back occurrence of `xxx → 0`), so completion adds nothing. Irreducible words avoiding `xxx`: `x, xx`. Basis `[e_1, x, xx]`, labels `["e_1", "x", "x*x"]`, **dim 3**. Structure constants identical to `build_monomial_algebra` (products are `0` or the concatenated path). `2L−1 = 5 ≤ D_default`.

### Fixture 2 — commutative square (non-monomial; the flagship)
`Q = Quiver([1,2,3,4], {"a": (1,2), "b": (2,4), "c": (1,3), "d": (3,4)})`, `relations=["a*b - c*d"]`. Arrow order `a<b<c<d`.
Relation `ab − cd`: both length 2, parallel `1→4`; lex `a<c` ⇒ `ab < cd` ⇒ `LM = cd`. Rule `cd → ab`. The single rule has no self-overlap (`suffix d ≠ prefix c`) and no inclusion, so completion adds nothing; `L = 2`, `2L−1 = 3 ≤ D_default`.
All paths of `kQ`: `e_1..e_4` (4), `a,b,c,d` (4), `ab, cd` (2) — total 10, no longer paths (both length-2 paths end at sink 4). Irreducible = all except `cd`. **Basis** `[e_1,e_2,e_3,e_4, a,b,c,d, ab]`, labels `["e_1","e_2","e_3","e_4","a","b","c","d","a*b"]`, **dim 9**. In the algebra `c·d = ab` (reduction `cd→ab`) — the commutativity is realized in the structure constants.
This algebra is `kA_2 ⊗ kA_2` (product of two `1→2` path algebras; the product quiver is exactly this commuting square). By Künneth, `HH^n(A⊗A) = ⊕_{i+j=n} HH^i(A)⊗HH^j(A)`; `kA_2` is directed hereditary with `HH^0 = k` (connected) and `HH^{≥1} = 0` in every characteristic (`A_2` is a tree, no oriented cycles). Hence **`HH^0 = 1`, `HH^1 = 0`, `HH^2 = 0`** over any field — the acceptance oracle. Dimensions are characteristic-independent (`−1` never vanishes; no coefficient collapses), so the same over `CC` and `GF(32003)`.

### Fixture 3 — completion that ADDS an S-polynomial (f.d., dim 4)
`Q = Quiver([1], {"x": (1,1), "y": (1,1)})` (two loops, order `x<y`), `relations=["y^2", "x*y", "y*x - x*x"]`.
Initial rules: `yy → 0`; `xy → 0`; and for `yx − xx` (both length 2, `xx < yx` since `x<y`) `LM = yx` ⇒ `yx → xx`. Start `R0 = {yy→0, xy→0, yx→xx}` is **not confluent**. The overlap of `w1=yy` (suffix `y`) with `w2=yx` (prefix `y`): `a=y, b=y, c=x`, ambiguity word `W = yyx`.
`S = tail(yy)⌢x − y⌢tail(yx) = (0) − y⌢(xx) = −yxx`. Reduce `−yxx` by `R0`: `yxx` has factor `yx` at position 0 → `xx⌢x = xxx`; `xxx` is irreducible under `R0`. Normal form `g = −xxx ≠ 0` ⇒ **add** rule `xxx → 0`.
With `R1 = {yy→0, xy→0, yx→xx, xxx→0}` every remaining ambiguity reduces to `0` (checked by hand): `(yy,yy)→yyy→0`; `(yy,yx)→ −yxx → −xxx → 0`; `(xy,yy)→xyy→0`; `(xy,yx)→xyx→0`; `(yx,xy)→yxy→xxy→0`; `(yx,xxx)→yxxx→0`; `(xxx,·)` self and cross all `→0`. So completion terminates at `R1`, one rule added. `L=3`, `2L−1 = 5 ≤ D_default` (= `max(8, 2·2+4) = 8`).
Irreducible words avoiding `yy, xy, yx, xxx`: length 1 `x, y`; length 2 `xx` only (`xy,yx,yy` forbidden); length 3 none (`xxx` forbidden, `xxy` has factor `xy`). **Basis** `[e_1, x, y, xx]`, labels `["e_1","x","y","x*x"]`, **dim 4**.

### Fixture 4 — two loops, no relations ⇒ `NotFiniteDimensionalError`
`Q = Quiver([1], {"a": (1,1), "b": (1,1)})` (a single vertex with a doubled arrow — the free algebra `k⟨a,b⟩`), `relations=[]`.
Completion of the empty rule set adds nothing (`L=0`, `2L−1 = −1 ≤ D`), so the degree-bound check passes; the certificate then enumerates irreducible words = **all** words in `a, b`. The automaton has a cycle at vertex 1 (self-loop `a`), so `irreducible_paths` raises `NotFiniteDimensionalError` naming the doubled-arrow growth (the loop). *Resolution note:* the literal Kronecker quiver `1 ⇉ 2` is finite-dimensional (path algebra dim 4); the intended infinite "double arrow, no relations" example is the two-loop quiver, whose two parallel loops generate arbitrarily long words `a, aa, ab, ba, …`. That is what is tested here. This exercises the Gröbner certificate directly with an empty rule set.

### Fixture 5 — bound exceeded ⇒ `AdmissibilityError` (distinct from Fixture 4)
`Q = Quiver([1], {"x": (1,1), "y": (1,1)})`, `relations=["x^2", "y^2", "y*x - x*y"]` (the quantum complete intersection with `q = 1`), built with an explicit tiny **`degree_bound=2`**.
Rules from the relations: `xx→0`, `yy→0`, and `yx − xy` (`xy < yx`) ⇒ `yx→xy`. All leading words have length `L = 2`, so certifying confluence needs ambiguity words up to length `2L−1 = 3` (e.g. `yyx`, `yxx`, `xxx`), which exceed `degree_bound = 2` and are never formed. The degree-bound check `2L−1 = 3 > 2` fires ⇒ **`AdmissibilityError`** telling the user to raise `degree_bound` to at least 3. This is *not* `NotFiniteDimensionalError`: the algebra is perfectly finite-dimensional; the *bound* is too small. With the default bound it certifies fine (Fixture 6).

### Fixture 6 — same algebra, default bound ⇒ certified, dim 4 (bonus; confluence + char-independence)
`relations=["x^2", "y^2", "y*x - x*y"]`, default `degree_bound`. Completion of `{xx→0, yy→0, yx→xy}` adds nothing: `(xx,xx)→xxx→0`, `(yy,yy)→yyy→0`, `(yy,yx)→ −yxy → 0`, `(yx,xx)→xyx→0` (via `yx→xy` then `xx→0`), etc. — all S-polynomials reduce to `0` by hand. Irreducible words avoiding `xx,yy,yx`: `x, y, xy`. **Basis** `[e_1, x, y, xy]`, labels `["e_1","x","y","x*y"]`, **dim 4**; in the algebra `y·x = x·y = xy`. Dimension is characteristic-independent (same over `CC`, `GF(2)`, `GF(32003)`).

---

### Task 1: `groebner` package skeleton, admissible `PathOrder`, step-event dataclasses

**Files:**
- Create: `src/quiverlab/groebner/__init__.py`, `src/quiverlab/groebner/order.py`, `src/quiverlab/groebner/events.py`, `tests/groebner/__init__.py` (empty), `tests/groebner/test_order.py`

**Interfaces:**
- Consumes: `Quiver`.
- Produces: `PathOrder` (frozen dataclass: `arrow_index: dict[str,int]`, `.key(word) -> tuple`, `.compare(u,v) -> int`, `.leading(comb) -> tuple|None`, `.sort_words(iterable) -> list`); factory `path_order(quiver) -> PathOrder` (arrow ranks = insertion order in `quiver.arrows`); plain dataclasses `Dispatch(route:str, reason:str, n_relations:int)` and `ReductionStep(word, rule_lead, before, after)` (inert trace hooks; formal trace wiring is Plan 07). All exported from `quiverlab.groebner`.

- [ ] **Step 1: Write the failing test**

`tests/groebner/__init__.py`: empty file.

`tests/groebner/test_order.py`:

```python
"""Length-lexicographic admissible path order (spec §5, component 3)."""
import itertools

from quiverlab.combinat import Quiver
from quiverlab.groebner import PathOrder, path_order, Dispatch, ReductionStep


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_arrow_ranks_follow_insertion_order():
    o = path_order(_square())
    assert o.arrow_index == {"a": 0, "b": 1, "c": 2, "d": 3}


def test_shorter_word_is_smaller():
    o = path_order(_square())
    assert o.compare(("a",), ("a", "b")) == -1
    assert o.compare(("a", "b"), ("a",)) == 1


def test_lex_breaks_ties_by_arrow_rank():
    o = path_order(_square())
    # equal length: a<c so a*b < c*d
    assert o.compare(("a", "b"), ("c", "d")) == -1
    assert o.compare(("c", "d"), ("a", "b")) == 1
    assert o.compare(("a", "b"), ("a", "b")) == 0


def test_leading_is_the_maximum_word():
    o = path_order(Quiver([1], {"x": (1, 1), "y": (1, 1)}))
    # y*x vs x*x with x<y: leading is y*x
    assert o.leading({("y", "x"): 1, ("x", "x"): 1}) == ("y", "x")
    # longer beats shorter regardless of lex
    assert o.leading({("x", "x", "x"): 1, ("y", "x"): 1}) == ("x", "x", "x")
    assert o.leading({}) is None


def test_two_sided_multiplicative_compatibility():
    """u < v  ==>  w u x < w v x for all composable contexts (checked combinatorially)."""
    o = path_order(Quiver([1], {"x": (1, 1), "y": (1, 1)}))
    letters = ["x", "y"]
    words = [tuple(p) for n in range(1, 4) for p in itertools.product(letters, repeat=n)]
    for u in words:
        for v in words:
            if o.compare(u, v) >= 0:
                continue
            for w in [()] + [(l,) for l in letters]:
                for x in [()] + [(l,) for l in letters]:
                    assert o.compare(w + u + x, w + v + x) < 0


def test_events_are_plain_dataclasses():
    d = Dispatch(route="groebner", reason="non-monomial", n_relations=1)
    assert (d.route, d.n_relations) == ("groebner", 1)
    s = ReductionStep(word=("c", "d"), rule_lead=("c", "d"),
                      before={("c", "d"): 1}, after={("a", "b"): 1})
    assert s.rule_lead == ("c", "d")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_order.py -q`
Expected: FAIL / error — `ModuleNotFoundError: No module named 'quiverlab.groebner'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/groebner/order.py`:

```python
"""Length-lexicographic admissible path order (spec §5, component 3).

Fix the total order on arrows = their INSERTION ORDER in Quiver.arrows: the
arrow at position 0 is smallest. For arrow-words (paths read LEFT TO RIGHT,
Assem-Simson-Skowronski) u, v:

    u < v  iff  len(u) < len(v),
             or len(u) == len(v) and, at the first index i where they differ,
                arrow_index[u[i]] < arrow_index[v[i]].

Empty words (vertices) have length 0 and are the minimum.

Admissibility -- this is a WELL-FOUNDED, TWO-SIDED MULTIPLICATIVE order:

 * Total: any two words compare by (length, then finite lex on ranks).
 * Well-founded: length is a nonnegative integer bounded below by 0; along any
   strictly descending chain the length is non-increasing, so it stabilizes,
   after which the words live in the finite set of words of that fixed length,
   on which lex is a well-order. Hence no infinite descending chain -- so
   reduction by any rule (which replaces a word by strictly smaller words)
   terminates.
 * Two-sided multiplicative: if u < v and w*u*x, w*v*x are both composable
   paths, then w*u*x < w*v*x. Indeed len(w*u*x) - len(w*v*x) = len(u) - len(v)
   <= 0; if the lengths differ the shorter wins directly; if len(u) == len(v)
   the common prefix w matches, the first differing index of u,v reappears
   (shifted by len(w)) as the first differing index of w*u*x, w*v*x with the
   SAME rank comparison, and the common suffix x is irrelevant.

This two-sided compatibility is exactly what makes leading words multiplicative,
so Bergman's Diamond Lemma applies: a reduction system all of whose ambiguities
resolve is confluent, and the irreducible words form a k-basis of kQ/I.

CONVENTION: because paths compose LEFT TO RIGHT, an overlap of leading words
w1, w2 is a word b that is a SUFFIX of w1 and a PREFIX of w2 (w1 = a*b,
w2 = b*c, ambiguity word a*b*c = w1*c = a*w2). Right-to-left treatments (most of
the Groebner/diamond-lemma literature) swap the roles of prefix and suffix; we
write every overlap, inclusion, S-polynomial and reduction in the left-to-right
convention.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class PathOrder:
    arrow_index: dict

    def key(self, word):
        """Total-order sort key: (length, ranks). Python tuple order == the order above."""
        return (len(word), tuple(self.arrow_index[a] for a in word))

    def compare(self, u, v):
        ku, kv = self.key(u), self.key(v)
        if ku < kv:
            return -1
        if ku > kv:
            return 1
        return 0

    def leading(self, comb):
        """The largest word appearing in the linear combination comb (word -> coeff),
        or None if comb is empty (the zero element). Callers keep comb free of
        zero coefficients, so every key is a genuine support word."""
        if not comb:
            return None
        return max(comb, key=self.key)

    def sort_words(self, words):
        """Ascending sort by the order (used for deterministic iteration)."""
        return sorted(words, key=self.key)


def path_order(quiver):
    """The admissible order whose arrow ranks are the insertion order of quiver.arrows."""
    return PathOrder(arrow_index={name: i for i, name in enumerate(quiver.arrows)})
```

`src/quiverlab/groebner/events.py`:

```python
"""Plain-dataclass step events emitted by the Groebner lowering, as trace hooks.

These are INERT: when a caller passes trace=[...] the lowering appends these
records; nothing else consumes them yet. The formal trace subsystem (typed
events, PDF/HTML/text renderers, eliding rules, golden-file tests) is Plan 07 --
this is only the emission boundary.
"""
from dataclasses import dataclass


@dataclass
class Dispatch:
    """Which lowering route Quiver.algebra took, and why."""
    route: str          # "monomial" | "groebner"
    reason: str
    n_relations: int


@dataclass
class ReductionStep:
    """One rewrite: the word occurrence reduced, the rule's leading word, and the
    linear combination (word -> domain element) before and after the step."""
    word: tuple
    rule_lead: tuple
    before: dict
    after: dict
```

`src/quiverlab/groebner/__init__.py`:

```python
"""quiverlab.groebner: noncommutative Groebner completion for path-algebra ideals
(spec §5, component 3). Emits the ReductionSystem consumed by Chouhy-Solotar
(Plan 04). Path composition is LEFT TO RIGHT throughout."""
from quiverlab.groebner.order import PathOrder, path_order  # noqa: F401
from quiverlab.groebner.events import Dispatch, ReductionStep  # noqa: F401

__all__ = ["PathOrder", "path_order", "Dispatch", "ReductionStep"]
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_order.py tests/test_no_floats.py -q`
Expected: all pass (order tests + the float-ban gate stays green).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/groebner tests/groebner
git commit -m "feat(groebner): admissible length-lex path order + step-event hooks

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Linear-combination arithmetic, reduction rules, and reduction to normal form

**Files:**
- Create: `src/quiverlab/groebner/reduction.py`, `tests/groebner/test_reduction.py`

**Interfaces:**
- Consumes: `PathOrder`, `ReductionStep`, `Domain`, `Quiver`.
- Produces (element = `dict[tuple[str,...], <domain element>]`, zero = `{}`, no zero coefficients kept):
  - `lc_add(comb, word, coeff, dom) -> None` (in place), `lc_sub(a, b, dom) -> dict`, `lc_scale_left/right` not needed.
  - `ReductionRule` (frozen dataclass): `.lead: tuple[str,...]`, `.tail: tuple[tuple[<domain element>, tuple[str,...]], ...]` (each word `< lead`), `.source`, `.target`.
  - `rule_from_comb(comb, order, quiver, dom) -> ReductionRule`.
  - `first_factor(word, rules) -> tuple[ReductionRule, int] | None` (leftmost occurrence of any rule's `lead`; ties broken by smallest position then rule order).
  - `reduce_comb(comb, rules, order, dom, trace=None) -> dict` (normal form; terminates by well-foundedness; deterministic — always rewrites the largest reducible word first).

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_reduction.py`:

```python
"""Reduction rules and normal forms over a Domain, left-to-right (spec §5 c.3)."""
from fractions import Fraction

from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import (
    ReductionRule, rule_from_comb, first_factor, reduce_comb, lc_sub,
)


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _c(dom, n):
    return dom.coerce(n)


def test_rule_from_monomial_comb_has_empty_tail():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    r = rule_from_comb({("x", "x", "x"): _c(dom, 1)}, o, Q, dom)
    assert r.lead == ("x", "x", "x")
    assert r.tail == ()
    assert (r.source, r.target) == (1, 1)


def test_rule_from_binomial_picks_leading_and_normalizes():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    # relation y*x - x*x  ==  {yx: 1, xx: -1};  x<y so leading is yx, tail = {xx: 1}
    comb = {("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}
    r = rule_from_comb(comb, o, Q, dom)
    assert r.lead == ("y", "x")
    assert r.tail == ((_c(dom, 1), ("x", "x")),)


def test_rule_normalizes_leading_coefficient_to_one():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    # 2*yx - 3*xx : leading yx (coeff 2), tail = (3/2) xx
    comb = {("y", "x"): _c(dom, 2), ("x", "x"): _c(dom, -3)}
    r = rule_from_comb(comb, o, Q, dom)
    assert r.lead == ("y", "x")
    assert r.tail == ((Fraction(3, 2), ("x", "x")),)


def test_first_factor_finds_leftmost_occurrence():
    Q = _two_loops()
    dom = QQ
    r = rule_from_comb({("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}, path_order(Q), Q, dom)
    hit = first_factor(("y", "x", "x"), [r])
    assert hit is not None
    rule, i = hit
    assert (rule.lead, i) == (("y", "x"), 0)
    assert first_factor(("x", "x"), [r]) is None


def test_reduce_single_monomial_rule():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    xxx = rule_from_comb({("x", "x", "x"): _c(dom, 1)}, o, Q, dom)  # xxx -> 0
    assert reduce_comb({("x", "x", "x", "x"): _c(dom, 1)}, [xxx], o, dom) == {}
    assert reduce_comb({("x", "x"): _c(dom, 1)}, [xxx], o, dom) == {("x", "x"): _c(dom, 1)}


def test_reduce_chains_rules_to_normal_form():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    yx = rule_from_comb({("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}, o, Q, dom)  # yx -> xx
    xxx = rule_from_comb({("x", "x", "x"): _c(dom, 1)}, o, Q, dom)                     # xxx -> 0
    # y*x*x  --yx-->  x*x*x  --xxx-->  0
    assert reduce_comb({("y", "x", "x"): _c(dom, 1)}, [yx, xxx], o, dom) == {}


def test_reduce_emits_trace_steps_when_requested():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    yx = rule_from_comb({("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}, o, Q, dom)
    trace = []
    reduce_comb({("y", "x"): _c(dom, 1)}, [yx], o, dom, trace=trace)
    assert len(trace) == 1
    assert trace[0].rule_lead == ("y", "x")


def test_lc_sub_cancels():
    dom = QQ
    a = {("a", "b"): _c(dom, 1)}
    b = {("a", "b"): _c(dom, 1)}
    assert lc_sub(a, b, dom) == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_reduction.py -q`
Expected: FAIL / error — `ModuleNotFoundError: No module named 'quiverlab.groebner.reduction'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/groebner/reduction.py`:

```python
"""Linear-combination arithmetic, reduction rules, and normal forms over a Domain.

An ELEMENT of kQ is a dict  word -> domain element  with no zero coefficients;
the zero element is {}. A word is a tuple of arrow names read LEFT TO RIGHT.
Reduction rewrites an occurrence of a rule's leading word W = u*lead*v into
u*tail*v; every tail word is strictly smaller than lead in the admissible order,
so reduction terminates (spec §5, component 3)."""
from dataclasses import dataclass

from quiverlab.groebner.events import ReductionStep


def lc_add(comb, word, coeff, dom):
    """comb[word] += coeff, in place; drop the entry if it becomes zero."""
    if dom.is_zero(coeff):
        return
    cur = comb.get(word)
    total = coeff if cur is None else dom.add(cur, coeff)
    if dom.is_zero(total):
        comb.pop(word, None)
    else:
        comb[word] = total


def lc_sub(a, b, dom):
    """Return a - b as a fresh element."""
    out = dict(a)
    for w, c in b.items():
        lc_add(out, w, dom.neg(c), dom)
    return out


@dataclass(frozen=True)
class ReductionRule:
    """lead -> tail: tail is a linear combination of words strictly SMALLER than
    lead, parallel to it (same source & target), with lead - sum(tail) in I.
    Stored as tail = ((coeff, word), ...) sorted descending by the order."""
    lead: tuple
    tail: tuple
    source: object
    target: object


def rule_from_comb(comb, order, quiver, dom):
    """Turn a nonzero element into a rule lead -> tail by dividing out the leading
    coefficient: lead = LM(comb), tail = -(1/LC)*(comb - LC*lead)."""
    lead = order.leading(comb)
    inv = dom.inv(comb[lead])
    tail_items = [(dom.neg(dom.mul(inv, c)), w) for w, c in comb.items() if w != lead]
    tail_items.sort(key=lambda cw: order.key(cw[1]), reverse=True)
    return ReductionRule(
        lead=lead,
        tail=tuple(tail_items),
        source=quiver.word_source(lead),
        target=quiver.word_target(lead),
    )


def first_factor(word, rules):
    """Leftmost occurrence of any rule's lead as a contiguous factor of word.
    Returns (rule, position) with the smallest position (ties: first such rule),
    or None if word is irreducible."""
    best = None
    for rule in rules:
        L = rule.lead
        n = len(L)
        for i in range(len(word) - n + 1):
            if word[i:i + n] == L:
                if best is None or i < best[1]:
                    best = (rule, i)
                break
    return best


def reduce_comb(comb, rules, order, dom, trace=None):
    """Normal form of comb under rules. Deterministic: at each step rewrite the
    LARGEST (in the order) reducible word. Terminates because every rewrite
    replaces a word by strictly smaller ones in a well-founded order; for a
    confluent system the result is independent of these choices."""
    work = {w: c for w, c in comb.items() if not dom.is_zero(c)}
    while True:
        target_word = None
        for w in sorted(work, key=order.key, reverse=True):
            if first_factor(w, rules) is not None:
                target_word = w
                break
        if target_word is None:
            return work
        rule, i = first_factor(target_word, rules)
        coeff = work.pop(target_word)
        u, v = target_word[:i], target_word[i + len(rule.lead):]
        if trace is not None:
            before = dict(work)
            before[target_word] = coeff
        for tc, tw in rule.tail:
            lc_add(work, u + tw + v, dom.mul(coeff, tc), dom)
        if trace is not None:
            trace.append(ReductionStep(word=target_word, rule_lead=rule.lead,
                                       before=before, after=dict(work)))
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_reduction.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/groebner/reduction.py tests/groebner/test_reduction.py
git commit -m "feat(groebner): linear-combination arithmetic, reduction rules, normal forms

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Overlap and inclusion ambiguity detection (left-to-right)

**Files:**
- Create: `src/quiverlab/groebner/overlap.py`, `tests/groebner/test_overlap.py`

**Interfaces:**
- Consumes: `ReductionRule`.
- Produces:
  - `Ambiguity` (frozen dataclass): `.kind: str` in `{"overlap","inclusion"}`, `.word: tuple[str,...]` (the ambiguity word containing both leads), `.left: ReductionRule`, `.right: ReductionRule`, `.a: tuple[str,...]`, `.c: tuple[str,...]`. Overlap: `left.lead = a⌢b`, `right.lead = b⌢c`, `word = a⌢b⌢c`. Inclusion: `left` is the outer (containing) rule, `right` the inner, `word = left.lead = a⌢right.lead⌢c`.
  - `overlaps(r1, r2, degree_bound) -> list[Ambiguity]` (all `b` = proper suffix of `r1.lead` = proper prefix of `r2.lead`, ambiguity word length `≤ degree_bound`).
  - `inclusions(r_outer, r_inner, degree_bound) -> list[Ambiguity]` (`r_inner.lead` a proper contiguous factor of `r_outer.lead`).
  - `all_ambiguities(rules, degree_bound) -> list[Ambiguity]` (all ordered pairs incl. a rule with itself).

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_overlap.py`:

```python
"""Overlap/inclusion ambiguities in the LEFT-TO-RIGHT convention (spec §5 c.3):
b is a SUFFIX of the first lead and a PREFIX of the second."""
from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import rule_from_comb
from quiverlab.groebner.overlap import overlaps, inclusions, all_ambiguities


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _rule(comb):
    Q = _two_loops()
    return rule_from_comb(comb, path_order(Q), Q, QQ)


def test_self_overlap_of_xyx():
    r = _rule({("x", "y", "x"): QQ.coerce(1), ("y",): QQ.coerce(1)})  # lead xyx (len3 > y)
    amb = overlaps(r, r, degree_bound=10)
    words = sorted(a.word for a in amb)
    # xyx suffix "x" = prefix "x" -> a=xy, c=yx -> word xyxyx (only overlap)
    assert words == [("x", "y", "x", "y", "x")]
    a0 = amb[0]
    assert a0.kind == "overlap"
    assert (a0.a, a0.c) == (("x", "y"), ("y", "x"))


def test_overlap_yy_with_yx():
    yy = _rule({("y", "y"): QQ.coerce(1)})               # yy -> 0
    yx = _rule({("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)})  # yx -> xx
    amb = overlaps(yy, yx, degree_bound=10)
    assert [a.word for a in amb] == [("y", "y", "x")]
    assert (amb[0].a, amb[0].c) == (("y",), ("x",))


def test_no_overlap_when_suffix_prefix_disagree():
    cd = _rule({("y", "x"): QQ.coerce(1)})   # lead yx: suffix "x" != prefix "y"
    assert overlaps(cd, cd, degree_bound=10) == []


def test_degree_bound_filters_long_ambiguities():
    r = _rule({("x", "y", "x"): QQ.coerce(1), ("y",): QQ.coerce(1)})
    assert overlaps(r, r, degree_bound=4) == []   # xyxyx has length 5 > 4
    assert len(overlaps(r, r, degree_bound=5)) == 1


def test_inclusion_detects_proper_factor():
    inner = _rule({("x", "x"): QQ.coerce(1)})               # xx -> 0
    outer = _rule({("x", "x", "x"): QQ.coerce(1)})          # xxx -> 0
    inc = inclusions(outer, inner, degree_bound=10)
    # xx occurs in xxx at positions 0 and 1
    assert sorted((a.a, a.c) for a in inc) == [((), ("x",)), (("x",), ())]
    assert all(a.kind == "inclusion" and a.word == ("x", "x", "x") for a in inc)


def test_all_ambiguities_covers_ordered_pairs_including_self():
    yy = _rule({("y", "y"): QQ.coerce(1)})
    yx = _rule({("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)})
    ambs = all_ambiguities([yy, yx], degree_bound=10)
    got = {(a.left.lead, a.right.lead, a.word) for a in ambs}
    # (yy,yy)->yyy ; (yy,yx)->yyx ; (yx,yy)? suffix x != prefix y -> none
    assert (("y", "y"), ("y", "y"), ("y", "y", "y")) in got
    assert (("y", "y"), ("y", "x"), ("y", "y", "x")) in got
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_overlap.py -q`
Expected: FAIL / error — `ModuleNotFoundError: No module named 'quiverlab.groebner.overlap'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/groebner/overlap.py`:

```python
"""Overlap and inclusion ambiguities between leading words, LEFT TO RIGHT.

For leading words w1, w2:
 * OVERLAP: a nonempty b that is a proper SUFFIX of w1 and a proper PREFIX of w2.
   Then w1 = a*b, w2 = b*c, ambiguity word W = a*b*c = w1*c = a*w2. (Right-to-left
   treatments swap prefix/suffix; we are left-to-right per the quiver convention.)
 * INCLUSION: w_inner is a proper contiguous factor of w_outer, w_outer = a*w_inner*c
   with (a, c) != ((), ()); ambiguity word W = w_outer.
Only ambiguity words of length <= degree_bound are produced (spec §5, component 3).
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Ambiguity:
    kind: str            # "overlap" | "inclusion"
    word: tuple          # the ambiguity word containing both leading words
    left: object         # overlap: rule with lead a*b ; inclusion: the OUTER rule
    right: object        # overlap: rule with lead b*c ; inclusion: the INNER rule
    a: tuple             # left/prefix context
    c: tuple             # right/suffix context


def overlaps(r1, r2, degree_bound):
    w1, w2 = r1.lead, r2.lead
    out = []
    # b nonempty, strictly shorter than both w1 and w2 (proper suffix / proper prefix)
    for blen in range(1, min(len(w1), len(w2))):
        b = w1[len(w1) - blen:]
        if b != w2[:blen]:
            continue
        a = w1[:len(w1) - blen]
        c = w2[blen:]
        word = a + b + c                      # == w1 + c
        if len(word) <= degree_bound:
            out.append(Ambiguity(kind="overlap", word=word, left=r1, right=r2, a=a, c=c))
    return out


def inclusions(r_outer, r_inner, degree_bound):
    W, w = r_outer.lead, r_inner.lead
    if len(w) >= len(W):                      # inner must be a PROPER factor
        return []
    if len(W) > degree_bound:
        return []
    out = []
    for i in range(len(W) - len(w) + 1):
        if W[i:i + len(w)] == w:
            a, c = W[:i], W[i + len(w):]
            if (a, c) == ((), ()):
                continue
            out.append(Ambiguity(kind="inclusion", word=W, left=r_outer, right=r_inner, a=a, c=c))
    return out


def all_ambiguities(rules, degree_bound):
    out = []
    for r1 in rules:
        for r2 in rules:
            out.extend(overlaps(r1, r2, degree_bound))
            if r1 is not r2:
                out.extend(inclusions(r1, r2, degree_bound))
    return out
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_overlap.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/groebner/overlap.py tests/groebner/test_overlap.py
git commit -m "feat(groebner): left-to-right overlap and inclusion ambiguity detection

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: S-polynomials and Buchberger–Mora completion with a degree bound

**Files:**
- Create: `src/quiverlab/groebner/complete.py`, `tests/groebner/test_complete.py`

**Interfaces:**
- Consumes: `Ambiguity`, `ReductionRule`, `reduce_comb`, `rule_from_comb`, `all_ambiguities`, `PathOrder`, `Domain`, `AdmissibilityError`.
- Produces:
  - `s_polynomial(amb, dom) -> dict` (the S-polynomial element, per the two reductions of the ambiguity word).
  - `complete(init_rules, order, quiver, dom, degree_bound, trace=None, max_rules=20000) -> list[ReductionRule]` (Buchberger–Mora; always terminates under a finite `degree_bound`; raises `AdmissibilityError` if `max_rules` is exceeded, naming the runaway growth).

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_complete.py`:

```python
"""S-polynomials + Buchberger-Mora completion (spec §5, component 3), left-to-right."""
from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import rule_from_comb
from quiverlab.groebner.overlap import overlaps
from quiverlab.groebner.complete import s_polynomial, complete


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _rules(combs, Q):
    o = path_order(Q)
    return [rule_from_comb(c, o, Q, QQ) for c in combs]


def test_s_polynomial_of_yy_and_yx():
    Q = _two_loops()
    yy, yx = _rules([{("y", "y"): QQ.coerce(1)},
                     {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)}], Q)
    amb = overlaps(yy, yx, degree_bound=10)[0]   # word yyx
    s = s_polynomial(amb, QQ)
    # tail(yy)=0 so front reduction is 0; back is y*tail(yx)=y*xx=yxx ; S = 0 - yxx
    assert s == {("y", "x", "x"): QQ.coerce(-1)}


def test_completion_fixture3_adds_xxx():
    """Fixture 3: {yy->0, xy->0, yx->xx} completes by ADDING xxx->0."""
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("y", "y"): QQ.coerce(1)},                                   # yy -> 0
        {("x", "y"): QQ.coerce(1)},                                   # xy -> 0
        {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)},        # yx -> xx
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    leads = sorted(r.lead for r in done)
    assert leads == [("x", "x", "x"), ("x", "y"), ("y", "x"), ("y", "y")]
    # the added rule is a genuine monomial rule xxx -> 0
    xxx = next(r for r in done if r.lead == ("x", "x", "x"))
    assert xxx.tail == ()


def test_completion_fixture6_quantum_ci_adds_nothing():
    """Fixture 6: {xx->0, yy->0, yx->xy} is already confluent."""
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("x", "x"): QQ.coerce(1)},                                   # xx -> 0
        {("y", "y"): QQ.coerce(1)},                                   # yy -> 0
        {("y", "x"): QQ.coerce(1), ("x", "y"): QQ.coerce(-1)},        # yx -> xy
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    assert sorted(r.lead for r in done) == [("x", "x"), ("y", "x"), ("y", "y")]


def test_completion_square_adds_nothing():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    o = path_order(Q)
    cd = rule_from_comb({("c", "d"): QQ.coerce(1), ("a", "b"): QQ.coerce(-1)}, o, Q, QQ)
    done = complete([cd], o, Q, QQ, degree_bound=8)
    assert [r.lead for r in done] == [("c", "d")]
    assert done[0].tail == ((QQ.coerce(1), ("a", "b")),)


def test_completion_terminates_and_is_confluent():
    """After completion, every ambiguity S-polynomial reduces to zero."""
    from quiverlab.groebner.overlap import all_ambiguities
    from quiverlab.groebner.reduction import reduce_comb
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("y", "y"): QQ.coerce(1)},
        {("x", "y"): QQ.coerce(1)},
        {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)},
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    for amb in all_ambiguities(done, degree_bound=8):
        assert reduce_comb(s_polynomial(amb, QQ), done, o, QQ) == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_complete.py -q`
Expected: FAIL / error — `ModuleNotFoundError: No module named 'quiverlab.groebner.complete'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/groebner/complete.py`:

```python
"""S-polynomials and Buchberger-Mora completion, LEFT TO RIGHT (spec §5, c.3).

For an ambiguity word W the S-polynomial reduces W two ways and subtracts:
 * OVERLAP (w1=a*b, w2=b*c, W=a*b*c): apply the w1 rule at the front (tail(w1)*c)
   minus apply the w2 rule at the back (a*tail(w2)).
 * INCLUSION (W=w_outer=a*w_inner*c): apply the outer rule (tail(w_outer)) minus
   apply the inner rule inside (a*tail(w_inner)*c).
Completion reduces each S-polynomial by the current rules; a nonzero normal form
yields a new rule. Only ambiguity words of length <= degree_bound are formed, and
there are finitely many words of length <= degree_bound over a finite arrow set,
so completion ALWAYS TERMINATES under a fixed bound."""
from quiverlab.errors import AdmissibilityError
from quiverlab.groebner.overlap import all_ambiguities
from quiverlab.groebner.reduction import lc_add, lc_sub, reduce_comb, rule_from_comb


def _apply_tail(rule, prefix, suffix, dom):
    """The element prefix * tail(rule) * suffix (each tail word wrapped in context)."""
    out = {}
    for tc, tw in rule.tail:
        lc_add(out, prefix + tw + suffix, tc, dom)
    return out


def s_polynomial(amb, dom):
    if amb.kind == "overlap":
        front = _apply_tail(amb.left, (), amb.c, dom)     # tail(w1) * c
        back = _apply_tail(amb.right, amb.a, (), dom)      # a * tail(w2)
        return lc_sub(front, back, dom)
    outer = _apply_tail(amb.left, (), (), dom)             # tail(w_outer)
    inner = _apply_tail(amb.right, amb.a, amb.c, dom)      # a * tail(w_inner) * c
    return lc_sub(outer, inner, dom)


def _is_factor(short, long):
    """True if `short` occurs as a contiguous factor of `long`."""
    n = len(short)
    return any(long[i:i + n] == short for i in range(len(long) - n + 1))


def _minimize_leads(rules):
    """Drop rules whose leading word contains ANOTHER rule's leading word as a
    proper factor (such a lead is already reducible, so the rule is redundant), and
    drop exact-duplicate leads (keep the first). The remaining leads form an
    antichain under the factor relation, each appearing once. This does not change
    the ideal, the set of irreducible words (forbidden = leads; a subsumed forbidden
    word is redundant), or confluence -- it is the standard reduced-basis cleanup,
    it keeps max-leading-length honest for the degree-bound check, and it keeps
    leading_words()/ambiguities() (the CS S-sequence seed) free of duplicates."""
    minimal, seen = [], set()
    for r in rules:
        if r.lead in seen:                          # exact-duplicate lead: keep the first
            continue
        if any(s.lead != r.lead and len(s.lead) < len(r.lead) and _is_factor(s.lead, r.lead)
               for s in rules):
            continue
        minimal.append(r)
        seen.add(r.lead)
    return minimal


def complete(init_rules, order, quiver, dom, degree_bound, trace=None, max_rules=20000):
    rules = list(init_rules)
    added = True
    while added:
        added = False
        for amb in all_ambiguities(rules, degree_bound):
            g = reduce_comb(s_polynomial(amb, dom), rules, order, dom, trace=trace)
            if g:                                          # nonzero normal form
                rules.append(rule_from_comb(g, order, quiver, dom))
                added = True
                if len(rules) > max_rules:
                    raise AdmissibilityError(
                        f"Groebner completion produced more than {max_rules} rules under "
                        f"degree_bound={degree_bound}: the reduction system is not stabilizing",
                        hint="the ideal may not be admissible; inspect the relations, or this "
                             "presentation is beyond the v1 certificate",
                    )
                break                                      # restart with the enlarged system
    return _minimize_leads(rules)
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_complete.py tests/test_no_floats.py -q`
Expected: all pass (Fixture 3 adds `xxx→0`; Fixtures 6 and 2 add nothing; post-completion confluence holds).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/groebner/complete.py tests/groebner/test_complete.py
git commit -m "feat(groebner): S-polynomials and Buchberger-Mora completion with degree bound

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Finiteness certificate — degree-bound check, irreducible basis, loud errors

**Files:**
- Create: `src/quiverlab/groebner/certificate.py`, `tests/groebner/test_certificate.py`

**Interfaces:**
- Consumes: completed `list[ReductionRule]`, `Quiver`, `core.monomial.irreducible_paths`, `AdmissibilityError`, `NotFiniteDimensionalError`.
- Produces:
  - `default_degree_bound(relations) -> int` (`max(8, 2*max_relation_length + 4)`).
  - `check_degree_bound(rules, degree_bound) -> None` (raise `AdmissibilityError` unless `2*L−1 ≤ degree_bound`, where `L = max leading-word length`).
  - `certified_irreducibles(quiver, rules) -> tuple[tuple[str,...], ...]` (nonempty irreducible words sorted `(len, word)` via the Plan-01 automaton; raises `NotFiniteDimensionalError` with an arrow cycle when infinite).

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_certificate.py`:

```python
"""Finiteness certificate: degree-bound check + irreducible-word enumeration
(spec §3.3, §5 component 3). Loud AdmissibilityError / NotFiniteDimensionalError."""
import pytest

from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import rule_from_comb
from quiverlab.groebner.complete import complete
from quiverlab.groebner.certificate import (
    default_degree_bound, check_degree_bound, certified_irreducibles,
)


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _rules(combs, Q):
    o = path_order(Q)
    return [rule_from_comb(c, o, Q, QQ) for c in combs]


def test_default_bound_covers_completion_growth():
    # Fixture 3 relations have max length 2 -> default 8, and completion reaches L=3 (need 5)
    class _R:
        max_length = 2
    assert default_degree_bound([_R(), _R()]) == 8


def test_check_degree_bound_passes_when_bound_large_enough():
    Q = _two_loops()
    rules = _rules([{("x", "x"): QQ.coerce(1)}], Q)   # L=2, need 2*2-1=3
    check_degree_bound(rules, 8)                       # no raise


def test_check_degree_bound_fixture5_raises_admissibility():
    """Fixture 5: quantum CI q=1 with degree_bound=2. L=2 needs bound >= 3."""
    Q = _two_loops()
    rules = _rules([
        {("x", "x"): QQ.coerce(1)},
        {("y", "y"): QQ.coerce(1)},
        {("y", "x"): QQ.coerce(1), ("x", "y"): QQ.coerce(-1)},
    ], Q)
    with pytest.raises(AdmissibilityError) as exc:
        check_degree_bound(rules, 2)
    assert "degree_bound=2" in str(exc.value)
    assert "at least 3" in str(exc.value)


def test_certified_irreducibles_fixture3():
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("y", "y"): QQ.coerce(1)},
        {("x", "y"): QQ.coerce(1)},
        {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)},
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    words = certified_irreducibles(Q, done)
    assert words == (("x",), ("y",), ("x", "x"))        # basis {x, y, xx}, dim 4 with vertex


def test_certified_irreducibles_fixture4_raises_not_finite():
    """Fixture 4: two loops, no relations -> infinite, names the loop cycle."""
    Q = Quiver([1], {"a": (1, 1), "b": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError) as exc:
        certified_irreducibles(Q, [])
    msg = str(exc.value)
    assert "infinite" in msg.lower()
    assert ("a" in msg or "b" in msg)                   # the doubled-arrow growth is named


def test_fixture5_and_fixture4_are_distinct_error_types():
    assert AdmissibilityError is not NotFiniteDimensionalError
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_certificate.py -q`
Expected: FAIL / error — `ModuleNotFoundError: No module named 'quiverlab.groebner.certificate'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/groebner/certificate.py`:

```python
"""Finiteness certificate for a completed reduction system (spec §3.3, §5 c.3).

After completion under a degree bound D, let L = max leading-word length. Every
overlap ambiguity word has length <= 2L-1 and every inclusion word length <= L,
so if 2L-1 <= D EVERY ambiguity was formed and reduced to 0 -- the system is
confluent and the irreducible words are a genuine k-basis of kQ/I. If 2L-1 > D
we refuse to certify (AdmissibilityError: raise the bound). Given a confluent
system we enumerate irreducible words with the Plan-01 automaton
(core.monomial.irreducible_paths): finitely many -> certified finite-dimensional;
a cycle -> NotFiniteDimensionalError naming the offending arrow cycle."""
from quiverlab.core.monomial import irreducible_paths
from quiverlab.errors import AdmissibilityError


def default_degree_bound(relations):
    """A CONSERVATIVE heuristic bound. If completion grows leading words past it,
    the certificate refuses to certify with a loud AdmissibilityError (never a wrong
    answer) telling the user to raise degree_bound; it never silently certifies an
    unconfluent system. Generous enough for every admissible v1 fixture."""
    m = max((r.max_length for r in relations), default=1)
    return max(8, 2 * m + 4)


def check_degree_bound(rules, degree_bound):
    L = max((len(r.lead) for r in rules), default=0)
    need = 2 * L - 1
    if need > degree_bound:
        raise AdmissibilityError(
            f"Groebner completion cannot be certified under degree_bound={degree_bound}: "
            f"leading words reach length {L}, so overlap ambiguities can have length up to "
            f"{need}, which exceeds the bound",
            hint=f"raise degree_bound to at least {need}",
        )


def certified_irreducibles(quiver, rules):
    """Nonempty irreducible words (avoiding every rule's leading word as a factor),
    sorted (len, word) exactly as core.monomial.irreducible_paths -- so the Groebner
    basis indexing agrees elementwise with the Plan-01 monomial route. Raises
    NotFiniteDimensionalError (naming an arrow cycle) when infinitely many exist."""
    forbidden = [r.lead for r in rules]
    return tuple(irreducible_paths(quiver, forbidden))
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_certificate.py tests/test_no_floats.py -q`
Expected: all pass (Fixture 5 → `AdmissibilityError`; Fixture 4 → `NotFiniteDimensionalError` naming a loop; Fixture 3 basis `{x, y, xx}`).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/groebner/certificate.py tests/groebner/test_certificate.py
git commit -m "feat(groebner): finiteness certificate with loud admissibility/finiteness errors

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: `ReductionSystem` — the object frozen for Chouhy–Solotar (Plan 04)

**Files:**
- Create: `src/quiverlab/groebner/system.py`, `tests/groebner/test_system.py`
- Modify: `src/quiverlab/groebner/__init__.py` (export `ReductionSystem`, `ReductionRule`, `Ambiguity`, `build_reduction_system`)

**Interfaces:**
- Consumes: everything above; `parse_relations` (only for a string-input convenience), `AdmissibilityError`, `NotFiniteDimensionalError`, `RelationError`, fields.
- Produces the **frozen** public shape (Plan 04 consumes it verbatim — see the frozen-interface block at the end of this plan):

```python
@dataclass(frozen=True)
class ReductionSystem:
    quiver: object                       # the Quiver
    domain: object                       # the exact Domain coefficients live in
    order: PathOrder                     # the admissible length-lex order
    rules: tuple                         # tuple[ReductionRule, ...], confluent, minimal leads
    irreducibles: tuple                  # tuple[tuple[str,...], ...], certified-finite nonempty
                                         #   irreducible words, sorted (len, word)
    degree_bound: int                    # the bound completion/certification used

    def leading_words(self) -> tuple: ...                 # tuple(r.lead for r in rules)
    def reduce(self, comb: dict) -> dict: ...             # normal form of a linear combination
    def normal_form(self, word: tuple) -> dict: ...       # reduce a single word
    def ambiguities(self) -> tuple: ...                   # tuple[Ambiguity, ...], the S-sequence
                                                          #   seed (all resolve to 0)
    is_confluent: bool = True            # always True for a built system
```

- `build_reduction_system(quiver, relations, field, degree_bound=None, trace=None) -> ReductionSystem` where `relations` is a list of parsed `Relation` objects (or strings; strings are parsed against `quiver`). Raises `AdmissibilityError` (length < 2, or bound exceeded), `RelationError` (coefficient vanishes in the field), `NotFiniteDimensionalError` (infinite basis).

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_system.py`:

```python
"""ReductionSystem: the object Chouhy-Solotar (Plan 04) consumes verbatim."""
import pytest
from dataclasses import fields as dc_fields

from quiverlab.combinat import Quiver
from quiverlab.fields import CC, QQ
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError, RelationError
from quiverlab.groebner import ReductionSystem, ReductionRule, build_reduction_system
from quiverlab.groebner.reduction import ReductionRule as _RR
from quiverlab.groebner.overlap import Ambiguity


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_frozen_public_shape():
    """Plan 04 depends on these field/method names -- freeze them."""
    names = {f.name for f in dc_fields(ReductionSystem)}
    assert {"quiver", "domain", "order", "rules", "irreducibles",
            "degree_bound", "is_confluent"} <= names
    for meth in ("leading_words", "reduce", "normal_form", "ambiguities"):
        assert callable(getattr(ReductionSystem, meth))
    assert ReductionRule is _RR


def test_build_square_system():
    rs = build_reduction_system(_square(), ["a*b - c*d"], CC)
    assert rs.leading_words() == (("c", "d"),)
    assert rs.irreducibles == (("a",), ("b",), ("c",), ("d",), ("a", "b"))
    assert rs.is_confluent is True


def test_reduce_and_normal_form_use_the_rules():
    rs = build_reduction_system(_square(), ["a*b - c*d"], CC)
    one = rs.domain.coerce(1)
    # c*d reduces to a*b
    assert rs.normal_form(("c", "d")) == {("a", "b"): one}
    assert rs.reduce({("c", "d"): one}) == {("a", "b"): one}


def test_ambiguities_of_completed_system_all_resolve():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rs = build_reduction_system(Q, ["y^2", "x*y", "y*x - x*x"], QQ)   # Fixture 3
    assert sorted(r.lead for r in rs.rules) == [
        ("x", "x", "x"), ("x", "y"), ("y", "x"), ("y", "y")]
    for amb in rs.ambiguities():
        assert isinstance(amb, Ambiguity)
        assert rs.reduce(_s(amb, rs)) == {}


def _s(amb, rs):
    from quiverlab.groebner.complete import s_polynomial
    return s_polynomial(amb, rs.domain)


def test_build_rejects_length_one_relation():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    with pytest.raises(AdmissibilityError):
        build_reduction_system(Q, ["a*b - c"], CC)   # c has length 1


def test_build_raises_admissibility_when_bound_too_small():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    with pytest.raises(AdmissibilityError):
        build_reduction_system(Q, ["x^2", "y^2", "y*x - x*y"], CC, degree_bound=2)  # Fixture 5


def test_build_raises_not_finite_dimensional():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError):
        build_reduction_system(Q, ["x*y - y*x"], CC)   # commutative k[x,y], infinite


def test_build_rejects_relation_vanishing_in_field():
    from quiverlab.fields import GF
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    # 2*x*y is a single monomial with coefficient 2; over GF(2) it vanishes.
    with pytest.raises(RelationError):
        build_reduction_system(Q, ["2*x*y"], GF(2))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_system.py -q`
Expected: FAIL / error — `ImportError: cannot import name 'ReductionSystem'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/groebner/system.py`:

```python
"""ReductionSystem: the certified confluent reduction system for kQ/I, and the
object Chouhy-Solotar (Plan 04, arXiv:1406.2300) consumes verbatim (spec §5 c.3,
§6). Pairs (leading word, reduction) + the admissible order + the ambiguity
S-sequence, satisfying CS's reduction-finiteness requirement (finitely many rules,
finitely many irreducible words, every ambiguity resolves). LEFT TO RIGHT."""
from dataclasses import dataclass

from quiverlab.errors import RelationError
from quiverlab.groebner.order import PathOrder, path_order
from quiverlab.groebner.reduction import ReductionRule, reduce_comb, rule_from_comb
from quiverlab.groebner.overlap import all_ambiguities
from quiverlab.groebner.complete import complete
from quiverlab.groebner.certificate import (
    default_degree_bound, check_degree_bound, certified_irreducibles,
)


@dataclass(frozen=True)
class ReductionSystem:
    quiver: object
    domain: object
    order: PathOrder
    rules: tuple
    irreducibles: tuple
    degree_bound: int
    is_confluent: bool = True

    def leading_words(self):
        return tuple(r.lead for r in self.rules)

    def reduce(self, comb):
        """Normal form of a linear combination (word -> domain element)."""
        return reduce_comb(comb, self.rules, self.order, self.domain)

    def normal_form(self, word):
        """Reduce a single word to its normal form (a linear combination)."""
        return reduce_comb({word: self.domain.one()}, self.rules, self.order, self.domain)

    def ambiguities(self):
        """All overlap/inclusion ambiguities of the completed system (each resolves
        to 0); the seed of the Chouhy-Solotar S-sequence (spec §6)."""
        return tuple(all_ambiguities(self.rules, self.degree_bound))


def _relations(quiver, relations):
    if relations and isinstance(relations[0], str):
        from quiverlab.combinat.relations import parse_relations
        return parse_relations(list(relations), quiver)
    return list(relations)


def build_reduction_system(quiver, relations, field, degree_bound=None, trace=None):
    from quiverlab.errors import AdmissibilityError
    rels = _relations(quiver, relations)
    for rel in rels:
        if rel.min_length < 2:
            raise AdmissibilityError(
                f"relation {rel!r} has a path of length {rel.min_length}: the ideal is not "
                "inside the square of the arrow ideal",
                hint="admissible relations use paths of length >= 2",
            )
    # Build the domain from the coefficient entries (0 and 1 always included so the
    # domain matches the monomial route's on monomial inputs).
    raw = [field.parse_entry(0), field.parse_entry(1)]
    for rel in rels:
        for c, _w in rel.terms:
            raw.append(field.parse_entry(c))
    dom = field.make_domain(raw)
    order = path_order(quiver)
    init_rules = []
    for rel in rels:
        comb = {}
        for c, w in rel.terms:
            cc = dom.coerce(field.parse_entry(c))
            if not dom.is_zero(cc):
                comb[w] = cc
        if not comb:
            raise RelationError(
                f"relation {rel!r} vanishes in {dom.name}, so it says 0 = 0",
                hint="a relation must have a nonzero coefficient in the chosen field; "
                     "drop it or change the field",
            )
        init_rules.append(rule_from_comb(comb, order, quiver, dom))
    if degree_bound is None:
        degree_bound = default_degree_bound(rels)
    rules = complete(init_rules, order, quiver, dom, degree_bound, trace=trace)
    check_degree_bound(rules, degree_bound)                 # AdmissibilityError if too small
    irreducibles = certified_irreducibles(quiver, rules)    # NotFiniteDimensionalError if infinite
    return ReductionSystem(quiver=quiver, domain=dom, order=order, rules=tuple(rules),
                           irreducibles=irreducibles, degree_bound=degree_bound)
```

Then update `src/quiverlab/groebner/__init__.py`:

```python
"""quiverlab.groebner: noncommutative Groebner completion for path-algebra ideals
(spec §5, component 3). Emits the ReductionSystem consumed by Chouhy-Solotar
(Plan 04). Path composition is LEFT TO RIGHT throughout."""
from quiverlab.groebner.order import PathOrder, path_order  # noqa: F401
from quiverlab.groebner.events import Dispatch, ReductionStep  # noqa: F401
from quiverlab.groebner.reduction import ReductionRule  # noqa: F401
from quiverlab.groebner.overlap import Ambiguity  # noqa: F401
from quiverlab.groebner.system import ReductionSystem, build_reduction_system  # noqa: F401

__all__ = [
    "PathOrder", "path_order", "Dispatch", "ReductionStep",
    "ReductionRule", "Ambiguity", "ReductionSystem", "build_reduction_system",
]
```

- [ ] **Step 4: Run the suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_system.py tests/test_no_floats.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/groebner/system.py src/quiverlab/groebner/__init__.py tests/groebner/test_system.py
git commit -m "feat(groebner): ReductionSystem object frozen for Chouhy-Solotar (Plan 04)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: General `kQ/I` lowering to `Algebra` + `Quiver.algebra` dispatch

**Files:**
- Create: `src/quiverlab/groebner/lower.py`, `tests/groebner/test_lower.py`
- Modify: `src/quiverlab/combinat/quiver.py` (route non-monomial relations through the Gröbner lowering; emit a `Dispatch` event when `trace` is given; add optional `degree_bound=None, trace=None` kwargs)
- Modify: `tests/core/test_monomial.py` (retarget the superseded Plan-01 gateway test `test_nonmonomial_waits_for_plan03`: after this task, `Q.algebra(relations=["a*b - c"])` routes through the Gröbner engine and raises `AdmissibilityError` — `c` has length 1 — not `NotImplementedError`)

**Interfaces:**
- Consumes: `ReductionSystem`/`build_reduction_system`, `Algebra`, `Quiver`, `Dispatch`, fields.
- Produces:
  - `groebner_algebra(quiver, relations, field, degree_bound=None, trace=None) -> Algebra` — build the reduction system, take its certified irreducible-path basis (trivial vertices first, then irreducibles sorted `(len, word)` — **identical basis order and labels to `core.monomial.build_monomial_algebra`**), fill structure constants `T[i][j] = normal_form(word_i ⌢ word_j)` in irreducible-word coordinates, `unit = Σ e_v`. Attaches `_quiver`, `_relations`. Not pre-unit-adapted (the bar path calls `.unit_adapted()` lazily, exactly as in Plan 01).
  - `Quiver.algebra(relations=(), field=None, degree_bound=None, trace=None)` — monomial inputs → Plan-01 `build_monomial_algebra` (unchanged); otherwise → `groebner_algebra`. On monomial inputs the two routes produce **elementwise-equal** `Algebra` objects (asserted).

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_lower.py`:

```python
"""General kQ/I lowering + Quiver.algebra dispatch (spec §3.3, §5 components 3, 4)."""
from quiverlab.combinat import Quiver
from quiverlab.core.monomial import build_monomial_algebra
from quiverlab.combinat.relations import parse_relations
from quiverlab.fields import CC
from quiverlab.groebner import Dispatch
from quiverlab.groebner.lower import groebner_algebra


def _loop():
    return Quiver([1], {"x": (1, 1)})


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def _same_algebra(A, B):
    return (A.dim == B.dim and A.basis_labels == B.basis_labels
            and A.T == B.T and A.unit == B.unit)


def test_fixture1_monomial_route_equivalence():
    """x^3: the Groebner route and the Plan-01 monomial route agree elementwise."""
    Q = _loop()
    rels = parse_relations(["x^3"], Q)
    A_mono = build_monomial_algebra(Q, rels, CC)
    A_grob = groebner_algebra(Q, rels, CC)
    assert _same_algebra(A_mono, A_grob)
    assert A_grob.dim == 3
    assert A_grob.basis_labels == ["e_1", "x", "x*x"]


def test_fixture2_square_dim_nine():
    Q = _square()
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    assert A.dim == 9
    assert A.basis_labels == ["e_1", "e_2", "e_3", "e_4", "a", "b", "c", "d", "a*b"]


def test_square_realizes_commutativity_in_structure_constants():
    """In the algebra, c*d = a*b (the reduction cd -> ab)."""
    Q = _square()
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    labels = A.basis_labels

    def e(name):
        v = [A.domain.zero()] * A.dim
        v[labels.index(name)] = A.domain.one()
        return v

    cd = A.multiply(e("c"), e("d"))
    ab = A.multiply(e("a"), e("b"))
    assert cd == ab
    assert cd[labels.index("a*b")] == A.domain.one()


def test_dispatch_records_route():
    Q = _square()
    trace = []
    Q.algebra(relations=["a*b - c*d"], field=CC, trace=trace)
    routes = [ev.route for ev in trace if isinstance(ev, Dispatch)]
    assert routes == ["groebner"]

    trace2 = []
    Q2 = _loop()
    Q2.algebra(relations=["x^3"], field=CC, trace=trace2)
    assert [ev.route for ev in trace2 if isinstance(ev, Dispatch)] == ["monomial"]


def test_monomial_dispatch_still_uses_plan01_path():
    """A monomial input through Quiver.algebra equals the direct monomial build."""
    Q = _loop()
    A_dispatch = Q.algebra(relations=["x^3"], field=CC)
    A_mono = build_monomial_algebra(Q, parse_relations(["x^3"], Q), CC)
    assert _same_algebra(A_dispatch, A_mono)


def test_fixture6_quantum_ci_dim_four():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    A = Q.algebra(relations=["x^2", "y^2", "y*x - x*y"], field=CC)
    assert A.dim == 4
    assert A.basis_labels == ["e_1", "x", "y", "x*y"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_lower.py -q`
Expected: FAIL / error — `ModuleNotFoundError: No module named 'quiverlab.groebner.lower'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/groebner/lower.py`:

```python
"""Lower a general (non-monomial) kQ/I to a Plan-01 structure-constant Algebra
(spec §3.3, §5 components 3-4). The basis is the certified irreducible-path basis
from the reduction system; structure constants are products reduced to normal form.

The basis is ordered EXACTLY as core.monomial.build_monomial_algebra: trivial
vertex idempotents first (in quiver.vertices order), then irreducible words sorted
(len, word). On a monomial input the reduction system's forbidden words equal the
relation words and reduction sends every reducible product to 0, so the resulting
Algebra is elementwise-equal to the monomial route (a required cross-check)."""
from quiverlab.core.algebra import Algebra
from quiverlab.groebner.system import build_reduction_system


def groebner_algebra(quiver, relations, field, degree_bound=None, trace=None):
    rs = build_reduction_system(quiver, relations, field,
                                degree_bound=degree_bound, trace=trace)
    dom = rs.domain
    zero, one = dom.zero(), dom.one()

    basis = [("e", v) for v in quiver.vertices] + [("p", w) for w in rs.irreducibles]
    index = {b: i for i, b in enumerate(basis)}
    m = len(basis)

    def src(b):
        return b[1] if b[0] == "e" else quiver.word_source(b[1])

    def tgt(b):
        return b[1] if b[0] == "e" else quiver.word_target(b[1])

    def product_vector(bi, bj):
        """Coordinate vector of b_i * b_j (concatenate, then reduce to normal form)."""
        vec = [zero] * m
        if tgt(bi) != src(bj):
            return vec
        if bi[0] == "e":                       # e * y = y
            vec[index[bj]] = one
            return vec
        if bj[0] == "e":                       # x * e = x
            vec[index[bi]] = one
            return vec
        nf = rs.reduce({bi[1] + bj[1]: one})   # normal form of the concatenated word
        for word, coeff in nf.items():
            vec[index[("p", word)]] = coeff
        return vec

    T = [[product_vector(bi, bj) for bj in basis] for bi in basis]
    unit = [zero] * m
    for v in quiver.vertices:
        unit[index[("e", v)]] = one
    labels = [f"e_{b[1]}" if b[0] == "e" else "*".join(b[1]) for b in basis]
    return Algebra(dom, T, unit, basis_labels=labels, _quiver=quiver,
                   _relations=list(relations))
```

Now modify `src/quiverlab/combinat/quiver.py` — replace the `algebra` method body:

```python
    def algebra(self, relations=(), field=None, degree_bound=None, trace=None):
        """Build kQ/I over the field (default CC). Monomial presentations route
        through the Plan-01 monomial path; general (non-monomial) relations route
        through the Groebner engine (Plan 03). Paths compose LEFT TO RIGHT.

        degree_bound: cap on ambiguity length during completion (default: adaptive).
        trace: optional list; step events (Dispatch, ReductionStep) are appended
        (inert hooks -- formal trace rendering is Plan 07)."""
        from quiverlab.combinat.relations import parse_relations

        if field is None:
            from quiverlab.fields import CC
            field = CC
        rels = parse_relations(list(relations), self)
        if all(r.is_monomial for r in rels):
            if trace is not None:
                from quiverlab.groebner.events import Dispatch
                trace.append(Dispatch(route="monomial",
                                      reason="every relation is a single monomial",
                                      n_relations=len(rels)))
            from quiverlab.core.monomial import build_monomial_algebra
            return build_monomial_algebra(self, rels, field)
        if trace is not None:
            from quiverlab.groebner.events import Dispatch
            trace.append(Dispatch(route="groebner",
                                  reason="at least one relation is non-monomial",
                                  n_relations=len(rels)))
        from quiverlab.groebner.lower import groebner_algebra
        return groebner_algebra(self, rels, field, degree_bound=degree_bound, trace=trace)
```

(The old `raise NotImplementedError(...)` block is removed.)

- [ ] **Step 4: Verify the new test passes and observe the one expected Plan-01 regression**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_lower.py tests/test_no_floats.py -q`
Expected: all pass (Fixture 1 both-route equivalence; Fixture 2 dim 9 with `c*d = a*b`; Fixture 6 dim 4; dispatch events).

Now run the **full** suite:
Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
Expected: **exactly one** failure — `tests/core/test_monomial.py::test_nonmonomial_waits_for_plan03` — with `DID NOT RAISE <class 'NotImplementedError'>` (it now raises `AdmissibilityError`, because the dispatch routes `a*b - c` through the Gröbner engine and the length-1 summand `c` is not admissible). This is the intended, correct behavior change; the next step retargets the stale test. If any OTHER test fails, stop and fix the implementation before proceeding.

- [ ] **Step 5: Retarget the superseded Plan-01 gateway test**

In `tests/core/test_monomial.py`, replace the `test_nonmonomial_waits_for_plan03` function

```python
def test_nonmonomial_waits_for_plan03():
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    with pytest.raises(NotImplementedError):
        Q.algebra(relations=["a*b - c"])
```

with (same gateway intent — a non-monomial relation is admitted by the engine and admissibility-checked; `AdmissibilityError` is already imported at the top of this file):

```python
def test_nonmonomial_routes_through_groebner_and_checks_admissibility():
    # Plan 03: non-monomial relations now lower through the Groebner engine. A
    # summand of length 1 (here `c`) is not admissible -- the ideal is not inside
    # the square of the arrow ideal -- so the gateway raises AdmissibilityError.
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    with pytest.raises(AdmissibilityError):
        Q.algebra(relations=["a*b - c"])
```

Then run the full suite again:
Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
Expected: all pass (the suite is now genuinely green — the retarget removed the one expected regression).

- [ ] **Step 6: Commit**

```bash
git add src/quiverlab/groebner/lower.py src/quiverlab/combinat/quiver.py tests/groebner/test_lower.py tests/core/test_monomial.py
git commit -m "feat(groebner): general kQ/I lowering to Algebra + Quiver.algebra dispatch

Routes non-monomial relations through the Groebner engine and retargets the
Plan-01 gateway test (a*b - c now raises AdmissibilityError, not NotImplementedError).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Property tests — unique normal forms, monomial-route equivalence, characteristic independence

**Files:**
- Create: `tests/groebner/test_properties.py`

**Interfaces:**
- Consumes: everything above.
- Produces: no source changes — property/invariance guarantees over the built system.

Properties asserted:
1. **Unique normal forms after completion** — reducing an element by applying rules in *random* legal orders always yields the same normal form as the deterministic reducer (Church–Rosser of the completed, confluent system).
2. **Monomial-route equivalence (property form)** — for a family of monomial presentations, `groebner_algebra` and `build_monomial_algebra` agree elementwise.
3. **Characteristic independence on Fixture 2** — the commutative square has the same basis dimension and the same `HH^0..HH^1` dims over exact ℂ (`CC`) and `GF(32003)`.

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_properties.py`:

```python
"""Property/invariance tests for the Groebner engine (spec §5 c.3; §8 ring 4)."""
import itertools
import random

from quiverlab.combinat import Quiver
from quiverlab.combinat.relations import parse_relations
from quiverlab.core.monomial import build_monomial_algebra
from quiverlab.fields import CC, GF, QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import first_factor, lc_add
from quiverlab.groebner.system import build_reduction_system
from quiverlab.groebner.lower import groebner_algebra


def _random_reduce(comb, rules, order, dom, rng):
    """Reduce by applying rules at RANDOMLY chosen reducible occurrences. For a
    confluent system this must reach the same normal form as reduce_comb."""
    work = {w: c for w, c in comb.items() if not dom.is_zero(c)}
    while True:
        reducible = [w for w in work if first_factor(w, rules) is not None]
        if not reducible:
            return work
        w = rng.choice(reducible)
        rule, i = first_factor(w, rules)
        coeff = work.pop(w)
        u, v = w[:i], w[i + len(rule.lead):]
        for tc, tw in rule.tail:
            lc_add(work, u + tw + v, dom.mul(coeff, tc), dom)


def test_unique_normal_form_random_orders_agree():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rs = build_reduction_system(Q, ["x^2", "y^2", "y*x - x*y"], QQ)   # Fixture 6, confluent
    rng = random.Random(20260718)
    letters = ["x", "y"]
    words = [tuple(p) for n in range(1, 6) for p in itertools.product(letters, repeat=n)]
    for w in words:
        deterministic = rs.reduce({w: QQ.coerce(1)})
        for _ in range(5):
            assert _random_reduce({w: QQ.coerce(1)}, rs.rules, rs.order, rs.domain, rng) \
                == deterministic


def test_monomial_route_equivalence_family():
    cases = [
        (Quiver([1], {"x": (1, 1)}), ["x^2"]),
        (Quiver([1], {"x": (1, 1)}), ["x^4"]),
        (Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}), ["a*b"]),
        (Quiver([1], {"x": (1, 1), "y": (1, 1)}), ["x^2", "x*y", "y*x"]),
    ]
    for Q, rels in cases:
        parsed = parse_relations(rels, Q)
        A_mono = build_monomial_algebra(Q, parsed, CC)
        A_grob = groebner_algebra(Q, parsed, CC)
        assert A_mono.dim == A_grob.dim
        assert A_mono.basis_labels == A_grob.basis_labels
        assert A_mono.T == A_grob.T
        assert A_mono.unit == A_grob.unit


def test_characteristic_independence_square():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    dims = {}
    for field, key in [(CC, "CC"), (GF(32003), "p")]:
        A = Q.algebra(relations=["a*b - c*d"], field=field)
        assert A.dim == 9
        dims[key] = A.hochschild_cohomology(1).dims
    assert dims["CC"] == dims["p"]
    assert dims["CC"] == [1, 0]         # HH^0 = 1 (center), HH^1 = 0 (kA_2 x kA_2, Kunneth)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_properties.py -q`
Expected: PASS (all machinery already exists). If any property fails, it exposes a real confluence/equivalence/characteristic bug in Tasks 1–7 — fix the offending module (do not weaken the property). This task exists to lock the invariants; there is no new source.

- [ ] **Step 3: (no implementation — invariance lock)**

If Step 2 failed, apply the systematic-debugging skill to the implicated module and re-run. Otherwise proceed.

- [ ] **Step 4: Run the full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add tests/groebner/test_properties.py
git commit -m "test(groebner): unique normal forms, monomial equivalence, characteristic independence

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: End-to-end acceptance + frozen-interface statement for Plan 04

**Files:**
- Create: `tests/groebner/test_acceptance.py`
- Modify: `README.md` (add the general-`kQ/I` example that the acceptance test executes verbatim)

**Interfaces:**
- Consumes: the public surface — `Quiver`, `Quiver.algebra`, `Algebra.hochschild_cohomology`, `CC`, `GF`, `build_reduction_system`, `ReductionSystem`.
- Produces: the Plan-03 exit criterion (README example runs in CI; HH dims asserted through the existing bar path over the Gröbner-built `Algebra`) and the frozen `ReductionSystem` interface statement below.

- [ ] **Step 1: Write the failing test**

`tests/groebner/test_acceptance.py`:

```python
"""Plan 03 acceptance: a general (non-monomial) kQ/I built by Quiver.algebra, with
Hochschild cohomology computed through the existing Plan-01 bar path, exact and
characteristic-independent. Also asserts the ReductionSystem shape Plan 04 consumes."""
from dataclasses import fields as dc_fields

from quiverlab.combinat import Quiver
from quiverlab.fields import CC, GF
from quiverlab.groebner import ReductionSystem, build_reduction_system


def _commutative_square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_readme_general_kqi_hochschild():
    """The README example: commutative square kQ/(a*b - c*d), a genuine non-monomial
    presentation, lowered through Groebner completion and fed to the bar complex."""
    Q = _commutative_square()
    A = Q.algebra(relations=["a*b - c*d"], field=CC)      # non-monomial -> Groebner route
    assert A.dim == 9
    # HH^0 = center = 1 ; HH^1 = 0 ; HH^2 = 0  (kA_2 x kA_2 by Kunneth; both factors
    # are trees, HH^{>=1}(kA_2)=0). HH^2 via the bar oracle on a 9-dim algebra is a
    # ~4608x576 exact rank: ~1.5s over GF(32003), ~17s over exact CC -- the sparse
    # differential lets rref terminate early; well within the max_cells guard.
    hh = A.hochschild_cohomology(2)
    assert hh.dims == [1, 0, 0]

    # exact over any field, same dimensions:
    Ap = Q.algebra(relations=["a*b - c*d"], field=GF(32003))
    assert Ap.dim == 9
    assert Ap.hochschild_cohomology(2).dims == hh.dims


def test_reduction_system_is_the_plan04_contract():
    rs = build_reduction_system(_commutative_square(), ["a*b - c*d"], CC)
    # exact frozen field set + methods Chouhy-Solotar (Plan 04) depends on
    assert {f.name for f in dc_fields(ReductionSystem)} == {
        "quiver", "domain", "order", "rules", "irreducibles", "degree_bound", "is_confluent"}
    assert rs.leading_words() == (("c", "d"),)
    assert rs.irreducibles == (("a",), ("b",), ("c",), ("d",), ("a", "b"))
    # every ambiguity of the completed system resolves to zero (CS reduction-finiteness)
    from quiverlab.groebner.complete import s_polynomial
    for amb in rs.ambiguities():
        assert rs.reduce(s_polynomial(amb, rs.domain)) == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/groebner/test_acceptance.py -q`
Expected: PASS if Tasks 1–8 are complete (this task is the assembly/acceptance). If it fails, the failure localizes the missing piece; fix it before proceeding. (Written test-first so the acceptance criterion is explicit.)

- [ ] **Step 3: Update the README**

Append to `README.md`:

```markdown
## General quivers with relations (kQ/I)

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3, 4],
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
A = Q.algebra(relations=["a*b - c*d"], field=CC)   # commutative square, exact
print(A.dim)                                        # 9
print(A.hochschild_cohomology(1))                   # HH^0 = 1  HH^1 = 0
```

Non-monomial relations are completed with an exact noncommutative Gröbner
(Buchberger–Mora overlap) engine and certified finite-dimensional; a
non-admissible or infinite presentation fails loudly with `AdmissibilityError`
or `NotFiniteDimensionalError`, never a hang.
```

- [ ] **Step 4: Run the full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
Expected: all pass (Plan-01 suite + all `tests/groebner/`). Confirm the float gate `tests/test_no_floats.py` is green.

- [ ] **Step 5: Commit**

```bash
git add tests/groebner/test_acceptance.py README.md
git commit -m "test(groebner): end-to-end acceptance for general kQ/I + HH through the bar path

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Frozen interface for Plan 04 (Chouhy–Solotar)

Plan 04 consumes the following **verbatim** (do not rename fields, change types, or alter method signatures without a coordinated interface bump). Element representation everywhere: a linear combination is `dict[tuple[str, ...], <domain element>]` (word → nonzero coefficient), the zero element is `{}`; words are arrow-name tuples read left-to-right.

```python
# quiverlab.groebner.order
@dataclass(frozen=True)
class PathOrder:
    arrow_index: dict                       # arrow name -> rank (insertion order)
    def key(self, word) -> tuple: ...        # (len, ranks)
    def compare(self, u, v) -> int: ...      # -1 / 0 / 1
    def leading(self, comb) -> tuple | None  # max support word of a linear combination
    def sort_words(self, words) -> list      # ascending

def path_order(quiver) -> PathOrder

# quiverlab.groebner.reduction
@dataclass(frozen=True)
class ReductionRule:
    lead: tuple                              # leading word (a path, left-to-right)
    tail: tuple                              # ((coeff, word), ...), every word < lead
    source: object                           # = quiver.word_source(lead)
    target: object                           # = quiver.word_target(lead)

# quiverlab.groebner.overlap
@dataclass(frozen=True)
class Ambiguity:
    kind: str                                # "overlap" | "inclusion"
    word: tuple                              # ambiguity word containing both leads
    left: ReductionRule                      # overlap: lead a*b ; inclusion: OUTER rule
    right: ReductionRule                     # overlap: lead b*c ; inclusion: INNER rule
    a: tuple                                 # prefix context
    c: tuple                                 # suffix context

# quiverlab.groebner.system
@dataclass(frozen=True)
class ReductionSystem:
    quiver: object
    domain: object                           # the exact Domain
    order: PathOrder
    rules: tuple                             # tuple[ReductionRule, ...], confluent
    irreducibles: tuple                      # tuple[tuple[str, ...], ...], certified finite,
                                             #   nonempty irreducible words, sorted (len, word)
    degree_bound: int
    is_confluent: bool = True

    def leading_words(self) -> tuple         # tuple(r.lead for r in rules)
    def reduce(self, comb: dict) -> dict     # normal form of a linear combination
    def normal_form(self, word: tuple) -> dict   # reduce a single word
    def ambiguities(self) -> tuple           # tuple[Ambiguity, ...]; the CS S-sequence seed

def build_reduction_system(quiver, relations, field,
                           degree_bound=None, trace=None) -> ReductionSystem
```

The Chouhy–Solotar S-sequence `S_1, S_2, …` seeds from `rules` (leads) and `ambiguities()`; the CS differentials call `reduce`/`normal_form` over `domain`; the irreducible-path basis (trivial vertices + `irreducibles`) matches the `Algebra` that `groebner_algebra` emits, so CS classes transport onto the same basis the bar path uses. `is_confluent is True` certifies CS's reduction-finiteness requirement (finitely many rules, finitely many irreducibles, every ambiguity resolves to `0`).

**Two flags for Plan 04.** (i) Because completion minimizes leads to an antichain under the factor relation, `ambiguities()` yields ONLY overlap (2-)ambiguities — there are no inclusion ambiguities left. `S_1` is the set of `rules` (leading words), `S_2` is `ambiguities()`; the higher CS chains `S_n` (`n ≥ 3`) must be constructed by Plan 04 from `rules` + `order` (iterated overlaps), not read off this object. (ii) Rule leads are minimal, but rule **tails are NOT inter-reduced** — a `tail` word may be reducible by another rule. If CS needs tails in normal form, call `rs.reduce({...})` on them; do not assume `ReductionRule.tail` is already irreducible.

## Boundary notes

- **Trace (Plan 07).** `Dispatch` and `ReductionStep` are inert plain dataclasses appended to a caller-supplied `trace` list. Formal trace integration (typed event taxonomy, PDF/HTML/text renderers, eliding rules, golden-file tests, `verbose=True` default) is Plan 07; this plan only fixes the emission points.
- **General exact coefficients.** The Plan-01 relation grammar accepts rational coefficients; those are coerced into the chosen `Domain` (so completion/reduction are coefficient-domain-generic and work over `CC`/`GF(q)` unchanged). Relation strings with non-rational exact coefficients (e.g. `"i"`, `"E(3)"`) are a later grammar extension (families such as `QuantumCI`, Plan 06) and need no change to the Gröbner engine — only the parser.
- **Performance.** `complete` recomputes all ambiguities each pass and `reduce` rescans words; correctness-first, mirroring `fields/linalg.py`. Fast/incremental completion is out of scope for v1 (the certificate, not speed, is the v1 flagship).
