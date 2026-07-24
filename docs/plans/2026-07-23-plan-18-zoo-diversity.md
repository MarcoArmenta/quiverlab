# Plan 18 — Zoo diversity audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** the *standing* zoo (`families/zoo.py` + `zoo_catalog.json` + the batch
scan surface) contains mixed-length-tip (straddling) monomial and multi-vertex
presentations — the two shapes whose absence hid both 2026-07-22 bugs — with
diversity gates so curation can never silently regress (backlog Tier-1 item 5).

**Architecture:** today all 121 catalog records are one-vertex 2-generator
non-monomial algebras, and `build_from_record` hard-wires `Quiver([1], loops)` —
multi-vertex records are inexpressible. (1) Extend the record schema with optional
`"vertices"` (list) + `"arrows"` (list of `[name, s, t]`, list order = the index
space of `rules` words) — legacy records keep the loop path byte-for-byte.
(2) Curate five new standing records appended to `zoo_catalog.json` under a fresh
provenance block: `straddle_xx_yy_xyx` / `straddle_xx_yy_yxy` (monomial, mixed
tips {2,3} — the Plan-12 Bardzell-straddle shape; the catalog currently has ZERO
monomial records), `line_abc_cde` (6-vertex line quiver, tips `abc`/`cde` with the
straddling overlap `abcde` — the Plan-12+13 double witness), `comm_square`
(`kQ/(ab−cd)`, non-monomial multi-vertex), `cn_3_2` (`kZ₃/rad²`, periodic
multi-vertex). (3) The batch surface follows through the ONE declared touch-point
`_reduction_system_args` (+ the `reduction_system` builder gaining optional
`arrows`/`vertices` args; legacy specs byte-unchanged, mirroring the
`max_transient_bytes` back-compat law). (4) Gates: catalog well-formedness admits
the new keys; a diversity test pins "≥ 2 monomial mixed-tip records, ≥ 3
multi-vertex records"; each new record is live-certified (recorded `dim` == built
`A.dim`; HH via the minimal/corner engine ≡ bar over two primes) so the entries
are load-bearing, not dead JSON.

**Tech Stack:** JSON catalog + relation-string rendering (`_render_terms`);
engines already handle every new record (Plans 12/13/16). Tests in
`tests/families/` + `tests/batch/` → deep bucket.

## Global Constraints

- No float literals in `src/` (AST gate); catalog coefficients are exact ints.
- Python is always `.venv/bin/python`; tests via
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q …`.
- Legacy catalog records and legacy batch specs stay **byte-identical** (the
  existing zoo goldens — `open_33_0` depth-16 HH, `open2_33_712` — must pass
  unchanged).
- The read-only bank is untouched; `zoo_catalog.json` is quiverlab's own curated
  file (provenance blocks record every curation).
- Oracles live, never hardcoded: bar (`engine.hh_engine.hochschild_homology_dims`)
  vs minimal/corner (`minimal_homology_dims`) on every new record.
- Conventional commits; green at every commit; merge/push only when Marco asks.

---

### Task 1: Multi-vertex record schema in `build_from_record`

**Files:**
- Modify: `src/quiverlab/families/zoo.py`
- Test: `tests/families/test_zoo.py` (append)

**Interfaces:**
- Consumes: `Quiver(vertices, {name: (s, t)})`, `Q.algebra(relations=[...], field=…)`.
- Produces: `build_from_record(rec, field=None)` accepting optional
  `rec["vertices"]` + `rec["arrows"]` (list of `[name, s, t]`; `rules` words index
  the arrow list). Task 2's new catalog records and Task 3's batch builder rely on
  this schema.

- [ ] **Step 1: Write the failing test** — append to `tests/families/test_zoo.py`:

```python
def test_build_from_record_multivertex_schema():
    """Records may carry vertices + arrows (list order = rules index space);
    legacy records (no arrows) keep the one-vertex loop path."""
    from quiverlab.families.zoo import build_from_record
    rec = {"name": "cn_3_2_inline", "ngen": 3, "dim": 6,
           "vertices": [1, 2, 3],
           "arrows": [["a", 1, 2], ["b", 2, 3], ["c", 3, 1]],
           "rules": [[[0, 1], []], [[1, 2], []], [[2, 0], []]]}
    A = build_from_record(rec)
    assert A.dim == 6                                  # kZ_3/rad^2
    assert A.zoo_name == "cn_3_2_inline"
    nonmono = {"name": "sq_inline", "ngen": 4, "dim": 9,
               "vertices": [1, 2, 3, 4],
               "arrows": [["a", 1, 2], ["b", 2, 4], ["c", 1, 3], ["d", 3, 4]],
               "rules": [[[2, 3], [[1, [0, 1]]]]]}     # cd -> ab  (ab - cd = 0)
    B = build_from_record(nonmono)
    assert B.dim == 9                                  # commutative square
```

- [ ] **Step 2: Run to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/families/test_zoo.py::test_build_from_record_multivertex_schema -p no:cacheprovider`
Expected: FAIL (legacy path builds loops on one vertex → wrong dim / build error)

- [ ] **Step 3: Implement** — in `src/quiverlab/families/zoo.py`, generalize the
  word rendering to a names tuple and branch on `"arrows"`:

```python
def _word(idxs, gens=_GENS):
    return "*".join(gens[i] for i in idxs)


def _render_terms(terms, gens=_GENS):
    # (unchanged body, but every _word(...) call becomes _word(..., gens))
    ...


def build_from_record(rec, field=None):
    arrows = rec.get("arrows")
    if arrows is not None:
        # multi-vertex record (Plan 18): arrows = [[name, s, t], ...] in rules
        # index order; vertices listed explicitly
        gens = tuple(name for name, _s, _t in arrows)
        Q = Quiver(list(rec["vertices"]), {name: (s, t) for name, s, t in arrows})
    else:
        gens = _GENS[:rec["ngen"]]
        Q = Quiver([1], {g: (1, 1) for g in gens})
    rels = []
    for lead, tail in rec["rules"]:
        if not tail:
            rels.append(_word(lead, gens))                      # lead -> 0 (monomial)
        else:
            rels.append(_render_terms([(1, lead)] + [(-c, w) for c, w in tail], gens))
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("han_conjecture", "chouhy_solotar")
    A.zoo_name = rec["name"]
    return A
```

- [ ] **Step 4: Run the new test + the whole zoo suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/families/test_zoo.py -p no:cacheprovider`
Expected: PASS (legacy goldens byte-unchanged + the new schema test)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/families/zoo.py tests/families/test_zoo.py
git commit -m "feat(families): zoo records may carry vertices+arrows -- multi-vertex reduction systems (Plan 18)"
```

### Task 2: Curate the five standing records + diversity gates

**Files:**
- Modify: `src/quiverlab/families/zoo_catalog.json` (append provenance block + 5 records)
- Modify: `tests/families/test_zoo.py` (append gates; extend well-formedness)

**Interfaces:**
- Consumes: Task 1's schema.
- Produces: catalog records `straddle_xx_yy_xyx` (dim 6), `straddle_xx_yy_yxy`
  (dim 6), `cn_3_2` (dim 6), `comm_square` (dim 9), `line_abc_cde` (dim 16); the
  diversity-gate test names Task 4's docs reference.

- [ ] **Step 1: Append to `zoo_catalog.json`** (via a python script to keep JSON
  exact): a provenance dict
  `{"_provenance": {"source": "Plan-18 diversity audit (Plan-12/13 battery fixtures promoted)", "date": "2026-07-23", "curation": "mixed-length-tip monomial (straddling) + multi-vertex standing records; schema adds optional vertices/arrows"}}`
  and the records:

```json
{"name": "straddle_xx_yy_xyx", "ngen": 2, "dim": 6,
 "rules": [[[0, 0], []], [[1, 1], []], [[0, 1, 0], []]]}
{"name": "straddle_xx_yy_yxy", "ngen": 2, "dim": 6,
 "rules": [[[1, 1], []], [[0, 0], []], [[1, 0, 1], []]]}
{"name": "cn_3_2", "ngen": 3, "dim": 6, "vertices": [1, 2, 3],
 "arrows": [["a", 1, 2], ["b", 2, 3], ["c", 3, 1]],
 "rules": [[[0, 1], []], [[1, 2], []], [[2, 0], []]]}
{"name": "comm_square", "ngen": 4, "dim": 9, "vertices": [1, 2, 3, 4],
 "arrows": [["a", 1, 2], ["b", 2, 4], ["c", 1, 3], ["d", 3, 4]],
 "rules": [[[2, 3], [[1, [0, 1]]]]]}
{"name": "line_abc_cde", "ngen": 5, "dim": 16, "vertices": [1, 2, 3, 4, 5, 6],
 "arrows": [["a", 1, 2], ["b", 2, 3], ["c", 3, 4], ["d", 4, 5], ["e", 5, 6]],
 "rules": [[[0, 1, 2], []], [[2, 3, 4], []]]}
```

- [ ] **Step 2: Extend the well-formedness test and add the gates** — in
  `tests/families/test_zoo.py`, replace `assert rec["ngen"] in (2, 3)` with the
  schema-aware check, and append:

```python
# in test_catalog_is_bundled_and_well_formed, replace the ngen assert with:
#        if "arrows" in rec:
#            assert rec["ngen"] == len(rec["arrows"])
#            assert {v for _n, s, t in rec["arrows"] for v in (s, t)} <= set(rec["vertices"])
#        else:
#            assert rec["ngen"] in (2, 3)


def _is_monomial(rec):
    return all(not tail for _lead, tail in rec["rules"])


def _tip_lengths(rec):
    return {len(lead) for lead, _tail in rec["rules"]}


def test_zoo_diversity_gates():
    """The standing zoo must keep the two shapes whose absence hid the 2026-07-22
    bugs: mixed-length-tip MONOMIAL records (Bardzell straddling) and multi-vertex
    records. Curation that drops them fails here."""
    cat = load_catalog()
    straddle_mono = [r for r in cat if _is_monomial(r) and len(_tip_lengths(r)) > 1]
    multivertex = [r for r in cat if "arrows" in r]
    assert len(straddle_mono) >= 2, "zoo lost its straddling monomial records"
    assert len(multivertex) >= 3, "zoo lost its multi-vertex records"


NEW_RECORDS = ("straddle_xx_yy_xyx", "straddle_xx_yy_yxy", "cn_3_2",
               "comm_square", "line_abc_cde")


@pytest.mark.parametrize("name", NEW_RECORDS)
def test_new_records_are_live_certified(name):
    """Each Plan-18 record: recorded dim == built dim, and HH_0..3 via the
    minimal/corner engine == the normalized bar complex over two primes.
    Live oracles -- the records are load-bearing, not dead JSON."""
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    from quiverlab.engine.resolutions_minimal import minimal_homology_dims
    from quiverlab.families.zoo import build_from_record
    rec = next(r for r in load_catalog() if r.get("name") == name)
    for p in (32003, 2):
        A = build_from_record(rec, field=GF(p))
        assert A.dim == rec["dim"], f"{name}: catalog dim {rec['dim']} != built {A.dim}"
        E = to_engine(A.unit_adapted())
        mh = minimal_homology_dims(E, 3, primes=(p,))[p]
        bh = hochschild_homology_dims(E, 3, primes=(p,))[p]
        assert mh == bh[:len(mh)], f"{name} p={p}: {mh} != {bh}"
```

- [ ] **Step 3: Run the zoo suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/families/test_zoo.py -p no:cacheprovider`
Expected: PASS — including the legacy `open_33_0` / `open2_33_712` goldens
(byte-compat) and `test_zoo_yields_algebras_up_to_dim_max_in_order` now sweeping
the dim-6/9 newcomers.

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/families/zoo_catalog.json tests/families/test_zoo.py
git commit -m "feat(families): five standing diversity records + zoo diversity gates (Plan 18)"
```

### Task 3: Batch surface follows

**Files:**
- Modify: `src/quiverlab/batch/builders.py` (`_reduction_system` gains optional
  `vertices`/`arrows`)
- Modify: `src/quiverlab/batch/scan.py` (`_reduction_system_args` appends them
  when present)
- Test: `tests/batch/test_open_zoo_broaden.py` (append)

**Interfaces:**
- Consumes: Task 1's `build_from_record` schema; Task 2's `cn_3_2` record.
- Produces: `reduction_system` specs carrying quiver data; `analyze(spec)` routes
  multi-vertex records through the (Plan-13 corner) minimal engine.

- [ ] **Step 1: Write the failing test** — append to
  `tests/batch/test_open_zoo_broaden.py`:

```python
def test_multivertex_record_flows_through_specs_and_analyze():
    """A multi-vertex catalog record survives the spec adapter and analyze():
    the batch scan surface serves the Plan-18 diversity records."""
    from quiverlab.batch.scan import analyze, open_zoo_to_specs
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_minimal import minimal_homology_dims
    from quiverlab.families.zoo import build_from_record, load_catalog
    rec = next(r for r in load_catalog() if r.get("name") == "cn_3_2")
    specs = open_zoo_to_specs([rec], primes=(32003,))
    assert len(specs) == 1 and specs[0]["builder"] == "reduction_system"
    out = analyze(specs[0])
    E = to_engine(build_from_record(rec, field=GF(32003)).unit_adapted())
    ref = minimal_homology_dims(E, specs[0]["N"], primes=(32003,))[32003]
    got = out["HH"][32003] if isinstance(out["HH"], dict) else out["HH"]
    assert list(got)[:len(ref)] == ref[:len(list(got))]
```

(Adapt the final `out` unpacking to `analyze`'s actual record shape — read it in
`src/quiverlab/batch/scan.py::analyze` before finalizing; the assertion's content
— batch HH ≡ direct minimal-engine HH — must not weaken.)

- [ ] **Step 2: Run to verify it fails** (spec adapter drops the quiver data →
  builder builds the wrong algebra or crashes)

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/batch/test_open_zoo_broaden.py -p no:cacheprovider`

- [ ] **Step 3: Implement** — `builders.py`:

```python
def _reduction_system(ngen, rules, name, arrows=None, vertices=None, field=None):
    """Open-zone builder: materialise A = k<g>/I from a confluent reduction system,
    routed through quiverlab's zoo reduction-system -> Algebra path.  Multi-vertex
    records (Plan 18) carry arrows/vertices; legacy args are byte-unchanged."""
    rec = {"ngen": ngen, "rules": rules, "name": name, "dim": None}
    if arrows is not None:
        rec["arrows"] = arrows
        rec["vertices"] = vertices
    return _build_from_record(rec, field=field)
```

`scan.py::_reduction_system_args`:

```python
def _reduction_system_args(entry):
    """Catalog entry -> the ``reduction_system`` builder's args.

    Matches ``builders._reduction_system(ngen, rules, name[, arrows, vertices])``;
    legacy (local) entries keep the 3-arg shape byte-for-byte (spec back-compat)."""
    args = [entry["ngen"], entry["rules"], entry["name"]]
    if "arrows" in entry:
        args += [entry["arrows"], entry["vertices"]]
    return args
```

- [ ] **Step 4: Run the batch suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q tests/batch/ -p no:cacheprovider`
Expected: PASS (legacy spec goldens byte-unchanged + the new flow test)

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/batch/builders.py src/quiverlab/batch/scan.py tests/batch/test_open_zoo_broaden.py
git commit -m "feat(batch): reduction_system specs carry quiver data -- multi-vertex scan surface (Plan 18)"
```

### Task 4: Docs + suites + backlog

**Files:**
- Modify: `CLAUDE.md` (status paragraph gains Plan 18)
- Modify: `docs/plans/ROADMAP.md` (DELIVERED row 18)
- Modify: `docs/plans/DEEPER-ENGINES-BACKLOG.md` (item 5 tick — done at branch
  start; also tick/annotate the Tier-2 "Han's-conjecture batch campaigns" note
  that the multi-vertex scan surface now exists)
- Modify: `docs/internals/` families/zoo chapter if one exists (check
  `grep -rn "zoo" docs/internals/`)

- [ ] **Step 1: Doc edits** (facts per Architecture block)

- [ ] **Step 2: Full suites**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m deep -p no:cacheprovider`
Expected: PASS (~30 min)
Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast -p no:cacheprovider`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md docs/plans/ROADMAP.md docs/plans/DEEPER-ENGINES-BACKLOG.md docs/internals/
git commit -m "docs: Plan-18 status -- standing zoo diversity delivered"
```

## Validation matrix

1. Multi-vertex schema builds `kZ₃/rad²` (dim 6) and `kQ/(ab−cd)` (dim 9) from
   records; legacy records byte-unchanged (existing zoo goldens).
2. Five standing records live-certified: catalog dim == built dim; minimal ≡ bar
   HH over 2 primes each.
3. Diversity gates: ≥ 2 straddling-monomial + ≥ 3 multi-vertex records — curation
   regressions fail loudly.
4. Batch: multi-vertex record → spec → `analyze` ≡ direct minimal-engine HH;
   legacy specs byte-identical.
5. Full deep + fast suites green.

## Status

- [ ] Executed (fill in on completion)
