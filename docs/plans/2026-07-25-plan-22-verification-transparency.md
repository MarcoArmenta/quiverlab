# Plan 22 — Verification transparency

Date: 2026-07-25. Branch: `plan-22-verification-transparency`.
Backlog: Tier 1a item ("display in the repo that everything is unit tested, and
say HOW"). This is a **docs + audit** plan — no `src/` change, no new algebra.

## Goal (Marco's words, honored exactly)

Display in the repo that everything is unit tested, and say HOW things are tested:
**both against QPA and against theory from the literature** — we construct many
examples the literature has already computed, or theorems we know give a result,
and test against them. Highest mathematical rigour and accuracy. QPA does not
implement everything quiverlab has; we compare to it whenever we can. **Never
overclaim — audit before asserting.**

## Method

Before writing a single sentence, the test tree was read and the suite counted
with the project's own tools (`pytest --collect-only -m <bucket>` in the worktree
venv, `.[dev,fast,docs,web]`). Every fact on the new page is traced to a named
test file that was actually opened and read.

## Audit findings (headline numbers)

- **1377 tests collected** (worktree venv with `[dev,fast,docs,web]`), partitioned
  by `tests/conftest.py` into disjoint buckets:
  - `fast` 504 (every OS × Python CI cell; includes the 161-test `webapp` suite),
  - `deep` 868 (one Linux cell, run **twice** — numba and pure-Python
    `QUIVERLAB_NO_NUMBA=1`),
  - `qpa` 5 (live GAP/QPA; weekly CI job, skipped locally),
  - `slow` 0 (none individually marked at present).
- **Two oracle classes are in force and were verified in the source:**
  1. *Theory / literature pins on constructed examples* — the normalized bar
     complex as the base oracle; multi-prime cross-checks over `{32003, 2, 3, 5}`;
     the read-only hanlab bank's byte-level closed forms; named literature-value
     pins (Happel 1989; Buchweitz–Green–Madsen–Solberg 2005; Bergh–Erdmann 2008;
     Loday); self-certifying internal identities (`d∘d=0`, CS order condition,
     Leibniz as the cup's sign arbiter, chain-map / roundtrip gates); closed-form
     family pins; second-model oracles (Connes λ-complex, relative-Tor Betti,
     socle-criterion Nakayama form); Bardzell chain-count pins.
  2. *Cross-engine + QPA agreement* — bar ≡ minimal ≡ Bardzell ≡ CS degreewise
     wherever two engines overlap; numba vs pure parity; live GAP/QPA recomputation
     of Hochschild dims (`HH^n = Ext^n_{A^e}(A,A)`) and module self-Ext.
- **Honest gaps found and scoped on the page (not hidden):**
  - QPA has **no native HH function** and **no cup/cap/bracket, cyclic homology,
    CS resolution, deep degrees, or Frobenius/Nakayama**; live QPA covers only
    Hochschild dims and module *self*-Ext, over **QQ or prime GF(p) only**. Where
    QPA cannot reach, the page names the theory oracle that covers that ground
    (bar complex; bank byte-level closed forms; cross-engine agreement;
    self-certifying identities; second models).
  - `complexity` is a lower-bound estimate off local / single-vertex inputs;
    `is_symmetric` off GF(p) uses a Schwartz–Zippel inner-automorphism sweep,
    **loud when inconclusive** — never a silent wrong answer.
  - `webapp/` and `docs/gui/` are non-algebraic glue: they delegate all
    mathematics to `import quiverlab` and are unit-tested for API / schema /
    isolation / artifacts, not oracle-validated mathematics. The float-ban AST
    gate exempts them by design (it scans `src/` only).

## Deliverables

1. `docs/verification.md` — the methodology (the two oracle classes in Marco's
   framing), a subsystem → oracles → test-file table with audited facts, the
   marker/bucket scheme with counts, the CI matrix, the honest-scope section, and
   the standing rule.
2. `mkdocs.yml` nav — `verification.md` wired in right after "Under the hood".
3. `README.md` — a concise "How quiverlab is verified" section after the web
   interface / usage sections and before "Status", linking to the full page.
4. This plan doc; `ROADMAP.md` + `DEEPER-ENGINES-BACKLOG.md` mark Tier 1a done;
   `CLAUDE.md` status line updated.

## Standing rule established by this plan

**Every future plan adds its new oracles to `docs/verification.md` as part of its
acceptance** — the same way each plan already updates the internals chapters. The
page is the single living record of how each shipped feature is verified.

## Acceptance

- `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast`
  green.
- `DISABLE_MKDOCS_2_WARNING=true .venv/bin/mkdocs build --strict` exits 0
  (the new page is in the nav; no `[i][j]`-style autoref-breaking sequences).
- No `src/` change; the float-ban gate and every existing release/docs test stay
  green.
