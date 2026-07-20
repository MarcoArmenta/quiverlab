# Development, release, and the QPA cross-check

## Continuous integration

- **`ci.yml`** — a fast test suite on every OS × Python 3.10–3.13 cell, plus two
  deep Linux legs running the full suite on the numba and pure-Python engine paths
  (`QUIVERLAB_NO_NUMBA=1`), plus a float-ban lint. Tests are bucketed by the
  `fast` / `deep` / `qpa` markers.
- **`qpa.yml`** — a Linux-only job that installs GAP + QPA (`passagemath-gap[qpa]`)
  and runs the cross-check suite (`-m qpa`), mandatory there.
- **`docs.yml`** — builds this site `--strict` (executing the tutorials) and deploys
  to GitHub Pages.
- **`paper.yml`** — compiles the JOSS paper draft to PDF.
- **`release.yml`** — on a `v*` tag: build, `twine check`, and publish to PyPI via
  OIDC trusted publishing (no API token).

## Releasing (semver 0.x → 1.0 at JOSS acceptance)

1. Bump `version` in `pyproject.toml` (and `__version__`); update `CHANGELOG.md`.
2. Commit, then `git tag vX.Y.Z && git push --tags`.
3. `release.yml` builds and publishes; the tag must equal the pyproject version.

## The QPA cross-check (the mathematics)

QPA has no Hochschild cohomology function, so `quiverlab[qpa]` assembles it via the
**enveloping algebra**. For `A = kQ/I` with enveloping algebra `A^e = A^{op} ⊗ A`,
Hochschild cohomology is `A^e`-module Ext of `A` with itself — `Ext^n_{A^e}(A, A)`:

$$ HH^n(A) \;=\; \mathrm{Ext}^n_{A^e}(A, A). $$

`A.crosscheck("hochschild", n)` scripts, in GAP: build the quiver and `PathAlgebra`
over the same field, form `EnvelopingAlgebra(A)`, present `A` as a right `A^e`-module,
take a minimal projective resolution, and read `dim Ext^k` for `k = 0..n`. It then
compares these to quiverlab's own `hochschild_cohomology(n).dims` and fails loudly on
any disagreement (both use QPA's `ExtAlgebraGenerators(-, n)[1]` dimension series).
`A.crosscheck("module_ext", M, n)` does the analogous check for module **self-Ext**
`Ext^*(M, M)`; distinct-module `Ext(M, N)` (needing `ExtOverAlgebra` + iterated
syzygies) is a flagged post-v1 extension. This is an independent-implementation oracle
(spec §8 ring 3), not a dependency of the core.
