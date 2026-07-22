# Changelog

All notable changes to quiverlab are documented here. This project adheres to
[Semantic Versioning](https://semver.org) (0.x during battle-testing; 1.0 at JOSS
acceptance).

## [Unreleased]
### Added
- **In-browser GUI on the docs landing page (Plan 10).** Draw a quiver, type
  relations, pick CC/GF(p^n), and compute HH ranges, Cartan/Coxeter, gl.dim and
  the center — no installation; the page runs the repo's own wheel under
  Pyodide in a Web Worker. Worked-steps report printable to PDF; TikZ/JSON
  downloads; copy-paste Python reproduction. GUI requests use the Plan-09
  schema (`kind: "quiver"`), and `docs/gui/runner.py` is the execution
  semantics a future server tier reuses.
- **Live wait estimates in the GUI.** While configuring, the page shows the
  algebra's dimension and an honest coarse estimate (a few seconds / under a
  minute / a few minutes / could be long), predicts when a request would hit
  the engine's cell cap, and surfaces relation/field errors before Compute;
  while computing, elapsed time plus a re-scaling estimate. Calibrated to the
  visitor's machine at engine start.
- Optional `[qpa]` backend: `A.crosscheck(...)` (independent QPA recomputation).
- GitHub Actions CI (matrix + engine-path legs), docs site, JOSS paper draft.
- Modernized packaging (PEP 639 SPDX license), README, community files.
