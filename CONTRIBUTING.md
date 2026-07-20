# Contributing to quiverlab

Thank you for your interest! quiverlab is a pure-Python library for exact
computation with finite-dimensional algebras (quivers with relations and their
Hochschild theory).

## Development setup

```bash
git clone https://github.com/MarcoArmenta/quiverlab
cd quiverlab
python3.12 -m venv .venv && . .venv/bin/activate
pip install -e ".[dev,fast]"          # add ,docs for the docs site; ,qpa on macOS/Linux
```

## Running tests

The suite is split by marker (see `pyproject.toml`):

```bash
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 pytest -m fast     # quick, cross-platform
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 pytest             # full (adds -m deep)
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 pytest -m qpa      # needs [qpa] + GAP (macOS/Linux)
QUIVERLAB_NO_NUMBA=1 pytest                              # the pure-Python engine path
```

Always throttle threads with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2`.

## Non-negotiables

- **Exact arithmetic only.** No float or complex literals anywhere under `src/`;
  the AST gate `tests/test_no_floats.py` enforces this and must stay green.
- **Loud failure over silent approximation.** Errors carry a fix-it hint.
- **Tests first.** Add a failing test, make it pass, keep the full suite green.
- **Conventional commits** (`feat:`, `fix:`, `test:`, `docs:`, `ci:`, `build:`).

## Reporting problems / asking for help

Open a GitHub issue. For the web GUI there is also an in-app feedback form.
