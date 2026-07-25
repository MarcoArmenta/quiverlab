"""quiverlab.hpc -- the HPC/container batch tier + offline-app CLI (Plan 28).

The compute CLI lives in the WHEEL (`quiverlab-hpc` console script /
`python -m quiverlab.hpc`); the container image is a thin convenience wrapper.

Import policy (pinned by tests/hpc/test_import_boundary.py): importing
``quiverlab.hpc`` pulls in ONLY the public ``import quiverlab`` surface and the
stdlib -- never fastapi / uvicorn / jinja2 / pydantic / the ``webapp`` package.
The single sanctioned engine-internal reach (``quiverlab.engine.deepen``) is
isolated in ``deepen_runner.py`` and imported lazily there; the offline ``gui``
verb imports ``webapp`` lazily, inside its handler only.

Submodules:
  * ``resources``      -- exact-int host detection (cores / mem / GPUs).
  * ``spec``           -- stdlib request parse+validate + the compute dispatch.
  * ``deepen_runner``  -- the sanctioned public wrapper over ``engine.deepen``.
  * ``report``         -- result.json -> LaTeX/HTML/text via ``quiverlab.trace``.
  * ``cli`` / ``__main__`` -- the verbs.
"""
from quiverlab.hpc.spec import (  # noqa: F401
    ComputeError, SpecError, RESULT_SCHEMA, parse_request, run,
)

__all__ = ["ComputeError", "SpecError", "RESULT_SCHEMA", "parse_request", "run"]
