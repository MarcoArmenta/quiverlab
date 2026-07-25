"""HTTP-boundary adapter over the wheel's HPC spec core (Plan 28).

The compute dispatch that used to live here was promoted into the wheel
(``quiverlab.hpc.spec``) so the container/CLI tier and the webapp share ONE
implementation. This module now:

  * keeps the pydantic validation at the HTTP boundary (the app builds a
    ``ComputeRequest`` before calling in);
  * DELEGATES the actual compute to :func:`quiverlab.hpc.spec.run`, passing the
    request as ``req.model_dump(by_alias=True)`` and returning the byte-identical
    result dict (pinned by ``tests/webapp/test_runner_delegation.py``: the result
    dicts AND the Plan-25 ``canonical_key`` are unchanged across a fixture
    corpus);
  * preserves the public surface other webapp modules import -- ``RunError``,
    ``run_spec``, ``build_algebra`` -- so ``app.py`` / ``instant.py`` /
    ``worker.py`` / ``bigjobs.py`` are untouched.

Execution semantics (no user-code exec; result-block shapes mirroring
``docs/gui/runner.py``; ``verbose`` pinned off; worked-steps written with an
explicit ``out_dir``) are unchanged -- they now live in the spec core."""
from __future__ import annotations

from typing import Callable

from quiverlab.hpc.spec import ComputeError, SpecError
from quiverlab.hpc.spec import build_algebra as _core_build_algebra
from quiverlab.hpc.spec import run as _core_run

from webapp.server.schema import ComputeRequest


class RunError(Exception):
    """A refusal the caller (worker/sync endpoint) turns into an honest error
    result. ``error_type`` is the library exception class name (verbatim) or one
    of the runner's own tags (``CatalogError``, ``ResultTooLarge``, ...).

    This is the identity other webapp modules catch; the spec core raises its own
    ``ComputeError`` which is translated to this class at the boundary below."""

    def __init__(self, error_type: str, message: str):
        super().__init__(f"{error_type}: {message}")
        self.error_type = error_type
        self.message = message


def build_algebra(spec):
    """Build a library ``Algebra`` from a validated ``AlgebraSpec`` (a pydantic
    model or a plain dict). Delegates to the spec core; a spec-core
    ``ComputeError`` (its FieldError / CatalogError tags) is re-raised as a
    ``RunError``, while a raw ``QuiverlabError`` (bad relation, non-prime field)
    propagates for ``app.py`` to tag by class name -- exactly as before."""
    payload = spec.model_dump() if hasattr(spec, "model_dump") else spec
    try:
        return _core_build_algebra(payload)
    except ComputeError as exc:
        raise RunError(exc.error_type, exc.message)


def run_spec(req: ComputeRequest, artifact_dir,
             progress_cb: Callable[[dict], None] | None = None,
             result_max_bytes: int | None = None) -> dict:
    """Run a validated request against the public quiverlab API and write
    ``result.json`` + artifacts into ``artifact_dir``. Returns the result dict,
    byte-identical to the pre-Plan-28 runner (the webapp never sets an ``hpc``
    block or a ``result_schema``, so the spec core takes its default path)."""
    try:
        return _core_run(req.model_dump(by_alias=True), artifact_dir,
                         progress_cb=progress_cb, result_max_bytes=result_max_bytes)
    except ComputeError as exc:
        raise RunError(exc.error_type, exc.message)
    except SpecError as exc:                       # defensive: pydantic already validated
        raise RunError(type(exc).__name__, str(exc))
