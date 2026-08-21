"""Wave 4 (v1.0.1, 2026-08-21): the worked-steps filename must fold the engine /
provenance, not just ``(repr(algebra), kind, top)``.

Observed directly: ``hochschild_cohomology(4, engine="auto")`` and
``engine="bar"`` on ``k[x]/(x^3)`` over GF(5) both wrote the SAME
``quiverlab_traces/HHc_<hash>.html`` -- the second silently overwriting the first
despite a different ``engine:`` line (hanlab fast rank vs the normalized bar
complex) in the body. The stem now folds the ``HHTable`` provenance (engine,
coefficient module, dims) and a digest of the deterministic JSON record, so two
reports whose bodies differ never share a filename, while identical computations
still map to one file (deterministic, idempotent).
"""
import pathlib

import pytest

import quiverlab
from quiverlab.fields import GF
import quiverlab.trace.writer as W


def _run(A, engine, monkeypatch, tmp_path):
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    monkeypatch.chdir(tmp_path)
    # verbose=True forces the worked-steps report to be written (the report is
    # otherwise gated by the global quiverlab.verbose, off under pytest).
    A.hochschild_cohomology(4, engine=engine, verbose=True)


def test_engine_provenance_distinguishes_trace_filename(tmp_path, monkeypatch):
    A = quiverlab.truncated_polynomial(3, field=GF(5))
    _run(A, "auto", monkeypatch, tmp_path)   # fast GF(p) hanlab engine
    _run(A, "bar", monkeypatch, tmp_path)    # normalized bar oracle
    out = tmp_path / "quiverlab_traces"
    htmls = sorted(out.glob("HHc_*.html"))
    # Two DISTINCT reports -- the auto run must not have been clobbered.
    assert len(htmls) == 2, [p.name for p in htmls]
    bodies = [p.read_text(encoding="utf-8") for p in htmls]
    assert bodies[0] != bodies[1]
    joined = "\n".join(bodies)
    assert "hanlab" in joined and "normalized bar complex" in joined
    # the JSON sidecars are likewise distinct and both present.
    assert len(sorted(out.glob("HHc_*.json"))) == 2


def test_identical_computation_reuses_one_filename(tmp_path, monkeypatch):
    # Determinism / idempotency: the SAME computation must map to ONE file, not a
    # new name each run (guards the events digest against nondeterminism).
    A = quiverlab.truncated_polynomial(3, field=GF(5))
    _run(A, "bar", monkeypatch, tmp_path)
    _run(A, "bar", monkeypatch, tmp_path)
    out = tmp_path / "quiverlab_traces"
    assert len(sorted(out.glob("HHc_*.html"))) == 1
    assert len(sorted(out.glob("HHc_*.json"))) == 1
