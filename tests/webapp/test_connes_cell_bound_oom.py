"""Regression: a single over-cap bar product must NOT OOM/SIGKILL the whole
request (the hotfix-report-reps-oom bug).

The GF(p) ``connes_b`` path silently ignored its ``max_cells`` parameter (honoured
on the generic-Domain branch and by cup/cap/bracket) and built the DENSE
degree-(top+1) bar boundary matrix. On QuantumCI ``connes_b:0..7`` that is ``b_8``,
an 8748x26244 int64 array (~1.8 GB) whose rref/nullspace copies peaked ~6 GB and
SIGKILLed the memory-capped worker (exit 137) -- reproducibly, with pdf=true on the
curated qci-q2 example. It was NOT the worked-steps report path: ``render_html`` is
~17 MB and the peak is entirely in the ``connes_b`` compute (confirmed: the identical
request minus the report peaks the same 6.3 GB; the pdf correlation was an uncapped
CLI-vs-capped worker confound).

The fix: (1) ``connes_b`` honours ``max_cells`` and refuses LOUDLY before allocating
(over GF(p) it has no Chouhy-Solotar route, so -- like an explicit bar engine -- it
raises rather than silently OOMing); (2) the server runner degrades per-item on
``DepthLimitError`` into an honest error block, PARITY with the Pyodide runner
(``docs/gui/runner.py``), so every other invariant + the worked-steps report still
ship. Unmarked (contract/infrastructure), and cheap: the guard trips before any
matrix is built, so this never allocates gigabytes in CI.
"""
import pytest

from quiverlab.errors import DepthLimitError
from webapp.server.runner import build_algebra, run_spec
from webapp.server.schema import ComputeRequest


def _req(compute, pdf):
    return ComputeRequest.model_validate({
        "schema": 2,
        "algebra": {"kind": "family", "family": "QuantumCI",
                    "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": compute,
        "artifacts": {"pdf": pdf, "tikz": False}})


def test_connes_over_cap_raises_before_allocating():
    """The library guard trips BEFORE any dense matrix is built (cheap + loud);
    a feasible top still computes. connes_b:0..6 needs b_7 = 25.5M cells > 4M."""
    A = build_algebra(_req(["dimension"], False).algebra)
    cb = A.connes_differentials(4)                  # feasible top: computes
    assert set(cb.matrices) == {0, 1, 2, 3}
    with pytest.raises(DepthLimitError) as ei:      # over-cap top: refuses loudly
        A.connes_differentials(6)
    assert "max_cells" in str(ei.value)


def test_run_spec_completes_when_one_product_is_over_cap(tmp_path):
    """The whole request completes with pdf=true: connes_b becomes an honest
    per-item error block, every other invariant + the worked-steps report ship."""
    res = run_spec(_req(["hh_cohomology:0..2", "connes_b:0..6", "dimension"], True),
                   tmp_path)
    r = res["results"]
    # (a) completion -- all kinds present, no exception, no OOM.
    assert set(r) == {"hh_cohomology", "connes_b", "dimension"}
    # connes_b is the honest refusal, not a computed table.
    assert r["connes_b"]["error"]["type"] == "DepthLimitError"
    assert "max_cells" in r["connes_b"]["error"]["message"]
    # the other invariants stand, exact.
    assert r["dimension"]["value"] == 4
    assert isinstance(r["hh_cohomology"]["dims"], list)
    # (b) the report is written, states the honest refusal, and keeps the HH table.
    html = (tmp_path / "trace_steps.html").read_text(encoding="utf-8")
    assert "DepthLimitError" in html and "not computed" in html
    assert "Hochschild" in html
    assert (tmp_path / "trace.json").exists()


def test_over_cap_product_never_perturbs_the_other_results(tmp_path):
    """(c) The other results are byte-identical with and without the over-cap
    product riding along -- the refusal is isolated to its own block."""
    a = run_spec(_req(["hh_cohomology:0..3"], False), tmp_path / "a")
    b = run_spec(_req(["hh_cohomology:0..3", "connes_b:0..6"], False), tmp_path / "b")
    assert a["results"]["hh_cohomology"] == b["results"]["hh_cohomology"]
    assert b["results"]["connes_b"]["error"]["type"] == "DepthLimitError"
