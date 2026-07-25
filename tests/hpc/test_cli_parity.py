"""CLI == public-API parity (Plan 28). A tiny config run through the CLI must
produce exactly what the direct ``import quiverlab`` calls compute, and the
spec-core ``run`` must agree with the CLI (minus the CLI-only ``result_schema``
envelope field)."""
import json
from pathlib import Path

import quiverlab as ql
from quiverlab.hpc import cli
from quiverlab.hpc.spec import RESULT_SCHEMA
from quiverlab.hpc.spec import run as spec_run

_LOOP_CFG = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_cohomology:0..3", "cartan", "dimension"],
    "artifacts": {"pdf": False, "tikz": False},
}


def _write(tmp_path, cfg, name="c.json"):
    p = tmp_path / name
    p.write_text(json.dumps(cfg), encoding="utf-8")
    return p


def test_cli_run_matches_direct_api(tmp_path):
    out = tmp_path / "result.json"
    rc = cli.main(["run", str(_write(tmp_path, _LOOP_CFG)), "-o", str(out)])
    assert rc == 0
    result = json.loads(out.read_text(encoding="utf-8"))

    A = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x*x"], field=ql.GF(2))
    assert result["results"]["hh_cohomology"]["dims"] == \
        list(A.hochschild_cohomology(3, verbose=False).dims)
    mat = A.cartan_matrix()
    expected = [list(r) for r in (mat.tolist() if hasattr(mat, "tolist") else mat)]
    assert result["results"]["cartan"]["matrix"] == expected == [[3]]
    assert result["results"]["dimension"]["value"] == A.dim
    assert result["result_schema"] == RESULT_SCHEMA


def test_cli_equals_spec_core_run(tmp_path):
    out = tmp_path / "result.json"
    cli.main(["run", str(_write(tmp_path, _LOOP_CFG)), "-o", str(out)])
    cli_result = json.loads(out.read_text(encoding="utf-8"))

    core = spec_run(_LOOP_CFG, tmp_path / "core")
    # The CLI adds result_schema; strip it for the comparison.
    cli_result.pop("result_schema", None)
    assert json.dumps(cli_result, sort_keys=True, default=str) == \
        json.dumps(core, sort_keys=True, default=str)


def test_container_ci_tiny_fixture_if_present(tmp_path):
    """If the container sibling's tiny fixture is present, the CLI runs it and the
    same document validates as a webapp ComputeRequest (cross-sibling contract)."""
    fixture = Path(__file__).resolve().parents[2] / "container" / "ci-tiny.yaml"
    if not fixture.exists():
        import pytest
        pytest.skip("container/ci-tiny.yaml not present in this worktree")
    out = tmp_path / "tiny.json"
    rc = cli.main(["run", str(fixture), "-o", str(out)])
    assert rc == 0
    result = json.loads(out.read_text(encoding="utf-8"))
    assert result["results"]["hh_cohomology"]["dims"][0] >= 1
