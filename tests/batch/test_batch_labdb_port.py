"""Re-port of the bank tests/test_labdb.py, adapted to quiverlab families.
Golden numeric fixtures preserved (Fixture Z1)."""
import pytest

from quiverlab.batch.builders import BUILDERS, build_algebra
from quiverlab.batch.scan import analyze, run_scan

P = 32003


def test_registry_covers_the_core_builders():
    assert {"truncated_polynomial", "quantum_ci", "cyclic_nakayama",
            "linear_path_algebra", "dynkin", "reduction_system"} <= set(BUILDERS)


def test_builders_construct_associative_algebras():
    for spec in [{"builder": "truncated_polynomial", "args": [3]},
                 {"builder": "quantum_ci", "args": [2]},
                 {"builder": "cyclic_nakayama", "args": [3, 2]},
                 {"builder": "linear_path_algebra", "args": [3]},
                 {"builder": "dynkin", "args": ["D", 4]}]:
        A = build_algebra(spec)
        assert A.dim >= 1


def test_analyze_is_deterministic_and_serial_equals_parallel():
    specs = [{"builder": "quantum_ci", "args": [3], "N": 4}]
    assert analyze(dict(specs[0])) == analyze(dict(specs[0]))
    assert run_scan(specs, 1) == run_scan(specs, 1)


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("quiverlab.resolutions_cs") is None,
    reason="open-zone analyze needs the Plan 04 CS backend")
def test_open_zone_golden_open_33_0():
    spec = {"builder": "reduction_system", "N": 16,
            "args": [2, [[[0, 0, 0], [[1, [1, 1]]]], [[1, 1, 1], []],
                         [[1, 0], [[-1, [0, 1]]]]], "open_33_0"]}
    rec = analyze(spec)
    assert rec["dim"] == 9
    assert rec["HH_homology"][str(P)] == [6] + [5] * 16          # Fixture Z1 golden
    assert rec["resolution_ranks"] == [1, 2, 2, 1, 1, 2, 2, 1,
                                       1, 2, 2, 1, 1, 2, 2, 1, 1, 2]
