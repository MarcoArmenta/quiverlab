"""Plan 69 / R33 -- the `barcode` module-side compute kind end-to-end through the
server runner (quiverlab.hpc.spec delegation): a schema-2 A_5 filtration request returns
the interval barcode {[1,5],[2,2],[4,4]} with the Escolar-Hiraoka / Gabriel citations.
Unmarked (webapp dir is extras-gated; oracle markers live in tests/modules)."""
import pathlib

from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

# The A_5 filtration H_0 module (dims (1,2,1,2,1)) as Plan-26 per-arrow BLOCK maps.
_A5_BARCODE = {
    "schema": 2,
    "algebra": {"kind": "quiver", "vertices": [1, 2, 3, 4, 5],
                "arrows": {"a": [1, 2], "b": [2, 3], "c": [3, 4], "d": [4, 5]},
                "relations": [], "field": {"kind": "QQ"}},
    "module": {"dims": {"1": 1, "2": 2, "3": 1, "4": 2, "5": 1},
               "maps": {"a": [[1], [0]], "b": [[1, 1]], "c": [[1], [0]], "d": [[1, 1]]}},
    "compute": ["barcode"],
    "artifacts": {"pdf": False, "tikz": False},
}


def test_a5_barcode_block(tmp_path):
    res = run_spec(ComputeRequest.model_validate(_A5_BARCODE), tmp_path)
    block = res["results"]["barcode"]
    assert block["kind"] == "persistence" and block["n"] == 5
    assert block["field_robust"] is True and block["diagram"] is None
    got = sorted((b["birth"], b["death"], b["multiplicity"]) for b in block["bars"])
    assert got == [(1, 5, 1), (2, 2, 1), (4, 4, 1)]
    assert any(b["death"] == 5 and b["essential"] for b in block["bars"])
    assert "escolar_hiraoka" in block["references"]
    assert any(c[0] == "EscolarHiraoka2016" for c in block["citations"])


def test_barcode_refusal_is_typed_block_not_500(tmp_path):
    # A non-A_n / non-CL quiver (kD4 star) -> the barcode block is an {"error": ...}
    # entry (a clean typed refusal), never a 500.
    body = dict(_A5_BARCODE)
    body = {**_A5_BARCODE,
            "algebra": {"kind": "quiver", "vertices": [0, 1, 2, 3],
                        "arrows": {"a": [1, 0], "b": [2, 0], "c": [3, 0]},
                        "relations": [], "field": {"kind": "QQ"}},
            "module": {"builtin": {"kind": "simple", "vertex": 0}}}
    res = run_spec(ComputeRequest.model_validate(body), tmp_path)
    block = res["results"]["barcode"]
    assert "error" in block and "escolar_hiraoka" in block["references"]
