"""Plan 71 / R12 -- the ``hh_lie_module`` HH kind end-to-end through the server runner
(quiverlab.hpc.spec via webapp.server.runner). Unmarked (the extras-gated tests/webapp
dir; the Plan-32 ruling keeps oracle_* markers out of webapp/gui/hpc)."""
import json

from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest


def _block(body, tmp_path):
    return run_spec(ComputeRequest.model_validate(body), tmp_path)["results"]["hh_lie_module"]


_KRON_QQ = {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": [1, 2],
                        "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                        "field": {"kind": "QQ"}},
            "compute": ["hh_lie_module:0..2"],
            "artifacts": {"pdf": False, "tikz": False}}

_WITT_GF3 = {"schema": 2,
             "algebra": {"kind": "quiver", "vertices": [1],
                         "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                         "field": {"kind": "GF", "p": 3, "n": 1}},
             "compute": ["hh_lie_module:0..2"],
             "artifacts": {"pdf": False, "tikz": False}}


def test_kronecker_block_shape(tmp_path):
    b = _block(_KRON_QQ, tmp_path)
    assert b["kind"] == "hh_lie_module" and b["hh_dims"] == [1, 3, 0] and b["hh1_dim"] == 3
    assert b["module_axiom_ok"] is True and b["inner_acts_zero"] is True
    deg1 = next(e for e in b["summands"] if e["n"] == 1)
    assert len(deg1["parts"]) == 1 and deg1["parts"][0]["dim"] == 3        # irreducible L(2)
    # a degree-1 weight table that is a symmetric sl2-string (PATTERN, not literal ints)
    w1 = next(e for e in b["weights"] if e["n"] == 1)
    vals = sorted(int(pair[0][0]) for pair in w1["weights"])
    assert len(vals) == 3 and vals[0] == -vals[2] and vals[1] == 0
    assert all(pair[1] == 1 for pair in w1["weights"])                     # each mult 1
    assert "alsolotar_toupie" in b["references"]                          # ALS toupie anchor
    assert any("toupie" in str(pair[1]).lower() for pair in b["citations"])


def test_char_p_weight_gate(tmp_path):
    b = _block(_WITT_GF3, tmp_path)
    assert b["weights"] is None and b["char0_note"]        # weight/torus char-0 gated
    assert b["module_axiom_ok"] is True                    # field-general part still computed
    # summands are NOT asserted None -- decompose's guard governs them independently
    assert b["summands"] is not None


def test_default_top_and_citations(tmp_path):
    """A bare 'hh_lie_module' (no range) defaults to top = 2; citations resolve."""
    body = dict(_KRON_QQ, compute=["hh_lie_module"])
    b = _block(body, tmp_path)
    assert b["top"] == 2 and b["hh_dims"] == [1, 3, 0]
    assert b["citations"] and all(len(pair) == 2 for pair in b["citations"])
