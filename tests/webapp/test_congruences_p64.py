"""The congruences algebra-level compute kind (Plan 64 / R26): served by hpc.spec, mirrored
by the Pyodide twin (docs/gui/runner.py), both runners byte-identical via the shared builder
tautilting.congruence.congruences_block. kA2 -> |L|=5 (pentagon N5), |Con|=5, #wide=5 (M3),
complete; the wild 2-Kronecker trips status='budget' with every lattice invariant omitted.
Unmarked (extras-gated dir), per the Plan-32 cross-runner ruling."""
import json
import pathlib
import tempfile


def _congruences_request(budget=512):
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                    "relations": [], "field": {"kind": "QQ"}},
        "compute": [f"congruences:{budget}"],
        "artifacts": {"pdf": False, "tikz": False},
    }


def test_congruences_block_shape(tmp_path):
    from quiverlab.hpc.spec import parse_request, run
    req = parse_request(_congruences_request(512))
    out = run(req, tmp_path)
    b = out["results"]["congruences"]
    assert b["complete"] and b["status"] == "complete" and b["n"] == 2
    assert b["lattice"]["size"] == 5 and b["lattice"]["join_irreducibles"] == 3
    assert b["lattice"]["num_covers"] == 5
    assert b["lattice"]["is_semidistributive"] is True
    assert b["lattice"]["is_modular"] is False and b["lattice"]["is_distributive"] is False
    assert b["congruences"]["size"] == 5 and b["congruences"]["num_join_irreducibles"] == 3
    assert b["congruences"]["is_distributive"] is True
    # the forcing "V": one brick below two (one shared lower endpoint, two distinct uppers)
    rels = b["congruences"]["forcing_order"]["relations"]
    lowers = {r[0] for r in rels}
    uppers = {r[1] for r in rels}
    assert len(rels) == 2 and len(lowers) == 1 and len(uppers) == 2
    assert b["wide"]["size"] == 5 and b["wide"]["constructions_agree"] is True
    assert sorted(len(lab) for lab in b["wide"]["labels"]) == [0, 1, 1, 1, 2]
    assert "dirrt_lattice_torsion" in b["references"]
    assert b["citations"]                                # DIRRT/BCZ/Enomoto/... resolved


def _load_gui_runner():
    import importlib.util
    p = str(pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p64", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_twin_parity():
    # run the same request through the server (hpc.spec) and the Pyodide twin
    # (docs/gui/runner.py); assert the congruences block is byte-identical.
    from quiverlab.hpc.spec import parse_request, run
    req_dict = _congruences_request(512)
    with tempfile.TemporaryDirectory() as d:
        server = run(parse_request(req_dict), d)["results"]["congruences"]
    gr = _load_gui_runner()
    assert json.loads(gr.run_build(json.dumps(req_dict)))["ok"]
    twin = json.loads(gr.compute_one("congruences:512"))["block"]
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(twin, sort_keys=True, default=str))


def test_estimator_budget_is_not_a_degree():
    # the congruences budget (512) must NOT be read as homological degree 512 by the tier
    # classifier (like tau_tilting/ar_quiver). _max_degree must ignore the pair budget.
    from webapp.server.estimator import _max_degree
    from webapp.server.schema import ComputeRequest
    req = ComputeRequest.model_validate(_congruences_request(512))
    assert _max_degree(req) == 0


def test_tau_tilting_infinite_status(tmp_path):
    # 2-Kronecker with a small budget -> complete False, status "budget", note set,
    # lattice/congruences/wide all null, no crash (no partial-lattice lie).
    from quiverlab.hpc.spec import parse_request, run
    # GF(32003) (char > dim 4, in P45 scope), budget 8: the 2-Kronecker is
    # tau-tilting-infinite, so the graph never closes -- the honest refusal at a cheap budget
    # over a fast field (the test_wild_budget_status precedent; QQ budget 40 exceeds 120s).
    req = parse_request({
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2],
                    "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                    "field": {"kind": "GF", "p": 32003, "n": 1}},
        "compute": ["congruences:8"],
        "artifacts": {"pdf": False, "tikz": False}})
    b = run(req, tmp_path)["results"]["congruences"]
    assert b["complete"] is False and b["status"] == "budget"
    assert b["lattice"] is None and b["congruences"] is None and b["wide"] is None
    assert b["note"]
