"""The tau_cluster algebra-level compute kind (Plan 66 / R29): served by hpc.spec, mirrored by
the Pyodide twin (tests/gui/test_tau_cluster_runner_twin.py), both byte-identical via the shared
builder tautilting.cluster_morphism.tau_cluster_block. kA3 -> #objects=14, face vector
(14,49,49,14) [classifying space], g-fan sphere (1,9,21,14), 6 generators, 4 atom + 2 comm,
abelianization Z^3. Unmarked (extras-gated dir, Plan-32 ruling). QQ-scope."""
from quiverlab.hpc.spec import parse_request
from quiverlab.hpc.spec import run as spec_run


def _req(vertices, arrows, relations, budget=512):
    return {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": vertices, "arrows": arrows,
                        "relations": relations, "field": {"kind": "QQ"}},
            "compute": ["tau_cluster:%d" % budget],
            "artifacts": {"pdf": False, "tikz": False}}


_KA3 = _req([1, 2, 3], {"a": [1, 2], "b": [2, 3]}, [])


def test_tau_cluster_block_shape(tmp_path):
    out = spec_run(parse_request(_KA3), tmp_path)
    b = out["results"]["tau_cluster"]
    assert b["complete"]
    assert b["category"]["object_count"] == 14
    # H1: face_vector = classifying space (f_0 = #wide); g_fan_face_vector = the SPHERE
    assert b["category"]["face_vector"] == [14, 49, 49, 14]
    assert b["category"]["g_fan_face_vector"] == [1, 9, 21, 14]
    assert b["category"]["face_vector"][0] == b["category"]["object_count"]   # f_0 == #wide
    assert sum(b["category"]["face_vector"]) == b["category"]["morphism_count"]  # == 126
    assert b["category"]["euler_characteristic"] == 0                          # NOT (-1)^3
    assert b["picture_group"]["num_generators"] == 6
    assert (b["picture_group"]["num_atom"], b["picture_group"]["num_commutation"]) == (4, 2)
    assert b["picture_group"]["abelianization_rank"] == 3                      # H2: Z^3 (6-3), not Z^2
    assert "hanson_igusa" in b["references"] and "igusa_todorov_cat0" in b["references"]


def test_tau_tilting_infinite_status(tmp_path):
    # 2-Kronecker, small budget -> complete False, status budget/error, category+picture_group
    # null, note set, no crash. Over GF(32003) budget=8 (char > dim, in P45 scope) to dodge the
    # >120s QQ BFS (the P64 precedent -- adjust-to-reality); the refusal is IDENTICAL.
    body = {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2],
                        "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                        "field": {"kind": "GF", "p": 32003, "n": 1}},
            "compute": ["tau_cluster:8"],
            "artifacts": {"pdf": False, "tikz": False}}
    out = spec_run(parse_request(body), tmp_path)
    b = out["results"]["tau_cluster"]
    assert b["complete"] is False and b["status"] in ("budget", "error")
    assert b["category"] is None and b["picture_group"] is None and b["note"]


def test_block_carries_citations_and_kind(tmp_path):
    out = spec_run(parse_request(_KA3), tmp_path)
    b = out["results"]["tau_cluster"]
    assert b["kind"] == "tau_cluster"
    # the block names the P66 sources in `references` (snake registry keys)...
    for k in ("buan_marsh_wide", "hanson_igusa", "igusa_todorov_weyman", "igusa_todorov_cat0"):
        assert k in b["references"], k
    # ...and every reference resolves to a (bibtex_key, human) citation pair.
    assert b["citations"] and all(len(pair) == 2 for pair in b["citations"])
    assert len(b["citations"]) == len(b["references"])
