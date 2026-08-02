"""End-to-end regression: Marco's dim-7-over-CC desktop session (2026-08-02).

He drew Q = (1 --a--> 2, loop b: 2->2) with relation ``b*b*b`` (dim 7 over CC),
attached a no-code module and a second module N, and checked the full compute
menu. The request died with ``ComputeError: FieldError: cannot read PythonMPQ
MPQ(1,1) as an exact scalar`` -- the CS-native cup/cap products fed sympy's own
QQ-domain elements back through the exact-scalar reader, a FALSE refusal.

This reconstructs that session (faithful to the saved ``fielderror.html`` state)
and runs it through ``run_spec``. It must COMPLETE with no FieldError and no error
blocks. ``bracket`` is intentionally omitted: over CC it is a SEPARATE, by-design
honest refusal (Plan 35: the Gerstenhaber bracket is served over GF(p) only), and
it sits after cup/cap in the compute order -- Marco never reached it.
"""
from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

ALG = {"kind": "quiver", "vertices": [1, 2],
       "arrows": {"a": [1, 2], "b": [2, 2]},
       "relations": ["b*b*b"], "field": {"kind": "CC"}}

# The main module M (right): the saved GUI state.
MODULE = {"dims": {"1": 1, "2": 3},
          "maps": {"a": [[1], [1], [0]],
                   "b": [[0, 1, 0], [0, 0, 1], [0, 0, 0]]},
          "side": "right"}

# The second argument N. Ext reads it as a RIGHT module (target x source blocks);
# Tor reads the SAME N as a LEFT module (source x target -- the transpose).
_EXT_A = [[1, 0], [1, 1], [0, 0], [0, 0]]
_EXT_B = [[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]]
EXT_TARGET = {"dims": {"1": 2, "2": 4}, "maps": {"a": _EXT_A, "b": _EXT_B}, "side": "right"}
_T = lambda m: [list(r) for r in zip(*m)]
TOR_TARGET = {"dims": {"1": 2, "2": 4},
              "maps": {"a": _T(_EXT_A), "b": _T(_EXT_B)}, "side": "left"}

# Marco's checked kinds, in GUI order, minus the GF(p)-only bracket. Small degree
# tops keep the regression cheap; the false refusal was degree-independent.
COMPUTE = [
    "hh_cohomology:0..2", "hh_homology:0..2", "cup:0..2", "cap:0..2", "connes_b:0..2",
    "cartan", "coxeter_polynomial", "global_dimension", "center",
    "dimension_vector", "rad_top_soc", "tau", "tau_minus",
    "projective_dimension", "injective_dimension",
    "projective_resolution:0..3", "injective_resolution:0..3", "decompose",
    "ext:0..2", "tor:0..2",
]


def test_marco_dim7_cc_session_completes(tmp_path):
    req = ComputeRequest.model_validate({
        "schema": 2, "algebra": ALG, "compute": COMPUTE,
        "module": MODULE, "ext_target": EXT_TARGET, "tor_target": TOR_TARGET,
    })
    result = run_spec(req, tmp_path)
    results = result["results"]

    # Every requested kind is present ...
    expected = {c.split(":")[0] for c in COMPUTE}
    assert expected <= set(results)

    # ... and none of them is an honest error block (the FieldError is gone, and
    # nothing else silently degraded).
    errs = {k: v["error"] for k, v in results.items()
            if isinstance(v, dict) and "error" in v}
    assert errs == {}, f"unexpected error blocks: {errs}"

    # The products that actually broke computed real tables over CC.
    assert results["cup"]["tables"]
    assert results["cap"]["tables"]
