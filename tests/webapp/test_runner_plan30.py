"""Plan 30 -- webapp runner dispatch: `decompose`, `tor`, the AR-translate input
certificate on tau/tau^-, and the resolution `summands` field.

The dispatch lives in the SPEC CORE (`quiverlab.hpc.spec`); the runner delegates.
`decompose` uses the Plan-30 Part-A engine (`quiverlab.modules.decompose`); `tor`
uses the Plan-29 engine (`quiverlab.modules.tor`). Both are imported LAZILY, so
this branch is green whether or not those engines have landed:

  * decompose tests skip if the engine is absent (it is present in this worktree);
  * the tor test asserts the honest "engine unavailable" refusal when tor is
    absent, and the real dims when it is present -- green in BOTH states.
"""
import pathlib
import tempfile

import pytest

from webapp.server.runner import RunError, run_spec
from webapp.server.schema import ComputeRequest
from webapp.server.security import is_safe_error_type

try:
    import quiverlab.modules.decompose as _decompose_mod  # noqa: F401
    _HAVE_DECOMPOSE = True
except ImportError:
    _HAVE_DECOMPOSE = False

try:
    import quiverlab.modules.tor as _tor_mod  # noqa: F401
    _HAVE_TOR = True
except ImportError:
    _HAVE_TOR = False

_needs_decompose = pytest.mark.skipif(
    not _HAVE_DECOMPOSE, reason="Plan-30 Part-A decompose engine not present yet")

# A_2 quiver 1 --a--> 2 over GF(5): char 5 > small dims, so the decomposition
# engine can certify (its trace-form radical is reliable when char > dim).
_A2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}}
# k[x]/(x^3) over GF(2): char 2 <= dim, where certification is LOUD (unreliable).
_LOOP = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
         "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}}
_M2 = {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}
# S_1 (+) S_2 over A_2: dims {1:1, 2:1}, a acts by the zero block -> decomposable.
_S1S2 = {"dims": {"1": 1, "2": 1}, "maps": {"a": [[0]]}}


def _run(algebra, compute, module=None, tor_target=None):
    body = {"schema": 2, "algebra": algebra, "compute": compute,
            "artifacts": {"pdf": False, "tikz": False}}
    if module is not None:
        body["module"] = module
    if tor_target is not None:
        body["tor_target"] = tor_target
    with tempfile.TemporaryDirectory() as d:
        return run_spec(ComputeRequest.model_validate(body), pathlib.Path(d))


# --------------------------------------------------------------------------- #
# Worked-steps bundle for module computations (Marco #5; Part-C trace hooks)
# --------------------------------------------------------------------------- #

def test_pdf_module_request_produces_worked_steps_bundle(tmp_path):
    # A pdf-requesting module computation auto-emits the exhaustive .tex/.pdf bundle
    # (like HH), and the .tex names the resolution steps verbatim -- the acceptance
    # bar (an algebraist can replay by hand). Byte-stability is unaffected: this
    # activates only on artifacts.pdf=True, which no delegation golden uses.
    body = {"schema": 2, "algebra": _LOOP,
            "compute": ["projective_resolution:0..3"],
            "artifacts": {"pdf": True, "tikz": False}, "module": _M2}
    r = run_spec(ComputeRequest.model_validate(body), tmp_path)
    assert r["meta"]["pdf"] in ("trace.pdf",
                                "PDF toolchain (pdflatex/tectonic) not found -- "
                                "worked steps in trace_steps.html")
    tex = tmp_path / "trace.tex"
    assert tex.exists(), "the downloadable .tex must be written (Plan 30 C1)"
    txt = tex.read_text(encoding="utf-8")
    assert "projective" in txt and "P_" in txt


# --------------------------------------------------------------------------- #
# Resolution `summands` field (Marco #3): P_1^2 (+) P_3 as LaTeX, per term
# --------------------------------------------------------------------------- #

def test_projective_resolution_has_latex_summands():
    b = _run(_LOOP, ["projective_resolution:0..3"], _M2)["results"]["projective_resolution"]
    assert b["summands"] == ["P_{1}", "P_{1}", "P_{1}", "P_{1}"]
    # the raw fields stay for backward-compat (only the RENDERING drops columns)
    assert "betti" in b and "terms" in b and len(b["terms"]) == len(b["summands"])


def test_injective_resolution_has_latex_summands():
    b = _run(_LOOP, ["injective_resolution:0..2"], _M2)["results"]["injective_resolution"]
    assert b["summands"] == ["I_{1}", "I_{1}", "I_{1}"]


def test_multiplicity_and_ordering_in_summands_latex():
    # P_1 (+) P_1 (+) P_2 renders as "P_{1}^{2} \\oplus P_{2}" (sorted, exponent).
    from quiverlab.hpc.spec import _summands_latex
    assert _summands_latex([1, 2, 1], "P") == "P_{1}^{2} \\oplus P_{2}"
    assert _summands_latex([3, 1], "I") == "I_{1} \\oplus I_{3}"
    assert _summands_latex([], "P") == "0"


# --------------------------------------------------------------------------- #
# decompose kind
# --------------------------------------------------------------------------- #

@_needs_decompose
def test_decompose_splits_a_direct_sum():
    b = _run(_A2, ["decompose"], _S1S2)["results"]["decompose"]
    assert b["kind"] == "decompose" and b["side"] == "right"
    assert b["iso_classes"] == 2
    dvs = sorted(tuple(sorted(s["dim_vector"].items())) for s in b["summands"])
    assert dvs == [(("1", 0), ("2", 1)), (("1", 1), ("2", 0))]
    assert all(s["multiplicity"] == 1 and s["indecomposable"] for s in b["summands"])
    assert b["references"] == ["assem_book"] and len(b["citations"][0]) == 2


@_needs_decompose
def test_decompose_reports_multiplicities():
    # S_2 (+) S_2 : one iso class, multiplicity 2. dim[1]=0, so a's block is 2x0
    # (empty) -- omit the map; the runner fills the zero action.
    b = _run(_A2, ["decompose"], {"dims": {"2": 2}})["results"]["decompose"]
    assert b["iso_classes"] == 1 and b["summands"][0]["multiplicity"] == 2


@_needs_decompose
def test_decompose_is_a_loud_typed_refusal_when_uncertifiable():
    # char 2 <= dim 2 (k[x]/x^3): the engine cannot certify and refuses LOUDLY ->
    # a client-safe QuiverlabError-typed refusal (a clean 4xx, never a 500).
    with pytest.raises(RunError) as exc:
        _run(_LOOP, ["decompose"], _M2)
    assert is_safe_error_type(exc.value.error_type), exc.value.error_type


# --------------------------------------------------------------------------- #
# tau / tau^- input certificate (Marco #1)
# --------------------------------------------------------------------------- #

@_needs_decompose
def test_tau_of_indecomposable_input_is_certified():
    b = _run(_A2, ["tau"], {"builtin": {"kind": "simple", "vertex": 2}})["results"]["tau"]
    assert b["indecomposable"] is True
    assert "decomposition" not in b and "note_key" not in b


@_needs_decompose
def test_tau_of_decomposable_input_carries_decomposition_and_note():
    b = _run(_A2, ["tau"], _S1S2)["results"]["tau"]
    assert b["indecomposable"] is False
    assert b["note_key"] == "mod.tau_additive"
    assert len(b["decomposition"]) == 2
    assert all(s["indecomposable"] for s in b["decomposition"])


@_needs_decompose
def test_tau_minus_is_also_certified():
    b = _run(_A2, ["tau_minus"], _S1S2)["results"]["tau_minus"]
    assert b["indecomposable"] is False and b["note_key"] == "mod.tau_additive"


def test_tau_omits_certificate_when_uncertifiable_byte_stable():
    # char 2 <= dim 2: certification is LOUD, so the tau block OMITS the certificate
    # (no false claim) and stays byte-identical to the pre-Plan-30 shape -- tau
    # itself is still computed. This is what keeps the module_basic golden green.
    b = _run(_LOOP, ["tau"], _M2)["results"]["tau"]
    assert "indecomposable" not in b and "decomposition" not in b
    assert set(b) == {"kind", "side", "is_zero", "dimvec", "dim",
                      "references", "citations"}


# --------------------------------------------------------------------------- #
# tor kind: honest either way (engine present or absent)
# --------------------------------------------------------------------------- #

def test_tor_dispatch_is_honest_about_the_engine():
    compute, module = ["tor:0..2"], {"dims": {"1": 1, "2": 1}, "maps": {"a": [[1]]}}
    tor_target = {"dims": {"1": 1}}                 # side defaults to left
    if _HAVE_TOR:
        b = _run(_A2, compute, module, tor_target=tor_target)["results"]["tor"]
        assert b["kind"] == "tor" and b["top"] == 2
        assert isinstance(b["dims"], list) and "target" in b
        assert b["references"] == ["minimal_resolution", "module_ext"]
    else:
        with pytest.raises(RunError) as exc:
            _run(_A2, compute, module, tor_target=tor_target)
        assert is_safe_error_type(exc.value.error_type), exc.value.error_type
        assert "Tor" in exc.value.message and "Plan 29" in exc.value.message
