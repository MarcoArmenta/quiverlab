"""Shared derived-fingerprint block for the three runners (Plan 43 / Task 6).

One scalar kind on the algebra block (schema v1). Built once here so the HPC
runner (``hpc/spec.py``), the Pyodide twin (``docs/gui/runner.py``) and the webapp
runner return byte-identical blocks. The two-algebra compare panel is deferred to
P50 (it needs a second-algebra request field -- a schema change)."""
from quiverlab.derived.fingerprint import derived_fingerprint

_REFERENCES = ["happel_triangulated", "rickard_derived", "lenzing_delapena_spectral"]

_SCOPE = ("a derived-invariant fingerprint; equal values are a necessary condition "
          "for derived equivalence, not a proof")


def _fmt_dims(dims):
    return "[" + ", ".join(str(d) for d in dims) + "]"


def _fp_latex(fp):
    """A compact LaTeX rendering of the fingerprint tuple. A field captured as an
    error renders its message in text; no field is ever silently dropped."""
    lines = []

    def _row(label, val):
        if isinstance(val, dict) and "error" in val:
            lines.append(rf"\text{{{label}: (unavailable -- {val['error']})}}")
        else:
            lines.append(f"{label} &= {val}")

    cox = fp.get("coxeter_polynomial")
    if isinstance(cox, dict):
        _row(r"\text{Coxeter polynomial}", cox)
    else:
        lines.append(rf"\text{{Coxeter polynomial }} p(t) &= {cox}")
    _row(r"\det C", fp.get("cartan_det"))
    smith = fp.get("cartan_smith")
    if isinstance(smith, dict):
        _row(r"\text{Cartan Smith factors}", smith)
    else:
        lines.append(rf"\text{{Cartan Smith factors}} &= {_fmt_dims(smith)}")
    lines.append(rf"\dim HH^\bullet &= {_fmt_dims(fp['hh_cohomology_dims'])}")
    lines.append(rf"\dim HH_\bullet &= {_fmt_dims(fp['hh_homology_dims'])}")
    lines.append(rf"\dim HC_\bullet &= {_fmt_dims(fp['cyclic_dims'])}")
    lines.append(rf"\dim Z(A) &= {fp['center_dim']}")
    lines.append(rf"\text{{gl.dim}} &: \text{{{fp['gl_dim']}}}")
    return r"\begin{aligned}" + r" \\ ".join(lines) + r"\end{aligned}"


def derived_fingerprint_block(A, top=4):
    fp = derived_fingerprint(A, top)
    return {"kind": "derived_fingerprint", "top": top, "fingerprint": fp,
            "latex": _fp_latex(fp), "scope": _SCOPE,
            "references": list(_REFERENCES)}   # citations added by the caller
