"""Shared derived-fingerprint + silting blocks for the three runners (Plan 43 / Task 6,
Plan 67 / Task 5).

Scalar kinds on the algebra block (schema v1). Built once here so the HPC runner
(``hpc/spec.py``), the Pyodide twin (``docs/gui/runner.py``) and the webapp runner
return byte-identical blocks. The two-algebra compare panel is deferred to P50 (it needs
a second-algebra request field -- a schema change)."""
from quiverlab.derived.fingerprint import compare_fingerprints, derived_fingerprint

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


# --------------------------------------------------------------------------- #
# derived_compare: the two-algebra panel deferred by P43. It fingerprints BOTH
# algebras and compares them with the library's HONEST verdict -- "distinguished
# by <field>" / "not distinguished by these invariants", NEVER an equivalence
# claim (equal fingerprints are a necessary, not sufficient, condition). Shared
# by both runners (byte-identical twins); ``references`` present, citations added
# by the caller, exactly like every other block here.
# --------------------------------------------------------------------------- #

_COMPARE_SCOPE = ("equal fingerprints are a NECESSARY condition for derived "
                  "equivalence, not a proof; a difference in any field is a rigorous "
                  "obstruction to derived equivalence, but agreement never certifies it")


def _verdict_text(cmp):
    """The library's honest verdict as display text: 'distinguished by <fields>' when
    any field differs, else the verbatim 'not distinguished by these invariants'. Never
    an equivalence claim."""
    if cmp["distinguished_by"]:
        return "distinguished by " + ", ".join(cmp["distinguished_by"])
    return "not distinguished by these invariants"


def derived_compare_block(A, B, top=4):
    """The ``derived_compare`` block: the derived fingerprints of ``A`` and ``B`` side
    by side plus the honest comparison verdict (P43 ``compare_fingerprints``). Never
    asserts (in)equivalence -- a distinguishing field is an obstruction, agreement is
    only a necessary condition."""
    fa = derived_fingerprint(A, top)
    fb = derived_fingerprint(B, top)
    cmp = compare_fingerprints(fa, fb)
    return {"kind": "derived_compare", "top": top,
            "fingerprint_a": fa, "fingerprint_b": fb,
            "distinguished_by": cmp["distinguished_by"],
            "incomparable_fields": cmp["incomparable_fields"],
            "verdict": cmp["verdict"],           # the library's raw verdict string
            "verdict_text": _verdict_text(cmp),  # display form: "distinguished by ..."
            "scope": _COMPARE_SCOPE,
            "references": list(_REFERENCES)}     # citations added by the caller


# --------------------------------------------------------------------------- #
# silting: the algebra-level Plan-67 block (verifier verdict on the regular object +
# single-mutation neighbours + a bounded-radius exploration + the co-t-structure record).
# Schema v1 scalar kind (radius, budget), sized on A.dim; citations added by the caller.
# --------------------------------------------------------------------------- #
_SILTING_REFERENCES = ["aihara_iyama_silting", "oppermann_silting_quivers",
                       "jorgensen_cotstructures"]

_SILTING_SCOPE = (
    "silting objects in K^b(proj A) (Aihara-Iyama): presilting is DECIDED on the exact "
    "positive window Hom_{D^b}(T, T[n>0]) = 0; generation is three-valued -- True on the "
    "tilting / 2-term / local classes, else 'unknown' (det(g_proj) = +-1 in the projective "
    "K0 basis is necessary but NOT sufficient; thick subcategories are not K0-classified, "
    "Krah phantom). The silting quiver can be INFINITE (kA2 already is) and transitivity is "
    "proven only for local / hereditary / canonical (AI Thm 1.2), so the exploration is "
    "bounded-radius with a loud status -- certified complete ONLY for local, never a general "
    "enumeration. The co-t-structure is a documentation record (coheart add(T)), not a "
    "computed subcategory.")


def _silting_key_str(key):
    """A byte-stable string for a shift-insensitive silting key (a frozenset of summand
    fingerprints): sorted so it is deterministic across runs (a frozenset repr is not)."""
    return "|".join(sorted(str(part) for part in key))


def _summand_dimvecs(cx):
    """Per-degree dim-vectors of a perfect-complex summand, JSON-stable (string keys)."""
    return {str(d): {str(v): int(m) for v, m in cx.term(d).dimension_vector().items() if m}
            for d in cx.degrees()}


def silting_block(A, radius=3, budget=64):
    """The ``silting`` block for algebra ``A`` (Plan 67): the verifier verdict on the
    regular object ``A = (+)_v P_v``, its single-mutation neighbours (left + right), a
    bounded-radius exploration (loud status; complete only for local), and the
    co-t-structure record. ``references`` are citation KEYS; the caller adds ``citations``."""
    from quiverlab.modules.complexes import ChainComplex
    from quiverlab.derived.silting import (is_silting_object, silting_neighbors,
                                           bounded_silting_exploration, co_t_structure_of)
    verts = list(A.quiver.vertices)
    n = len(verts)
    regular = [ChainComplex.stalk(A.projective(v), 0) for v in verts]
    rep = is_silting_object(regular)
    block = {
        "kind": "silting", "n": n,
        "regular": {"is_silting": rep.is_silting, "is_presilting": rep.is_presilting,
                    "k0_basis": rep.k0_basis,
                    "generation_certified_by": rep.generation_certified_by,
                    "window": list(rep.window),
                    "g_matrix": [list(row) for row in rep.g_matrix], "det": rep.det},
        "scope": _SILTING_SCOPE,
        "references": list(_SILTING_REFERENCES),   # citation KEYS; caller adds "citations"
    }
    neighbors = []
    for direction in ("left", "right"):
        for (i, mut) in silting_neighbors(regular, direction):
            if mut is None:                        # mutation degenerates -- never a silent skip
                neighbors.append({"summand": i, "direction": direction,
                                  "is_silting": None, "summand_dimvecs": None})
                continue
            neighbors.append({
                "summand": i, "direction": direction,
                "is_silting": is_silting_object(mut).is_silting,
                "summand_dimvecs": [_summand_dimvecs(c) for c in mut]})
    block["neighbors"] = neighbors
    ex = bounded_silting_exploration(A, radius=radius, budget=budget)
    block["exploration"] = {
        "status": ex.status, "radius": ex.radius, "finite_class": ex.finite_class,
        "vertices": [{"key": _silting_key_str(v["key"]), "is_initial": v["is_initial"],
                      "silting": v["silting"]} for v in ex.vertices],
        "arrows": [{"from": i, "to": j, "direction": lab["direction"],
                    "summand": lab["summand"]} for (i, j), lab in sorted(ex.arrows.items())],
    }
    ct = co_t_structure_of(regular)
    block["co_t_structure"] = {
        "bounded": ct["bounded"], "coheart_is_addT": ct["coheart_is_addT"],
        "aisle": ct["aisle"], "coaisle": ct["coaisle"], "references": list(ct["references"]),
        "coheart": [{"dimvecs": {str(d): {str(v): int(m) for v, m in dv.items() if m}
                                 for d, dv in c["dimvecs"].items()}}
                    for c in ct["coheart"]],
    }
    return block
