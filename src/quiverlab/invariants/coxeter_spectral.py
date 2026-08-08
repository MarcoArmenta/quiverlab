"""Plan 58 (record R20): the certified Coxeter *spectral report* on an Algebra.

Assembles the shipped, exact spectral primitives -- ``coxeter_polynomial`` /
``coxeter_matrix`` (``invariants/cartan.py``), ``spectral_radius`` /
``mahler_measure`` / ``is_cyclotomic_product`` and the Plan-58
``certify_real_algebraic`` / ``cyclotomic_factorization`` / ``off_circle_root_count``
helpers (``invariants/spectral.py``) -- into a single report:

  * the exact ``ZZ[x]`` cyclotomic factorization with ``Phi_n`` labels;
  * a ``cyclotomic`` / ``quasi_unipotent`` verdict and the finite ``coxeter_order``
    (``Phi^m = I``, verified by EXACT matrix power) or ``None`` with an honest reason;
  * the exact count of roots outside the unit circle;
  * the spectral radius ``rho`` and Mahler measure ``M`` as CERTIFIED ALGEBRAIC
    NUMBERS (minimal polynomial + rational isolating interval + root index) -- never a
    float -- OR a loud per-field refusal on a complex-dominant spectrum;
  * the class-conditional Lehmer note (documentation only).

P58 adds NO new root-finding or spectral machinery: it is a certification +
labelling + assembly layer. See docs/plans/2026-08-07-plan-58-coxeter-spectral.md.

THE COMPLEX-DOMINANT GATE (§D3, load-bearing). ``rho`` / ``M`` are always non-negative
reals, EXACTLY computed by the shipped ``spectral_radius`` / ``mahler_measure`` (sound
for complex off-circle roots). But ``sympy.minimal_polynomial`` HANGS (~121 s, measured
on a dim-5 algebra) -- it does NOT raise -- on the complex-modulus form
``sqrt(CRootOf*CRootOf)`` that a complex-dominant spectrum yields, so an
exception-based guard never fires and the offline tier (no wall-clock kill) hangs on a
SMALL quiver. Instead we decide real-dominance up front, deterministically and cheaply
(~0.4 ms), with the shipped predicate ``_real_roots_suffice(_noncyclotomic_part(chi))``,
and only then call ``certify_real_algebraic``. All-cyclotomic short-circuits to
``rho = M = 1``. The outside-circle count is reported in BOTH branches."""

import math

import sympy as sp

from quiverlab.errors import QuiverlabError

_COMPLEX_DOMINANT_MSG = (
    "complex-dominant spectrum: ρ = |complex root|; "
    "certified-algebraic-number extraction is out of scope; "
    "outside-circle count reported")

_REFERENCES = ["dlPena2014mahler", "dlPena2013cyclotomic",
               "dlPenaTakane1990spectral", "lenzing_delapena_spectral", "assem_book"]

_SCOPE = (
    "exact spectral analysis of the Coxeter polynomial χ = charpoly(-C^{-T}C); "
    "ρ and M are certified algebraic numbers (minimal polynomial + rational "
    "isolating interval). The finite/tame/wild reading is a theorem on the hereditary "
    "locus (see form_type).")

# The Lehmer dichotomy is a CLASS-CONDITIONAL THEOREM (accessible algebras), not a
# computed verdict, and the general Lehmer problem is open -- quiverlab emits no verdict
# on it. mu_0 is referenced by its MINIMAL POLYNOMIAL (never a float; §D4).
_LEHMER_NOTE = (
    "de la Peña (arXiv:1310.1910): for an accessible algebra, M(χ) is either "
    "1 or ≥ μ₀, where μ₀ is the Mahler measure of Lehmer's "
    "polynomial -- the algebraic number whose minimal polynomial is "
    "t^10 + t^9 - t^7 - t^6 - t^5 - t^4 - t^3 + t + 1 (the smallest known Salem number). "
    "quiverlab does not decide accessibility and does not rule on the open general "
    "Lehmer problem; this documents the class-conditional dichotomy only.")


def coxeter_spectral(A) -> dict:
    """The certified Coxeter spectral report dict (no ``kind``/``references``/
    ``citations``/``scope``/``lehmer_class_note`` -- the raw analysis; see
    :func:`coxeter_spectral_block` for the runner block). Any field that raises on
    this input (singular / absent Cartan -> ``coxeter_polynomial``; complex-dominant
    or otherwise uncertifiable ``rho`` / ``M`` -> §D3) is captured as
    ``{"error": <the loud message>}`` per field -- never a crash, never a silent
    default (the Plan-30 tau-block / P43 fingerprint precedent)."""
    from quiverlab.invariants.spectral import (
        spectral_radius, mahler_measure, certify_real_algebraic,
        cyclotomic_factorization, off_circle_root_count,
        _noncyclotomic_part, _real_roots_suffice)
    from quiverlab.engine.coxeter_spectrum import is_cyclotomic_product

    report = {}
    # 1. Coxeter polynomial (guarded: singular / absent Cartan refuses loudly).
    try:
        chi = A.coxeter_polynomial().as_expr()
    except QuiverlabError as exc:
        msg = str(exc)
        for key in ("coxeter_polynomial", "coxeter_latex", "cyclotomic",
                    "factorization", "quasi_unipotent", "coxeter_order",
                    "coxeter_order_reason", "outside_unit_circle_count",
                    "spectral_radius", "mahler_measure"):
            report[key] = {"error": msg}
        return report
    report["coxeter_polynomial"] = str(chi)
    report["coxeter_latex"] = sp.latex(chi)

    # 2. cyclotomic verdict + Phi_n-labelled factorization.
    cyclotomic = bool(is_cyclotomic_product(chi))
    report["cyclotomic"] = cyclotomic
    report["quasi_unipotent"] = cyclotomic          # eigenvalue-level; n&s from factors
    report["factorization"] = cyclotomic_factorization(chi)

    # 3. outside-unit-circle root count (always reported, exact).
    report["outside_unit_circle_count"] = off_circle_root_count(chi)

    # 4. rho + M via the deterministic real-dominant gate (§D3).
    q = _noncyclotomic_part(chi)
    if q is None or _real_roots_suffice(q):
        # all-cyclotomic (rho = M = 1) or real-dominant: minimal_polynomial is
        # sub-second; certify. (certify_real_algebraic is guarded per field.)
        for key, alpha in (("spectral_radius", spectral_radius(chi)),
                           ("mahler_measure", mahler_measure(chi))):
            try:
                report[key] = certify_real_algebraic(alpha)
            except QuiverlabError as exc:
                report[key] = {"error": str(exc)}
    else:
        report["spectral_radius"] = {"error": _COMPLEX_DOMINANT_MSG}
        report["mahler_measure"] = {"error": _COMPLEX_DOMINANT_MSG}

    # 5. finite Coxeter order Phi^m = I (§D1): lcm of the recognised indices, VERIFIED
    #    by an exact matrix power (never is_diagonalizable).
    if not cyclotomic:
        report["coxeter_order"] = None
        report["coxeter_order_reason"] = (
            "not quasi-unipotent: some eigenvalue has modulus > 1 (ρ > 1), so the "
            "Coxeter transformation has infinite order")
    else:
        ns = [f["cyclotomic_index"] for f in report["factorization"]
              if isinstance(f.get("cyclotomic_index"), int)]
        m = math.lcm(*ns) if ns else 1
        try:
            Phi = sp.Matrix(A.coxeter_matrix())
            if Phi**m == sp.eye(Phi.rows):
                report["coxeter_order"] = m
                report["coxeter_order_reason"] = (
                    "Φ^%d = I (finite order %d; minimal by construction as the lcm "
                    "of the cyclotomic indices)" % (m, m))
            else:
                report["coxeter_order"] = None
                report["coxeter_order_reason"] = (
                    "quasi-unipotent but the Coxeter transformation has a nontrivial "
                    "Jordan block ⇒ infinite order (typical of tame/affine type)")
        except QuiverlabError as exc:
            report["coxeter_order"] = {"error": str(exc)}
            report["coxeter_order_reason"] = {"error": str(exc)}
    return report


def coxeter_spectral_block(A) -> dict:
    """The no-code ``coxeter_spectral`` compute block (schema v1): the
    :func:`coxeter_spectral` report plus ``{"kind", "references", "scope",
    "lehmer_class_note", "latex"}``. SHARED by both runners so their blocks are
    byte-identical; each runner adds ``citations`` from ``references``."""
    block = dict(coxeter_spectral(A))
    block["kind"] = "coxeter_spectral"
    block["references"] = list(_REFERENCES)
    block["scope"] = _SCOPE
    block["lehmer_class_note"] = _LEHMER_NOTE
    block["latex"] = _summary_latex(block)
    return block


def _summary_latex(report) -> str:
    """A compact ``\\begin{aligned}`` summary of the report (exact latex forms only;
    empty when the Coxeter polynomial itself refused)."""
    cp = report.get("coxeter_polynomial")
    if isinstance(cp, dict):                         # coxeter_polynomial refused
        return ""
    lines = [r"\chi(t) &= %s" % report.get("coxeter_latex", "")]
    lines.append(r"\text{cyclotomic} &= \text{%s}"
                 % ("yes" if report.get("cyclotomic") else "no"))
    order = report.get("coxeter_order")
    if isinstance(order, int):
        lines.append(r"\text{Coxeter order} &= %d" % order)
    elif order is None:
        lines.append(r"\text{Coxeter order} &= \infty")
    lines.append(r"\#\{\,|z|>1\,\} &= %s" % report.get("outside_unit_circle_count"))
    for label, key in ((r"\rho", "spectral_radius"), ("M", "mahler_measure")):
        v = report.get(key)
        if isinstance(v, dict) and "latex" in v:
            lines.append(r"%s &= %s" % (label, v["latex"]))
    return r"\begin{aligned} %s \end{aligned}" % (r" \\ ".join(lines))
