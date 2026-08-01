"""The unified worked-steps step-event taxonomy (spec §3.8).

This module is the SINGLE import surface for trace events. It RE-EXPORTS the
event dataclasses shipped by Plan 03 (groebner) and Plan 04 (Chouhy-Solotar)
from their home modules -- it never redefines them -- and DEFINES the one new
event, RankStep, which carries a numeric differential matrix over a field plus
its rank (the bar/fast engines' notion of a differential; the CS
DifferentialEvent carries symbolic bimodule terms instead).

Reconciliation notes (see the plan's taxonomy table):
  * The requested `Differential` is Plan 04's DifferentialEvent.
  * AmbiguityEvent (Plan 04) is folded in.
  * Dispatch (Plan 03) is reused for BOTH the construction route (monomial vs
    groebner) AND the engine/resolution choice (bar/fast/bardzell/chouhy-solotar);
    same dataclass, distinguished by the value of `route`.
"""
from dataclasses import dataclass

from quiverlab.groebner.events import Dispatch, ReductionStep  # noqa: F401
from quiverlab.resolutions_cs.trace import (  # noqa: F401
    AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep,
)


@dataclass
class RankStep:
    """One rank computation of a differential matrix over a stated field.

    `matrix` is a list[list[str]] of domain-element string renderings (kept small),
    or None when `elided` is True (matrix larger than the elision threshold -- only
    shape + rank are retained). `side` is "cochain" (d^n : C^n -> C^{n+1}) or
    "chain" (b_n : C_n -> C_{n-1}). `nrows`/`ncols` are the matrix dimensions
    (= dim of the target/source cochain space)."""
    degree: int
    side: str
    nrows: int
    ncols: int
    rank: int
    field: str
    matrix: object = None
    elided: bool = False
    note: str = ""


# --------------------------------------------------------------------------- #
# Plan 30, Part C: module-computation worked-step events. These carry the actual
# ALGEBRA a module computation manipulates -- projective/injective resolution
# terms and their differentials as matrices, the Hom/tensor collapse of Ext/Tor
# with the per-degree rank bookkeeping, and free-form narrative lines (the chosen
# projective-cover generators, the D/Tr steps of an AR translate). They mirror the
# RankStep contract: a matrix is a list[list[str]] of domain-element renderings,
# or None with `elided=True` + a stated-shape `note` when it exceeds the threshold.
# --------------------------------------------------------------------------- #

@dataclass
class ModuleTerm:
    """One term of a module (co)resolution: the direct sum of indecomposable
    projectives (P_v) or injectives (I_v). `summands` is the vertex list WITH
    repetition (e.g. [1, 1, 3] renders as ``P_1^{2} + P_3``); `sym` is "P" or "I";
    `kind` is "projective" or "injective"; `dim` is the k-dimension; `dimvec` maps
    vertex -> multiplicity (a dict with str/int keys, or None)."""
    degree: int
    kind: str
    sym: str
    summands: object
    dim: int
    dimvec: object = None


@dataclass
class ModuleDifferential:
    """A differential of a module (co)resolution rendered AS A MATRIX over the
    stated field. `dom_summands`/`cod_summands` are the source/target summand
    vertex lists (with repetition); `cod_is_module` is True for the augmentation
    d_0: Q_0 -> M (whose target is the module M itself, not a projective term).
    `sym` is "P"/"I", `kind` "projective"/"injective". `matrix` is list[list[str]]
    or None (elided). `symbol` is the differential's TeX name (e.g. "d_1").

    `dom_is_module` / `cod_is_module` (Plan 34 adds the domain flag) render that
    endpoint as the module ``M`` itself rather than a direct sum of summands -- used for
    the self-maps of a module (an arrow action ``rho_M(a): M -> M`` or a splitting
    endomorphism), keeping the map declaration narrow. Both default False (byte-stable);
    the text/HTML renderers read only the summand lists, so a flagged endpoint falls back
    to its ``*_summands`` rendering there (still honest: M = (+)_v M e_v)."""
    degree: int
    kind: str
    sym: str
    symbol: str
    dom_summands: object
    cod_summands: object
    nrows: int
    ncols: int
    field: str
    cod_is_module: bool = False
    dom_is_module: bool = False
    mod_name: str = "M"
    matrix: object = None
    elided: bool = False
    note: str = ""


@dataclass
class ExtDegree:
    """One degree of an Ext/Tor computation over a minimal resolution: the
    collapsed Hom (Ext) / tensor (Tor) space dimension at degree n, the connecting
    map (delta^n for Ext / d_n for Tor) as a matrix, its rank and the neighbouring
    rank, and the resulting Ext^n / Tor_n dimension. `op` is "Ext" or "Tor"."""
    degree: int
    op: str
    space_dim: int
    rank_here: int
    rank_prev: int
    result_dim: int
    nrows: int
    ncols: int
    field: str
    matrix: object = None
    elided: bool = False
    note: str = ""


@dataclass
class StepNote:
    """A free-form narrative worked-step line (a projective-cover generator choice,
    the D/Tr steps of an AR translate, ...). `text` is the headline; `detail` is an
    optional indented continuation. Rendered verbatim (escaped per format).

    `heading` (Plan 34) marks a line that OPENS a new worked step: the LaTeX renderer
    turns it into a numbered ``\\paragraph{Step N. ...}`` run-in heading (homework
    style); the text/HTML renderers, which read only `text`/`detail`, show it as a
    plain paragraph (the flag defaults False, so every pre-Plan-34 StepNote renders
    exactly as before -- byte-stable)."""
    text: str
    detail: str = ""
    heading: bool = False


@dataclass
class ResultDims:
    """The AUTHORITATIVE final (co)homology dimensions the engine actually returned
    (``HHTable.dims`` / ``.kind``), recorded into the trace as the Result source (Plan 34
    fix). The three human renderers show THIS as the "Result" line (with the correct
    ``HH^`` / ``HH_`` variance) rather than re-deriving from the per-degree events, so a
    trace can never mislabel or misreport what the engine computed; ``render_text.derive_dims``
    stays a CROSS-CHECK (writer.py raises on drift). ``note`` carries an honest one-liner
    when the engine records no per-degree worked steps (e.g. the fast GF(p) engine), so
    the report still has a Result line without fabricating steps.

    ``kind`` is ``"HH^"`` (cohomology) or ``"HH_"`` (homology); ``dims`` is the list of
    integer dimensions in degree order; ``note`` defaults to ``""``."""
    kind: str
    dims: object
    note: str = ""


# --------------------------------------------------------------------------- #
# Plan 35: the HH-product worked-step event. One ``ProductStep`` is a single
# bidegree's block of a product chapter -- a structure-constant table
# (cup/cap/bracket) rendered as its nonzero-product EQUATION LINES, or one induced
# Connes differential B rendered as a MATRIX. It mirrors the other events'
# exact-string contract: every coefficient is a domain-element string, never a float.
# --------------------------------------------------------------------------- #

@dataclass
class ProductBasis:
    """Plan 35 UNIT 2: the EXPLICIT (co)cycle representatives of an HH-product chapter
    -- per (side, degree) the basis classes (each a labeled term-sum + a sparse
    coordinate vector), the ordered (co)chain enumeration the vectors index into, and
    the annihilating differential -- as the JSON-safe ``{side: {str(degree): ...}}``
    shape of ``quiverlab.hochschild.products.blocks()`` (captured at table-build time by
    ``quiverlab.hochschild.basis_reps``). Emitted ONCE per chapter, before the
    ``ProductStep`` tables; the renderer turns it into per-degree sub-sections (ordered
    basis -> explicit classes -> differential + a one-line verification sentence) that
    the structure-constant tables then reference.

    All three payload fields default ``None`` -- a legacy product object without the
    explicit-reps fields emits no ``ProductBasis`` (or one whose fields are ``None``),
    and the renderers fall back to the naming-only legend (tolerance)."""
    kind: str
    basis_classes: object = None
    chain_basis: object = None
    differentials: object = None


@dataclass
class ProductStep:
    """One block of an HH-product worked-steps chapter (cup / cap / bracket /
    connes_b). ``kind`` is the product kind. ``degrees`` is the bidegree tuple --
    ``(p, q)`` for cup/bracket, ``(p, n)`` for cap, ``(n,)`` for one Connes
    differential. ``heading`` is the TeX map label (e.g.
    ``HH^{p} \\otimes HH^{q} \\to HH^{p+q}``), typeset by the renderer. ``lines`` is
    the list of TeX equation lines spelling out the nonzero products in the recorded
    class bases (cup/cap/bracket); ``matrix`` is the induced-B matrix
    (``list[list[str]]`` of exact coefficient strings) for connes_b, else ``None``.
    ``note`` is an optional one-liner (a vanishing bidegree, the induced rank)."""
    kind: str
    degrees: object
    heading: str = ""
    lines: object = ()
    matrix: object = None
    note: str = ""


__all__ = [
    "Dispatch", "ReductionStep", "AmbiguityEvent", "ResolutionTerm",
    "DifferentialEvent", "LiftStep", "RankStep",
    "ModuleTerm", "ModuleDifferential", "ExtDegree", "StepNote", "ResultDims",
    "ProductStep", "ProductBasis", "ALL_EVENTS",
]

# The complete tuple of trace event types.  Renderers validate their input
# against this and refuse loudly on anything else -- a foreign object in an
# event stream is a caller bug (e.g. an unpacked (events, result) tuple) and
# silently skipping it would drop worked steps from the report.
ALL_EVENTS = (Dispatch, ReductionStep, AmbiguityEvent, ResolutionTerm,
              DifferentialEvent, LiftStep, RankStep,
              ModuleTerm, ModuleDifferential, ExtDegree, StepNote, ResultDims,
              ProductStep, ProductBasis)
