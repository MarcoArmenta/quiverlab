"""The homological string-algebra test (Plan 59 / R34).

Suarez-Alvarez (Alg. Rep. Theory 26 (2023) 1759-1772, arXiv:2105.02948): among
REPRESENTATION-FINITE algebras, A is a STRING algebra iff the middle term of EVERY
extension of indecomposable modules has at most two indecomposable summands. The
quantifier is over all Ext^1(X, Y) classes, NOT just AR sequences.

Honest semi-decision (see the Plan-59 design note); verdict in
{"not_string", "string_over_ground_field", "inconclusive"}:
  * REFUTE ("not_string", k-bar-sound): a >= 3-summand middle E. "String" is a property
    of (Q, I) and is field-independent, and a >= 3 decomposition base-changes up to
    k-bar, so a ground-field witness certifies non-string. Only the terminal split step
    is char-robust; the REFUTE PIPELINE (knit / almost_split / decompose) needs char 0
    or char > dim, so refute runs over QQ (small primes only on tiny algebras). Checked
    classes: the split (2), every AR sequence, and a basis of each Ext^1.
  * CONFIRM ("string_over_ground_field", INCONCLUSIVE w.r.t. the k-bar theorem): finite
    ground field only -- exhaustive scalar-class enumeration of every Ext^1(X, Y), each
    middle certified <= 2 (needs char > dim E or char 0). Over k-bar there may be more
    indecomposables / non-rational classes, so this does NOT prove the theorem's
    "string"; the definitive string answer is the syntactic is_string.
  * else INCONCLUSIVE, loud reason.

CONTRACT (one thing): RAISES iff a >= 3 witness is found AND is_string(A) is True (one
engine is a bug -- k-bar-sound). Otherwise RETURNS the verdict; a
"string_over_ground_field" with is_string False records a kbar_gap_note (never raised).
Refuses loudly for rep-infinite / self-injective / presentation-less A (out of scope --
the theorem is stated for rep-finite algebras).
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field as _dc_field

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.recognizers import is_string
from quiverlab.modules.ar import (_combine_matrices, _ext1_data,
                                   almost_split_sequence, knit_ar_quiver)
from quiverlab.modules.decompose import decompose
from quiverlab.modules.morphism import ModuleHom
from quiverlab.modules.ses import ShortExactSequence
from quiverlab.modules.yoneda import baer_extension

_REFS = ("suarez_alvarez", "butler_ringel", "assem_book")


@dataclass(frozen=True)
class StringHomologyVerdict:
    verdict: str                 # not_string | string_over_ground_field | inconclusive
    is_string: bool              # the P38 syntactic arbiter (definitive k-bar "string")
    witness: "dict | None"       # unified schema (see _mk_witness)
    kbar_gap_note: "str | None"  # confirm-vs-is_string mismatch (never raised)
    checked: dict
    reason: "str | None"
    references: tuple = _dc_field(default=_REFS)


# --------------------------------------------------------------------------- #
# scope / field helpers
# --------------------------------------------------------------------------- #
def _field_order(dom):
    """Finite-field size (p for GF(p), p^n for GF(p^n)) or None for an infinite field."""
    q = getattr(dom, "q", None)
    if q is None:
        q = getattr(dom, "p", None)
    return q


def _require_in_scope(A, knit_budget_modules, knit_budget_dim):
    if A.quiver is None:
        raise QuiverlabError(
            "homological_string_test: needs a quiver-presented algebra",
            hint="the extension/middle-term test walks the knitted indecomposables; "
                 "structure-constant algebras carry no quiver")
    ar = knit_ar_quiver(A, budget_modules=knit_budget_modules,
                        budget_dim=knit_budget_dim)
    if not ar.is_complete or ar.status != "complete":
        raise QuiverlabError(
            "homological_string_test: A is not representation-finite in scope "
            f"(AR knit status={ar.status!r}); the Suarez-Alvarez characterization is "
            "stated for rep-finite algebras",
            hint="self-injective and rep-infinite inputs are refused (knit "
                 "unsupported/budget) -- this test does not extend past them")
    return ar


# --------------------------------------------------------------------------- #
# middle terms of extensions
# --------------------------------------------------------------------------- #
def _ext1_pack(A, X, Y):
    """(cocycle_mats, terms, dmats, width) for Ext^1(X, Y). ``len(cocycle_mats)`` is
    dim Ext^1(X, Y) (the cohomology class basis); ``width`` = # columns of ``P_1``."""
    cocycle_mats, _cob, terms, dmats = _ext1_data(A, X, Y)
    width = len(dmats[1][0]) if (len(dmats) > 1 and dmats[1] and dmats[1][0]) else 0
    return cocycle_mats, terms, dmats, width


def _middle_from_pack(X, Y, coeffs, pack, dom):
    """The middle E of the Ext^1(X, Y) class with coordinates ``coeffs`` over the cocycle
    basis (the almost_split_sequence idiom). Self-certifies: ``baer_extension`` refuses a
    non-cocycle, the Yoneda sequence asserts exactness, and a ShortExactSequence wrapper
    re-certifies 0 -> Y -> E -> X -> 0."""
    cocycle_mats, terms, dmats, width = pack
    f = _combine_matrices(cocycle_mats, coeffs, Y.dim, width, dom)
    seq = baer_extension(X, Y, f, terms, dmats)          # 0 -> Y -> E -> X -> 0
    seq.assert_exact()
    E = seq.modules[1]
    ShortExactSequence(
        ModuleHom(seq.modules[0], E, seq.maps[0], check=False),
        ModuleHom(E, seq.modules[2], seq.maps[1], check=False))   # check=True: exactness
    return E


def _middle_of_class(A, X, Y, coeffs):
    """The middle E of the Ext^1(X, Y) class ``coeffs`` (fresh pack)."""
    return _middle_from_pack(X, Y, coeffs, _ext1_pack(A, X, Y), A.domain)


# --------------------------------------------------------------------------- #
# counting summands
# --------------------------------------------------------------------------- #
def _decompose_count(E):
    """(count, dimvecs): the EXACT number of indecomposable summands WITH multiplicity,
    and their dim vectors (repeated by multiplicity). May raise QuiverlabError when
    ``decompose`` cannot certify (char <= dim E, no split)."""
    parts = decompose(E)
    count = 0
    dvs = []
    for mod, mult in parts:
        dv = mod.dimension_vector()
        dvs.extend([dv] * mult)
        count += mult
    return count, dvs


def _bounded_split_count(E, budget=512):
    """A LOWER BOUND on the number of indecomposable summands via the bounded Fitting
    split (``decompose``'s ``_try_split``), never raising. Used only when ``decompose``
    refuses on char grounds -- a >= 3 split found here still base-changes up to k-bar, so
    it is a valid REFUTE witness in any characteristic."""
    from quiverlab.modules.decompose import _try_split
    from quiverlab.modules.hom import hom_space
    stack = [E]
    leaves = []
    while stack:
        M = stack.pop()
        if M.dim == 0:
            continue
        H = hom_space(M, M)
        split = _try_split(M, H, budget)
        if split is None:
            leaves.append(M)
        else:
            stack.extend(split)
    return len(leaves), [m.dimension_vector() for m in leaves]


def _summands_at_least(E, k):
    """(count, dimvecs) for the REFUTE side: exact via ``decompose``; on a char-blocked
    refusal, a bounded Fitting-split LOWER BOUND. A ``count >= k`` is a valid witness in
    any characteristic (a found split base-changes up)."""
    try:
        return _decompose_count(E)
    except QuiverlabError:
        return _bounded_split_count(E)


def _confirm_summands(E):
    """For the CONFIRM side, which needs a CERTIFIED bound:
    ("ok", count) | ("witness", count, dimvecs) | ("undecidable",). ``decompose``
    refusing (char <= dim E) makes the class char-undecidable -- never a guess."""
    try:
        count, dvs = _decompose_count(E)
    except QuiverlabError:
        return ("undecidable",)
    if count >= 3:
        return ("witness", count, dvs)
    return ("ok", count)


# --------------------------------------------------------------------------- #
# witness schema
# --------------------------------------------------------------------------- #
def _fmt_coeffs(coeffs):
    return [str(c) for c in coeffs]


def _mk_witness(X, Y, klass, coeffs, E, cnt, dvs):
    """Unified witness schema across BOTH refute sub-routes. ``klass`` in
    {"ar_socle","basis"}; ``coeffs`` is None on the AR route (the almost-split socle
    class is not a coordinate vector), the unit/combination vector on the basis route."""
    return {"X": X.dimension_vector(), "Y": Y.dimension_vector(),
            "class": klass, "coeffs": (None if coeffs is None else _fmt_coeffs(coeffs)),
            "E_dimvec": E.dimension_vector(),
            "summand_count": cnt, "summand_dimvecs": list(dvs)}


# --------------------------------------------------------------------------- #
# the exhaustive finite-ground-field confirm
# --------------------------------------------------------------------------- #
def _scalar_class_reps(d, dom, p):
    """Representatives of P^{d-1}(GF(p)): every nonzero vector up to a nonzero scalar
    (leading nonzero coordinate = 1). Middles of scalar-multiple cocycles are isomorphic
    extensions, so one rep per scalar class is exhaustive."""
    for vec in itertools.product(range(p), repeat=d):
        lead = next((x for x in vec if x != 0), None)
        if lead != 1:                            # skip zero (lead is None) and non-unit leads
            continue
        yield [dom.coerce(x) for x in vec]


def _exhaustive_confirm(A, indec, budget_classes):
    """Finite prime-field exhaustive middle check.
    ("all_le_2", n) | ("witness", witness) | ("undecidable", pair|None) | ("budget", n).
    """
    from quiverlab.fields.primefield import PrimeField
    dom = A.domain
    if not isinstance(dom, PrimeField):
        return ("undecidable", None)             # non-prime finite field: not enumerated here
    p = dom.p
    n_checked = 0
    for X in indec:
        for Y in indec:
            pack = _ext1_pack(A, X, Y)
            d = len(pack[0])
            if d == 0:
                continue
            n_reps = (p ** d - 1) // (p - 1)
            if n_checked + n_reps > budget_classes:
                return ("budget", n_checked)
            for coeffs in _scalar_class_reps(d, dom, p):
                n_checked += 1
                E = _middle_from_pack(X, Y, coeffs, pack, dom)
                res = _confirm_summands(E)
                if res[0] == "undecidable":
                    return ("undecidable",
                            (X.dimension_vector(), Y.dimension_vector()))
                if res[0] == "witness":
                    return ("witness",
                            _mk_witness(X, Y, "basis", coeffs, E, res[1], res[2]))
    return ("all_le_2", n_checked)


# --------------------------------------------------------------------------- #
# the public test
# --------------------------------------------------------------------------- #
def _raise_contradiction(witness):
    raise QuiverlabError(
        "homological_string_test: found an extension of indecomposables with a "
        ">= 3-summand middle, so A is NOT a string algebra, yet is_string(A) is True -- "
        "one engine is a bug",
        hint=f"witness middle dim vector {witness['E_dimvec']} "
             f"({witness['summand_count']} summands)")


def homological_string_test(A, budget_classes=4096, budget_pairs=4096,
                            knit_budget_modules=128, knit_budget_dim=32):
    """Three-valued homological string test; see the module docstring for the contract.

    ``knit_budget_modules`` / ``knit_budget_dim`` cap the AR knit that decides rep-
    finiteness. A rep-INFINITE algebra produces unboundedly large indecomposables, so a
    modest ``knit_budget_dim`` trips the knit (status="budget") and the test refuses
    quickly; a rep-finite algebra whose indecomposables all fit under the caps knits
    complete. The defaults cover the standard small rep-finite zoo; a genuine rep-finite
    algebra with an indecomposable of dim > ``knit_budget_dim`` is refused (raise the cap
    for it) rather than silently mis-classified.
    """
    ar = _require_in_scope(A, knit_budget_modules, knit_budget_dim)
    syn = is_string(A)
    dom = A.domain
    q = _field_order(dom)
    indec = [v["module"] for v in ar.vertices]
    route = []
    n_classes = 0
    witness = None

    # -- REFUTE (a): every AR-sequence middle (0 -> tauY -> E -> Y -> 0) --
    route.append("ar")
    for Y in indec:
        try:
            ses = almost_split_sequence(Y)       # raises on a projective end -> skip
        except QuiverlabError:
            continue
        cnt, dvs = _summands_at_least(ses.M, 3)
        if cnt >= 3:
            witness = _mk_witness(Y, ses.L, "ar_socle", None, ses.M, cnt, dvs)
            break

    # -- REFUTE (b): a basis of each Ext^1(X, Y), middle 0 -> Y -> E -> X -> 0 --
    if witness is None:
        route.append("basis")
        for X in indec:
            if witness is not None or n_classes > budget_classes:
                break
            for Y in indec:
                pack = _ext1_pack(A, X, Y)
                d = len(pack[0])
                for j in range(d):
                    coeffs = [dom.one() if i == j else dom.zero() for i in range(d)]
                    n_classes += 1
                    E = _middle_from_pack(X, Y, coeffs, pack, dom)
                    cnt, dvs = _summands_at_least(E, 3)
                    if cnt >= 3:
                        witness = _mk_witness(X, Y, "basis", coeffs, E, cnt, dvs)
                        break
                if witness is not None or n_classes > budget_classes:
                    break

    checked = {"n_indec": len(indec), "route": list(route),
               "field": getattr(dom, "name", str(dom)), "classes": n_classes}

    # -- THE one loud disagreement: a k-bar-sound non-string witness vs is_string True. --
    if witness is not None:
        if syn is True:
            _raise_contradiction(witness)
        return StringHomologyVerdict("not_string", syn, witness, None, checked, None)

    if n_classes > budget_classes:
        return StringHomologyVerdict(
            "inconclusive", syn, None, None, checked,
            "budget: the basis refute pass exceeded budget_classes before closing")

    # -- CONFIRM (finite ground field only): exhaustive scalar-class enumeration. --
    if q is None:
        return StringHomologyVerdict(
            "inconclusive", syn, None, None, checked,
            "infinite ground field: the quantifier over ALL Ext^1 classes cannot be "
            "exhausted; checked split/AR/basis only (no >= 3 witness). Defers to is_string.")

    status, data = _exhaustive_confirm(A, indec, budget_classes)
    if status == "witness":
        checked["route"].append("exhaustive")
        if syn is True:
            _raise_contradiction(data)
        return StringHomologyVerdict("not_string", syn, data, None, checked, None)
    if status == "all_le_2":
        checked["route"].append("exhaustive")
        checked["classes"] = data
        gap = (None if syn else
               "string_over_ground_field but is_string(A) is False: the finite ground "
               f"field GF({q}) did not witness the non-stringness the k-bar theorem "
               "guarantees (more indecomposables / non-rational classes over k-bar)")
        return StringHomologyVerdict(
            "string_over_ground_field", syn, None, gap, checked, None)
    if status == "undecidable":
        pair = "" if data is None else f" of Ext^1{data}"
        return StringHomologyVerdict(
            "inconclusive", syn, None, None, checked,
            f"char-undecidable: a middle{pair} could not be certified (char <= dim, or "
            "the finite field is non-prime); defers to is_string")
    # budget
    return StringHomologyVerdict(
        "inconclusive", syn, None, None, checked,
        f"budget: the exhaustive confirm exceeded budget_classes at {data} classes")


def string_homological_block(A):
    """The ``string_homological`` algebra-only compute kind (Task 5)."""
    try:
        v = homological_string_test(A)
    except QuiverlabError as exc:
        return {"error": str(exc), "references": list(_REFS)}
    return {"verdict": v.verdict, "is_string": v.is_string, "witness": v.witness,
            "kbar_gap_note": v.kbar_gap_note, "checked": v.checked, "reason": v.reason,
            "references": list(v.references)}
