"""Certified arrow removal / addition Hochschild reductions (Plan 72, R6).

Cibils–Lanzilotta–Marcos–Solotar, ``1812.07655`` / Proc. Amer. Math. Soc. 148
(2020) 2421–2432. An **inert arrow** (CLMS **Def. 3.1**) is an arrow of ``Q`` that
appears in **no** minimal relation of ``I``. Deleting a set ``D`` of inert arrows
gives ``B = A ∖ D = kQ'/I'`` (``Q' = Q ∖ D``, ``I' = I`` verbatim) and CLMS prove a
**clean homology isomorphism** ``HH_n(A) ≅ HH_n(B)`` for ``n ≥ 2`` (**Thm 3.2**) —
the certificate asserts this range ONLY. ``HH_0 = A/[A,A]`` is **provably**
unchanged (an inert arrow never lies on an oriented cycle of a finite-dimensional
``kQ/I``, so no cyclic-path class moves); ``HH_1`` is reported informationally. By
contrast **cohomology carries an Ext correction** (**Thm 4.2**) nonzero for
``n ≥ 2``, distinct from the ``n = 0, 1`` center/disconnection deltas.

The dual **arrow addition** (**Thm 3.5 / 3.6**) builds ``A = B_F = T_B(N)``,
finite-dimensional **iff** adding ``F`` creates no relative cycle.

**Shared substrate (P73 seam).** ``inert_arrows`` / ``remove_arrows`` /
``add_arrows`` are the reduction primitives P73 (Han recognizer + Jacobi–Zariski)
consumes — no Han logic here. Whichever plan's constructor lands first OWNS it;
P72's ``remove_arrows`` is the **inert special case** of P73's more general
``arrow_removal_subalgebra(A, F)`` (arbitrary ``F``, induced ideal ``I_B = ker π``).
The cross-consistency test ``remove_arrows ≡ arrow_removal_subalgebra`` on inert
``F`` lives in whichever plan merges second (DD-B2).
"""
from dataclasses import dataclass, field
from fractions import Fraction

from quiverlab.errors import NotFiniteDimensionalError, QuiverlabError


# --------------------------------------------------------------------------- #
# inert arrows (Def. 3.1)
# --------------------------------------------------------------------------- #
def _relation_words(A):
    """The set of arrow names appearing in ANY term of ANY relation of ``I``
    (monomial AND binomial/non-monomial — every ``(coeff, word)`` is scanned)."""
    used = set()
    for rel in A.relations:
        for _coeff, word in rel.terms:
            used.update(word)
    return used


def _require_presented(A, what):
    if A.quiver is None or A.relations is None:
        raise QuiverlabError(
            f"{what} needs a quiver presentation (to read Def. 3.1 off the "
            f"relations); this algebra is structure-constant only",
            hint="present A via Quiver(...).algebra(...)")


def inert_arrows(A):
    """The arrows of ``Q`` that appear in no minimal relation of ``I`` (CLMS Def.
    3.1), read off ``A.relations``. Structure-constant-only ⇒ loud."""
    _require_presented(A, "inert_arrows")
    used = _relation_words(A)
    return [a for a in A.quiver.arrows if a not in used]


def inert_certificate(A):
    """``{arrow: witness}`` for the arrows deleted: an inert arrow's witness is
    'appears in no relation'; a non-inert arrow's witness names a relation it
    appears in (used for the loud refusal message)."""
    used_by = {}
    for rel in A.relations:
        for _coeff, word in rel.terms:
            for a in word:
                used_by.setdefault(a, repr(rel))
    cert = {}
    for a in A.quiver.arrows:
        cert[a] = used_by.get(a)              # None => inert
    return cert


# --------------------------------------------------------------------------- #
# relation -> grammar string (I' = I verbatim)
# --------------------------------------------------------------------------- #
def _coeff_sign_mag(coeff):
    if isinstance(coeff, Fraction):
        sign = -1 if coeff < 0 else 1
        mm = abs(coeff)
        mag = str(mm.numerator) if mm.denominator == 1 else f"{mm.numerator}/{mm.denominator}"
        return sign, mag
    return 1, str(int(coeff))


def _relation_to_grammar(rel):
    """A ``Relation`` back to a ``combinat.relations`` grammar string (coefficient
    BEFORE the arrows, signs folded into ' + ' / ' - '). Mirrors
    ``families.trivial_extension._relation_string`` so ``I' = I`` rebuilds verbatim
    over QQ and GF(p)."""
    parts = []
    for idx, (coeff, word) in enumerate(rel.terms):
        path = "*".join(word)
        sign, mag = _coeff_sign_mag(coeff)
        if idx == 0:
            sep = "" if sign > 0 else "-"
        else:
            sep = " + " if sign > 0 else " - "
        parts.append(f"{sep}{path}" if mag == "1" else f"{sep}{mag}*{path}")
    return "".join(parts)


# --------------------------------------------------------------------------- #
# deletion / addition
# --------------------------------------------------------------------------- #
def remove_arrows(A, arrows):
    """``B = A ∖ arrows`` with ``Q' = Q ∖ arrows`` and ``I' = I`` verbatim. Loud
    unless every named arrow is inert (Def. 3.1), naming the witnessing relation."""
    from quiverlab.combinat.quiver import Quiver
    _require_presented(A, "remove_arrows")
    arrows = list(arrows)
    cert = inert_certificate(A)
    for a in arrows:
        if a not in A.quiver.arrows:
            raise QuiverlabError(f"remove_arrows: {a!r} is not an arrow of Q",
                                 hint=f"arrows are {sorted(A.quiver.arrows)}")
        if cert.get(a) is not None:
            raise QuiverlabError(
                f"remove_arrows: arrow {a!r} appears in relation {cert[a]} — not "
                f"inert; deletion would change the ideal",
                hint="only inert arrows (Def. 3.1) may be removed")
    drop = set(arrows)
    new_arrows = {n: e for n, e in A.quiver.arrows.items() if n not in drop}
    rel_strings = [_relation_to_grammar(rel) for rel in A.relations]
    Q2 = Quiver(list(A.quiver.vertices), new_arrows)
    return Q2.algebra(relations=rel_strings, field=A.domain)


def _vertex_idems(A):
    labels = A.basis_labels or []
    out = {}
    for v in A.quiver.vertices:
        lab = f"e_{v}"
        if lab in labels:
            out[v] = A._basis_vec(labels.index(lab))
    return out


def _corner_nonzero(A, idems, i, j):
    """Is ``e_i A e_j ≠ 0`` (a nonzero A-path from ``i`` to ``j``)?"""
    dom = A.domain
    ei, ej = idems[i], idems[j]
    for k in range(A.dim):
        x = A.multiply(ei, A.multiply(A._basis_vec(k), ej))
        if any(not dom.is_zero(v) for v in x):
            return True
    return False


def add_arrows(A, new_arrows):
    """The dual ``A_F = T_A(N)`` (Thm 3.5 / 3.6): add the arrows ``new_arrows``
    (``{name: (s, t)}``) to ``Q``, keeping ``I`` verbatim. Loud refusal when ``F``
    creates a **relative cycle** (Thm 3.6: a directed cycle in the link graph — a
    link ``a → b`` exists when ``e_{t(a)} A e_{s(b)} ≠ 0`` — makes ``A_F``
    infinite-dimensional)."""
    from quiverlab.combinat.quiver import Quiver
    _require_presented(A, "add_arrows")
    new_arrows = dict(new_arrows)
    idems = _vertex_idems(A)
    # relative-cycle gate: DFS for a directed cycle in the link graph on new arrows.
    names = list(new_arrows)
    adj = {a: [] for a in names}
    for a in names:
        _sa, ta = new_arrows[a]
        for b in names:
            sb, _tb = new_arrows[b]
            if _corner_nonzero(A, idems, ta, sb):
                adj[a].append(b)
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {a: WHITE for a in names}
    cyclic = [False]

    def dfs(u):
        color[u] = GRAY
        for w in adj[u]:
            if color[w] == GRAY:
                cyclic[0] = True
            elif color[w] == WHITE:
                dfs(w)
        color[u] = BLACK

    for a in names:
        if color[a] == WHITE:
            dfs(a)
    if cyclic[0]:
        raise QuiverlabError(
            "add_arrows: adding these arrows creates a relative cycle — A_F would "
            "be infinite-dimensional (Thm 3.6)",
            hint="the new arrows link through A-paths into an oriented cycle")
    merged = dict(A.quiver.arrows)
    merged.update(new_arrows)
    rel_strings = [_relation_to_grammar(rel) for rel in A.relations]
    Q2 = Quiver(list(A.quiver.vertices), merged)
    try:
        return Q2.algebra(relations=rel_strings, field=A.domain)
    except NotFiniteDimensionalError as exc:
        raise QuiverlabError(
            "add_arrows: A_F is infinite-dimensional (a relative cycle survived)",
            hint="Thm 3.6 finiteness fails for this arrow set") from exc


# --------------------------------------------------------------------------- #
# the certified HH reduction
# --------------------------------------------------------------------------- #
@dataclass
class ArrowRemovalReport:
    algebra: object                    # A
    reduced: object                    # B = A ∖ D
    removed: list                      # D (the inert arrows deleted)
    inert_certificate: dict            # {arrow: "appears in no relation"} witnesses
    hom_A: list                        # HH_n(A)
    hom_B: list                        # HH_n(B)
    hom_iso_from: int                  # 2 (HH_n iso for n >= 2, Thm 3.2 -- the ONLY asserted range)
    hom_agrees: bool                   # HH_n(A) == HH_n(B) for n >= 2
    hom_low_agrees: list               # [HH_0 same, HH_1 same]  (HH_0 provably True)
    coh_A: list = None                 # HH^n(A)
    coh_B: list = None                 # HH^n(B)
    coh_correction: list = None        # degree-indexed (0..top): HH^n(A)-HH^n(B) for n>=2, 0 for n<2
    coh_low_delta: list = None         # [HH^0(A)-HH^0(B), HH^1(A)-HH^1(B)]  (center/disconnection)
    status: str = "complete"           # "complete" | "unsupported" | "error"
    note: str = ""


def arrow_removal(A, arrows=None, top=6, *, side="both"):
    """The certified HH reduction along inert-arrow deletion (CLMS Thm 3.2 / 4.2).
    ``arrows=None`` auto-detects and removes ALL inert arrows (the maximal
    reduction). Reports the clean ``HH_{≥2}`` homology isomorphism (Thm 3.2 — the
    ONLY asserted range; ``HH_0`` provably invariant, ``HH_{0,1}`` reported) and the
    honest cohomology split: ``coh_correction`` (``n ≥ 2``, the Thm 4.2 Ext term)
    separated from ``coh_low_delta`` (``n = 0, 1``, center/disconnection)."""
    _require_presented(A, "arrow_removal")
    if arrows is None:
        arrows = inert_arrows(A)
    arrows = list(arrows)
    cert_full = inert_certificate(A)
    removed_cert = {a: "appears in no relation" for a in arrows}
    B = remove_arrows(A, arrows)
    rep = ArrowRemovalReport(
        algebra=A, reduced=B, removed=arrows, inert_certificate=removed_cert,
        hom_A=[], hom_B=[], hom_iso_from=2, hom_agrees=False, hom_low_agrees=[])
    if side in ("both", "homology"):
        rep.hom_A = A.hochschild_homology(top, verbose=False).dims
        rep.hom_B = B.hochschild_homology(top, verbose=False).dims
        rep.hom_agrees = rep.hom_A[2:] == rep.hom_B[2:]
        rep.hom_low_agrees = [rep.hom_A[0] == rep.hom_B[0], rep.hom_A[1] == rep.hom_B[1]]
    if side in ("both", "cohomology"):
        rep.coh_A = A.hochschild_cohomology(top, verbose=False).dims
        rep.coh_B = B.hochschild_cohomology(top, verbose=False).dims
        rep.coh_correction = [0, 0] + [rep.coh_A[n] - rep.coh_B[n]
                                       for n in range(2, top + 1)]
        rep.coh_low_delta = [rep.coh_A[0] - rep.coh_B[0], rep.coh_A[1] - rep.coh_B[1]]
    return rep


# --------------------------------------------------------------------------- #
# GUI / webapp block builder (algebra-only, budget-carrying kind)
# --------------------------------------------------------------------------- #
_ARROW_REFERENCES = ("clms_arrow_removal", "han_conjecture")


def arrow_removal_block(A, top=6):
    """The ``arrow_removal`` compute-kind block (all three tiers). Auto-detects the
    inert arrows, builds ``B = A ∖ (inert)``, and renders the certificate. Loud
    refusals become a clean ``status``/``error`` block (never a 500)."""
    block = {"kind": "arrow_removal", "references": list(_ARROW_REFERENCES)}
    try:
        rep = arrow_removal(A, top=top)
    except QuiverlabError as exc:
        block.update(status="unsupported", error=str(exc),
                     note="arrow-removal reduction unavailable for this input")
        return block
    block.update(
        status=rep.status,
        algebra=repr(A).splitlines()[0],
        reduced=repr(rep.reduced).splitlines()[0],
        removed=rep.removed,
        inert_certificate=rep.inert_certificate,
        hom_A=rep.hom_A, hom_B=rep.hom_B,
        hom_iso_from=rep.hom_iso_from, hom_agrees=rep.hom_agrees,
        hom_low_agrees=rep.hom_low_agrees,
        coh_A=rep.coh_A, coh_B=rep.coh_B,
        coh_correction=rep.coh_correction, coh_low_delta=rep.coh_low_delta,
        note=rep.note,
    )
    return block
