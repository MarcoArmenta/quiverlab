"""Stable-category fractional Calabi-Yau dimension of self-injective algebras
(Plan 53 / R24).

In the stable module category ``mod-bar A`` of a self-injective ``A`` (triangulated,
Ivanov-Volkov arXiv:1212.2619 sec 1.3), the suspension is ``Sigma = Omega^{-1}`` (Heller's
cosyzygy) and the Serre functor is ``S = Omega . nu`` (``nu`` = the Nakayama functor
``D Hom_A(-, A)``). Following Herschend-Iyama-Oppermann / Keller, ``mod-bar A`` is
**fractionally Calabi-Yau of dimension m/ell** iff ``S^ell = Sigma^m``, i.e. the least
positive ``ell`` for which ``S^ell X ~ Sigma^m X`` for all objects ``X`` (some ``m``).
Because ``nu`` and ``Omega`` commute, ``S^ell = Omega^ell nu^ell``; fixing the sign to
match the ``ell=1`` Ivanov-Volkov integer criterion ``Omega^{n+1} ~ nu^{-1}`` gives the
certified identity, per generator ``X`` (stable iso, projective summands dropped):

    S^ell(X) ~ Sigma^{+m}(X)   <=>   Omega^{m+ell}( nu^ell(X) ) ~ X.

**Honest scope -- TIER 3 (weak-on-generators).** The certificate is checked OBJECT-WISE
on a generating set (the non-projective simples + the ``Omega^i nu^j``-orbit reps met in
the search), NOT as a natural isomorphism of functors. This is a NECESSARY condition for
the tier-2 weak-CY property (object-wise on ALL objects) and hence for the tier-1 STRONG
CY property (``S ~ Sigma^n`` as functors); "as functors" is reserved for the strong tier
and appears NOWHERE in what this module certifies. Object-wise agreement on a generating
set need not propagate to the whole category, so the search can UNDER-report BOTH the
numerator ``m`` (defined only mod the Sigma-period) AND the denominator ``ell`` (a
spurious small ``ell`` may pass on the generators while failing on a non-generator). The
tier-1 functorial certificate ``Omega^{n+1}_{A^e}(A) ~ (A^v)_phi`` (Ivanov-Volkov Thm 1.8,
via the minimal ``A^e`` engine) is the named successor on the verification page.

**The numerator ``m`` is a residue mod the Sigma-period ``p``** (``Sigma^p X ~ X`` on a
periodic stable category): ``S^ell X ~ Sigma^m X`` implies ``~ Sigma^{m+kp} X`` for all
``k``. The canonical representative is "smallest ``ell >= 1``, then smallest ``m >= 0``",
with ``p`` reported in ``sigma_period`` when the search determines it.

Leans entirely on shipped module operations -- ``nu``/``nu^-`` (``modules/ar.py``),
``Omega``/``Omega^-`` (``modules/resolution.py``), ``decompose`` (stable normal form),
``is_isomorphic`` (the loud-when-undecidable certificate) -- no new math engine. The
``decompose`` char caveat is load-bearing (loud over ``char <= dim`` on GF(p)); batteries
run over ``QQ`` / ``GF(32003)``."""
from __future__ import annotations

from dataclasses import dataclass
from math import gcd

from quiverlab.errors import QuiverlabError
from quiverlab.modules.decompose import decompose
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.homdims import _is_projective
from quiverlab.modules.morphism import direct_sum

_CHECKED_ON = ("simples + Omega/nu-orbit reps (object-wise on a generating set; "
               "anchored by literature + one (m,ell) for all generators; can "
               "under-report m AND ell)")
_TIER = "weak-on-generators"


@dataclass
class FractionalCY:
    """The stable-category fractional Calabi-Yau dimension ``(m, ell)`` of a
    self-injective algebra (Plan 53 / R24), certified at the **weak-on-generators tier**
    (object-wise on the simples + orbit reps -- a NECESSARY condition for weak/strong CY,
    never a functor iso; the ``tier`` field says so). ``cy_dimension`` is the reduced
    ``m/ell``; ``weakly_n_cy`` is the Ivanov-Volkov integer ``n = m`` when ``ell == 1``
    (else ``None``); ``sigma_period`` is the order of ``Sigma`` on the generators (the
    numerator's mod-ambiguity, ``None`` when the search did not determine it)."""
    m: "int | None"
    ell: "int | None"
    cy_dimension: "str | None"
    weakly_n_cy: "int | None"
    sigma_period: "int | None"
    status: str                       # "certified" | "budget" | "shift-trivial"
    tier: str                         # ALWAYS "weak-on-generators" (tier-3 certificate)
    checked_on: str
    certificate: list                 # per-generator [(name, "Omega^{m+ell} nu^ell ~ id")]

    def __repr__(self):
        if self.status == "certified":
            base = (f"stable CY dimension {self.cy_dimension} (weak-on-generators)")
            if self.ell == 1:
                base += f"; weakly {self.weakly_n_cy}-CY"
            return base
        if self.status == "shift-trivial":
            return ("the suspension Omega^{-1} is trivial (A is radical-square-zero "
                    "self-injective local); CY dimension degenerate, reported (0, 1)")
        return ("not certified within the (m,ell) search window (self-injective but the "
                "generators were not S-periodic in budget -- likely representation-"
                "infinite)")


def _stable_nf(X):
    """The stable normal form of ``X``: the direct sum of the NON-projective
    indecomposable summands of ``decompose(X)`` (projective summands dropped). A module
    all of whose summands are projective is ``0`` in the stable category, returned as a
    genuine dim-0 module (``Omega`` of the projective). Loud on the ``decompose`` char
    caveat (``char <= dim`` over GF(p)), unchanged."""
    if X.dim == 0:
        return X
    parts = []
    for S, mult in decompose(X):           # loud on char caveat
        if _is_projective(S):
            continue
        parts += [S] * mult
    if not parts:                          # X is projective => 0 in mod-bar A
        return X.syzygy()                  # Omega(projective) = 0 (a genuine dim-0 module)
    D, _, _ = direct_sum(*parts)
    return D


# a one-word alias the tests import as the module's stable normal form
_nf = _stable_nf


def _generators(A):
    """The generating set: the NON-projective simple modules (they generate ``mod-bar A``
    as a thick triangulated subcategory -- every module has a finite radical filtration by
    simples). A projective simple is dropped (it is ``0`` in the stable category)."""
    gens = []
    for v in A.quiver.vertices:
        S = A.simple(v)
        if not _is_projective(S):
            gens.append((v, S))
    return gens


def _Sfun(X):
    """The Serre functor value ``S X = Omega(nu X)`` in ``mod-bar A`` (stable normal form
    applied). ``nu`` and ``Omega`` are computed literally (the commuting ``nu Omega ~
    Omega nu`` is used only in the docstring derivation, never relied on in code)."""
    return _stable_nf(X.nakayama().syzygy())


def _sigma_period(gens, dim_budget, m_window):
    """The order ``p`` of ``Sigma = Omega^{-1}`` on the generators = the order of
    ``Omega`` on them (same order): least ``p >= 1`` with ``Omega^p g ~ g`` for every
    generator, searched to ``m_window``; ``None`` if not found in budget."""
    cur = [_stable_nf(S.syzygy()) for _v, S in gens]       # Omega^1
    for p in range(1, m_window + 1):
        if all(is_isomorphic(cur[i], gens[i][1]) for i in range(len(gens))):
            return p
        cur = [_stable_nf(x.syzygy()) for x in cur]
        if any(x.dim > dim_budget for x in cur):
            return None
    return None


def fractional_calabi_yau(A, *, ell_max=4, m_window=8, dim_budget=4096):
    """The stable-category fractional Calabi-Yau dimension ``(m, ell)`` of a
    self-injective algebra ``A`` (Plan 53 / R24; Ivanov-Volkov arXiv:1212.2619).

    Certified at the **weak-on-generators tier** (object-wise on the non-projective
    simples + the ``Omega/nu``-orbit reps -- a NECESSARY condition for weak/strong CY,
    NEVER a functor isomorphism; see the module docstring). Bounded ``(m, ell)`` search
    with the derived sign ``S^ell ~ Sigma^{+m}`` (``Omega^m(S^ell G) ~ G``, ``m >= 0``
    only -- the canonical representative for the mod-``p`` ambiguity). Returns a
    :class:`FractionalCY`:
      * ``status="certified"`` with ``(m, ell)``, reduced ``m/ell``, ``weakly_n_cy = m``
        iff ``ell == 1``, and ``sigma_period`` when the search found it;
      * ``status="shift-trivial"`` = ``(0, 1)`` when ``Sigma = Omega^{-1} ~ id`` on every
        generator (radical-square-zero self-injective local);
      * ``status="budget"`` when no ``(m, ell)`` is found in ``ell_max``/``m_window`` or a
        module iterate exceeds ``dim_budget`` (honest -- likely rep-infinite, whose
        modules need not be ``S``-periodic).

    RAISES ``QuiverlabError`` for non-self-injective ``A`` (``nu`` is not an
    autoequivalence and ``mod-bar A`` is not triangulated otherwise). The
    ``is_isomorphic`` loud refusal propagates unchanged (never a silent ``budget`` where
    an iso was undecidable); the ``decompose`` char caveat propagates too."""
    from quiverlab.modules.ext import is_selfinjective
    if A.quiver is None or not is_selfinjective(A):
        raise QuiverlabError(
            "fractional CY of the stable category is defined for self-injective A only; "
            "nu is not an autoequivalence and mod-bar A is not triangulated otherwise",
            hint="for a finite-global-dimension algebra see the derived-category "
                 "Calabi-Yau successor (the stable category is trivial there)")
    gens = _generators(A)
    if not gens:
        # every simple is projective => A is semisimple => mod-bar A is trivial (0).
        return FractionalCY(m=0, ell=1, cy_dimension="0/1", weakly_n_cy=0,
                            sigma_period=1, status="shift-trivial", tier=_TIER,
                            checked_on=_CHECKED_ON, certificate=[])

    # (3) shift-trivial: Sigma = Omega^{-1} ~ id on every generator.
    if all(is_isomorphic(_stable_nf(S.cosyzygy()), S) for _v, S in gens):
        cert = [(str(v), "Omega^{-1} ~ id (Sigma trivial)") for v, _S in gens]
        return FractionalCY(m=0, ell=1, cy_dimension="0/1", weakly_n_cy=0,
                            sigma_period=1, status="shift-trivial", tier=_TIER,
                            checked_on=_CHECKED_ON, certificate=cert)

    # (4) search (smallest ell, then smallest m >= 0) for S^ell(G) ~ Sigma^m(G).
    G = [S for _v, S in gens]
    for ell in range(1, ell_max + 1):
        T = list(G)
        for _ in range(ell):                        # T = S^ell(G)
            T = [_Sfun(t) for t in T]
            if any(t.dim > dim_budget for t in T):
                return _budget()
        cur = list(T)                               # Omega^0(T)
        for m in range(0, m_window + 1):
            if all(is_isomorphic(cur[i], G[i]) for i in range(len(G))):
                return _certified(A, gens, m, ell, dim_budget, m_window)
            cur = [_stable_nf(x.syzygy()) for x in cur]   # Omega^{m+1}(T)
            if any(x.dim > dim_budget for x in cur):
                return _budget()
    return _budget()


def _certified(A, gens, m, ell, dim_budget, m_window):
    g = gcd(m, ell) or 1
    cy = f"{m // g}/{ell // g}"
    p = _sigma_period(gens, dim_budget, m_window)
    cert = [(str(v), f"Omega^{{{m + ell}}} nu^{{{ell}}} ~ id") for v, _S in gens]
    return FractionalCY(m=m, ell=ell, cy_dimension=cy,
                        weakly_n_cy=(m if ell == 1 else None), sigma_period=p,
                        status="certified", tier=_TIER, checked_on=_CHECKED_ON,
                        certificate=cert)


def _budget():
    return FractionalCY(m=None, ell=None, cy_dimension=None, weakly_n_cy=None,
                        sigma_period=None, status="budget", tier=_TIER,
                        checked_on=_CHECKED_ON, certificate=[])


def is_fractionally_calabi_yau(A, **kw):
    """True iff :func:`fractional_calabi_yau` CERTIFIES a fractional CY dimension
    (``status == "certified"``) for the self-injective algebra ``A`` (Plan 53). The
    degenerate shift-trivial case and the budget cap are NOT reported as certified."""
    return fractional_calabi_yau(A, **kw).status == "certified"


# ---------------------------------------------------------------------------
# The `fractional_cy` no-code block (Plan 53 Task E) -- ONE shared library builder,
# so the server (hpc.spec) and Pyodide (docs/gui/runner) twins are byte-identical by
# construction. Algebra-level scalar kind (routes through _dispatch, NOT
# _dispatch_module; schema stays v1). A non-self-injective input is caught into
# {"error": ...} (a clean typed 4xx in the webapp, never a 500 -- the Plan-26 precedent).
# Each runner adds the resolved `citations`.
# ---------------------------------------------------------------------------
def fractional_cy_block(A):
    """The ``fractional_cy`` block for algebra ``A`` (Plan 53 / R24): the stable-category
    fractional Calabi-Yau dimension, certified at the weak-on-generators tier. The ``tier``
    field is ALWAYS ``"weak-on-generators"`` and is rendered beside the value so the
    displayed result never overstates the certificate. A non-self-injective input is
    caught into ``{"error": ...}`` (the loud not-self-injective message), never raised out
    of the block. Returns the block minus ``citations`` (each runner resolves
    ``references`` to citation pairs)."""
    refs = ["ivanov_volkov", "erdmann_skowronski_scy", "assem_book"]
    try:
        fcy = fractional_calabi_yau(A)
    except QuiverlabError as exc:
        return {"kind": "fractional_cy", "error": str(exc)}
    return {
        "kind": "fractional_cy",
        "cy_dimension": fcy.cy_dimension,
        "m": fcy.m, "ell": fcy.ell,
        "weakly_n_cy": fcy.weakly_n_cy,
        "sigma_period": fcy.sigma_period,
        "status": fcy.status,
        "tier": fcy.tier,
        "checked_on": fcy.checked_on,
        "certificate": [list(c) for c in fcy.certificate],
        "text": repr(fcy),
        "references": refs,
    }
