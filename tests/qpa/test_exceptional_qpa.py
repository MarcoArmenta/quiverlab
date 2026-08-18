"""QPA (GAP) probe for the exceptional-sequence surface (Plan 65 / R27+R28).

QPA 1.37 has NO exceptional-sequence surface: no braid mutation, no classical
exceptional enumeration, no tau-exceptional / Jasso reduction. NamesGVars() carries
no matching name. So the SEQUENCE surface is theory/self-cert-oracled (the count
`n!*h^n/|W|` classical + `n!*#sTt` tau, the materialisation cross-check, the braid-
orbit transitivity certificate) and this probe is an honest fail-if-appears guard.

The exceptional MODULES, however, ARE crosscheckable: an exceptional module is a
rigid BRICK (End = k, Ext^1(M,M) = 0), and QPA's ExtAlgebraGenerators recomputes
Ext^*(M,M) -- an input-level anchor on the Dynkin zoo (every kA_n indecomposable is
a rigid brick, live-verified `#exceptional == #indec`).

qpa-marked (dir): skips locally without GAP, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Names QPA might plausibly use for an exceptional-sequence surface. NOTE: we do NOT
# scan for "braid" -- GAP core ships BraidGroup (Artin groups), an unrelated hit.
_EXC_NAMES = (
    "ExceptionalSequence",
    "ExceptionalSequences",
    "IsExceptionalSequence",
    "CompleteExceptionalSequences",
    "ExceptionalModules",
    "TauExceptionalSequence",
    "TauExceptionalSequences",
    "TauPerpendicularReduction",
)


def test_qpa_has_no_exceptional_sequence_surface():
    lg = session.libgap_handle()

    # (1) Scan NamesGVars() FIRST -- IsBoundGlobal("Foo") REGISTERS "Foo" into
    #     NamesGVars(), so scanning after the queries would echo them back and
    #     falsely "find" a surface (the Plan-35 products-probe precedent).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    exc_like = sorted(n for n in gvar_names
                      if "exceptionalseq" in n.lower() or "tauexceptional" in n.lower())

    # (2) None of the named exceptional-sequence entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _EXC_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not exc_like:
        pytest.skip(
            "QPA 1.37 exposes no exceptional-sequence surface (no ExceptionalSequence* / "
            "TauExceptional*; NamesGVars has no matching name). Covering oracle: the Plan-65 "
            "classical closed-form n!*h^n/|W| + tau signed-count n!*#sTt (materialised) "
            "batteries + the braid-orbit transitivity certificate; the exceptional MODULES "
            "are anchored below via ExtAlgebraGenerators.")

    # If QPA ever grows an exceptional-sequence surface, FAIL loudly so a real
    # crosscheck is wired in and this honest-scope skip is revisited.
    pytest.fail(
        "QPA now exposes an exceptional-sequence surface -- wire a real crosscheck against "
        f"the Plan-65 classical/tau enumerators. Found present={present}, names={exc_like}.")


def test_exceptional_modules_are_bricks_no_self_ext_via_qpa():
    """Input-level anchor (MINOR b): every kA_3 indecomposable is an exceptional MODULE = a
    rigid BRICK, and QPA independently confirms BOTH halves PER MODULE -- ``dim End(M) = 1``
    (QPA ``HomOverAlgebra``, the Plan-37 bridge) and ``dim Ext^1(M, M) = 0`` (QPA
    ``ExtOverAlgebra``, the Plan-27 idiom). The exceptional SEQUENCE surface has no QPA
    verb (the guard above), but the exceptional-module criterion does."""
    from quiverlab.families.basic import linear_path_algebra
    from quiverlab.fields import QQ
    from quiverlab.modules.exceptional import is_exceptional_module
    from quiverlab.qpa import scripts
    from quiverlab.qpa.crosscheck import _graded, crosscheck_hom_glue
    A = linear_path_algebra(3, field=QQ)
    for rec in A.ar_quiver().vertices:                    # every kA_3 indecomposable
        M = rec["module"]
        assert is_exceptional_module(A, M) is True        # our BRICK criterion (M3)
        rep = crosscheck_hom_glue(A, M, M)                # dim End(M) = dim Hom(M,M) vs QPA
        rep.assert_agree()
        assert rep.ours["hom_dim"] == 1                   # brick: End_A(M) = k (QPA-confirmed)
        dv, arr = _graded(A, M)                           # rigidity: dim Ext^1(M,M) = 0 via QPA
        script = (scripts.quiver_and_algebra_script(A) + "\n"
                  + scripts.module_decl(A, dv, arr, "MM")
                  + "\nLength(ExtOverAlgebra(MM, MM)[2]);")
        assert int(session.run(script)) == 0 == A.ext(M, M, 1)
