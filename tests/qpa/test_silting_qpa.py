"""QPA has NO silting surface (Plan 67 honest scope). QPA 1.37 ships no silting-object /
silting-mutation / silting-quiver surface -- so there is nothing to crosscheck the P67
verifier / mutation / exploration against. This probe SKIPS with an honest message, and
FAILS LOUDLY if a silting surface ever appears in QPA (so the verification-page
honest-scope entry stays truthful) -- the Plan-35 no-surface precedent.

qpa-marked: skips locally, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import pytest

from quiverlab.qpa import session

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")

# Every name QPA might plausibly use for a silting surface.
_SILTING_NAMES = (
    "SiltingObjects",
    "SiltingMutation",
    "SiltingQuiver",
    "SiltingComplexes",
    "IsSiltingComplex",
    "SiltingMutationComplex",
    "MutationOfSiltingObject",
)


def test_qpa_has_no_silting_surface():
    lg = session.libgap_handle()

    # (1) Scan the global name table FIRST -- for no name mentioning "silting". This must
    #     precede any IsBoundGlobal query: in GAP, IsBoundGlobal("Foo") REGISTERS "Foo"
    #     into NamesGVars() (as an unbound known name), so scanning after the queries would
    #     echo them back and falsely "find" a surface (the Plan-35 precedent).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    silting_like = sorted(n for n in gvar_names if "silting" in n.lower())

    # (2) None of the named silting entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")')) for name in _SILTING_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not silting_like:
        pytest.skip(
            "QPA 1.37 exposes no silting surface (no SiltingObjects / SiltingMutation / "
            "SiltingQuiver; NamesGVars has no 'silting' name). Covering oracles for P67: "
            "the self-cert batteries (presilting window, mutant re-verification, "
            "involution) + the AI/Oppermann literature pins + the P45 2-term cross-check."
        )

    # If a future QPA ever grows a silting surface, FAIL loudly so this honest-scope skip
    # is revisited and a real crosscheck is wired in + the verification page updated.
    pytest.fail(
        "QPA now exposes a silting surface -- wire a live crosscheck against the Plan-67 "
        "verifier / mutation / exploration and update the verification page's honest-scope "
        f"entry. Found: present={present}, silting_like_names={silting_like}.")
