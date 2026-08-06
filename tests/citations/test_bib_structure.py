"""Structural gate on references.bib: every entry closes its own braces.

Merge-train line-dedup has twice eaten a lone ``}`` closer, silently merging
two BibTeX entries into one (Voigt1977 swallowed DlabRingel1989; earlier,
armenta_coxeter_calculus swallowed IgusaTodorov2005) -- the swallowed key
vanishes and the swallowing key renders the WRONG bibliographic data. The
formatter tolerates it, so nothing failed loudly until a frozen golden
happened to cite the damaged key. This gate makes the damage class a direct
test failure (contract & infrastructure, unmarked per the Plan-32 ruling).
"""
import pathlib
import re

_BIB = (pathlib.Path(__file__).resolve().parents[2]
        / "src" / "quiverlab" / "citations" / "references.bib")


def _entries():
    text = _BIB.read_text(encoding="utf-8")
    return [("@" + chunk) for chunk in re.split(r"(?m)^@", text)[1:]]


def test_every_entry_is_brace_balanced():
    bad = [e.splitlines()[0] for e in _entries()
           if e.count("{") != e.count("}")]
    assert not bad, f"unbalanced bib entries (a merge ate a closer?): {bad}"


def test_no_entry_header_inside_another_entry():
    # a swallowed closer shows up as a second `@type{key,` line inside one chunk
    bad = [e.splitlines()[0] for e in _entries()
           if len(re.findall(r"(?m)^@\w+\{", e)) != 1]
    assert not bad, f"nested entry headers (a merge ate a closer?): {bad}"
