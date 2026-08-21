"""The READ-ONLY HanLab bank directory resolves from ``QUIVERLAB_BANK_DIR`` first
(v1.0.1 wave 2 item 2c).

The byte-level bank oracle (``test_battery_bank_oracle.py``, 12 tests marked
``oracle_literature`` + ``oracle_crossengine``) was pinned to a single laptop's
absolute path, so it skipped on CI and on every collaborator checkout while STILL
counting toward the audited oracle totals.  ``_resolve_bank`` now tries the
``QUIVERLAB_BANK_DIR`` env var, then a repo-relative sibling, then the original
absolute path -- so the oracle is reachable off Marco's machine, and his checkout
keeps working unchanged.

These are contract tests (UNMARKED, per the Plan-32 ruling): they must NOT carry an
oracle-class marker, or they would perturb the audited class counts the release gate
pins.  The bank need not be present for them to run -- they exercise the resolver,
not the oracle.
"""
import pathlib

from tests.resolutions_cs import test_battery_bank_oracle as bank


def _make_fake_bank(root):
    """A directory carrying the bank marker file ``hanlab/resolutions_cs.py``."""
    (root / "hanlab").mkdir(parents=True)
    (root / "hanlab" / "resolutions_cs.py").write_text("# stub bank\n", encoding="utf-8")
    return root


def test_env_var_is_honoured(tmp_path, monkeypatch):
    """QUIVERLAB_BANK_DIR wins over the sibling and the hardcoded absolute path."""
    fake = _make_fake_bank(tmp_path / "bank")
    monkeypatch.setenv("QUIVERLAB_BANK_DIR", str(fake))
    assert bank._resolve_bank() == fake


def test_env_var_takes_precedence_over_the_real_bank(tmp_path, monkeypatch):
    """Even when the hardcoded/sibling bank exists on this machine, an explicit
    QUIVERLAB_BANK_DIR pointing at a valid bank is used instead."""
    fake = _make_fake_bank(tmp_path / "elsewhere")
    monkeypatch.setenv("QUIVERLAB_BANK_DIR", str(fake))
    resolved = bank._resolve_bank()
    assert resolved == fake
    assert (resolved / "hanlab" / "resolutions_cs.py").exists()


def test_env_var_without_the_marker_falls_through(monkeypatch):
    """A QUIVERLAB_BANK_DIR that does not carry ``hanlab/resolutions_cs.py`` is not
    selected; resolution falls through to the sibling / absolute path."""
    monkeypatch.setenv("QUIVERLAB_BANK_DIR", "/definitely/not/a/bank/zzz123")
    resolved = bank._resolve_bank()
    assert resolved != pathlib.Path("/definitely/not/a/bank/zzz123")


def test_absent_env_still_yields_a_path_for_the_skipif(monkeypatch):
    """With no env var set, the resolver still returns a Path (the last candidate
    when nothing is reachable), so the module-level ``skipif`` can fire -- collection
    never errors, it SKIPS."""
    monkeypatch.delenv("QUIVERLAB_BANK_DIR", raising=False)
    assert isinstance(bank._resolve_bank(), pathlib.Path)
