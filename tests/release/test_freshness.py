"""Plan 08 freshness gate as a test. `check()` returns only PREREQUISITE drift
(Plans 03-07 surfaces + docs sources); the deprecated-license line is an
informational Task-3 TODO handled by `license_todo()` and never a STOP, so this
test asserts prerequisite freshness directly and is agnostic to whether the license
has been SPDX-fixed yet."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent.parent / "scripts"))
import release_freshness as rf  # noqa: E402


def test_prerequisite_surfaces_present():
    drift = rf.check()
    assert drift == [], "Plan 08 prerequisites drifted:\n" + "\n".join(drift)
