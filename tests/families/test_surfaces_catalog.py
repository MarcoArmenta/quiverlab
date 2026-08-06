"""Surface constructors are catalog-listed + top-level exported, skipped in the webapp
family form (non-scalar args, the zoo precedent) (Plan 48)."""
import quiverlab as ql
from quiverlab.families.discover import families


def test_surface_constructors_catalogued_and_exported():
    names = families().names()
    for nm in ("fan_triangulation", "annulus_triangulation", "jacobian_of"):
        assert nm in names
        assert getattr(ql, nm, None) is not None          # top-level export contract


def test_webapp_form_skips_surface_inputs():
    from webapp.server.catalog import _iter_families
    listed = {name for name, _ in _iter_families()}
    assert "fan_triangulation" not in listed and "jacobian_of" not in listed
    assert "annulus_triangulation" not in listed
