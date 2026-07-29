import pytest

from webapp.server.catalog import build_catalog, family_names, validate_family, CatalogError


def test_catalog_lists_families():
    cat = build_catalog()
    names = {f["name"] for f in cat["families"]}
    assert "QuantumCI" in names
    assert cat["fields"]["GF"]["needs"] == ["p"]


def test_catalog_excludes_zoo():
    # zoo yields a batch of algebras, not one buildable family (amendment pt 2).
    cat = build_catalog()
    assert "zoo" not in {f["name"] for f in cat["families"]}


def test_validate_known_family_ok():
    validate_family("QuantumCI", {"q": 1})  # no raise


def test_validate_unknown_family_raises():
    with pytest.raises(CatalogError):
        validate_family("no_such_family", {})


def test_validate_unknown_param_raises():
    with pytest.raises(CatalogError):
        validate_family("QuantumCI", {"nope": 1})


# --------------------------------------------------------------------------- #
# Curated form metadata (Marco 2026-07-28: explain every parameter, prefill
# defaults). Every family must carry a bilingual summary; every parameter a
# bilingual help line; and the PREFILL the form derives (default, else example,
# empty omitted) must actually BUILD through the shared spec core -- the form
# never offers a family nobody can compute with (TrivialExtension was listed
# but unbuildable until the wrapper-family fix).
# --------------------------------------------------------------------------- #
def _prefill_params(fam: dict) -> dict:
    out = {}
    for p in fam["params"]:
        v = p.get("default")
        if v is None:
            v = p.get("example")
        if v is None or v is False:
            continue
        out[p["name"]] = v
    return out


def test_every_family_has_bilingual_summary_and_param_help():
    for fam in build_catalog()["families"]:
        assert set(fam.get("summary", {})) == {"en", "es"}, fam["name"]
        for p in fam["params"]:
            assert set(p.get("help", {})) == {"en", "es"}, (fam["name"], p["name"])


def test_every_family_prefill_builds():
    from quiverlab.hpc.spec import build_algebra
    for fam in build_catalog()["families"]:
        params = _prefill_params(fam)
        A = build_algebra({"kind": "family", "family": fam["name"],
                           "params": params,
                           "field": {"kind": "GF", "p": 5, "n": 1}})
        assert A.dim > 0, (fam["name"], params)
