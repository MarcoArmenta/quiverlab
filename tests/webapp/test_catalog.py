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
