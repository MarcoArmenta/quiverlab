"""families() discoverability (spec §3.4)."""
from quiverlab import families


def test_catalog_lists_every_family_with_signature_and_route():
    listing = families()
    names = set(listing.names())
    assert {
        "NakayamaAlgebra", "PathAlgebra", "TruncatedPathAlgebra", "RadicalSquareZero",
        "IncidenceAlgebra", "QuantumCI", "ExteriorAlgebra", "PreprojectiveAlgebra",
        "TrivialExtension", "TensorProduct", "zoo",
    } <= names
    info = listing.by_name("QuantumCI")
    assert info.route == "general"
    assert "q=" in info.signature
    assert "quantum_ci" in info.citations
    s = str(listing)
    assert "NakayamaAlgebra" in s and "monomial" in s


def test_star_import_exposes_the_catalog():
    ns = {}
    exec("from quiverlab import *", ns)
    for name in ("NakayamaAlgebra", "QuantumCI", "zoo", "families", "bibliography"):
        assert name in ns
