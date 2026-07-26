"""Registry: quiverlab citation keys -> the papers behind each algorithm and family.
Annotations live HERE (the web /literature page and quiverlab.bibliography() consume
them). Loud failure on unknown keys (spec §3.9)."""
import difflib
import pathlib
import re
from dataclasses import dataclass, field

from quiverlab.errors import CitationError

_BIB = pathlib.Path(__file__).with_name("references.bib")


@dataclass(frozen=True)
class Reference:
    key: str            # the public quiverlab key (e.g. "bardzell")
    bibtex_key: str     # the @-entry id in references.bib (e.g. "Bardzell1997")
    kind: str           # "algorithm" | "family" | "field" | "foundation"
    title: str
    annotation: str     # one sentence: what it underpins
    tags: tuple = field(default_factory=tuple)


def _r(key, bibtex_key, kind, title, annotation, *tags):
    return Reference(key, bibtex_key, kind, title, annotation, tuple(tags))


REGISTRY: dict = {r.key: r for r in [
    _r("bardzell", "Bardzell1997", "algorithm",
       "Alternating syzygies of monomial algebras",
       "The minimal projective bimodule resolution for monomial algebras "
       "(quiverlab's Bardzell engine and the truncated/radical-square-zero families).",
       "resolution", "monomial"),
    _r("chouhy_solotar", "ChouhySolotar2015", "algorithm",
       "Projective resolutions via ambiguities",
       "The general kQ/I bimodule resolution from a reduction system -- quiverlab's "
       "Chouhy-Solotar engine for non-monomial algebras.",
       "resolution", "general"),
    _r("bracket_liftings", "NegronWitherspoon2016", "algorithm",
       "Gerstenhaber bracket via homotopy liftings",
       "The Gerstenhaber bracket computed directly on a non-bar resolution "
       "(with Volkov2019), transported onto bar representatives.",
       "bracket"),
    _r("bracket_liftings_volkov", "Volkov2019", "algorithm",
       "Gerstenhaber bracket on an arbitrary resolution",
       "A bracket formula valid on any projective bimodule resolution "
       "(companion to Negron-Witherspoon).",
       "bracket"),
    _r("minimal_resolution", "GSZ2001", "algorithm",
       "Minimal projective resolutions",
       "The Green-Solberg-Zacharia minimal module resolution algorithm "
       "(quiverlab's minimal engine and module Ext).",
       "resolution", "module"),
    _r("module_ext", "GSZ2001", "algorithm",
       "Module Ext via minimal resolutions",
       "Module-level Ext^n over minimal resolutions (Plan 05 module engine).",
       "module"),
    _r("bar", "Hochschild1945", "foundation",
       "Hochschild cohomology via the (normalized) bar complex",
       "Hochschild's original definition of the cohomology of an associative algebra; "
       "the normalized bar complex is quiverlab's HH^*/HH_* oracle in any characteristic.",
       "resolution", "bar"),
    _r("happel_question", "Happel1989", "foundation",
       "Happel's question",
       "Whether finite global dimension is equivalent to eventual vanishing of HH^n "
       "-- the motivating question for the hereditary and truncated families.",
       "conjecture"),
    _r("quantum_ci", "BGMS2005", "family",
       "Quantum complete intersections",
       "The algebra k<x,y>/(x^2, y^2, xy + q yx): finite Hochschild cohomology with "
       "infinite global dimension (the QuantumCI family).",
       "family"),
    _r("qci_hh_oracle", "BerghErdmann2008", "family",
       "Hochschild (co)homology of quantum complete intersections",
       "Explicit HH^* / HH_* of quantum complete intersections -- the literature "
       "oracle QuantumCI results are checked against.",
       "family", "oracle"),
    _r("tensor_product", "CartanEilenberg1956", "family",
       "Kunneth formula for Hochschild (co)homology",
       "The Kunneth isomorphism HH^n(A(x)B) = (+)_{i+j=n} HH^i(A)(x)HH^j(B) that "
       "makes HH multiplicative on tensor factors -- the anchor for TensorProduct(A, B).",
       "family"),
    _r("hodge", "GerstenhaberSchack1987", "algorithm",
       "Hodge (lambda) decomposition",
       "The eigenspace splitting HH^n = (+) HH^{n,(i)} of commutative/tensor and "
       "incidence-algebra pieces.",
       "decomposition"),
    _r("cyclic", "Connes1985", "algorithm",
       "Cyclic homology",
       "Connes' B-operator and the SBI sequence -- quiverlab's cyclic homology.",
       "cyclic"),
    _r("cup", "Gerstenhaber1963", "algorithm",
       "Cup product on Hochschild cohomology",
       "The associative cup product on HH^* (Gerstenhaber-algebra structure).",
       "product"),
    _r("bracket", "Gerstenhaber1963", "algorithm",
       "Gerstenhaber bracket",
       "The graded Lie bracket making HH^* a Gerstenhaber algebra.",
       "bracket"),
    _r("gerstenhaber", "Gerstenhaber1963", "foundation",
       "Cohomology structure of an associative ring",
       "The definitional source of the cup product and Gerstenhaber bracket.",
       "foundation"),
    _r("conway", "Luebeck_ConwayPolynomials", "field",
       "Conway polynomials for finite fields",
       "Lubeck's Conway-polynomial tables fixing canonical generators of GF(p^n).",
       "field"),
    _r("finite_fields", "Luebeck_ConwayPolynomials", "field",
       "Finite field arithmetic",
       "Deterministic cross-compatible GF(q) arithmetic via Conway polynomials.",
       "field"),
    _r("path_algebra", "ASS2006", "family",
       "Bound quiver algebras kQ/I",
       "The path-algebra / bound-quiver formalism for PathAlgebra and the catalog.",
       "family"),
    _r("nakayama", "ASS2006", "family",
       "Nakayama (serial) algebras",
       "Serial algebras by Kupisch series -- the NakayamaAlgebra family.",
       "family"),
    _r("incidence", "ASS2006", "family",
       "Incidence algebras of posets",
       "The incidence algebra kP realized as a bound quiver -- the IncidenceAlgebra family.",
       "family"),
    _r("preprojective", "ASS2006", "family",
       "Preprojective algebras",
       "The preprojective algebra of a Dynkin quiver -- the PreprojectiveAlgebra family.",
       "family"),
    _r("assem_book", "ASS2006", "foundation",
       "Elements of the Representation Theory of Associative Algebras",
       "The standard reference for bound quivers and the representation theory quiverlab implements.",
       "book"),
    _r("han_conjecture", "Han2006", "foundation",
       "Han's conjecture",
       "Finite global dimension iff finite Hochschild homology dimension -- the "
       "conjecture the zoo scans probe.",
       "conjecture"),
    _r("priddy", "Priddy1970", "algorithm",
       "Koszul resolutions",
       "The Priddy PBW / G-quadratic certifier: a quadratic Gröbner basis (all "
       "reduction tips length 2) proves the algebra is Koszul (Plan 27).",
       "koszul"),
    _r("froberg_koszul", "Froberg1999", "algorithm",
       "Koszul algebras",
       "The Hilbert-series Koszulity criterion P(t)*C_A(-t) = I -- quiverlab's "
       "Fröberg numeric Koszulity falsifier (Plan 27).",
       "koszul"),
    _r("polishchuk_positselski", "PolishchukPositselski2005", "foundation",
       "Quadratic Algebras",
       "The quadratic-dual conventions (A^! = kQ^op/R^perp) behind quiverlab's "
       "Koszul dual and the E(A) = (A^!)^op cross-check (Plan 27).",
       "koszul"),
    # --- Plan 29: literature-oracle batteries ---
    _r("happel_trace", "Happel1997", "foundation",
       "The trace of the Coxeter matrix and Hochschild cohomology",
       "Happel's trace identity tr(Coxeter) = sum (-1)^i dim HH^i for finite "
       "global dimension -- the Hochschild/Coxeter cross-invariant consistency "
       "oracle (Plan 29).",
       "hochschild", "coxeter", "oracle"),
    _r("keller_cyclic_invariance", "Keller1998cyclic", "foundation",
       "Invariance and localization for cyclic homology of DG algebras",
       "Derived invariance of cyclic homology (with HH^*/HH_*): reflection-"
       "equivalent orientations of one graph share HH^*/HH_*/HC_* -- the derived-"
       "invariance oracle scheme.",
       "derived", "cyclic", "oracle"),
    _r("rickard_derived", "Rickard1989", "foundation",
       "Morita theory for derived categories",
       "Derived-equivalent algebras (e.g. a Brauer tree and its Brauer star) "
       "share Hochschild and cyclic homology -- the derived-invariance oracle.",
       "derived", "oracle"),
    _r("lenzing_meltzer_ruan", "LenzingMeltzerRuan2022nakayama", "family",
       "Nakayama algebras and Fuchsian singularities",
       "Exact Coxeter polynomials of the uniserial Nakayama algebras N_n(r) -- "
       "the spectral oracle for the Nakayama family (Plan 29).",
       "spectral", "nakayama", "oracle"),
    _r("lenzing_delapena_spectral", "LenzingdlPena2008spectral", "foundation",
       "Spectral analysis of finite dimensional algebras and singularities",
       "The Dynkin / extended-Dynkin / canonical Coxeter-polynomial tables and "
       "Happel's trace identity -- the spectral-invariant oracle.",
       "spectral", "coxeter", "oracle"),
    _r("delapena_mahler", "dlPena2014mahler", "foundation",
       "On the Mahler measure of the Coxeter polynomial of an algebra",
       "The wild star [2,3,7] realizes Lehmer's polynomial -- the "
       "spectral_radius / mahler_measure oracle (Plan 29).",
       "spectral", "oracle"),
    _r("redondo_roman_2014", "RedondoRoman2014", "family",
       "Hochschild cohomology of triangular string algebras and its ring structure",
       "HH^* of the triangular string algebras A_n (with the degree-(2m+1) "
       "revival) and its trivial positive-degree cup product -- the string-"
       "algebra HH oracle.",
       "hochschild", "oracle"),
    _r("taillefer_taft", "taillefer2001taft", "family",
       "Cyclic homology of the Taft algebras and of their Auslander algebras",
       "HH_* and HC_* of the cyclic Nakayama (Taft) algebras in characteristic "
       "zero -- the cyclic-homology oracle.",
       "cyclic", "nakayama", "oracle"),
    _r("cibils_radsq", "cibils1998radsq", "family",
       "Hochschild cohomology algebra of radical square zero algebras",
       "The parallel-path cochain complex for kQ/J^2 (with the characteristic-2 "
       "doubling of k[x]/(x^2)) -- the radical-square-zero HH oracle.",
       "hochschild", "oracle"),
    _r("cibils_incidence", "cibils1989incidence", "family",
       "Cohomology of incidence algebras and simplicial complexes",
       "HH^n of an incidence algebra equals the simplicial cohomology of the "
       "poset order complex -- the incidence-vs-nerve HH oracle.",
       "hochschild", "incidence", "oracle"),
    _r("redondo_incidence", "redondo2008incidence", "family",
       "Hochschild cohomology via incidence algebras",
       "The simplicial-cohomology identification of HH^* underpinning the "
       "incidence-vs-nerve oracle (with Cibils 1989).",
       "hochschild", "incidence", "oracle"),
    _r("cmrs_split", "cibilsmarcosredondosolotar2003", "foundation",
       "Cohomology of split algebras and of trivial extensions",
       "HH^1(T(A)) is never zero -- Z(A) is always a summand -- for every "
       "finite-dimensional A; the trivial-extension HH^1 oracle.",
       "hochschild", "oracle"),
    _r("crs_trivial_ext_hh1", "cibilsredondosaorin2004", "foundation",
       "The first cohomology group of the trivial extension of a monomial algebra",
       "The HH^1(T(A)) decomposition and Example 2.20 (the Z_5 cycle) -- the "
       "trivial-extension first-cohomology oracle.",
       "hochschild", "oracle"),
    _r("xhj_truncated", "xuhanjiang2007truncated", "foundation",
       "Hochschild cohomology of truncated quiver algebras",
       "For a truncated algebra kQ/R^N, dim HH^* is finite iff Q is acyclic -- "
       "the truncated finiteness boolean oracle.",
       "hochschild", "oracle"),
    _r("cibils_acyclic", "cibils1986nocycles", "foundation",
       "Hochschild homology of an algebra whose quiver has no oriented cycles",
       "An acyclic quiver has HH_n = 0 for n >= 1 -- the acyclic Hochschild-"
       "homology vanishing oracle.",
       "hochschild", "oracle"),
    _r("skowronski_yamagata", "SkowronskiYamagata2011", "foundation",
       "Frobenius Algebras I: Basic Representation Theory",
       "Symmetric and Frobenius representation theory, incl. the symmetric "
       "Nakayama criterion n | (L-1) -- the is_symmetric regression oracle.",
       "frobenius", "symmetric"),
    _r("schremmer_wpl", "Schremmer2025wpl", "foundation",
       "Weighted projective lines and Hochschild cohomology",
       "HH^* of the canonical algebras (dim HH^2 = t-3, after Happel LNM 1404) "
       "-- the canonical-algebra Hochschild oracle.",
       "hochschild", "oracle"),
]}


def all_keys() -> tuple:
    return tuple(REGISTRY)


def references_bib_path() -> pathlib.Path:
    return _BIB


def reference(key: str) -> Reference:
    try:
        return REGISTRY[key]
    except KeyError:
        near = difflib.get_close_matches(key, REGISTRY, n=3)
        hint = f"did you mean {near}?" if near else f"known keys: {sorted(REGISTRY)}"
        raise CitationError(f"unknown citation key {key!r}", hint=hint) from None


def bibtex(key: str) -> str:
    ref = reference(key)
    text = _BIB.read_text(encoding="utf-8")
    m = re.search(r"(@\w+\{" + re.escape(ref.bibtex_key) + r",.*?\n\})", text, re.S)
    if m is None:
        raise CitationError(
            f"{key!r} maps to {ref.bibtex_key!r} but that entry is not in references.bib",
            hint="references.bib and the registry are out of sync")
    return m.group(1)
