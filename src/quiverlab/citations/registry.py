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
    _r("weibel_homological", "Weibel1994", "foundation",
       "An Introduction to Homological Algebra",
       "The section-5.4 spectral-sequence page formulas and the strong-convergence "
       "theorem -- quiverlab's spectral-sequence engine (Plan 42).",
       "spectral", "resolution"),
    _r("barakat_homalg", "BarakatLangeHegermann2011", "algorithm",
       "Spectral filtrations via generalized morphisms",
       "The Grothendieck spectral sequence over general module categories (the homalg "
       "framing) -- quiverlab's Grothendieck / Cartan-Eilenberg change-of-rings preset.",
       "spectral"),
    _r("cartan_eilenberg", "CartanEilenberg1956", "foundation",
       "Homological Algebra",
       "The change-of-rings spectral sequence -- quiverlab's Cartan-Eilenberg preset "
       "(Plan 42).",
       "spectral"),
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
    _r("ars_book", "ARS1995", "foundation",
       "Representation Theory of Artin Algebras",
       "The Auslander-Reiten theory reference: almost-split sequences, irreducible "
       "maps, the AR quiver, and the Nakayama functor -- the ground truth for Plan 41.",
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
    _r("happel_trivial_extension", "Happel1988", "foundation",
       "Triangulated Categories in the Representation Theory of Finite Dimensional Algebras",
       "The trivial extension T(A) = A |x D(A) is symmetric for every "
       "finite-dimensional A; the repetitive-algebra framework -- the anchor for "
       "the certified double-quiver TrivialExtension presentation (Plan 31).",
       "frobenius", "symmetric"),
    _r("happel_triangulated", "Happel1988", "foundation",
       "Triangulated Categories in the Representation Theory of Finite Dimensional Algebras",
       "The derived-category reference: the Serre functor / AR triangles of "
       "D^b(mod A) exist iff gl.dim < infinity, tau_{D^b} = nu[-1] -- the ground "
       "truth for the Plan-43 derived surface.",
       "derived", "triangulated"),
    _r("schremmer_wpl", "Schremmer2025wpl", "foundation",
       "Weighted projective lines and Hochschild cohomology",
       "HH^* of the canonical algebras (dim HH^2 = t-3, after Happel LNM 1404) "
       "-- the canonical-algebra Hochschild oracle.",
       "hochschild", "oracle"),
    # --- Plan 40: C6 homological-dimensions family ---
    _r("igusa_todorov", "IgusaTodorov2005", "algorithm",
       "On the finitistic global dimension conjecture for Artin algebras",
       "The Igusa-Todorov functions phi/psi on the finite K0 + syzygy operator "
       "(phi = pd for finite projective dimension) -- quiverlab's Plan-40 IT engine.",
       "homdim", "finitistic"),
    _r("barrios_mata", "BarriosMataRama2020", "family",
       "The Igusa-Todorov phi function for truncated path algebras",
       "Closed forms for the phi/psi-dimension of truncated path algebras kQ/J^k "
       "(via kQ/J^2) -- the Plan-40 Igusa-Todorov literature oracle.",
       "homdim", "oracle"),
    _r("gelinas_delooping", "Gelinas2022", "foundation",
       "The depth, the delooping level and the finitistic dimension",
       "The delooping level dell(A) as an upper bound for the finitistic dimension "
       "-- the deferred Plan-40 Task-F invariant (honest-scope note).",
       "homdim", "finitistic"),
    # --- Plan 46: C5 gentle / string subsystem ---
    _r("butler_ringel", "ButlerRingel1987", "algorithm",
       "Auslander-Reiten sequences for string algebras",
       "Butler-Ringel: the string/band module classification and the hook/cohook "
       "description of the AR translate -- the ground truth for the string subsystem.",
       "modules"),
    _r("avella_geiss", "AvellaAlaminosGeiss2008", "algorithm",
       "Combinatorial derived invariants for gentle algebras",
       "The AG-invariant: a multiset of (n,m) pairs from permitted/forbidden threads; "
       "a DERIVED invariant, provably NOT complete.", "invariants"),
    _r("schroll_brauer", "Schroll2018", "family",
       "Brauer graph algebras (survey)",
       "The presentation of a Brauer graph algebra from a ribbon graph + multiplicities; "
       "the dimension and symmetric structure.", "families"),
    _r("wald_waschbusch", "WaldWaschbusch1985", "foundation",
       "Tame biserial algebras",
       "Biserial / special-biserial structure underlying string and Brauer graph "
       "algebras.", "families"),
    _r("bongartz_tilting", "Bongartz1981", "foundation",
       "Tilted algebras",
       "Bongartz's count criterion for tilting modules (# non-iso indecomposable "
       "summands = # vertices, given pd<=1 and self-Ext vanishing) and the Bongartz "
       "completion of a partial tilting module (Plan 44 / C7).",
       "tilting"),
    _r("derksen_weyman_zelevinsky", "DWZ2008", "family",
       "Quivers with potentials and their representations I",
       "The Jacobian algebra kQ/(cyclic derivatives) of a quiver with potential (Q, W) "
       "-- quiverlab's JacobianAlgebra constructor (Plan 44 / C7).",
       "family", "jacobian"),
    _r("labardini", "LabardiniFragoso2009", "family",
       "Quivers with potentials associated to triangulated surfaces",
       "Surface quivers with potentials whose Jacobian algebras are finite-dimensional "
       "(the framework for the Plan-44 Jacobian constructor; surface QPs land in P48).",
       "family", "jacobian"),
    _r("fomin_shapiro_thurston", "FominShapiroThurston2008", "foundation",
       "Cluster algebras and triangulated surfaces I",
       "The arc/triangulation combinatorics: the ideal-arc count n=6g-6+3(b+p)+Sum k_i, "
       "the admissibility exclusion list, and flip<->mutation -- the ground truth for the "
       "Plan-48 surface subsystem.",
       "families", "surfaces"),
    _r("fomin_zelevinsky_ca1", "FominZelevinsky2002", "foundation",
       "Cluster algebras I: Foundations",
       "The skew-symmetric matrix mutation mu_k that surface flip is certified against "
       "(Plan 48).",
       "families", "surfaces"),
    _r("abcp", "ABCP2010", "family",
       "Gentle algebras arising from surface triangulations",
       "For an unpunctured surface with boundary the Jacobian Jac(Q(T),W(T)) is a GENTLE "
       "algebra -- the Plan-48 v1 certifiability theorem (with Labardini 2009).",
       "families", "surfaces"),
    _r("hughes_waschbusche", "HughesWaschbusch1983", "family",
       "Trivial extensions of tilted algebras",
       "The repetitive algebra hat(A) and its connecting D(A) bimodule -- the source of "
       "quiverlab's finite repetitive-algebra slices (Plan 44 / C7).",
       "family", "repetitive"),
    _r("kac_canonical", "Kac1980", "foundation",
       "Infinite root systems, representations of graphs and invariant theory",
       "Kac's roots, Schur roots, and the canonical decomposition of a dimension "
       "vector -- the ground truth for Plan 49's canonical_decomposition.", "geometry"),
    _r("schofield_general_reps", "Schofield1992", "foundation",
       "General representations of quivers",
       "Schofield's generic hom/ext and the general-representation identities "
       "hom - ext = <a,b> underlying the canonical decomposition.", "geometry"),
    _r("derksen_weyman_canonical", "DerksenWeyman2002", "foundation",
       "On the canonical decomposition of quiver representations",
       "The Derksen-Weyman recursive algorithm for the canonical decomposition "
       "(Plan 49 ships the Dynkin case; Euclidean/wild is the named deferral).", "geometry"),
    _r("zwara_degenerations", "Zwara2000", "foundation",
       "Degenerations of finite-dimensional modules are given by extensions",
       "Zwara: the degeneration order equals the extension order for Artin algebras "
       "-- half of Plan 49's degeneration_order theorem.", "geometry"),
    _r("bongartz_degenerations", "Bongartz1996", "foundation",
       "On degenerations and extensions of finite dimensional modules",
       "Bongartz: degeneration = hom order for representation-finite algebras -- the "
       "computable form Plan 49's degeneration_order uses.", "geometry"),
    _r("voigt_rigidity", "Voigt1977", "foundation",
       "Induzierte Darstellungen ... (Voigt's lemma)",
       "Voigt's lemma: Ext^1(M,M) = 0 => the orbit of M is open (rigid => open orbit); "
       "the codim = dim Ext^1(M,M) equality on hereditary algebras.", "geometry"),
    # --- Plan 47: quasi-hereditary algebras + recollements ---
    _r("dlab_ringel", "DlabRingel1989", "foundation",
       "Quasi-hereditary algebras",
       "The definition of quasi-hereditary algebras, standard modules Delta(i), and the "
       "quasi-heredity test used in Plan 47; qh => finite gl.dim.",
       "quasihereditary"),
    _r("ringel_dual", "Ringel1991", "foundation",
       "Good filtrations and the characteristic tilting module",
       "The characteristic tilting module and the Ringel dual R(A) = End_A(T)^op "
       "(Plan 47).",
       "quasihereditary", "tilting"),
    _r("cps", "CPS1988", "foundation",
       "Finite-dimensional algebras and highest weight categories",
       "Highest weight categories and the idempotent recollement (eAe, A/AeA) of "
       "Plan 47.",
       "quasihereditary", "recollement"),
    _r("bbd", "BBD1982", "foundation",
       "Faisceaux pervers",
       "The origin of recollement and the six-functor formalism (Plan 47).",
       "recollement"),
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
