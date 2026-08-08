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
    _r("oke_koszul", "Oke2021", "algorithm",
       "Bracket structure on Hochschild cohomology of Koszul quiver algebras "
       "using homotopy liftings",
       "The homotopy-lifting bracket formulation quiverlab's native CS bracket "
       "implements (Plan 51); source of the section-7 Koszul-quiver worked tables.",
       "bracket"),
    _r("witherspoon_gsm204", "WitherspoonGSM204", "foundation",
       "Hochschild Cohomology for Algebras",
       "The textbook derivation of the homotopy-lifting Gerstenhaber bracket "
       "(the expository anchor for quiverlab's native CS bracket).",
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
    _r("bv_tradler", "Tradler2008", "algorithm",
       "The BV algebra on Hochschild cohomology from infinity inner products",
       "The symmetric-algebra BV operator Delta on HH^* induced by a symmetric, "
       "invariant, nondegenerate inner product -- quiverlab's symmetric (nu-inner) "
       "BV route (Plan 54).",
       "bv"),
    _r("bv_lzz", "LambreZhouZimmermann2016", "algorithm",
       "BV structure on HH^* for Frobenius algebras with semisimple Nakayama automorphism",
       "The Frobenius-algebra BV operator under the semisimple-nu hypothesis -- "
       "quiverlab's semisimple-nu / twisted BV route (Plan 54).",
       "bv"),
    _r("bv_volkov", "Volkov2016bv", "algorithm",
       "BV-differential on Hochschild cohomology of Frobenius algebras",
       "The BV differential under ord(nu) coprime to char k (over GF(p): the "
       "squarefree-minpoly / p does not divide ord(nu) phrasing recorded in "
       "quiverlab's provenance) -- Plan 54.",
       "bv"),
    _r("bv_biklz", "BianItagakiKouLyuZhou2026", "algorithm",
       "BV algebra on HH^* of self-injective Nakayama algebras",
       "The self-injective Nakayama BV structure; Sec 3.2 (e=1) gives the explicit "
       "k[x]/(x^N) Delta values that back quiverlab's char-sensitivity + Delta-rank "
       "literature oracle (Plan 54).",
       "bv"),
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
    # Plan 57 (R21 + R37): Liu degree theory, the Chaio-Liu radical of the module
    # category, CMMS radical-square-zero, and Ringel's directing/trichotomy reference.
    _r("liu_degrees", "Liu1992degrees", "foundation",
       "Degrees of irreducible maps and the shapes of Auslander-Reiten quivers",
       "Liu's left/right degrees of irreducible morphisms and their control of the "
       "AR-quiver shape -- the R21 degree theory Plan 57 computes on the knit.", "ar"),
    _r("liu_semistable", "Liu1993semistable", "foundation",
       "Semi-stable components of an Auslander-Reiten quiver",
       "Liu's component classification (sectional paths, semistable/directed "
       "components) underlying Plan 57's partition and rep-directed recognizer.", "ar"),
    _r("chaio_liu_radical", "ChaioLiu2013", "foundation",
       "A note on the radical of a module category",
       "Chaio-Liu: rep-finiteness through the infinite radical; the nilpotency of "
       "rad(mod A) in the rep-finite case is the maximal depth of composites -- the "
       "theorem behind Plan 57's nilpotency index and rad^inf=0 gate.", "ar"),
    _r("cmms_radsq", "CMMS1994radsq", "foundation",
       "Module categories with infinite radical square zero are of finite type",
       "CMMS: (rad^inf)^2 = 0 implies representation-finite -- the class oracle "
       "justifying Plan 57's finite-nilpotency-index rep-finiteness certificate "
       "(documented, not a per-instance decider off the rep-finite domain).", "ar"),
    _r("ringel_tame", "Ringel1984tame", "foundation",
       "Tame Algebras and Integral Quadratic Forms",
       "Ringel LNM 1099: directing modules, the postprojective/regular/preinjective "
       "trichotomy -- the R21 component-invariant reference.", "book"),
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
    # --- Plan 58: certified Coxeter spectral analysis (record R20) ---
    _r("dlPena2014mahler", "dlPena2014mahler", "foundation",
       "On the Mahler measure of the Coxeter polynomial of an algebra",
       "de la Pena's class-conditional Mahler-measure dichotomy for accessible "
       "algebras (M = 1 or M >= mu_0, mu_0 = Lehmer's number) -- the Plan-58 "
       "Lehmer-class documentation note (shares the dlPena2014mahler eprint with "
       "the Plan-29 delapena_mahler oracle key).",
       "spectral", "coxeter", "oracle"),
    _r("dlPena2013cyclotomic", "dlPena2013cyclotomic", "foundation",
       "Algebras whose Coxeter polynomials are products of cyclotomic polynomials",
       "de la Pena's separation of periodic (finite-order) Coxeter transformations "
       "from merely cyclotomic-type (quasi-unipotent) ones -- the Plan-58 "
       "quasi-unipotent / finite-order verdict.",
       "spectral", "coxeter"),
    _r("dlPenaTakane1990spectral", "dlPenaTakane1990spectral", "foundation",
       "Spectral properties of Coxeter transformations and applications",
       "de la Pena-Takane: reality of the dominant Coxeter eigenvalue on the "
       "wild-hereditary locus -- the Plan-58 real-dominant spectral-radius "
       "foundation (Arch. Math. 55 (1990) 120-134).",
       "spectral", "coxeter"),
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
    _r("chaparro_schroll_solotar", "ChaparroSchrollSolotar2020", "foundation",
       "On the Lie algebra structure of the first Hochschild cohomology of gentle "
       "and Brauer graph algebras",
       "Determines HH^1(A, M) with different coefficients for gentle algebras via "
       "ribbon-graph combinatorics -- the Plan-52 gentle HH^1-with-coefficients "
       "provenance (numeric pin BLOCKED-until-transcribed).",
       "hochschild", "coefficients"),
    _r("lindell_rubio_relative", "LindellRubio2024", "foundation",
       "On the first relative Hochschild cohomology and the contracted fundamental group",
       "The Lie-algebra structure of the first RELATIVE (vertex-relative, E = kQ_0) "
       "Hochschild cohomology, with radical-square-zero computations -- the Plan-52 "
       "relative HH provenance (numeric pin BLOCKED-until-transcribed).",
       "hochschild", "relative"),
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
    # --- Plan 53: phidim/psidim as algebra invariants + LIT + fractional CY ---
    _r("fernandes_lanzilotta_mendoza", "FernandesLanzilottaMendoza2015", "foundation",
       "The Phi-dimension: a new homological measure",
       "phidim(A) = sup phi(M) as an algebra invariant + derived-equivalence invariance "
       "of its finiteness -- the Plan-53 phidim/psidim ground truth.",
       "homdim", "finitistic"),
    _r("bravo_lanzilotta_mendoza_vivero", "BravoLanzilottaMendozaVivero2021", "foundation",
       "Generalised Igusa-Todorov functions and Lat-Igusa-Todorov algebras",
       "The LIT algebras + the proof-carrying finitistic bound psi_D(V) + n + 1 -- the "
       "Plan-53 LIT finitistic certificate.",
       "homdim", "finitistic"),
    _r("ivanov_volkov", "IvanovVolkov2012", "family",
       "Stable Calabi-Yau dimension of self-injective algebras of finite type",
       "The Serre functor S = Omega.nu, suspension Sigma = Omega^{-1}, and the criterion "
       "Omega^{n+1} ~ nu^{-1} -- the Plan-53 fractional-CY ground truth (Table 1 deferred).",
       "modules", "oracle"),
    _r("erdmann_skowronski_scy", "ErdmannSkowronski2006", "foundation",
       "The stable Calabi-Yau dimension of tame symmetric algebras",
       "Introduced the stable Calabi-Yau dimension of a self-injective algebra as the "
       "weak CY dimension of mod-bar A -- the Plan-53 fractional-CY foundation.",
       "modules"),
    _r("geiss_leclerc_schroer", "GeissLeclercSchroer2006", "family",
       "Rigid modules over preprojective algebras",
       "The stable category of a Dynkin preprojective algebra is 2-Calabi-Yau -- the "
       "Plan-53 Pi(Delta) fractional-CY (2,1) literature oracle.",
       "modules", "oracle"),
    # --- Plan 46: C5 gentle / string subsystem ---
    _r("butler_ringel", "ButlerRingel1987", "algorithm",
       "Auslander-Reiten sequences for string algebras",
       "Butler-Ringel: the string/band module classification and the hook/cohook "
       "description of the AR translate -- the ground truth for the string subsystem.",
       "modules"),
    # --- Plan 59: R34 homological string test + R35 toupie algebras ---
    _r("suarez_alvarez", "SuarezAlvarez2023", "algorithm",
       "A simple homological characterization of string algebras of finite rep. type",
       "Suarez-Alvarez: among rep-finite algebras, string <=> the middle term of EVERY "
       "extension of indecomposables has <= 2 summands (all Ext^1 classes, not just AR "
       "sequences) -- the homological string test.", "modules"),
    _r("alsolotar_toupie", "ArtensteinLanzilottaSolotar2020", "family",
       "Hochschild cohomology of toupie algebras",
       "Artenstein-Lanzilotta-Solotar: toupie = unique source/sink + a parallel "
       "branches; a-Kronecker HH^* = [1, a^2-1, 0, ..]; HH^1 contains sl_a (char 0), "
       "a = # direct source->sink arrows.", "families", "hochschild"),
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
    # --- Plan 68: R32 skew-gentle algebras ---
    _r("he_zhou_zhu", "HeZhouZhu2020", "family",
       "A geometric model for the module category of a skew-gentle algebra",
       "He-Zhou-Zhu: the skew-gentle triple (Q,Sp,I), the idempotent specialization "
       "eps^2=eps, the split/doubled quiver, and support tau-tilting via the orbifold "
       "model -- the primary skew-gentle source.", "families", "modules"),
    _r("chen_skew_gentle", "Chen2022skewgentle", "family",
       "A characteristic free approach to skew-gentle algebras",
       "Chen: the char-free idempotent split construction (sec 3), the selfinjective "
       "classification (Sp=empty & gentle selfinjective), the K-theory/dim relation, and "
       "Gorensteinness -- the split constructor's binding source, valid in char 2.",
       "families"),
    _r("amiot_skew_gentle", "Amiot2021skewgentle", "algorithm",
       "Indecomposable objects in the derived category of a skew-gentle algebra via orbifolds",
       "Amiot: skew-gentle as Z2-skew-group of a gentle algebra; the orbifold/double-cover "
       "geometric model backing the geometric-vs-engine tau-tilting cross-check.", "modules"),
    _r("garcia_lavoue", "GarciaLavoue2026", "algorithm",
       "Brick-finite skew-gentle algebras are representation-finite",
       "Garcia-Lavoue Thm 3.1 (char != 2): brick-finite <=> rep-finite for skew-gentle -- "
       "composed with DIJ (brick-finite <=> tau-tilting-finite) gives the rep-type "
       "certificate.", "modules"),
    _r("geiss_delapena", "GeissDeLaPena1999", "foundation",
       "Auslander-Reiten components for clans",
       "Geiss-de la Pena: the original skew-gentle / clan definition (char != 2).",
       "families"),
    _r("crawley_boevey_clans", "CrawleyBoevey1989", "foundation",
       "Functorial filtrations II: clans and the Gelfand problem",
       "Crawley-Boevey: the classification of indecomposables for clans, underlying the "
       "special-string re-gluing (the +/- forms).", "modules"),
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
    _r("air_tau_tilting", "AIR2014", "foundation",
       "tau-tilting theory",
       "Adachi-Iyama-Reiten: support tau-tilting pairs, mutation, the exchange graph, "
       "g-vectors, and the bijection with functorially finite torsion classes -- the "
       "ground truth for the Plan-45 / C4 tau-tilting engine.",
       "tau-tilting"),
    _r("demonet_iyama_jasso", "DIJ2019", "foundation",
       "tau-tilting finite algebras, bricks, and g-vectors",
       "DIJ: tau-tilting-finiteness <=> finite g-fan <=> finitely many bricks; the "
       "counting identities and the wall-and-chamber / g-vector fan (Plan 45 / C4).",
       "tau-tilting"),
    _r("king_stability", "King1994", "foundation",
       "Moduli of representations of finite-dimensional algebras",
       "King's theta-(semi)stability and GIT walls -- the wall-and-chamber structure the "
       "Plan-45 / C4 engine draws.",
       "tau-tilting", "stability"),
    # --- Plan 63: wall-and-chamber structure via bricks (R25) ---
    _r("brustle_smith_treffinger", "BST2019", "foundation",
       "Wall and Chamber Structure for finite-dimensional Algebras",
       "Brustle-Smith-Treffinger: the wall D(M) = {theta : M theta-semistable} (Def 3.1-3.3), "
       "chambers <-> support tau-tilting pairs (Thm 1.2 / Cor 3.29), one wall = many facets "
       "(Rem 3.19) -- the ground truth for Plan 63. NB walls=D(brick) is DIJ, not BST.",
       "tau-tilting", "stability"),
    _r("asai_semibricks", "Asai2020", "foundation",
       "Semibricks",
       "Asai: semibricks <-> functorially finite torsion classes <-> support tau-tilting "
       "modules (Thm 1.3 / Prop 1.6) -- the brick/semibrick INDEXING that labels walls and "
       "chambers (Plan 63). Contains no g-vector/fan/wall/chamber content -- do not cite for "
       "the geometry.",
       "tau-tilting", "stability"),
    _r("kaipel_treffinger", "KT2023", "foundation",
       "Wall-and-chamber structures for finite-dimensional algebras and tau-tilting theory",
       "Kaipel-Treffinger: the definition + torsion-class/tau-tilting relationship, with the "
       "worked kA2 (Ex 13: D(P1) a ray) and cyclic rad^2-Nakayama N3^2 (Ex 15: 14 chambers) "
       "examples -- Plan 63 literature oracles.",
       "tau-tilting", "stability"),
    # --- Plan 56: pi1(Q,I) + strongly simply connected (coverings) ---
    _r("assem_delapena", "AssemDelaPena1996", "foundation",
       "The fundamental groups of a triangular algebra",
       "The fundamental group pi1(Q,I) of a presentation and the Hom(pi1, k+) embedding "
       "into HH^1 for triangular algebras -- the Plan-56 pi1 + Hurewicz foundation.",
       "coverings", "pi1"),
    _r("martinez_villa_delapena", "MartinezVillaDelaPena1983", "foundation",
       "The universal cover of a quiver with relations",
       "The original homotopy relation ~_I on walks (minimal relations glue parallel "
       "paths) that quiverlab's pi1^ab computes by exact linear algebra on I/(rad.I + I.rad).",
       "coverings", "pi1"),
    _r("le_meur_pi1", "LeMeur2005", "foundation",
       "The fundamental group of a triangular algebra without double bypasses",
       "Thm 1.1: a char-0 triangular algebra with no double bypasses has a privileged "
       "presentation whose pi1 surjects onto every other -- the Plan-56 no-bypass True route.",
       "coverings", "pi1"),
    _r("crs_hurewicz", "CibilsRedondoSolotar2010hurewicz", "foundation",
       "Fundamental group of Schurian categories and the Hurewicz isomorphism",
       "For a Schurian category the Hurewicz map Hom(pi1, k+) -> HH^1 is an isomorphism "
       "-- the equality case of the Plan-56 Hom(pi1,k+) <= dim HH^1 cross-check.",
       "coverings", "pi1"),
    _r("crs_gradings", "CibilsRedondoSolotar2010gradings", "foundation",
       "Connected gradings and the fundamental group",
       "The intrinsic pi1 (inverse limit over connected gradings) and the oracle "
       "pi1(k[x]/(x^p)) = Z x C_p in char p -- documents WHY the presentation group (Z) "
       "is not the algebra invariant (Plan-56 intrinsic refusal).",
       "coverings", "pi1"),
    _r("crs_intrinsic", "CibilsRedondoSolotar2009intrinsic", "foundation",
       "The intrinsic fundamental group of a linear category",
       "The intrinsic fundamental group of a linear category -- the not-bounded-computable "
       "invariant quiverlab refuses loudly in favour of the presentation pi1 (Plan 56).",
       "coverings", "pi1"),
    _r("briggs_ryd_tori", "BriggsRubioyDegrassi2023", "foundation",
       "Maximal tori in HH^1 and the fundamental group",
       "Every maximal torus of HH^1(A) is dual to some fundamental group of A -- the "
       "maximal-torus refinement of the Plan-56 Hurewicz cross-check.",
       "coverings", "pi1"),
    _r("skowronski_ssc", "Skowronski1993", "foundation",
       "Simply connected algebras and Hochschild cohomologies",
       "A triangular algebra is strongly simply connected iff every full convex "
       "subcategory satisfies the separation condition -- the Plan-56 R16 recognizer "
       "(the P62 tame/wild gate).",
       "coverings", "ssc"),
    _r("act_left_right", "ACT2004", "foundation",
       "The left and the right parts of a module category",
       "Assem-Coelho-Trepode: L_A / R_A via predecessor/successor closure of pd<=1 / id<=1, the "
       "Ext-injective criterion tau^{-1}X notin L_A, and the support algebra A_lambda = End of "
       "the projectives in L_A -- the ground truth for Plan 55.", "recognizer"),
    _r("aclv_supports", "ACLV2011", "foundation",
       "Algebras determined by their supports",
       "Assem-Castonguay-Lanzilotta-Vargas: A_lambda / A_rho are products of tilted algebras for "
       "ada algebras (quasi-tilted in general), D L_A = R_{A^op}, and the laura complement "
       "ind A minus (L_A u R_A) -- non-empty even for ada; feeds P60 (tilted) and P61 (laura/ada).",
       "recognizer"),
    _r("organising_module_category", "AACV2021", "foundation",
       "Organising the module category",
       "Alvares-Assem-Castonguay-Vargas survey of the left/right parts, supports, and the "
       "quasi-tilted/laura/ada organisation of mod A -- Plan 55's secondary reference.", "survey"),
    # Plan 60 (R17): the tilted-algebra recognizer via the Liu-Skowronski faithful section.
    _r("liu_tilted_1993", "Liu1993", "foundation",
       "Tilted algebras and generalized standard Auslander-Reiten components",
       "Liu (independently with Skowronski): A is tilted iff Gamma_A has a faithful "
       "generalized-standard component with a section -- the criterion Plan 60 searches for. "
       "Venue verified (Arch. Math. 61 (1993) 12-19).",
       "recognizer"),
    _r("liu_another_2014", "Liu2014", "foundation",
       "Another characterization of tilted algebras (arXiv:1409.2054)",
       "Liu: A is tilted iff Gamma_A contains a FAITHFUL CUT Delta with Hom(X, tau Y)=0 "
       "(Thm 2.6); the cut is a finite/local object (weakly convex), the documented "
       "rep-infinite extension path. Ringel's slice theorem (Thm 1.9(2)) is the "
       "reconstruction certificate Plan 60 uses.",
       "recognizer"),
    _r("happel_ringel_tilted", "HappelRingel1982", "foundation",
       "Tilted algebras",
       "Happel-Ringel: the origin of tilted algebras A = End_H(T) (H hereditary, T tilting); "
       "tilted => gl.dim <= 2 -- the theorem gate Plan 60 uses to refute high-gl.dim algebras.",
       "recognizer"),
    _r("hrs_quasitilted", "HRS1996", "foundation",
       "Tilting in Abelian Categories and Quasitilted Algebras",
       "Happel-Reiten-Smalo: quasi-tilted = (QT1) gl.dim <= 2 AND (QT2) every indec pd<=1 or "
       "id<=1; QT2 alone => gl.dim <= 3. The definitional ground truth for Plan 61's "
       "quasi-tilted rung.", "recognizer"),
    _r("coelho_lanzilotta_weakly_shod", "CL2003weaklyshod", "foundation",
       "Weakly shod algebras",
       "Coelho-Lanzilotta: shod = every indec pd<=1 or id<=1 (=> gl.dim <= 3); weakly shod = "
       "bounded irreducible-morphism paths from an injective to a projective (WSA). Plan 61's "
       "shod + weakly-shod rungs.", "recognizer"),
    _r("assem_coelho_laura", "AC2003laura", "foundation",
       "Two-sided gluings of tilted algebras",
       "Assem-Coelho: laura = ind A minus (L_A u R_A) is finite. Plan 61's laura rung (trivially "
       "true in representation-finite scope; the finite complement is the reported datum).",
       "recognizer"),
    _r("smith_almost_laura", "Smith2007almostlaura", "foundation",
       "Almost laura algebras",
       "Smith: the almost-laura generalisation of laura algebras -- context for the laura "
       "landscape; Plan 61 cites it for the class, not for the elementary rep-finite triviality.",
       "recognizer"),
    _r("bft_quasitilted_quiver", "BFT2017quasitilted", "foundation",
       "On the quiver with relations of a quasitilted algebra and applications",
       "Bordino-Fernandez-Trepode: the quiver-with-relations structure of quasitilted algebras "
       "-- a secondary reference on Plan 61's quasi-tilted rung.", "recognizer"),
    _r("bongartz_criterion", "Bongartz1984", "foundation",
       "A criterion for finite representation type",
       "Bongartz: a simply connected algebra is representation-finite iff its Tits form is "
       "weakly positive (iff it has no critical convex subcategory) -- the Plan-62 rep-finite "
       "axis of the tame/wild trichotomy.",
       "tits", "tame_wild"),
    _r("bdps_tame_tits", "BrustleDlPSkowronski2011", "foundation",
       "Tame algebras and Tits quadratic forms",
       "Bruestle-de la Pena-Skowronski: a strongly simply connected algebra over an "
       "algebraically closed field is tame iff its Tits form is weakly nonnegative -- the "
       "Plan-62 tame axis (the theorem gated on the P56 strong-simple-connectivity certificate).",
       "tits", "tame_wild"),
    _r("kasjan_skowronski", "KasjanSkowronski2019", "foundation",
       "On the tame-wild dichotomy for strongly simply connected algebras",
       "Kasjan-Skowronski (arXiv:1905.06028): the tame/wild dichotomy statements Plan 62 "
       "consumes for strongly simply connected algebras (companion source to BdlPS).",
       "tits", "tame_wild"),
    _r("ovsienko_forms", "Ovsienko1978", "foundation",
       "Integral weakly positive forms",
       "Ovsienko: every positive root of a weakly positive unit form has coordinates <= 6 -- "
       "the cited complete decision behind Plan-62 weak positivity (the box-6 branch-and-bound).",
       "tits"),
    _r("vonhohne_wnn", "vonHohne1996", "foundation",
       "On weakly non-negative unit forms and tame algebras",
       "von Hohne: the classified hypercritical unit forms driving the Plan-62 weak-nonnegativity "
       "decision (the primary route; NOT bounded to <= 9 variables -- T_{2,3,7} is a 10-variable "
       "hypercritical form).",
       "tits", "tame_wild"),
    _r("delapena_banach26", "DelaPenaBanach26", "foundation",
       "Algebras with hypercritical Tits form",
       "de la Pena (Banach Center Publ. 26): the printed hypercritical Tits-form list -- the "
       "Plan-62 weak-nonnegativity cross-oracle for the classified data.",
       "tits", "tame_wild"),
    _r("bjp_quadratic_forms", "BarotJimenezDlP2019", "foundation",
       "Quadratic Forms: Combinatorics and Numerical Results",
       "Barot-Jimenez-Gonzalez-de la Pena: the combinatorics of critical / hypercritical unit "
       "forms and the search-box exposition underpinning Plan 62.",
       "tits", "tame_wild"),
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
