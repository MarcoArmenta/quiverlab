"""Byte-stability acceptance for the Plan-28 runner refactor.

``webapp/server/runner.py`` no longer holds the compute dispatch -- it DELEGATES
to the wheel's ``quiverlab.hpc.spec``. This test pins that the delegation is
BYTE-IDENTICAL to the pre-refactor runner across a fixture corpus covering
family / quiver / module / ext computes:

  (a) the returned result dict is byte-identical (``json.dumps(sort_keys=True)``);
  (b) the Plan-25 ``canonical_key`` is unchanged.

The goldens in ``_runner_goldens.json`` were frozen from the CURRENT runner
BEFORE the refactor (a regression fence: if the delegation ever perturbs a result
byte or a cache key, this fails).

Re-freezes are DELIBERATE feature changes, listed here so an accidental drift can
never hide behind one:

  * 2026-07-29 (``module_left_a2``): the ``projective_dimension`` /
    ``injective_dimension`` blocks gained the ``latex`` key. Without it the
    draw-page renderer typeset a literal "undefined" for both (Marco's
    example-a). The canonical keys are request-derived and are unchanged.
  * 2026-07-31 (all six, v0.1.0 release): the embedded ``quiverlab_version``
    stamp moved ``0.1.0.dev0`` -> ``0.1.0`` with the library bump. Gated
    re-freeze: each regenerated blob was asserted byte-identical to the old
    one under the version substitution alone, so no math moved. ``_V`` below
    and the frozen ``canonical_key`` entries deliberately KEEP ``0.1.0.dev0``
    -- the key pin is version-explicit and anchors the canonicalizer
    independently of the running library.
  * 2026-08-01 (``products_loop_gf2``): ADDED for Plan 35 -- the new
    cup/cap/bracket/connes_b compute kinds. One entry (``cup:0..2`` +
    ``connes_b:0..2`` over the ``x^3`` loop / GF(2)); the existing six entries
    are untouched.
  * 2026-08-01 (``products_loop_gf2`` re-freeze, explicit representatives): the
    product blocks gained the additive ``basis_classes`` / ``chain_basis`` /
    ``differentials`` fields (Marco: the class symbols must PRINT their explicit
    (co)cycle + coordinate vector + annihilating differential). Gated re-freeze:
    the regenerated blob was asserted byte-identical to the old one after DELETING
    those three keys everywhere, so no pre-existing content moved. The
    ``canonical_key`` is request-derived and UNCHANGED (no new request fields).
  * 2026-08-01 (``module_ext`` re-freeze, Plan 35 UNIT 2 -- resolution term bases):
    the ``projective_resolution`` block gained the additive ``term_basis`` field (the
    ordered concatenated path bases of each term's summands, so a reader can map each
    differential column to the basis vector it acts on). Gated re-freeze, scoped BY
    KEY: the regenerated blob was asserted byte-identical to the old one after
    deleting the ``term_basis`` key everywhere -- and ``term_basis`` is a NEW key on
    resolution blocks ONLY (the pre-existing overloaded ``differentials`` key was not
    touched), so no pre-existing content moved. The ``canonical_key`` is
    request-derived and UNCHANGED. (``products_loop_gf2`` and the other five entries
    are untouched -- only ``module_ext`` carries an actual resolution block.)
  * 2026-08-01 (``module_ext`` re-freeze, Plan 35 wave 3a -- explicit Ext
    representatives): the ``ext`` block gained the additive ``basis_classes`` /
    ``chain_basis`` / ``differentials`` fields (the per-degree cocycle representatives
    over the ordered Hom basis + the annihilating coboundary). Gated re-freeze, scoped
    BY KEY: the regenerated blob was asserted byte-identical to the old one after
    deleting those three keys from the EXT block ONLY (the ``projective_resolution``
    block's overloaded LIST-shaped ``differentials`` was NOT touched -- ext/tor ship the
    dict shape ``{str(degree): ...}``), so no pre-existing content moved. The
    ``canonical_key`` is request-derived and UNCHANGED. (The other six entries are
    untouched -- only ``module_ext`` computes ``ext``.)
  * 2026-08-01 (``module_ext`` re-freeze, Plan 35 wave 3c -- Yoneda exact sequences):
    the ``ext`` block gained the additive ``interpretation`` key (per class, the
    constructed + self-certified n-fold exact sequence 0 -> N -> Q -> ... -> M -> 0
    realizing it). Gated re-freeze, scoped BY KEY: the regenerated blob was asserted
    byte-identical to the old one after deleting the ``interpretation`` key from the EXT
    block ONLY, so no pre-existing content moved. ``canonical_key`` is request-derived
    and UNCHANGED. (Only ``module_ext`` computes ``ext``.)
  * 2026-08-01 (``family_hh_cartan`` + ``quiver_hh_homology`` re-freeze, Plan 35
    wave 3d -- explicit HH representatives): the plain ``hh_cohomology`` /
    ``hh_homology`` dims blocks gained the additive ``basis_classes`` / ``chain_basis``
    / ``differentials`` / ``inner_dims`` fields (the per-degree (co)cycle
    representatives + coordinate vectors + annihilating differential, so HH^0's centre,
    HH^1's derivations, HH^2's deformation cochain and HH_0's commutator residues can be
    read off). Gated re-freeze, scoped BY KEY: each regenerated blob was asserted
    byte-identical to the old one after deleting those four keys from the
    ``hh_cohomology`` / ``hh_homology`` block ONLY, so no pre-existing content moved.
    ``canonical_key`` is request-derived and UNCHANGED. (The other five entries carry no
    plain HH dims block and are untouched.)
  * 2026-08-03 (``family_hh_cartan`` + ``products_loop_gf2`` + ``module_ext``
    re-freeze, Marco report pass 2 -- tensor separators + Ext/Tor resolution):
    (a) bar-route cochain labels now separate tensor factors with `` (x) `` like
    the chain side, never ``|`` (``[x (x) y -> v]``; Marco: keep only ⊗ -- a
    ``|`` read as something other than a k-tensor); (b) the ``ext``/``tor``
    blocks gained the additive ``resolution`` key (the ⊕-decomposition + Betti
    numbers of the resolution of M the engine walked, shown BEFORE the data).
    Gated re-freeze: each regenerated blob was asserted byte-identical to the
    old one after mapping `` (x) `` back to ``|`` inside bracketed labels and
    deleting ``resolution`` from ext/tor blocks -- so nothing else moved.
    ``canonical_key`` is request-derived and UNCHANGED.
  * 2026-08-03 (``module_ext`` re-freeze, Marco report pass -- resolved provenance):
    the ``ext`` / ``tor`` blocks gained the additive ``resolved`` key ({module, side,
    resolution} -- WHICH module the engine resolved and by which resolution, so the
    report/GUI can state the objects before the numbers). Gated re-freeze, scoped BY
    KEY: the regenerated blob was asserted byte-identical to the old one after
    deleting ``resolved`` from the ext/tor blocks ONLY, so no pre-existing content
    moved. ``canonical_key`` is request-derived and UNCHANGED. (Only ``module_ext``
    computes ext; no golden computes tor.)
  * 2026-08-05 (``module_basic`` re-freeze, Plan 37 C1 -- Loewy series): the
    ``rad_top_soc`` block gained the additive ``series`` key (the Loewy / radical
    layers top-to-bottom, ``Module.loewy_layers()`` -- a list of str-keyed
    composition-factor multiplicity dicts, so the report and both GUIs render the
    stacked Loewy diagram). Gated re-freeze, scoped BY KEY: the regenerated blob was
    asserted byte-identical to the old one after DELETING the ``series`` key from the
    ``rad_top_soc`` block, so no pre-existing content moved. ``series`` is a NEW key
    on ``rad_top_soc`` blocks ONLY. The ``canonical_key`` is request-derived and
    UNCHANGED (no new request fields). Only ``module_basic`` carries a ``rad_top_soc``
    block; the other six entries are untouched.
  * 2026-08-05 (``ext_algebra_exterior_gf7`` + ``recognizers_gentle_a3`` ADDED,
    Plan 38): the two new algebra-block compute kinds ``ext_algebra`` (Yoneda
    Ext-algebra + three-valued Koszulity, on the Koszul poster child
    k<x,y>/(x^2, y^2, x*y+y*x) over GF(7)) and ``recognizers`` (the eight
    structural recognizers + Dynkin/Euclidean type + form type, on gentle A3 over
    GF(5)). Both are schema v1 (they act on the algebra block, no new request
    block), so the existing seven entries are UNTOUCHED (verified byte-identical
    on re-dump before appending). ``canonical_key`` is request-derived.
  * 2026-08-05 (``almost_split_a3_s2`` ADDED, Plan 41 -- AR completion): a NEW
    golden for the ``almost_split`` module compute kind (the almost-split sequence
    0 -> tau M -> E -> M -> 0 of S_2 over kA3 / GF(5)). Additive: a brand-new entry
    keyed ``almost_split_a3_s2``; the existing seven entries are untouched and were
    asserted byte-identical before it was added. The block carries ``tau`` (a full
    representation), ``middle.summands`` (E's Krull-Schmidt summands via the shared
    serializer), and ``references``/``citations`` (assem_book + the new ars_book).
  * 2026-07-31 (ALL entries re-freeze, Marco ADDENDUM 2 -- json_guide): the result
    envelope gained ONE additive top-level key, ``json_guide`` -- a per-computation list
    of ``{object, path, note}`` recipes for recovering every computed object from
    result.json (``quiverlab.trace.json_guide.build_json_guide``, self-validating). Gated
    re-freeze: each regenerated blob was asserted byte-identical to the old one after
    DELETING the ``json_guide`` top-level key, so no pre-existing content moved. The
    ``canonical_key`` is request-derived and UNCHANGED (no new request fields).
  * 2026-08-05 (``derived_fingerprint_a4`` ADDED, Plan 43): a NEW fixture for the
    ``derived_fingerprint`` scalar compute kind (schema v1, kA4 = path 1->2->3->4 over
    GF(5) -- the necessary-condition invariant tuple: Coxeter polynomial + Cartan
    det/Smith + HH/HC/centre + gl.dim). Pure addition: the eleven existing entries were
    verified byte-identical BEFORE the new one was appended (the delegation test passed
    on all eleven unchanged). Both runners share ``derived.block.derived_fingerprint_block``,
    so the Pyodide twin agrees (``test_derived_fingerprint_p43.test_twin_parity``); the
    cyclic-homology field is an honest per-field ``{error}`` (the generic (b,B) mixed
    complex has no CS route and blows up over GF(p) for this dim). ``canonical_key`` is
    request-derived.
  * 2026-08-05 (``strings_gentle_a3`` ADDED, Plan 46 C5): a NEW fixture for the
    ``strings`` algebra-only scalar compute kind (gentle A3/(ab) over GF(5) -- the
    recognizer verdicts + string census + band presence + rep-type + AG invariant).
    Pure addition: the eleven existing entries were verified byte-identical BEFORE the
    new one was appended (the delegation test passed on all eleven unchanged). Its
    ``canonical_key`` is request-derived; ``result_json`` was frozen from the server
    runner, and the Plan-46 cross-runner test (tests/gui/test_strings_runner_twin,
    tests/webapp/test_strings_block_p46) asserts the Pyodide twin agrees.
  * 2026-08-05 (``homological_profile_kA2`` ADDED, Plan 40 C6): a NEW fixture for the
    ``homological_profile`` scalar compute kind (kA2 over GF(7) -- global / finitistic
    / dominant / Gorenstein dimensions + Igusa-Todorov phi/psi of the sum of simples).
    Pure addition: the seven existing entries were verified byte-identical BEFORE the
    new one was appended (the delegation test passed on all seven unchanged). Its
    ``canonical_key`` is request-derived; the ``result_json`` was frozen from the
    server runner, and the Plan-40 cross-runner test asserts the Pyodide twin agrees.
  * 2026-08-05 (Plan 42, NEW golden ``ss_hochschild_dualnumbers``): the
    ``ss_hochschild`` compute kind -- the Hochschild ``(b, B)`` spectral-sequence block
    (E_inf page dims + abutment == HC + convergence prose) on ``k[x]/(x^2)`` over
    GF(5), schema v1. A pure ADDITION (a new golden key); every pre-existing entry is
    byte-identical (the parametrized byte-identity + canonical-key tests confirm it),
    and the ``canonical_key`` is request-derived (no new request fields -- an
    algebra-only range kind).
  * 2026-08-05 (``module_tilting_check_kA2`` ADDED, Plan 44 C7): a NEW fixture for the
    ``tilting_check`` module compute kind (kA2 over GF(5), module = the projective P1 --
    NOT tilting: one indecomposable summand for two vertices). Pure addition: every
    pre-existing entry was verified byte-identical BEFORE the new one was appended (the
    delegation test passed on all of them unchanged). Its ``canonical_key`` is
    request-derived; the ``result_json`` was frozen from the server runner, and the
    Plan-44 cross-runner tests assert the Pyodide twin agrees on the math subkeys.
  * 2026-08-05 (``orbit_geometry_kA3_s2`` ADDED, Plan 49 C8): a NEW fixture for the
    ``orbit_geometry`` module compute kind (kA3 over GF(32003), module = the simple
    S_2). Reports orbit dim + Voigt rigidity + honest codim + (hereditary Dynkin) the
    Kac canonical decomposition. Pure addition: every pre-existing entry was verified
    byte-identical BEFORE the new one was appended (the delegation test passed on all
    of them unchanged). Its ``canonical_key`` is request-derived; the ``result_json``
    was frozen from the server runner, and the Plan-49 cross-runner tests
    (``tests/webapp/test_orbit_geometry_p49.py`` / ``tests/gui/``) assert the Pyodide
    twin is byte-identical via the shared ``orbit_geometry_block`` builder.
  * 2026-08-05 (``orbit_geometry_kA3_s2`` re-freeze, P49 devil's-advocate round):
    the ``orbit_geometry`` block gained the additive ``canonical_of`` key (=
    "dimension_vector" -- the canonical decomposition names the GENERIC module
    of d, not M; Marco's name-the-object rule). Gated re-freeze: the regenerated
    blob is byte-identical to the old one after deleting ``canonical_of``.
  * 2026-08-05 (``tau_tilting_kA2`` ADDED, Plan 45 C4): a NEW fixture for the
    ``tau_tilting`` ALGEBRA-level compute kind (kA2 over GF(7), budget 512 -- the full
    run: 5 support τ-tilting pairs, complete, the n=2 fan with 5 chambers, the four-way
    counts all 5). Pure addition: every pre-existing entry was verified byte-identical
    BEFORE the new one was appended (the delegation test passed on all 15 unchanged). Its
    ``canonical_key`` is request-derived (the budget rides in the ``compute`` string, no
    new request field); the ``result_json`` was frozen from the server runner, and the
    Plan-45 cross-runner test asserts the Pyodide twin agrees byte-for-byte.
  * 2026-08-07 (``congruences_kA2`` ADDED, Plan 64 / R26): a NEW fixture for the
    ``congruences`` ALGEBRA-level compute kind (kA2 over QQ, budget 512 -- the full run:
    the pentagon N5 torsion lattice |L|=5, |Con(tors A)|=5, the forcing "V" on 3 bricks,
    #wide=5 = M3). Pure addition: every pre-existing entry was verified byte-identical
    BEFORE the new one was appended (the delegation test passed on all 18 unchanged). Its
    ``canonical_key`` is request-derived (the pair budget rides in the ``compute`` string,
    no new request field), so schema stays v1 and no existing golden re-freezes; the
    ``result_json`` was frozen from the server runner, and the Plan-64 cross-runner tests
    (``tests/webapp/test_congruences_p64.py`` / ``tests/gui/test_congruences_runner_twin.py``)
    assert the Pyodide twin is byte-identical via the shared ``congruences_block`` builder.
  * 2026-08-05 (ALL 18 result_json re-frozen, v0.2.0 bump at the P50 gate): the
    embedded ``quiverlab_version`` moved 0.1.0 -> 0.2.0. Gated re-freeze: every
    regenerated blob is byte-identical to its predecessor after mapping the
    version string back (asserted for all 18 before writing). ``canonical_key``
    values are UNCHANGED -- they are pinned against the frozen ``_V`` constant
    below, deliberately decoupled from the live version.
  * 2026-08-07 (``homological_profile_kA2`` re-freeze, Plan 53 R23): the
    ``homological_profile`` block gained FOUR strictly-additive keys -- ``phidim`` /
    ``psidim`` (the Igusa-Todorov phi/psi-dimensions as ALGEBRA invariants), ``phi_spectrum``
    (the phi-value set + gaps, Barrios-Mata-Rama), and ``lit`` (the Lat-Igusa-Todorov
    finitistic certificate) -- plus three new ``references`` (``fernandes_lanzilotta_mendoza``
    / ``barrios_mata`` / ``bravo_lanzilotta_mendoza_vivero``) with their resolved
    ``citations``, and an extended ``reproduce`` snippet (``A.phi_dim()``,
    ``A.psi_dim()``, ``A.phi_spectrum()``, ``A.finitistic_certificate()``). GATED
    re-freeze: EVERY pre-existing block entry (``global_dimension`` / ``finitistic`` /
    ``dominant`` / ``gorenstein`` / ``igusa_todorov``) was asserted byte-identical before
    writing -- kA2 is HEREDITARY (gl.dim 1), so its ``finitistic`` entry is UNTOUCHED by
    the Task-C LIT wiring (that flip only affects gl.dim=infinity LIT-family inputs). The
    ``canonical_key`` is request-derived and UNCHANGED (the keys are additive, no request
    field moved).
  * 2026-08-07 (``homological_profile_kxx3`` + ``fractional_cy_kxx3`` ADDED, Plan 53):
    two NEW fixtures over the self-injective ``k[x]/(x^3)`` (loop x, ``x*x*x`` = 0, GF(7)).
    ``homological_profile_kxx3`` is the CROSS-RUNNER home of Task C's ``None -> 0``
    ``finitistic`` flip (self-injective => ``finitistic.upper = 0``, ``note`` starts
    ``"LIT"``, ``lit.family == "self-injective"``); ``fractional_cy_kxx3`` is the new
    ``fractional_cy`` compute kind (stable CY dimension ``1/1``, weakly 1-CY, tier
    ``weak-on-generators``). Pure ADDITION: every pre-existing entry was verified
    byte-identical BEFORE the two were appended; both runners share the library builders
    (``homological_profile`` / ``fractional_cy_block``), so the Pyodide twin agrees
    (``tests/webapp/test_phidim_fcy_gui_p53.py``). ``canonical_key`` is request-derived.
  * 2026-08-07 (``fundamental_group_square_gf7`` + ``simply_connected_zito`` ADDED,
    Plan 56): two NEW fixtures for the algebra-scalar compute kinds
    ``fundamental_group`` (the commutative square WITH the commutativity relation over
    GF(7) -- pi1^ab = 0) and ``simply_connected`` (the Zito example over QQ -- simply
    connected True, strongly-simply-connected False at vertex 2). Both are schema v1
    (they act on the algebra block, no new request block), so every pre-existing entry
    is byte-identical (verified by re-dumping the goldens dict byte-for-byte before
    appending). Both runners share the library builders
    (``invariants.coverings_block.fundamental_group_block`` / ``simply_connected_block``),
    so the Pyodide twin agrees (``tests/webapp/test_coverings_exposure_p56.py::
    test_twin_parity``). ``canonical_key`` is request-derived.
  * 2026-08-07 (``coxeter_spectral_3kronecker_qq`` ADDED, Plan 58 / R20): a NEW fixture
    for the ``coxeter_spectral`` algebra-scalar compute kind (the 3-Kronecker over QQ --
    wild, non-cyclotomic: χ = t²-7t+1, ρ = M = (7+3√5)/2 as a CERTIFIED algebraic number
    with minpoly [1,-7,1] and rational isolating interval (6,7), one root outside the
    unit circle, Coxeter order None, the class-conditional Lehmer note). Pure ADDITION:
    the round-trip of the whole goldens file (indent=1) was asserted byte-identical to the
    original BEFORE appending, so every pre-existing entry is untouched. Schema v1 (acts on
    the algebra block, no new request block). Both runners share
    ``invariants.coxeter_spectral.coxeter_spectral_block``, so the Pyodide twin is
    byte-identical (``tests/webapp/test_coxeter_spectral_p58.py::test_twin_parity``).
    ``canonical_key`` is request-derived.
  * 2026-08-07 (``radical_filtration_kA3`` + ``ar_invariants_kA3`` ADDED, Plan 57):
    two NEW algebra-only fixtures over the hereditary ``kA_3`` (1->2->3, GF(32003)).
    ``radical_filtration_kA3`` pins the ``radical_filtration`` block (nilpotency index
    3, ``rad^inf = 0``, layer profile ``[9, 3]``); ``ar_invariants_kA3`` pins the
    ``ar_invariants`` block (representation-directed, Liu-degree table). Pure ADDITION:
    every pre-existing entry was verified content-identical BEFORE the two were
    appended (the generator reproduces the file bytes exactly, then re-dumps with the
    same settings). Both runners share the library builders
    (``modules.radical.radical_filtration_block`` /
    ``modules.ar_invariants.ar_invariants_block``), so the Pyodide twin agrees
    (``tests/webapp/test_radical_filtration_p57.py``, ``tests/gui/test_radical_runner_twin.py``).
    ``canonical_key`` is request-derived (the budget rides in the ``compute`` string).
  * 2026-08-07 (``left_right_parts_kA3`` ADDED, Plan 55 R15): a NEW fixture for the
    ``left_right_parts`` ALGEBRA-level compute kind (schema v1, kA3 = 1->2->3 over QQ,
    budget 256 -- the module-category atlas: L_A = R_A = ind A (6 indecomposables), empty
    complement, both support algebras = A). Pure ADDITION: every pre-existing entry was
    verified byte-identical BEFORE the new one was appended (the reserialization pre-check
    asserted the file round-trips unchanged, and the delegation test passed on all prior
    entries). The budget rides in the ``compute`` string (no new request field), so the
    ``canonical_key`` is request-derived; both runners share
    ``modules.left_right.left_right_parts_block``, so the Pyodide twin agrees
    (``tests/webapp/test_left_right_parts_p55.py`` / ``tests/gui/test_left_right_runner_twin.py``).
  * 2026-08-07 (``hh_cohomology_dual_kA2`` ADDED, Plan 52): a NEW fixture for the
    schema-3 Hochschild-with-COEFFICIENTS block (kA2 over GF(5), the dual bimodule
    D(A), ``hh_cohomology:0..3``). Pure addition: every pre-existing entry was
    verified byte-identical BEFORE the new one was appended (the delegation test
    passed on all 18 unchanged). Its ``canonical_key`` is request-derived (the
    ``coefficients`` block rides the request only when present -- a coefficients-LESS
    request keeps sending schema 2 and drops the key via model_dump, so every prior
    entry's key is unchanged); the ``result_json`` was frozen from the server runner,
    and the Plan-52 cross-runner test asserts the Pyodide twin agrees on the block.
  * 2026-08-07 (``string_homological_kD4`` + ``toupie_a_kronecker`` ADDED, Plan 59):
    two NEW algebra-only scalar kinds over QQ. ``string_homological_kD4`` is the
    homological string test (R34) on kD4 (subspace) -- ``verdict == "not_string"``,
    ``is_string == False``, a 3-summand-middle witness. ``toupie_a_kronecker`` is the
    toupie block (R35) on the 3-Kronecker (three parallel arrows 1->2) -- ``is_toupie``,
    ``branch_count == 3``, ``hh == [1, 8, 0, 0, 0]``, ``sl_a.dim == 8``. Pure ADDITION:
    every pre-existing entry was verified byte-identical BEFORE the two were appended
    (the append is textual, existing bytes untouched); both runners share the library
    builders (``string_homological_block`` / ``toupie_block``), so the Pyodide twin
    agrees (``tests/gui/test_recognizer_runner_twin_p59.py``). ``canonical_key`` is
    request-derived (schema-1 algebra-only, no ``module`` block).
  * 2026-08-07 (``tilted_check_kA3`` ADDED, Plan 60 R17): a NEW fixture for the
    ``tilted_check`` ALGEBRA-level compute kind (schema v1, kA3 = 1->2->3 over QQ,
    ``compute == ["tilted_check:256"]``). The block is ``verdict == "tilted"``,
    ``reason == "hereditary"``, ``hereditary_type == "A_3"`` with the certified
    projective slice + the self-referential reconstruction. Pure ADDITION: every
    pre-existing entry was verified byte-identical BEFORE the append (added/removed/
    changed check == {tilted_check_kA3}/{}/[]); both runners share the library builder
    ``modules.tilted.tilted_check_block``, so the Pyodide twin agrees
    (``tests/webapp/test_tilted_check_p60.py``). ``canonical_key`` is request-derived
    (schema-1 algebra-only, no ``module`` block; the budget carries in the compute
    string, adding no request field).
  * 2026-08-07 (``recognizer_ladder_kA3`` ADDED, Plan 61 R18): a NEW fixture for the
    ``recognizer_ladder`` ALGEBRA-level compute kind (schema v1, kA3 = 1->2->3 over QQ,
    budget 256 -- the quasi-tilted/shod/weakly-shod/laura/ada ladder). Hereditary kA3 =>
    all five verdicts ``True``, empty laura complement, gl.dim 1; the ada/HH^1 block
    reports ``hh1_dim == 0`` with ``applicable == False`` and ``verdict == null`` because
    QQ is NOT algebraically closed (ACLV Theorem B's hypothesis unmet -- no SC verdict off
    CC). Pure ADDITION: every pre-existing entry was verified byte-identical BEFORE the
    append (textual, existing bytes untouched); both runners share the library builder
    (``modules.recognizers_ladder.recognizer_ladder_block``), so the Pyodide twin agrees
    (``tests/webapp/test_recognizer_ladder_p61.py`` /
    ``tests/gui/test_recognizer_ladder_twin.py``). ``canonical_key`` is request-derived
    (schema-1 algebra-only, budget rides in the ``compute`` string, no new request
    field).
  * 2026-08-08 (``recognizer_ladder_kA3`` re-freeze, Plan 61 fix round -- honest gl.dim):
    the ``recognizer_ladder`` block now carries ``gldim_exact`` alongside ``gldim`` (the
    house style of the ``global_dimension`` / ``homological_profile`` blocks), so an
    UNRESOLVED global dimension is presented as a certified lower bound, never a definite
    value (the fix for ``NakayamaAlgebra(kupisch=[3,3,2])``, whose ``global_dimension``
    returns ``value=32, exact=False``). kA3 is exact, so its block gains only
    ``"gldim_exact": true``. Gated re-freeze: the regenerated blob was asserted
    byte-identical to the old one after ADDING the single ``gldim_exact`` key (no other
    field changed), and every OTHER golden entry was verified byte-identical BEFORE the
    write. Both runners share the builder, so the twin stays byte-identical; the
    ``canonical_key`` is unchanged (request-derived, no request-shape change).
  * 2026-08-08 (``tame_wild_a5_cc`` ADDED, Plan 62 / R19): a NEW algebra-only scalar
    kind. A5 (linear 1->2->3->4->5) over CC -- the combinatorial Tits form + the
    rep-finite/tame/wild verdict gated on the P56 certificate: ``rep_type ==
    "rep-finite"``, ``weakly_positive == true``, ``is_unit_form == true``. Pure
    ADDITION: every pre-existing entry was verified byte-identical BEFORE the append
    (60 delegation asserts green first); both runners share the library builder
    (``invariants.tits_block.tame_wild_block``), so the Pyodide twin agrees
    (``tests/webapp/test_tame_wild_exposure_p62.py::test_twin_parity``).
    ``canonical_key`` is request-derived (schema-1 algebra-only, no ``module``
    block).
  * 2026-08-08 (``tame_wild_a5_cc`` RE-FREEZE, Plan 62 fix round -- honesty fixes):
    the tame/wild block gained a ``certified`` field (``"rep_infinite"`` when
    Bongartz certifies representation-infinite but rep_type stays None; ``null``
    here since A5 is rep-finite), and the ``reason`` / ``scope_note`` no longer say
    "over an algebraically closed field" for a char-0 field (the verdict is read by
    BASE CHANGE to the algebraic closure -- field-independent in char 0). Gated
    re-freeze, scoped BY KEY: the regenerated blob differs from the old one ONLY in
    ``results.tame_wild.{certified (added null),reason,scope_note}`` (verified by a
    structured diff); ``body`` and ``canonical_key`` are byte-identical, and every
    OTHER golden entry was asserted unchanged before the write.
  * 2026-08-07 (``bv_operator_kxx3`` ADDED, Plan 54 R2): a NEW fixture for the
    ``bv_operator`` HH-PRODUCT compute kind (the Batalin-Vilkovisky operator Delta)
    over the self-injective ``k[x]/(x^3)`` (loop x, ``x*x*x`` = 0, GF(7)),
    ``bv_operator:0..2``. k[x]/(x^3) is symmetric, so the block is served on the
    Tradler route (``hypothesis == "symmetric (Tradler AIF 2008)"``, ``matrices`` +
    ``ranks`` + ``bracket_check.agrees``). ``bv_operator`` is a member of
    ``PRODUCT_KINDS`` (spec.py), so its ``.blocks()`` IS the block (kind/top/hh_dims/
    matrices/ranks/hypothesis/nakayama/basis/window/references + resolved citations).
    Pure ADDITION: every pre-existing entry was verified byte-identical BEFORE the
    new one was appended (the goldens dict was re-dumped with the same settings, so
    prior bytes are untouched). Both runners share ``Algebra.bv_operator`` (the
    server ``spec._dispatch`` product branch and the Pyodide twin
    ``docs/gui/runner.py``), so the cross-runner block test agrees
    (``tests/webapp/test_bv_block_p54.py``). ``canonical_key`` is request-derived
    (the degree rides in the ``compute`` string; no new request field, so every
    prior entry's key is unchanged).
  * 2026-08-08 (``skew_gentle_arrow_sp2`` ADDED, Plan 68 / R32): a NEW algebra-only
    scalar-kind fixture -- the skew-gentle triple ``Q = 1 --a--> 2, Sp = {2}`` over QQ,
    expressed as the non-scalar constructor ``family: SkewGentleAlgebra`` (flattened
    triple params ``vertices``/``arrows``/``relations``/``special``, the
    ``BrauerGraphAlgebra`` precedent). The block pins the recognizer verdict, the split
    shape (3 vertices, dim 5), the HZZ Lemma 1.5 dim law (split == associated gentle ==
    5), rank ``|Q_0|+|Sp| = 3``, and the rep-type certificate (representation-finite).
    Pure ADDITION: the whole goldens file was asserted byte-identical to
    ``json.dumps(indent=1)+"\\n"`` BEFORE the single key was appended (order-preserving),
    so every pre-existing entry is untouched. Both runners share the library builder
    (``skewgentle.block.skew_gentle_block``) and the client accepts
    ``family: SkewGentleAlgebra`` narrowly, so the Pyodide twin is byte-identical
    (``tests/gui/test_skew_gentle_runner_twin.py``). ``canonical_key`` is request-derived
    (schema-1, no ``module`` block).
  * 2026-08-08 (``wall_chamber_kA2`` ADDED, Plan 63): a NEW algebra-only budget kind
    over QQ -- the wall-and-chamber structure via bricks on kA2 (``wall_chamber:512``),
    5 chambers / 3 walls, complete, render ``fan2d``, ``D(P1)`` a ray. Pure ADDITION:
    every pre-existing entry was verified byte-identical BEFORE the entry was appended
    (the goldens JSON round-trips through ``json.dumps(indent=1)`` byte-for-byte, so the
    new key is appended and no existing bytes move); both runners share the library
    builder (``tautilting.wallchamber.wall_chamber_structure``), so the Pyodide twin
    agrees (``tests/gui/test_wall_chamber_runner_twin.py``). ``canonical_key`` is
    request-derived (schema-1 algebra-only, the budget rides in the ``compute`` string).
  * 2026-08-08 (``silting_local`` ADDED, Plan 67 R30): a NEW fixture for the ``silting``
    ALGEBRA-level compute kind (schema v1, local k[x]/(x^2) = one vertex with a loop x and
    relation x*x over GF(32003), ``compute == ["silting:2,64"]``). The block's regular
    verdict is silting (in fact tilting), and the bounded exploration is
    ``status == "complete"`` / ``finite_class == "local"`` (a single vertex mod shift,
    Thm 2.26). Pure ADDITION: every pre-existing entry was verified byte-identical BEFORE
    the append (the JSON is re-dumped ``indent=1, sort_keys=False`` so stored key order is
    preserved). Both runners share the library builder ``derived.block.silting_block``, so
    the Pyodide twin agrees (``tests/webapp/test_silting_p67.py`` /
    ``tests/gui/test_silting_runner_twin.py``). ``canonical_key`` is request-derived
    (schema-1 algebra-only; the radius,budget rides in the ``compute`` string, no new
    request field).
  * 2026-08-08 (``exceptional_sequences_kA3`` ADDED, Plan 65 / R27+R28): a NEW
    algebra-only budget kind. kA3 (linear 1->2->3) over GF(7), ``exceptional_sequences
    :512`` -- the shared block dispatches BOTH halves: classical hereditary
    (``count == 16``, ``closed_form_count == 16 = n! h^n / |W|``, ``transitive ==
    true``, Dynkin ``A_3``) and tau-exceptional (``signed_count == 84 = 3! * 14``,
    ``stt_count == 14``). Pure ADDITION: every pre-existing entry was verified
    byte-identical BEFORE the append (the generator round-trips the file bytes, then
    re-dumps with the same ``indent=1`` settings; all 33 prior delegation asserts
    green first). Both runners share the library builder
    (``tautilting.exceptional.exceptional_sequences_block``), so the Pyodide twin
    agrees (``tests/gui/test_exceptional_runner_twin_p65.py``). ``canonical_key`` is
    request-derived (schema-1 algebra-only, the budget rides in the ``compute``
    string, no ``module`` block).
  * 2026-08-08 (``hh1_lie_kronecker`` ADDED, Plan 70 / R11): a NEW algebra-only DIM-budget
    kind. kK2 (the Kronecker quiver 1 => 2, no relations) over QQ, ``compute ==
    ["hh1_lie"]`` -- the shared block reports HH^1 = Der/Inn = sl2: ``dim == 3``,
    ``solvable == false``, ``perfect == true``, ``radical_dim == 0``, ``sl2_count == 1``,
    ``levi_type == "A1"``, ``toral_rank == 1``, with the QQ ``base_change_note``. Pure
    ADDITION: all 38 pre-existing entries were verified byte-identical BEFORE the append
    (the generator round-trips the file bytes, then re-dumps ``indent=1`` order-preserving).
    Both runners share the library builder (``invariants.hh1_lie.hh1_lie_block``), so the
    Pyodide twin agrees (``tests/gui/test_hh1_lie_runner_twin_p70.py``). ``canonical_key``
    is request-derived (schema-1 algebra-only, no ``module`` block).
  * 2026-08-08 (``hh_lie_module_kronecker`` ADDED, Plan 71 / R12): a NEW top-carrying HH
    kind (the ``bracket`` precedent). kK2 over QQ, ``compute == ["hh_lie_module:0..2"]``
    -- the shared block reports HH• as a graded Lie module over HH^1: ``hh_dims ==
    [1, 3, 0]``, ``hh1_dim == 3``, ``module_axiom_ok``/``inner_acts_zero`` true, the
    degree-1 indecomposable-summand entry (one part ``dim == 3``, the sl2-adjoint L(2)),
    and the char-0 weight table (a symmetric sl2-string). Pure ADDITION: all pre-existing
    entries were verified byte-identical BEFORE the append (the generator round-trips the
    file bytes, then re-dumps ``indent=1`` order-preserving). Both runners share the library
    builder (``hochschild.lie_module.hh_lie_module_block``), so the Pyodide twin agrees
    (``tests/gui/test_hh_lie_module_twin_p71.py``). ``canonical_key`` is request-derived
    (schema-2 algebra-only, no ``module`` block).
  * 2026-08-08 (``split_extension_kA2`` + ``arrow_removal_P1`` ADDED, Plan 72 /
    R5+R6): two NEW algebra-only TOP-DEGREE budget kinds. ``split_extension_kA2`` =
    kA2 (1->2) over GF(7), ``split_extension:4`` -- the trivial-extension Hochschild
    LES (``assembled == direct == [3,1,1,1,1]``, ``agrees``/``exact`` true, the
    grading-derivation ``HH^1 != 0`` witness). ``arrow_removal_P1`` = a loop x (x^2=0)
    plus an inert bridge c:1->2 over GF(7), ``arrow_removal:4`` (``removed == ['c']``,
    ``HH_n(A) == HH_n(B) == [3,1,1,1,1]``, ``hom_agrees`` true, ``coh_low_delta[0] ==
    -2`` the disconnection effect). Pure ADDITION: every pre-existing entry was
    verified byte-identical BEFORE the append (the generator round-trips the file
    bytes, then re-dumps with the same ``indent=1`` settings). Both runners share the
    library builders (``split_extension.split_extension_block`` /
    ``arrow_removal.arrow_removal_block``), so the Pyodide twin agrees
    (``tests/gui/test_split_arrow_runner_twin.py``). ``canonical_key`` is
    request-derived (schema-1 algebra-only, the top-degree budget rides in the
    ``compute`` string, no ``module`` block).
  * 2026-08-08 (``barcode_a5``): ADDED for Plan 69 (R33 persistence/TDA bridge) -- the
    new ``barcode`` module-side compute kind. One entry: a schema-2 request drawing the
    forward line ``A_5`` (1->2->3->4->5, QQ) with the filtration ``H_0`` module (dims
    (1,2,1,2,1), Plan-26 per-arrow block maps) + ``compute: ["barcode"]``. The block is
    the interval barcode ``{[1,5] essential, [2,2], [4,4]}`` (LIVE-VERIFIED). Appended AFTER
    dev's ``hh1_lie_kronecker`` (parse dev's dict, append ``barcode_a5`` LAST, re-dump
    ``indent=1``, NEVER ``sort_keys``); all 39 prior entries confirmed byte-identical first.
    Both runners share the library core builder (``quiverlab.modules.barcode.barcode_block``),
    so the Pyodide twin agrees (``tests/gui/test_barcode_runner_twin_p69.py``). ``canonical_key``
    is request-derived (the ``module`` block canonicalizes through the Plan-25 key; no new
    top-level request field).
  * 2026-08-08 (``skew_group_hh_z2dual`` ADDED, Plan 74 / R8): a NEW fixture for the
    ``SkewGroupAlgebra`` construction family + the ``skew_group_hh`` ALGEBRA-level
    TOP-DEGREE budget kind. Z/2 on k[x]/(x^2) (σ: x ↦ −x, arrow scalar −1) over QQ,
    ``compute == ["skew_group_hh:3"]`` -- the Ştefan conjugacy-class decomposition
    ``dims == direct_dims == [1,1,1,1]``, ``agrees == true``, the per-class summands
    (identity ``inv [1,1,1,1]`` + σ-twisted ``inv [0,0,0,0]``). Pure ADDITION: every
    pre-existing entry was verified byte-identical BEFORE the append (the goldens JSON
    round-trips through ``json.dumps(indent=1)`` + newline byte-for-byte, so the new key
    is appended LAST and no existing bytes move). Both runners share the library builder
    (``hochschild.skew_group.skew_group_hh_block``) and the family builder
    (``families.skew_group.build_skew_group_from_params``), so the Pyodide twin agrees
    (``tests/gui/test_skew_group_runner_twin_p74.py``). ``canonical_key`` is
    request-derived (family params; the SkewGroupAlgebra ``generators`` list is
    order-normalized in the schema so two orderings collide, the top-degree budget rides
    in the ``compute`` string -- no new top-level request field).
  * 2026-08-08 (``tau_cluster_kA2`` ADDED, Plan 66 / R29): a NEW algebra-only PAIR-budget
    kind. kA2 (1 -> 2, no relations) over QQ, ``compute == ["tau_cluster:512"]`` -- the
    shared block reports the tau-cluster morphism category W(A): ``object_count == 5`` (=
    #wide), classifying-space ``face_vector == [5, 11, 5]`` (f_0 = #wide, H1), g-fan sphere
    ``g_fan_face_vector == [1, 5, 5]``, ``euler_characteristic == -1``, ``is_kpi1 == true``
    (hereditary Dynkin), and the picture group (3 generators, 1 atom relation,
    ``abelianization_rank == 2``). Pure ADDITION: all 40 pre-existing entries were verified
    byte-identical BEFORE the append (the generator round-trips the file bytes, then re-dumps
    ``indent=1`` order-preserving). Both runners share the library builder
    (``tautilting.cluster_morphism.tau_cluster_block``), so the Pyodide twin agrees
    (``tests/gui/test_tau_cluster_runner_twin.py``). ``canonical_key`` is request-derived
    (schema-1 algebra-only, the pair budget rides in the ``compute`` string, no ``module``
    block).
    top-level request field).  * 2026-08-10 (``deformations_radsq3`` ADDED, Plan 78 / R13): a NEW algebra-only DIM
    BUDGET kind (the ``hh1_lie`` precedent, but with the plan's OWN ``DEFORM_MAXDIM = 32``
    -- NOT P70's 48, because the cost law is different: the obstruction bracket tracks
    HH-RICHNESS x resolution size, not ``A.dim``). The 3-cycle ``rad^2 = 0`` algebra
    (1->2->3->1, all paths of length 2 zero) over QQ, ``compute == ["deformations"]`` --
    the block reports HH^2/HH^3, the obstruction verdict, the Maurer-Cartan description
    and the ``dg_lie is True`` certificate that ``rad^2 = 0`` guarantees unconditionally.
    The CHEAP fixture was chosen on purpose: the plan's obstructed benchmark point
    (``QuantumCI(-1)``, HH^2 = 5) costs 85 s, which does not belong in a golden -- it is
    pinned instead by a ``deep``-marked twin test. Pure ADDITION: the goldens file
    round-trips ``json.dumps(indent=1)`` byte-for-byte, verified BEFORE the append, so the
    new key is appended LAST and no existing bytes move (git: 44 insertions, 0 deletions).
    Both runners share the library builder
    (``hochschild.deformations.deformations_block``), so the Pyodide twin agrees
    (``tests/gui/test_deformations_runner_twin_p78.py``). ``canonical_key`` is
    request-derived (schema-2 algebra-only, the dim budget rides in the ``compute``
    string -- no new top-level request field).
  * 2026-08-08 (``han_transport_ex53`` ADDED, Plan 73 -- Han bounded-extension
    transport): a NEW fixture for the ``han_transport`` algebra-level CERTIFICATE
    kind carrying the request-level ``new_arrows`` field (CLMS Ex. 5.3, ``F = {a}``
    over GF(32003) -- ``transport == 'bounded'``, tensor-nilpotent index 2,
    right-``B``-projective, ``pd_{B^e}`` finite via ``gl.dim B = 2``). Pure ADDITION:
    the 42 existing entries were verified byte-identical BEFORE the append (indent=1,
    order-preserving). ``new_arrows`` is dropped from ``model_dump`` when absent, so
    every existing request's ``canonical_key`` is byte-unchanged; both runners share
    ``invariants.han.han_transport_block``, so the Pyodide twin agrees
    (``tests/gui/test_han_runner_twin_p73.py``).
  * 2026-08-10 (``incidence_cohomology_b3`` ADDED, Plan 75 -- HH^* of an incidence
    algebra via the ORDER COMPLEX): a NEW fixture for the ``incidence_cohomology``
    degree-range kind on the Boolean lattice ``B_3`` over QQ (``dims == [1,0,0,0]``,
    face vector ``(8,19,18,6)``, contractible by the global bound). Pure ADDITION: the
    48 existing entries were verified byte-identical BEFORE the append (indent=1,
    order-preserving; git reported 74 insertions, 0 deletions). The poset input mode
    emits the SHIPPED ``IncidenceAlgebra`` family with ``elements`` normalized to
    ``null`` when the covers name every element, so NO new request field and NO new
    algebra ``kind`` exists and every pre-P75 request keys byte-unchanged; both runners
    share ``hochschild.simplicial.incidence_cohomology_block``, so the Pyodide twin
    agrees (``tests/webapp/test_incidence_p75.py``,
    ``tests/gui/test_incidence_runner_twin.py``).
  * 2026-08-17 (``koszul_kx3`` ADDED, Plan 77 -- generalized Koszulity): a NEW fixture
    for the ``koszul`` algebra-level kind on ``k[x]/(x^3)`` over QQ (``n_homogeneous ==
    3``; internal generation degrees ``[0,1,3,4,6,7,9,10,12]`` = Berger's ``delta(n)``;
    K2 generators in degrees ``[1,2]``; ``almost_koszul`` refused with ``q = 1``).
    Pure ADDITION: the 49 existing entries were verified byte-identical BEFORE the
    append (indent=1, order-preserving; git reported 28 insertions, 0 deletions), and
    the whole delegation suite was run green BEFORE generating it. ``koszul`` adds NO
    request field -- it is a new compute-list STRING parsed by the generic ``name:0..N``
    grammar (the ``ext_algebra`` sibling), so every pre-P77 request keys byte-unchanged;
    both runners share ``modules.nkoszul.koszul_profile_block``, so the Pyodide twin
    agrees (``tests/webapp/test_koszul_p77.py``, ``tests/gui/test_koszul_runner_twin.py``).
    NB: this golden's ``citations`` payload depends on the Plan-77 bib entries carrying
    a TRAILING COMMA on their last field -- without it the parser drops that field (see
    the named backlog item; 20 pre-existing entries are still affected).
  * 2026-08-17 (``cluster_category_kA3`` ADDED, Plan 79 -- the Amiot-Keller cluster
    category): a NEW fixture for the ``cluster_category`` ALGEBRA-ONLY, BUDGET-carrying
    kind on ``kA3`` over QQ (``num_indec == 9`` = the almost-positive roots;
    ``num_cluster_tilting.count == 14`` = Catalan(4), certified natively;
    ``cluster_tilted.dim == 6`` for the ``mu_1`` mutation = ``kZ3/J^2``; the 2-CY AR
    formula holding on all 36 ordered pairs). Pure ADDITION: the 50 existing entries were
    verified byte-identical BEFORE the append (indent=1, order-preserving; git reported 32
    insertions, 0 deletions) and the whole delegation suite was run green first.
    ``cluster_category`` adds NO request field -- it is a compute-list STRING on the
    budget grammar (the ``wall_chamber`` precedent), parsed identically in all THREE
    grammar sites -- so every pre-P79 request keys byte-unchanged; both runners share
    ``cluster.category.cluster_category_block``, so the Pyodide twin agrees
    (``tests/webapp/test_cluster_category_kind_p79.py``,
    ``tests/gui/test_cluster_runner_twin_p79.py``)."""
import json
import pathlib

import pytest

from webapp.server.cache import canonical_key
from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

_GOLDENS = json.loads(
    (pathlib.Path(__file__).parent / "_runner_goldens.json").read_text(encoding="utf-8"))
# The library version the goldens were frozen under (kept explicit so the cache-key
# pin is independent of the running version -- the key must reproduce the frozen one).
_V = "0.1.0.dev0"


@pytest.mark.parametrize("name", sorted(_GOLDENS))
def test_result_dict_is_byte_identical(name, tmp_path):
    g = _GOLDENS[name]
    req = ComputeRequest.model_validate(g["body"])
    result = run_spec(req, tmp_path)
    got = json.dumps(result, sort_keys=True, default=str)
    assert got == g["result_json"], f"result dict drifted for {name!r}"


@pytest.mark.parametrize("name", sorted(_GOLDENS))
def test_canonical_key_is_unchanged(name):
    g = _GOLDENS[name]
    req = ComputeRequest.model_validate(g["body"])
    key = canonical_key(req.model_dump(by_alias=True), _V)
    assert key == g["canonical_key"], f"cache key drifted for {name!r}"
