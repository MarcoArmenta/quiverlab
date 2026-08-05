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
    ``canonical_key`` is request-derived and UNCHANGED (no new request fields)."""
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
