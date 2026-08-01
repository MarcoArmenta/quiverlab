"""Every curated request must be reachable from the GUI compose path: kinds in
the GUI's push order, and the seed cache key a fresh seeding run would write must
be a stable fixed point (Plan 35 §5).

``_GUI_ORDER`` mirrors ``webapp/static/gui/gui.js::buildRequest`` (the offline GUI
twin ``docs/gui/gui.js`` is byte-identical). Push order there:

  * hh_cohomology, hh_homology                                     (gui.js:638-639)
  * cup, cap, bracket, connes_b -- the Plan-35 products, inserted   (gui.js:642-645)
    right after hh_homology; Task 10 landed them with a do-not-reorder comment
  * cartan, coxeter_polynomial, global_dimension, center           (gui.js:646)
  * the module-panel kinds                                         (gui.js:652-665)

Two kinds that have no GUI checkbox but are valid seed kinds are slotted at their
natural position: ``dimension`` (an algebra scalar, beside the other scalars); the
resolution/dimension relative order follows the curated seed convention
(``projective_resolution`` before ``projective_dimension``). The test enforces that
every seed request is INTERNALLY CONSISTENT with this one order and carries the
products where they are feasible.

Products carried per example is a FEASIBILITY decision (Task-12 probe, box ~120s;
every trim/omission is recorded beside the entry in ``webapp/precomputed/
manifest.yaml``). The four tractable examples carry the full product surface
(cup/cap/bracket/connes_b); the two deep bar-blow-up examples (``_DEEP_NO_PRODUCTS``,
dim >= 220) carry NONE -- the bar/tt product route's ``to_engine`` setup alone is
~290s (kz20) and the degree-2 cochain basis is 10.5M cells (max_cells forces CS or
OOM), so no product finishes the probe box at any degree, and bracket/connes_b have
no CS route at all. See the manifest for the per-example reasons.

Key stability: ``container/seed_cache.py`` writes the cache key as
``canonical_key(ComputeRequest.model_validate(request.json).model_dump(by_alias=
True))`` -- it VALIDATES then DUMPS. A live GUI/API request travels the same
validate->dump path, so the two share one key iff the canonicalisation is a fixed
point (re-validating the dumped form and re-dumping is byte-stable). That fixed
point -- not raw-body equality -- is the seeding invariant (a raw ``module`` block
carries ``side`` inside ``builtin`` and omits null ``dims``/``maps``; model_dump
normalises both, so raw != dump for module requests, yet the key is still stable)."""
import json
import pathlib

import pytest

from webapp.server.cache import canonical_key, library_version
from webapp.server.schema import ComputeRequest, parse_compute_item

_EX = pathlib.Path(__file__).resolve().parents[2] / "webapp" / "precomputed" / "examples"
_GUI_ORDER = ["hh_cohomology", "hh_homology", "cup", "cap", "bracket",
              "connes_b", "cartan", "coxeter_polynomial", "global_dimension",
              "center", "dimension", "dimension_vector", "rad_top_soc", "tau",
              "tau_minus", "projective_resolution", "injective_resolution",
              "projective_dimension", "injective_dimension", "decompose",
              "ext", "tor"]

# The deep bar-blow-up examples that carry NO products (Task-12 feasibility probe).
_DEEP_NO_PRODUCTS = {"nakayama-kz20-deep", "nakayama-kz24-deep"}
_PRODUCTS = {"cup", "cap", "bracket", "connes_b"}


@pytest.mark.parametrize("bundle", sorted(p.name for p in _EX.iterdir() if p.is_dir()))
def test_compute_list_is_in_gui_order(bundle):
    body = json.loads((_EX / bundle / "request.json").read_text(encoding="utf-8"))
    kinds = [parse_compute_item(s).kind for s in body["compute"]]
    order = {k: i for i, k in enumerate(_GUI_ORDER)}
    assert kinds == sorted(kinds, key=order.__getitem__), \
        f"{bundle}: compute list must follow the GUI compose order"
    if bundle in _DEEP_NO_PRODUCTS:
        assert not (_PRODUCTS & set(kinds)), \
            f"{bundle}: deep bar-blow-up example carries no products (see manifest.yaml)"
    else:
        assert {"cup", "cap", "connes_b"} <= set(kinds), \
            f"{bundle}: Plan 35 -- tractable curated examples carry the products"


@pytest.mark.parametrize("bundle", sorted(p.name for p in _EX.iterdir() if p.is_dir()))
def test_seed_key_is_a_stable_fixed_point(bundle):
    body = json.loads((_EX / bundle / "request.json").read_text(encoding="utf-8"))
    # The seed key, exactly as container/seed_cache.py computes it.
    dump1 = ComputeRequest.model_validate(body).model_dump(by_alias=True)
    seed_key = canonical_key(dump1, library_version())
    # A live request re-travels validate->dump; the key must be a fixed point.
    dump2 = ComputeRequest.model_validate(dump1).model_dump(by_alias=True)
    assert canonical_key(dump2, library_version()) == seed_key, \
        f"{bundle}: the seed cache key must be a stable validate->dump fixed point"
