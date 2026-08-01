"""Two guarantees over the curated seed requests (Plan 35 §5):
  (1) each request's compute list is in the CURATED SEED CONVENTION order
      (``_GUI_ORDER``), with the Plan-35 PRODUCTS in exact gui.js relative order
      right after ``hh_homology``; and
  (2) the seed cache key a fresh seeding run writes is a stable validate->dump
      fixed point.

``_GUI_ORDER`` is the CURATED SEED CONVENTION -- the fixed order the stored
``request.json`` files use and that the ``validate->dump`` canonicalisation
preserves. It is NOT byte-for-byte the ``webapp/static/gui/gui.js::buildRequest``
push order. It agrees with gui.js on the part that matters for Task 12 -- the
products ``cup, cap, bracket, connes_b`` are placed in exact gui.js relative order
immediately after ``hh_homology`` (gui.js:642-645; Task 10 landed them with a
do-not-reorder comment) -- and that products-after-hh placement is what the order
test guarantees going forward. It DIVERGES from gui.js on two pre-existing points
(seed-convention artifacts, NOT to be "fixed" by reordering requests -- that would
re-key all six bundles and lose content):

  * **resolutions-first.** ``_GUI_ORDER`` lists ``projective_resolution`` /
    ``injective_resolution`` BEFORE ``projective_dimension`` /
    ``injective_dimension`` / ``decompose``. gui.js pushes the module panel in the
    opposite block order -- ``projective_dimension, injective_dimension, decompose``
    (gui.js:652-655) come first, THEN ``projective_resolution, injective_resolution``
    (gui.js:656-659). So in gui.js ``decompose`` PRECEDES the resolutions; here it
    follows them.
  * **bare ``dimension``.** The curated requests carry a ``dimension`` kind (the
    algebra's total dimension) slotted beside the other scalars. gui.js has no
    checkbox for it and never pushes it.

Consequence (honest, intended): the curated requests are NOT fully GUI-composable
-- a user cannot reproduce them by ticking boxes (``dimension`` is not offered, and
the resolution/dimension/decompose block would come out in gui.js order). This test
is therefore NOT a GUI-composability proof; it pins the seed convention + the
products' GUI-relative placement + the key fixed point. The order test enforces that
every seed request is INTERNALLY CONSISTENT with this one convention and carries the
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
# The CURATED SEED CONVENTION order -- NOT the literal gui.js push order (see the
# module docstring): products in exact gui.js relative order after hh_homology, but
# it diverges from gui.js by listing the resolutions BEFORE the dimensions/decompose
# and by carrying a bare ``dimension`` kind the GUI never pushes. Do not "align" it
# to gui.js by reordering the requests -- that re-keys all six bundles.
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
        f"{bundle}: compute list must follow the curated seed-convention order " \
        "(_GUI_ORDER; products in gui.js relative order, see module docstring)"
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
