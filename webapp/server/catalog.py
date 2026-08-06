"""Catalog built by introspecting ``quiverlab.families()``. Generated, not
hand-maintained, so it cannot diverge from the library's real offering.

The real ``families()`` (confirmed by probing the installed library at Plans
10-19 state, per the 2026-07-24 amendment) returns a ``FamilyListing`` of
``FamilyInfo`` records with fields ``name``, ``signature``, ``route``,
``citations``, ``summary`` -- NOT ``(name, builder)`` pairs. Each family name is
a top-level ``quiverlab`` export, so the builder is ``getattr(quiverlab,
info.name)`` and its parameters come from ``inspect.signature(builder)``. The
``zoo`` entry is excluded: it yields a batch of algebras, not a single buildable
family. This module imports only the PUBLIC ``quiverlab`` surface."""
from __future__ import annotations

import inspect
from functools import lru_cache

import quiverlab as ql


class CatalogError(ValueError):
    pass


# Families that are NOT parameterized webapp forms: ``zoo`` (a batch iterator) and the
# Plan-44 constructions (Algebra/Module/Potential arguments -- surfaced as drawable
# presets, not a scalar form).
_NON_FORM_FAMILIES = frozenset({
    "BrauerGraphAlgebra",   # P46: graph-structured input, preset-surfaced
    "zoo", "OnePointExtension", "repetitive_slice", "JacobianAlgebra",
    # P48 surfaces: take Triangulation/MarkedSurface/int args (not scalar bool/int/str
    # forms) -- surfaced as drawable presets, the produced gentle algebra flows through
    # every existing compute kind.
    "fan_triangulation", "annulus_triangulation", "jacobian_of",
})


def _iter_families():
    """Yield ``(name, builder)`` for each buildable v1 family.

    ``ql.families()`` yields ``FamilyInfo`` records; the builder is the
    top-level ``quiverlab`` export named ``info.name``. ``zoo`` is skipped (it
    is a batch iterator, not a single family). A missing export is drift -- we
    raise loudly rather than silently drop it."""
    for info in ql.families():
        name = info.name
        # zoo is a batch iterator; the Plan-44 constructions take Algebra/Module/Potential
        # arguments (not scalar bool/int/str), so they do not fit the parameterized webapp
        # FORM -- they are surfaced as drawable presets instead (gui_build_hook), and skipped
        # here so the form builder never introspects their non-scalar signatures.
        if name in _NON_FORM_FAMILIES:
            continue
        builder = getattr(ql, name, None)
        if builder is None:
            raise CatalogError(
                f"families() lists {name!r} but quiverlab has no such export; "
                "the family catalog has drifted from the library surface"
            )
        yield name, builder


def _params_of(builder) -> list[dict]:
    out = []
    try:
        sig = inspect.signature(builder)
    except (TypeError, ValueError):
        return out
    for pname, p in sig.parameters.items():
        if pname == "field":
            continue
        default = None if p.default is inspect.Parameter.empty else p.default
        # bool is an int subclass -- test it FIRST (amendment pt 2).
        if isinstance(default, bool):
            kind = "bool"
        elif isinstance(default, int):
            kind = "int"
        else:
            kind = "str"
        out.append({"name": pname, "kind": kind, "default": default})
    return out


# --------------------------------------------------------------------------- #
# Curated per-family metadata (Marco 2026-07-28: the form must EXPLAIN each
# parameter and prefill an editable value, never demand JSON from scratch).
# ``example`` is the prefill for a parameter with no signature default; every
# example below is VERIFIED to build via quiverlab.hpc.spec.build_algebra
# (tests/webapp/test_catalog.py pins that for every family). Four-language
# (en/es/fr/zh, the i18n LANGS): the form picks help[lang], falling back to
# English. Wrapper families (TensorProduct / TrivialExtension)
# take PathAlgebra type strings -- see spec.py::_WRAPPER_FAMILIES.
# --------------------------------------------------------------------------- #
_TYPE_STRING_HELP = {"en": "Dynkin type string for a path algebra, e.g. A3, A5, D4, E6",
                     "es": "tipo de Dynkin para un álgebra de caminos, p. ej. A3, A5, D4, E6",
                     "fr": "chaîne de type Dynkin pour une algèbre de chemins, p. ex. A3, A5, D4, E6",
                     "zh": "路代数的 Dynkin 型字符串，例如 A3、A5、D4、E6"}
_FAMILY_META: dict[str, dict] = {
    "ExteriorAlgebra": {
        "summary": {"en": "Exterior algebra Λ(kⁿ) on n generators.",
                    "es": "Álgebra exterior Λ(kⁿ) con n generadores.",
                    "fr": "Algèbre extérieure Λ(kⁿ) à n générateurs.",
                    "zh": "n 个生成元上的外代数 Λ(kⁿ)。"},
        "params": {"n": {"example": 3,
                         "help": {"en": "number of generators n ≥ 1 (dimension 2ⁿ)",
                                  "es": "número de generadores n ≥ 1 (dimensión 2ⁿ)",
                                  "fr": "nombre de générateurs n ≥ 1 (dimension 2ⁿ)",
                                  "zh": "生成元个数 n ≥ 1（维数 2ⁿ）"}}},
    },
    "IncidenceAlgebra": {
        "summary": {"en": "Incidence algebra of a finite poset, given by its cover relations.",
                    "es": "Álgebra de incidencia de un conjunto ordenado finito, dado por sus cubiertas.",
                    "fr": "Algèbre d'incidence d'un ensemble ordonné fini, donnée par ses relations de couverture.",
                    "zh": "有限偏序集的关联代数，由其覆盖关系给出。"},
        "params": {
            "poset_or_covers": {"example": [[1, 2], [2, 3]],
                                "help": {"en": "cover relations as pairs [a, b] meaning a < b is a cover",
                                         "es": "cubiertas como pares [a, b]: a < b es una cubierta",
                                         "fr": "relations de couverture sous forme de paires [a, b] signifiant que a < b est une couverture",
                                         "zh": "覆盖关系以数对 [a, b] 给出，表示 a < b 是一个覆盖"}},
            "elements": {"example": None,
                         "help": {"en": "optional: every poset element, e.g. [1, 2, 3] — needed only when some element appears in no cover",
                                  "es": "opcional: todos los elementos, p. ej. [1, 2, 3] — solo necesario si algún elemento no aparece en ninguna cubierta",
                                  "fr": "facultatif : tous les éléments de l'ensemble ordonné, p. ex. [1, 2, 3] — nécessaire seulement lorsqu'un élément n'apparaît dans aucune couverture",
                                  "zh": "可选：偏序集的所有元素，例如 [1, 2, 3]——仅当某个元素不出现在任何覆盖中时才需要"}},
        },
    },
    "NakayamaAlgebra": {
        "summary": {"en": "Nakayama algebra — either from a Kupisch series, or cyclic/linear with n vertices and Loewy length l.",
                    "es": "Álgebra de Nakayama — por serie de Kupisch, o cíclica/lineal con n vértices y longitud de Loewy l.",
                    "fr": "Algèbre de Nakayama — soit à partir d'une série de Kupisch, soit cyclique/linéaire avec n sommets et longueur de Loewy l.",
                    "zh": "Nakayama 代数——可由 Kupisch 序列给出，或为具有 n 个顶点、Loewy 长度 l 的循环/线性代数。"},
        "params": {
            "kupisch": {"example": None,
                        "help": {"en": "Kupisch series as a list, e.g. [3, 3, 2] — leave n and l empty when set",
                                 "es": "serie de Kupisch como lista, p. ej. [3, 3, 2] — deja n y l vacíos si la usas",
                                 "fr": "série de Kupisch sous forme de liste, p. ex. [3, 3, 2] — laissez n et l vides si vous l'utilisez",
                                 "zh": "以列表形式给出的 Kupisch 序列，例如 [3, 3, 2]——设置它时请将 n 与 l 留空"}},
            "n": {"example": 4, "help": {"en": "number of vertices (used with l; leave kupisch empty)",
                                         "es": "número de vértices (con l; deja kupisch vacío)",
                                         "fr": "nombre de sommets (utilisé avec l ; laissez kupisch vide)",
                                         "zh": "顶点个数（与 l 一起使用；将 kupisch 留空）"}},
            "l": {"example": 3, "help": {"en": "Loewy length (radical vanishing degree)",
                                         "es": "longitud de Loewy (grado en que se anula el radical)",
                                         "fr": "longueur de Loewy (degré d'annulation du radical)",
                                         "zh": "Loewy 长度（根基消失的次数）"}},
            "cyclic": {"help": {"en": "cyclic quiver Z̃ₙ instead of the linear Aₙ line",
                                "es": "carcaj cíclico Z̃ₙ en lugar de la línea Aₙ",
                                "fr": "carquois cyclique Z̃ₙ au lieu de la ligne linéaire Aₙ",
                                "zh": "循环箭图 Z̃ₙ，而非线性的 Aₙ 直线"}},
        },
    },
    "PathAlgebra": {
        "summary": {"en": "Path algebra kQ of a Dynkin-type quiver (no relations).",
                    "es": "Álgebra de caminos kQ de un carcaj de tipo Dynkin (sin relaciones).",
                    "fr": "Algèbre de chemins kQ d'un carquois de type Dynkin (sans relations).",
                    "zh": "Dynkin 型箭图的路代数 kQ（无关系）。"},
        "params": {
            "type_or_quiver": {"example": "A4", "help": _TYPE_STRING_HELP},
            "orientation": {"help": {"en": "arrow orientation: linear (default) or alternating",
                                     "es": "orientación de flechas: linear (por defecto) o alternating",
                                     "fr": "orientation des flèches : linear (par défaut) ou alternating",
                                     "zh": "箭头定向：linear（线性，默认）或 alternating（交错）"}},
        },
    },
    "PreprojectiveAlgebra": {
        "summary": {"en": "Preprojective algebra Π(Q) of a Dynkin type (finite-dimensional).",
                    "es": "Álgebra preproyectiva Π(Q) de un tipo Dynkin (de dimensión finita).",
                    "fr": "Algèbre préprojective Π(Q) d'un type Dynkin (de dimension finie).",
                    "zh": "Dynkin 型的预投射代数 Π(Q)（有限维）。"},
        "params": {
            "type_or_quiver": {"example": "A3", "help": _TYPE_STRING_HELP},
            "degree_bound": {"example": None,
                             "help": {"en": "optional truncation degree (chosen automatically for Dynkin types)",
                                      "es": "grado de truncamiento opcional (automático para tipos Dynkin)",
                                      "fr": "degré de troncature facultatif (choisi automatiquement pour les types Dynkin)",
                                      "zh": "可选的截断次数（对 Dynkin 型自动选择）"}},
        },
    },
    "QuantumCI": {
        "summary": {"en": "Quantum complete intersection k⟨x,y⟩/(xᵃ, yᵇ, xy + q·yx).",
                    "es": "Intersección completa cuántica k⟨x,y⟩/(xᵃ, yᵇ, xy + q·yx).",
                    "fr": "Intersection complète quantique k⟨x,y⟩/(xᵃ, yᵇ, xy + q·yx).",
                    "zh": "量子完全交 k⟨x,y⟩/(xᵃ, yᵇ, xy + q·yx)。"},
        "params": {
            "q": {"example": 1,
                  "help": {"en": "exact scalar: an integer, or a fraction as a string like \"1/2\"",
                           "es": "escalar exacto: un entero, o una fracción como cadena, p. ej. \"1/2\"",
                           "fr": "scalaire exact : un entier, ou une fraction sous forme de chaîne, p. ex. \"1/2\"",
                           "zh": "精确标量：一个整数，或以字符串形式给出的分数，例如 \"1/2\""}},
            "a": {"help": {"en": "nilpotency exponent of x (a ≥ 2)",
                           "es": "exponente de nilpotencia de x (a ≥ 2)",
                           "fr": "exposant de nilpotence de x (a ≥ 2)",
                           "zh": "x 的幂零指数 (a ≥ 2)"}},
            "b": {"help": {"en": "nilpotency exponent of y (b ≥ 2)",
                           "es": "exponente de nilpotencia de y (b ≥ 2)",
                           "fr": "exposant de nilpotence de y (b ≥ 2)",
                           "zh": "y 的幂零指数 (b ≥ 2)"}},
        },
    },
    "RadicalSquareZero": {
        "summary": {"en": "Radical-square-zero algebra kQ/J² of a Dynkin-type quiver.",
                    "es": "Álgebra con radical al cuadrado cero kQ/J² de un carcaj de tipo Dynkin.",
                    "fr": "Algèbre à radical carré nul kQ/J² d'un carquois de type Dynkin.",
                    "zh": "Dynkin 型箭图的根基平方为零代数 kQ/J²。"},
        "params": {"quiver": {"example": "A3", "help": _TYPE_STRING_HELP}},
    },
    "TensorProduct": {
        "summary": {"en": "Tensor product A ⊗ B of two path algebras, over the chosen field.",
                    "es": "Producto tensorial A ⊗ B de dos álgebras de caminos, sobre el cuerpo elegido.",
                    "fr": "Produit tensoriel A ⊗ B de deux algèbres de chemins, sur le corps choisi.",
                    "zh": "两个路代数的张量积 A ⊗ B，在所选域上。"},
        "params": {"A": {"example": "A2", "help": _TYPE_STRING_HELP},
                    "B": {"example": "A2", "help": _TYPE_STRING_HELP}},
    },
    "TruncatedPathAlgebra": {
        "summary": {"en": "Truncated path algebra kQ/J^r — all paths of length ≥ r vanish.",
                    "es": "Álgebra de caminos truncada kQ/J^r — todo camino de longitud ≥ r se anula.",
                    "fr": "Algèbre de chemins tronquée kQ/J^r — tout chemin de longueur ≥ r s'annule.",
                    "zh": "截断路代数 kQ/J^r——所有长度 ≥ r 的路都为零。"},
        "params": {
            "type_or_quiver": {"example": "A4", "help": _TYPE_STRING_HELP},
            "r": {"example": 2,
                  "help": {"en": "truncation length r ≥ 2 (paths of length ≥ r are zero)",
                           "es": "longitud de truncamiento r ≥ 2 (los caminos de longitud ≥ r son cero)",
                           "fr": "longueur de troncature r ≥ 2 (les chemins de longueur ≥ r sont nuls)",
                           "zh": "截断长度 r ≥ 2（长度 ≥ r 的路为零）"}},
            "orientation": {"help": {"en": "arrow orientation: linear (default) or alternating",
                                     "es": "orientación de flechas: linear (por defecto) o alternating",
                                     "fr": "orientation des flèches : linear (par défaut) ou alternating",
                                     "zh": "箭头定向：linear（线性，默认）或 alternating（交错）"}},
        },
    },
    "TrivialExtension": {
        "summary": {"en": "Trivial extension T(A) = A ⋉ D(A) of a path algebra — always symmetric.",
                    "es": "Extensión trivial T(A) = A ⋉ D(A) de un álgebra de caminos — siempre simétrica.",
                    "fr": "Extension triviale T(A) = A ⋉ D(A) d'une algèbre de chemins — toujours symétrique.",
                    "zh": "路代数的平凡扩张 T(A) = A ⋉ D(A)——总是对称的。"},
        "params": {"A": {"example": "A3", "help": _TYPE_STRING_HELP}},
    },
}


def _merge_meta(name: str, params: list[dict]) -> dict:
    meta = _FAMILY_META.get(name, {})
    pmeta = meta.get("params", {})
    for p in params:
        m = pmeta.get(p["name"], {})
        if "example" in m:
            p["example"] = m["example"]
        if "help" in m:
            p["help"] = m["help"]
        p["required"] = p["default"] is None and m.get("example") is not None
    out = {"name": name, "params": params, "fields": ["CC", "GF"]}
    if "summary" in meta:
        out["summary"] = meta["summary"]
    return out


@lru_cache(maxsize=1)
def build_catalog() -> dict:
    families = [_merge_meta(name, _params_of(builder))
                for name, builder in _iter_families()]
    return {
        "families": sorted(families, key=lambda f: f["name"]),
        "fields": {"CC": {"needs": []}, "GF": {"needs": ["p"], "optional": ["n"]}},
    }


@lru_cache(maxsize=1)
def family_names() -> frozenset[str]:
    return frozenset(f["name"] for f in build_catalog()["families"])


def validate_family(name: str, params: dict) -> None:
    if name not in family_names():
        raise CatalogError(f"unknown family {name!r}; see GET /api/catalog")
    known = {p["name"] for f in build_catalog()["families"]
             if f["name"] == name for p in f["params"]}
    unknown = set(params) - known
    if unknown:
        raise CatalogError(f"family {name!r} got unknown params {sorted(unknown)}")
