"""Generate one API-reference page per module + a literate-nav SUMMARY (mkdocstrings
recipe). Run by the gen-files plugin at build time; writes virtual files only."""
from pathlib import Path

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()
root = Path(__file__).parent.parent
src = root / "src"

for path in sorted(src.rglob("*.py")):
    module_path = path.relative_to(src).with_suffix("")
    doc_path = path.relative_to(src).with_suffix(".md")
    full_doc_path = Path("reference", doc_path)
    parts = tuple(module_path.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
        doc_path = doc_path.with_name("index.md")
        full_doc_path = full_doc_path.with_name("index.md")
    elif parts[-1] == "__main__":
        continue
    if not parts or parts[-1].startswith("_"):
        continue                                  # skip private modules (_kernels, etc.)
    nav[parts] = doc_path.as_posix()
    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        fd.write(f"::: {'.'.join(parts)}")
    mkdocs_gen_files.set_edit_path(full_doc_path, path.relative_to(root))

# Section landing page: listed FIRST in the SUMMARY so section-index promotes
# it — the sidebar's "API Reference" header becomes clickable and reference/
# stops 404ing.
with mkdocs_gen_files.open("reference/index.md", "w") as fd:
    fd.write(
        "# API Reference\n\n"
        "Every public function and class, generated from the docstrings.\n\n"
        "Start at the [`quiverlab` package](quiverlab/index.md) — its top-level\n"
        "exports (`Quiver`, `Algebra`, `GF`, `CC`, `zoo`, `families`,\n"
        "`bibliography`, `sweep`, …) are what the tutorials and the\n"
        "landing-page GUI use. The subpackages below mirror the source tree.\n"
    )

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
    nav_file.write("* [API Reference](index.md)\n")
    nav_file.writelines(nav.build_literate_nav())
