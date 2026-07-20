"""Generate the docs References page from the single packaged bibliography
(src/quiverlab/citations/references.bib) via the library's own accessor. Run by
the gen-files plugin; writes a virtual bibliography.md only."""
import mkdocs_gen_files

from quiverlab import bibliography              # Plan 06: grouped/annotated text
from quiverlab.citations import references_bib_path

_raw = references_bib_path().read_text(encoding="utf-8")
with mkdocs_gen_files.open("bibliography.md", "w") as fd:
    fd.write("# References\n\n")
    fd.write("The curated, verified bibliography that quiverlab cites, rendered from "
             "the single packaged source `src/quiverlab/citations/references.bib`.\n\n")
    fd.write(str(bibliography()))                # annotated, grouped by algorithm/family/field
    fd.write("\n\n## Raw BibTeX\n\n```bibtex\n" + _raw + "\n```\n")
