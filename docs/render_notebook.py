#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "ipykernel",
#     "jupytext",
#     "nbconvert",
#     "typer",
# ]
# ///

# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Execute a jupytext notebook and render it as a doxygen page.

Several of the sympy-based code generators are written in jupytext's `percent`
format: plain Python that the build runs as a script, and that a notebook
viewer opens as a notebook whose prose explains the derivation. This turns one
of those into a markdown page for the doxygen site, so the explanation is
published with the outputs -- matrices, operation counts, the generated C++ --
filled in by actually running it.

The pipeline is jupytext -> execute -> nbconvert markdown -> doxygen-flavoured
markdown. Only the last step is specific to us: nbconvert emits LaTeX between
`$` delimiters, which doxygen's markdown parser mangles before MathJax ever
sees it, so those are rewritten to doxygen's `\\f$ ... \\f$` form, which
doxygen passes through verbatim.

Anything this needs beyond its own dependencies -- everything the notebook
imports -- has to be in the environment it is run in, e.g. via `uv run --with`.
"""

import re
import sys
from pathlib import Path
from typing import Annotated

import jupytext
import typer
from nbconvert import MarkdownExporter
from nbconvert.preprocessors import ExecutePreprocessor

app = typer.Typer(add_completion=False)

FENCE = re.compile(r"^\s*```")
# nbconvert renders stream and text/plain output as indented blocks; they are
# code, not prose, and must not be touched.
INDENTED = re.compile(r"^ {4}\S")
DISPLAY_MATH = re.compile(r"\$\$(.+?)\$\$", re.DOTALL)
INLINE_MATH = re.compile(r"(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)", re.DOTALL)


def to_doxygen_math(markdown: str) -> str:
    """Rewrite `$...$` and `$$...$$` into doxygen's formula commands.

    Doxygen only leaves LaTeX alone inside `\\f$ ... \\f$` (inline) and
    `\\f[ ... \\f]` (display); everywhere else it eats the backslashes and
    underscores first, and then errors out on what is left.

    Formulas are matched per paragraph rather than per line, because prose
    wraps and an inline formula can straddle the break. A blank line still
    ends the search: pairing `$` across a paragraph boundary would join two
    unrelated formulas into one, which is worse than not matching at all.
    """

    result = []
    paragraph = []
    in_fence = False

    def flush_paragraph():
        if not paragraph:
            return
        chunk = "\n".join(paragraph)
        chunk = DISPLAY_MATH.sub(r"\\f[\1\\f]", chunk)
        chunk = INLINE_MATH.sub(r"\\f$\1\\f$", chunk)
        if "$" in chunk.replace("\\f$", ""):
            # An odd delimiter: whatever was meant as a formula is now half
            # rewritten. Doxygen will fail on it, so say why here first.
            print(
                f"warning: unpaired '$' in\n{chunk}\n"
                "put literal dollar signs in `code spans`",
                file=sys.stderr,
            )
        result.append(chunk)
        paragraph.clear()

    for line in markdown.split("\n"):
        if FENCE.match(line):
            flush_paragraph()
            in_fence = not in_fence
            result.append(line)
        elif in_fence or not line.strip() or INDENTED.match(line):
            flush_paragraph()
            result.append(line)
        else:
            paragraph.append(line)
    flush_paragraph()

    return "\n".join(result)


def add_page_id(markdown: str, page_id: str, note: str) -> str:
    """Label the leading `# Title` so other pages can `@ref` it.

    The label is also what doxygen names the output file, giving the page a URL
    that does not move when the build directory does -- but only when the
    labelled title really is the first thing in the file, which is why `note`
    goes underneath it rather than above.
    """

    lines = markdown.split("\n")
    for i, line in enumerate(lines):
        if line.startswith("# "):
            lines[i] = f"{line.rstrip()} {{#{page_id}}}"
            lines.insert(i + 1, f"\n<!-- {note} -->")
            return "\n".join(lines)
    raise ValueError("notebook has no leading '# Title' markdown heading")


@app.command()
def main(
    source: Annotated[
        Path,
        typer.Argument(help="jupytext notebook to execute, in percent format"),
    ],
    output: Annotated[Path, typer.Argument(help="markdown file to write for doxygen")],
    page_id: Annotated[
        str, typer.Option(help="doxygen page label to attach to the title")
    ],
    timeout: Annotated[
        int, typer.Option(help="seconds a single cell may run for")
    ] = 600,
) -> None:
    notebook = jupytext.read(source)

    # Run in the notebook's own directory: relative paths in a notebook are
    # written with respect to where it lives, not where the build happens to be.
    ExecutePreprocessor(timeout=timeout, kernel_name="python3").preprocess(
        notebook, {"metadata": {"path": str(source.parent)}}
    )

    markdown, resources = MarkdownExporter().from_notebook_node(notebook)
    markdown = add_page_id(
        to_doxygen_math(markdown),
        page_id,
        f"Generated from {source.name} by docs/render_notebook.py. "
        "Edit that notebook, not this file.",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(markdown)

    # Figures, if the notebook ever grows any. The docs build puts this
    # directory on doxygen's IMAGE_PATH so they resolve.
    for name, data in (resources.get("outputs") or {}).items():
        (output.parent / name).write_bytes(data)


if __name__ == "__main__":
    app()
