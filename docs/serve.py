#!/usr/bin/env python3

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "typer",
#   "rich",
#   "livereload",
# ]
# ///

import logging
import re
import shlex
from pathlib import Path
from typing import Annotated

import livereload
import subprocess
import typer
from rich.console import Console
from rich.logging import RichHandler


logger = logging.getLogger("docs.serve")

app = typer.Typer()
console = Console()


def build_docs(build_dir: Path) -> None:
    subprocess.run(["cmake", "--build", build_dir, "--target", "docs"], check=True)


def find_source_dir(build_dir: Path) -> Path:
    cmake_cache = build_dir / "CMakeCache.txt"
    if not cmake_cache.is_file():
        logger.error("CMake cache file %s is missing", cmake_cache)
        raise typer.Exit(1)

    source_dir: Path | None = None
    with cmake_cache.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("Acts_SOURCE_DIR:STATIC="):
                source_dir = Path(line.split("=", 1)[1].strip())
                break

    if source_dir is None or not source_dir.is_dir():
        logger.error("Could not determine source directory from CMake cache")
        raise typer.Exit(1)

    return source_dir


@app.command()
def main(
    build_dir: Annotated[Path, typer.Option(file_okay=False, exists=True)] = Path(
        "build"
    ),
    port: int = 8000,
    verbose: Annotated[bool, typer.Option("-v", "--verbose")] = False,
):

    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(RichHandler())

    build_dir = build_dir.resolve()
    logger.info("Using build directory: %s", build_dir)
    source_dir = find_source_dir(build_dir)
    logger.debug("Source directory: %s", source_dir)
    docs_dir = (source_dir / "docs").resolve()
    if not docs_dir.is_dir():
        logger.error("Docs directory %s is missing", docs_dir)
        raise typer.Exit(1)

    logger.info("Docs directory: %s", docs_dir)

    logger.info("Initial documentation build...")
    build_docs(build_dir)

    output_dir = build_dir / "docs/html"
    if not output_dir.exists():
        console.print("[red]Error: output directory does not exist.[/red]")
        raise typer.Exit(1)

    doxygen_inputs, file_patterns = parse_doxygen_lists(build_dir)
    logger.debug(doxygen_inputs)
    logger.debug(file_patterns)

    server = livereload.Server()

    def rebuild(changed=None):
        if changed is not None:
            logger.info("Changes detected: %s => rebuilding...", changed)
        else:
            logger.info("Rebuilding documentation...")
        build_docs(build_dir)

    for target in expand_watch_patterns(docs_dir, doxygen_inputs, file_patterns):
        logger.debug("Watching %s", target)
        server.watch(target, rebuild)

    extra_files = [
        docs_dir / "header.html",
        docs_dir / "Doxyfile.in",
        docs_dir / "DoxygenLayout.xml",
        source_dir / "CONTRIBUTING.md",
    ]

    for extra in extra_files:
        if extra.exists():
            server.watch(str(extra), rebuild)

    server.watch(str(docs_dir / "examples/**/*.cpp"), rebuild)
    server.watch(str(docs_dir / "examples/**/*.py"), rebuild)
    server.watch(str(docs_dir / "*.bib"), rebuild)
    server.watch(str(docs_dir / "*.js"), rebuild)
    server.watch(str(docs_dir / "*.css"), rebuild)

    server.serve(root=str(output_dir), port=port)


def parse_doxygen_lists(build_dir: Path) -> tuple[list[str], list[str]]:
    """Return the INPUT and FILE_PATTERNS entries from the generated Doxyfile."""
    doxyfile = (build_dir / "docs" / "Doxyfile").resolve()
    if not doxyfile.is_file():
        raise FileNotFoundError(f"Could not find Doxyfile at {doxyfile}")

    targets = {"INPUT": [], "FILE_PATTERNS": []}
    collecting: str | None = None

    def _rstrip_continuation(value: str) -> tuple[str, bool]:
        value = value.rstrip()
        if value.endswith("\\"):
            return value[:-1].rstrip(), True
        return value, False

    def _extend(target_list: list[str], chunk: str) -> None:
        if not chunk:
            return
        tokens = shlex.split(chunk, comments=False, posix=True)
        for token in tokens:
            cleaned = token.rstrip(",")
            if cleaned:
                target_list.append(cleaned)

    with doxyfile.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            if collecting:
                chunk, carry_on = _rstrip_continuation(line)
                _extend(targets[collecting], chunk)
                if not carry_on:
                    collecting = None
                continue

            match = re.match(r"^(INPUT|FILE_PATTERNS)\s*(\+?=)\s*(.*)$", line)
            if not match:
                continue

            tag = match.group(1)
            chunk, carry_on = _rstrip_continuation(match.group(3))
            _extend(targets[tag], chunk)
            if carry_on:
                collecting = tag

    return targets["INPUT"], targets["FILE_PATTERNS"]


def expand_watch_patterns(
    docs_dir: Path, inputs: list[str], file_patterns: list[str]
) -> list[str]:
    """Return file or glob patterns for each Doxygen INPUT entry."""
    patterns = [pattern for pattern in (file_patterns or []) if pattern]
    resolved_docs = docs_dir.resolve()
    watch_targets: dict[str, None] = {}

    for entry in inputs:
        if not entry:
            continue

        path = Path(entry)
        if not path.is_absolute():
            path = (resolved_docs / path).resolve()
        else:
            path = path.resolve()

        if path.is_file():
            watch_targets[str(path)] = None
            continue

        if path.is_dir():
            if patterns:
                for pattern in patterns:
                    watch_targets[str(path / "**" / pattern)] = None
            else:
                watch_targets[str(path / "**" / "*")] = None
            continue

        logger.warning("Doxygen input %s does not exist; skipping.", path)

    return list(watch_targets.keys())


if "__main__" == __name__:
    app()
