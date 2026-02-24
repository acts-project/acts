#!/usr/bin/env python
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "gcovr==8.6",
#     "lxml",
#     "rich",
#     "typer",
# ]
# ///
import multiprocessing as mp
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
import shlex

from typing import Annotated

import typer
from rich.console import Console

app = typer.Typer(add_completion=False)
console = Console()

# Regex patterns for excluding files from coverage (matched against file paths as-is)
EXCLUDE_PATTERNS = [
    r"/boost/",
    r"json\.hpp",
    "/dependencies/",
]

# Paths relative to source directory (resolved to absolute when source_dir is known)
EXCLUDE_PATHS = [
    "Tests/",
    "Python/",
    "dependencies/",
]


def _resolve_excludes(source_dir: Path) -> list[str]:
    """Return exclude patterns: EXCLUDE_PATTERNS as-is plus EXCLUDE_PATHS prefixed with source_dir."""
    source_prefix = re.escape(source_dir.as_posix()) + r"/"
    return list(EXCLUDE_PATTERNS) + [source_prefix + p for p in EXCLUDE_PATHS]


def locate_executable(name: str, hint: str) -> str:
    path = shutil.which(name)
    if not path:
        console.print(hint, style="red")
        raise typer.Exit(1)
    return path


def gcovr_version(gcovr_exe: str) -> tuple[int, int] | None:
    version_text = subprocess.check_output([gcovr_exe, "--version"], text=True).strip()
    match = re.match(r"gcovr (\d+)\.(\d+)", version_text)
    if not match:
        console.print(
            f"Unexpected gcovr version output: {version_text}",
            style="yellow",
        )
        return None
    return (int(match.group(1)), int(match.group(2)))


@app.command()
def generate(
    build_dir: Annotated[Path, typer.Argument(help="CMake build directory")],
    gcov: Annotated[
        str | None,
        typer.Option(help="Path to gcov executable"),
    ] = None,
    jobs: Annotated[
        int, typer.Option("--jobs", "-j", help="Number of parallel jobs")
    ] = mp.cpu_count(),
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
    filter_xml: Annotated[
        bool,
        typer.Option(
            "--filter/--no-filter", help="Filter the coverage XML after generation"
        ),
    ] = False,
    html: Annotated[
        bool,
        typer.Option("--html/--no-html", help="Also generate HTML coverage report"),
    ] = True,
    html_theme: Annotated[
        str,
        typer.Option("--html-theme", "-t", help="HTML report theme (e.g. github.blue)"),
    ] = "github.blue",
) -> None:
    """Generate SonarQube XML and optionally HTML coverage reports from a CMake build directory using gcovr."""
    build_dir = build_dir.resolve()
    if not build_dir.is_dir():
        console.print(f"Build directory not found: {build_dir}", style="red")
        raise typer.Exit(1)
    if not (build_dir / "CMakeCache.txt").exists():
        console.print(
            f"Build directory missing CMakeCache.txt: {build_dir}",
            style="red",
        )
        raise typer.Exit(1)

    gcov_exe = gcov or locate_executable(
        "gcov",
        "gcov not installed. Install GCC coverage tooling or pass --gcov.",
    )
    gcovr_exe = locate_executable(
        "gcovr",
        "gcovr not installed. Use 'uv run --script' or install gcovr.",
    )

    version = gcovr_version(gcovr_exe)
    if version is not None:
        console.print(f"Found gcovr version {version[0]}.{version[1]}")
        if version < (5, 0):
            console.print(
                "Consider upgrading to a newer gcovr version.", style="yellow"
            )
        elif version == (5, 1):
            console.print(
                "Version 5.1 does not support parallel processing of gcov data.",
                style="red",
            )
            raise typer.Exit(1)

    coverage_dir = build_dir / "coverage"
    coverage_dir.mkdir(exist_ok=True)

    coverage_xml_path = coverage_dir / "cov.xml"
    raw_xml_path = coverage_dir / "cov_raw.xml" if filter_xml else coverage_xml_path

    with tempfile.TemporaryDirectory() as gcov_obj_dir:
        base_args = _build_gcovr_common_args(
            build_dir, gcov_exe, gcovr_exe, jobs, verbose, gcov_obj_dir
        )
        gcovr_cmd = base_args + ["--sonarqube", str(raw_xml_path)]
        if html:
            html_dir = coverage_dir / "html"
            html_dir.mkdir(exist_ok=True)
            html_path = html_dir / "index.html"
            gcovr_cmd += [
                "--html-nested",
                str(html_path),
                "--html-theme",
                html_theme,
            ]

        console.print(f"$ {shlex.join(gcovr_cmd)}")
        subprocess.run(gcovr_cmd, cwd=build_dir, check=True)

    if html:
        console.print(f"HTML coverage report written to {coverage_dir / 'html'}")

    if filter_xml:
        source_dir = Path(__file__).resolve().parent.parent
        xml_excludes = _resolve_excludes(source_dir) + ["^" + re.escape(build_dir.name)]
        filter_coverage_xml(raw_xml_path, coverage_xml_path, xml_excludes)
        raw_xml_path.unlink()
        console.print(f"Removed raw coverage file {raw_xml_path}")


def filter_coverage_xml(
    input_path: Path, output_path: Path, excludes: list[str]
) -> None:
    from lxml import etree

    patterns = [re.compile(p) for p in excludes]

    tree = etree.parse(input_path)  # noqa: S320
    root = tree.getroot()

    removed = 0
    for file_elem in root.findall("file"):
        path = file_elem.get("path", "")
        if any(p.search(path) for p in patterns):
            root.remove(file_elem)
            removed += 1

    remaining = len(root.findall("file"))
    console.print(f"Removed {removed} file entries, {remaining} remaining")

    deduped = 0
    for file_elem in root.findall("file"):
        lines: dict[int, etree._Element] = {}
        duplicates: list[etree._Element] = []
        for line_elem in file_elem.findall("lineToCover"):
            line_num = int(line_elem.get("lineNumber"))
            if line_num not in lines:
                lines[line_num] = line_elem
            else:
                existing = lines[line_num]
                if line_elem.get("covered") == "true":
                    existing.set("covered", "true")
                for attr in ("branchesToCover", "coveredBranches"):
                    new_val = line_elem.get(attr)
                    if new_val is not None:
                        old_val = existing.get(attr)
                        if old_val is None or int(new_val) > int(old_val):
                            existing.set(attr, new_val)
                duplicates.append(line_elem)
                deduped += 1
        for dup in duplicates:
            file_elem.remove(dup)

    if deduped:
        console.print(f"Deduplicated {deduped} lineToCover entries")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output_path, xml_declaration=True, encoding="utf-8")
    console.print(f"Filtered coverage written to {output_path}")


@app.command()
def filter(
    input: Annotated[Path, typer.Argument(help="Input coverage XML file")],
    output: Annotated[Path, typer.Argument(help="Output filtered coverage XML file")],
    exclude: Annotated[
        list[str] | None,
        typer.Option("--exclude", "-e", help="Regex pattern to exclude file paths"),
    ] = None,
    build_dir_name: Annotated[
        str | None,
        typer.Option(help="Build directory name to exclude from coverage paths"),
    ] = None,
) -> None:
    """Filter a SonarQube coverage XML file by removing file entries matching exclude patterns."""
    if not input.exists():
        console.print(f"Input file not found: {input}", style="red")
        raise typer.Exit(1)

    source_dir = Path(__file__).resolve().parent.parent
    excludes = list(exclude) if exclude is not None else _resolve_excludes(source_dir)
    if build_dir_name is not None:
        excludes.append("^" + re.escape(build_dir_name))

    filter_coverage_xml(input, output, excludes)


def _build_gcovr_common_args(
    build_dir: Path,
    gcov_exe: str,
    gcovr_exe: str,
    jobs: int,
    verbose: bool,
    gcov_object_directory: str,
) -> list[str]:
    script_dir = Path(__file__).resolve().parent
    source_dir = script_dir.parent.resolve()

    version = gcovr_version(gcovr_exe)
    extra_flags: list[str] = []
    if version is not None:
        if version >= (6, 0):
            extra_flags.append("--exclude-noncode-lines")
    if verbose:
        extra_flags.append("--verbose")

    excludes: list[str] = []
    for pattern in _resolve_excludes(source_dir):
        excludes.extend(["-e", pattern])
    excludes.extend(["-e", f"{build_dir.as_posix()}/"])

    return (
        [gcovr_exe]
        + ["-r", str(source_dir)]
        + ["--gcov-executable", gcov_exe]
        + ["--gcov-object-directory", gcov_object_directory]
        + ["-j", str(jobs)]
        + ["--merge-mode-functions", "separate"]
        + ["--gcov-ignore-errors", "source_not_found"]
        + ["--gcov-ignore-parse-errors", "suspicious_hits.warn"]
        + excludes
        + extra_flags
    )


if __name__ == "__main__":
    app()
