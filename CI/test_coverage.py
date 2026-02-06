#!/usr/bin/env python
# /// script
# requires-python = ">=3.11"
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
from pathlib import Path
import shlex

from typing import Annotated

import typer
from rich.console import Console

app = typer.Typer(add_completion=False)
console = Console()

DEFAULT_EXCLUDES = [
    r"(^|/)Tests/",
    r"/boost/",
    r"json\.hpp",
    r"(^|/)Python/",
    r"(^|/)dependencies/",
]


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


def validate_coverage_xml(path: Path, schema_path: Path) -> None:
    if not path.exists():
        console.print(f"Coverage XML not found: {path}", style="red")
        raise typer.Exit(1)
    xmllint = shutil.which("xmllint")
    if not xmllint:
        console.print("xmllint not available for XML validation.", style="red")
        raise typer.Exit(1)
    if not schema_path.exists():
        console.print(f"Coverage XSD not found: {schema_path}", style="red")
        raise typer.Exit(1)
    cmd = [xmllint, "--noout", "--schema", str(schema_path), str(path)]
    console.print(f"$ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        error_output = exc.stderr or exc.stdout or "xmllint failed"
        console.print(error_output.strip(), style="red")
        raise typer.Exit(1) from exc


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
) -> None:
    """Generate a SonarQube coverage XML report from a CMake build directory using gcovr."""
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
    extra_flags: list[str] = []
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
        elif version >= (6, 0):
            extra_flags.append("--exclude-noncode-lines")

    if verbose:
        extra_flags.append("--verbose")

    script_dir = Path(__file__).resolve().parent
    source_dir = script_dir.parent.resolve()
    coverage_dir = build_dir / "coverage"
    coverage_dir.mkdir(exist_ok=True)

    source_dir_posix = source_dir.as_posix()
    excludes = [
        "-e",
        f"{source_dir_posix}/Tests/",
        "-e",
        r".*/boost/.*",
        "-e",
        r".*json\.hpp",
        "-e",
        f"{source_dir_posix}/Python/",
        "-e",
        f".*{build_dir.name}.*",
        "-e",
        ".*dependencies.*",
    ]
    gcovr = [gcovr_exe]

    coverage_xml_path = coverage_dir / "cov.xml"
    raw_xml_path = coverage_dir / "cov_raw.xml" if filter_xml else coverage_xml_path
    schema_path = script_dir / "sonar_coverage.xsd"
    gcovr_sonar_cmd = (
        gcovr
        + ["-r", str(source_dir)]
        + ["--gcov-executable", gcov_exe]
        + ["-j", str(jobs)]
        + ["--merge-mode-functions", "separate"]
        + excludes
        + extra_flags
        + ["--sonarqube", str(raw_xml_path)]
    )
    console.print(f"$ {shlex.join(gcovr_sonar_cmd)}")
    subprocess.run(gcovr_sonar_cmd, cwd=build_dir, check=True)

    if filter_xml:
        filter_coverage_xml(raw_xml_path, coverage_xml_path, DEFAULT_EXCLUDES)

    validate_coverage_xml(coverage_xml_path, schema_path)


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
) -> None:
    """Filter a SonarQube coverage XML file by removing file entries matching exclude patterns."""
    if not input.exists():
        console.print(f"Input file not found: {input}", style="red")
        raise typer.Exit(1)

    filter_coverage_xml(
        input, output, exclude if exclude is not None else DEFAULT_EXCLUDES
    )


if __name__ == "__main__":
    app()
