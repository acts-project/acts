#!/usr/bin/env python
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "gcovr==8.6",
#     "rich",
#     "typer",
# ]
# ///
import multiprocessing as mp
import re
import shutil
import subprocess
from pathlib import Path

import typer
from typing import Annotated
from rich.console import Console

app = typer.Typer(add_completion=False)
console = Console()


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
def main(
    build_dir: Annotated[Path, typer.Argument(help="CMake build directory")],
    gcov: Annotated[
        str | None,
        typer.Option(None, "--gcov", help="Path to gcov executable"),
    ] = None,
) -> None:
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
    ]
    gcovr = [gcovr_exe]

    coverage_xml_path = coverage_dir / "cov.xml"
    schema_path = script_dir / "sonar_coverage.xsd"
    gcovr_sonar_cmd = (
        gcovr
        + ["-r", str(source_dir)]
        + ["--gcov-executable", gcov_exe]
        + ["-j", str(mp.cpu_count())]
        + ["--merge-mode-functions", "separate"]
        + excludes
        + extra_flags
        + ["--sonarqube", str(coverage_xml_path)]
    )
    console.print(f"$ {' '.join(gcovr_sonar_cmd)}")
    subprocess.run(gcovr_sonar_cmd, cwd=build_dir, check=True)
    validate_coverage_xml(coverage_xml_path, schema_path)

    gcovr_cmd = (
        gcovr
        + ["-r", str(source_dir)]
        + ["-j", str(mp.cpu_count())]
        + ["--gcov-executable", gcov_exe]
        + ["--merge-mode-functions", "separate"]
        + excludes
        + extra_flags
    )
    console.print(f"$ {' '.join(gcovr_cmd)}")
    subprocess.run(gcovr_cmd, cwd=build_dir, check=True)


if __name__ == "__main__":
    app()
