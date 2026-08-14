#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
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
Run every code generator listed in codegen/manifest.json and write the results
into a directory laid out the way cmake/ActsCodegen.cmake looks them up, i.e.
<output>/<key> for each key in the manifest.

Pointing ACTS_CODEGEN_PREBUILT_DIR at the result -- or shipping it as
`prebuilt-codegen/` next to the top-level CMakeLists.txt, which is what the
release source archive does -- lets a build skip the generators, and with them
the uv, Python and package downloads they need.

This needs uv and nothing else: no CMake, no compiler, and none of the C++
dependencies a real ACTS configure requires.
"""

import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Annotated

import typer

app = typer.Typer(add_completion=False)

SOURCE_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = SOURCE_ROOT / "codegen" / "manifest.json"


def generate(key: str, unit: dict, output_dir: Path, uv: str) -> tuple[str, str]:
    """Run one generator and return (key, error message or empty string)."""

    destination = output_dir / key
    destination.parent.mkdir(parents=True, exist_ok=True)

    command = [
        uv,
        "run",
        "--quiet",
        "--python",
        unit["python_version"],
        "--no-project",
    ]
    if unit["isolated"]:
        command.append("--isolated")
    for requirement in unit["with_requirements"]:
        command += ["--with-requirements", str(SOURCE_ROOT / requirement)]
    for package in unit["with"]:
        command += ["--with", str(SOURCE_ROOT / package)]
    command += [str(SOURCE_ROOT / unit["script"]), str(destination)]

    # The CMake path scrubs the environment before invoking uv, because a
    # configure can happen inside an LCG/Spack shell whose PYTHONPATH,
    # VIRTUAL_ENV etc. would otherwise leak into the generator. This script only
    # ever runs in CI, which starts from a clean environment to begin with, so
    # there is nothing to guard against here: just inherit it.
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return key, result.stderr.strip() or f"exited with {result.returncode}"
    if not destination.exists():
        return key, "generator succeeded but wrote no output"
    return key, ""


@app.command()
def main(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Directory to write generated code into"),
    ] = Path("prebuilt-codegen"),
    jobs: Annotated[
        int,
        typer.Option("--jobs", "-j", help="Number of generators to run in parallel"),
    ] = 0,
    clean: Annotated[
        bool,
        typer.Option("--clean/--no-clean", help="Remove the output directory first"),
    ] = True,
) -> None:
    uv = shutil.which("uv")
    if uv is None:
        typer.echo("error: uv not found on PATH", err=True)
        raise typer.Exit(1)

    units = json.loads(MANIFEST.read_text())["units"]

    if clean and output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)

    typer.echo(f"Generating {len(units)} units into {output}")

    failures = []
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=jobs or min(len(units), (os.cpu_count() or 1))
    ) as pool:
        futures = [
            pool.submit(generate, key, unit, output, uv) for key, unit in units.items()
        ]
        for future in concurrent.futures.as_completed(futures):
            key, error = future.result()
            if error:
                failures.append((key, error))
                typer.echo(f"  FAILED  {key}", err=True)
            else:
                typer.echo(f"  ok      {key}")

    if failures:
        typer.echo("", err=True)
        for key, error in failures:
            typer.echo(f"{key}:\n{error}\n", err=True)
        raise typer.Exit(1)

    typer.echo(f"Wrote {len(units)} files to {output}")


if __name__ == "__main__":
    sys.exit(app())
