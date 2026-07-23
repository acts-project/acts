#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "packaging",
# ]
# ///
# Prints information about the Python versions Acts supports, derived from
# `requires-python` in the top-level pyproject.toml so that the wheel metadata,
# cibuildwheel target matrix and the floor used to compile the requirements.txt
# lockfiles all share one definition.
#
# Usage:
#   CI/supported_python_versions.py --floor    # lowest version, e.g. 3.10
#   CI/supported_python_versions.py --cibw     # CIBW_BUILD value, e.g. "cp311-* cp312-*"

import argparse
import sys
import tomllib
from pathlib import Path

from packaging.specifiers import SpecifierSet


def main() -> None:
    pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
    with pyproject.open("rb") as f:
        raw = tomllib.load(f)["project"]["requires-python"]
    spec = SpecifierSet(raw)

    parser = argparse.ArgumentParser()
    parser.add_argument("--floor", action="store_true")
    parser.add_argument("--cibw", action="store_true")
    args = parser.parse_args()

    all_versions = [f"3.{i}" for i in range(100)]
    versions = list(spec.filter(all_versions))

    if args.floor:
        print(versions[0])
        return

    if args.cibw:
        tags = ["cp" + v.replace(".", "") + "-*" for v in versions]
        print(" ".join(tags))
        return

    parser.print_help()
    sys.exit(1)


if __name__ == "__main__":
    main()
