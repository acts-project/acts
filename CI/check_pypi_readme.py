#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "hatch-fancy-pypi-readme",
#     "readme-renderer[md]",
# ]
# ///
"""
Keep CI/pypi_readme.md in sync with the readme that is published on PyPI.

The PyPI readme is assembled by hatch-fancy-pypi-readme from fragments of other
files (see pyproject.toml). Checking in the result makes changes to the PyPI
page visible in review.

Usage:
    check_pypi_readme.py            # Regenerate CI/pypi_readme.md
    check_pypi_readme.py --check    # Fail if CI/pypi_readme.md is out of date
"""

import argparse
import difflib
import subprocess
import sys
from pathlib import Path

import readme_renderer.markdown

REPO_ROOT = Path(__file__).resolve().parent.parent
REFERENCE = REPO_ROOT / "CI" / "pypi_readme.md"


def generate() -> str:
    result = subprocess.run(
        [sys.executable, "-m", "hatch_fancy_pypi_readme", "pyproject.toml"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    # Match end-of-file-fixer, so pre-commit hooks do not fight each other
    return result.stdout.rstrip("\n") + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    parser.add_argument(
        "--check", action="store_true", help="Only check, do not update the file"
    )
    args = parser.parse_args()

    readme = generate()

    if readme_renderer.markdown.render(readme, variant="GFM") is None:
        print("ERROR: Generated PyPI readme does not render")
        return 1

    current = REFERENCE.read_text() if REFERENCE.exists() else ""
    if readme == current:
        return 0

    rel = REFERENCE.relative_to(REPO_ROOT)
    if not args.check:
        REFERENCE.write_text(readme)
        print(f"Updated {rel}")
        return 1

    sys.stdout.writelines(
        difflib.unified_diff(
            current.splitlines(keepends=True),
            readme.splitlines(keepends=True),
            fromfile=str(rel),
            tofile="generated",
        )
    )
    print(f"\nERROR: {rel} is out of date, run CI/check_pypi_readme.py")
    return 1


if __name__ == "__main__":
    sys.exit(main())
