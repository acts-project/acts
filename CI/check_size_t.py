#!/usr/bin/env python3

from pathlib import Path
import os
import argparse
from fnmatch import fnmatch
import re
import sys

ex = re.compile(r"(\b(?<!std::)size_t)\b")

github = "GITHUB_ACTIONS" in os.environ


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument(
        "--fix", action="store_true", help="Attempt to fix any license issues found."
    )
    p.add_argument("--exclude", "-e", action="append", default=[])

    args = p.parse_args()

    # walk over all files
    exit = 0
    for root, _, files in os.walk("."):
        root = Path(root)
        for filename in files:
            # get the full path of the file
            filepath = root / filename
            if filepath.suffix not in (
                ".hpp",
                ".cpp",
                ".ipp",
                ".h",
                ".C",
                ".c",
                ".cu",
                ".cuh",
            ):
                continue

            if any([fnmatch(str(filepath), e) for e in args.exclude]):
                continue

            changed_lines = handle_file(filepath, fix=args.fix)
            if len(changed_lines) > 0:
                exit = 1
                print()
                print(filepath)
                for i, oline in changed_lines:
                    print(f"{i}: {oline}")

                    if github:
                        print(
                            f"::error file={filepath},line={i+1},title=Do not use C-style size_t::Replace size_t with std::size_t"
                        )

    return exit


def handle_file(file: Path, fix: bool) -> list[tuple[int, str]]:
    content = file.read_text()
    lines = content.splitlines()

    changed_lines = []

    for i, oline in enumerate(lines):
        line, n_subs = ex.subn(r"std::size_t", oline)
        lines[i] = line
        if n_subs > 0:
            changed_lines.append((i, oline))

    if fix and len(changed_lines) > 0:
        file.write_text("\n".join(lines) + "\n")

    return changed_lines


if "__main__" == __name__:
    sys.exit(main())
