#!/usr/bin/env python3
"""Require every ``[[deprecated]]`` to carry a Doxygen ``@deprecated`` command.

Part of the API-surface effort. A C++ ``[[deprecated]]`` attribute only warns at
*compile* time; it never reaches the rendered documentation. Doxygen has a
separate ``@deprecated`` (a.k.a. ``\\deprecated``) command that collects the
entity into the generated "Deprecated List" and shows a note on its page. The
two are independent, so a symbol can be deprecated in code yet look perfectly
current in the docs -- exactly the surface our docs CI publishes.

This check pairs them mechanically: for every declaration carrying a
``[[deprecated]]`` attribute in a scanned header, the Doxygen comment block that
documents it must contain an ``@deprecated`` / ``\\deprecated`` command. It is a
deliberately dependency-free, regex/heuristic tool (no libclang, no compiler) --
it reads the same comment text Doxygen does, so "the tag Doxygen would render"
and "the tag this check sees" are the same string.

Association heuristic (no AST): from the attribute, walk upward past any
declaration-prefix lines that belong to the same statement (e.g. a
``template<...>`` line or a leading return type), stopping at a statement
boundary (blank line, ``;``, ``{``/``}``, an access specifier). The contiguous
run of comment lines found there is the entity's doc block. A violation is a
``[[deprecated]]`` whose doc block is missing or lacks an ``@deprecated``
command (including the case of no doc block at all).

    CI/public_api/check_deprecated_docs.py Core/include Plugins Fatras/include Alignment/include
    CI/public_api/check_deprecated_docs.py --write-baseline CI/public_api/deprecated_docs_baseline.txt <roots...>
    CI/public_api/check_deprecated_docs.py --baseline CI/public_api/deprecated_docs_baseline.txt <roots...>
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

DEPRECATED_ATTR_RE = re.compile(r"\[\[\s*(?:gnu\s*::\s*)?deprecated\b")
# Doxygen command, either @deprecated or \deprecated.
DEPRECATED_CMD_RE = re.compile(r"[@\\]deprecated\b")
# A line that is part of a documentation comment (/// , //! , or a /* ... */
# block body/opener). Continuation lines of a block comment start with ``*``.
COMMENT_LINE_RE = re.compile(r"^\s*(///|//!|/\*|\*)")
# A non-comment line that still belongs to the *same* declaration and may sit
# between the doc block and the attribute (template heads, attribute specifiers).
DECL_PREFIX_RE = re.compile(r"^\s*(template\b|requires\b|\[\[|explicit\b|virtual\b)")
ACCESS_SPEC_RE = re.compile(r"^\s*(public|private|protected)\s*:")

HEADER_SUFFIXES = {".hpp", ".ipp", ".cuh"}


github = "GITHUB_ACTIONS" in os.environ


def iter_header_files(root: Path):
    for dirpath, _dirs, files in os.walk(root):
        for name in files:
            fp = Path(dirpath) / name
            if fp.suffix in HEADER_SUFFIXES:
                yield fp


def doc_block_for(lines: list[str], attr_idx: int) -> list[str] | None:
    """Return the doc-comment lines documenting the declaration at ``attr_idx``.

    Walks upward from the attribute line, skipping declaration-prefix lines of
    the same statement, until it reaches the comment block or a boundary. Returns
    the comment lines (top-to-bottom) or ``None`` if there is no doc block.
    """
    i = attr_idx - 1
    # Skip decl-prefix lines that precede the attribute within the same statement
    # (e.g. a template head). Stop at any statement boundary.
    while i >= 0:
        stripped = lines[i].strip()
        if COMMENT_LINE_RE.match(lines[i]):
            break  # reached the doc block
        if stripped == "":
            return None  # blank line -> no doc block directly above
        if ACCESS_SPEC_RE.match(lines[i]):
            return None
        if stripped.endswith((";", "{", "}", "};")):
            return None  # previous statement/scope -> our decl is undocumented
        if DECL_PREFIX_RE.match(lines[i]):
            i -= 1
            continue
        # Some other declaration text above the attribute on its own line
        # (rare); treat as part of this declaration and keep looking up.
        i -= 1
    if i < 0:
        return None
    # Collect the contiguous comment block upward.
    end = i
    while i >= 0 and COMMENT_LINE_RE.match(lines[i]):
        i -= 1
    return lines[i + 1 : end + 1]


def find_violations(repo: Path, roots: list[str]) -> list[dict]:
    violations: list[dict] = []
    for root in roots:
        root_path = repo / root
        if not root_path.exists():
            continue
        for fp in iter_header_files(root_path):
            rel = fp.relative_to(repo)
            try:
                lines = fp.read_text(errors="ignore").splitlines()
            except OSError:
                continue
            for idx, line in enumerate(lines):
                if not DEPRECATED_ATTR_RE.search(line):
                    continue
                block = doc_block_for(lines, idx)
                has_doc = block is not None and len(block) > 0
                has_tag = has_doc and any(DEPRECATED_CMD_RE.search(bl) for bl in block)
                if has_tag:
                    continue
                violations.append(
                    {
                        "file": rel.as_posix(),
                        "line": idx + 1,
                        "reason": (
                            "no-doc-block" if not has_doc else "no-deprecated-tag"
                        ),
                    }
                )
    violations.sort(key=lambda v: (v["file"], v["line"]))
    return violations


def key(v: dict) -> str:
    # File-scoped, line-independent so unrelated edits above a site do not churn
    # the baseline. A file may legitimately have several deprecated members, so
    # include the line to keep them distinct within a file.
    return f"{v['file']}:{v['line']}"


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("roots", nargs="+", help="directories to scan for headers")
    p.add_argument("--repo", default=".", help="repository root")
    p.add_argument(
        "--baseline", help="ratchet file: only fail on violations not listed here"
    )
    p.add_argument(
        "--write-baseline",
        help="write current violation keys to this file and exit 0",
    )
    args = p.parse_args()

    repo = Path(args.repo).resolve()
    violations = find_violations(repo, args.roots)

    if args.write_baseline:
        Path(args.write_baseline).write_text(
            "\n".join(sorted(key(v) for v in violations)) + "\n"
        )
        print(f"wrote {len(violations)} baseline entries to {args.write_baseline}")
        return 0

    baseline: set[str] = set()
    if args.baseline and Path(args.baseline).exists():
        baseline = {
            l.strip() for l in Path(args.baseline).read_text().splitlines() if l.strip()
        }

    new = [v for v in violations if key(v) not in baseline]

    def emit(v: dict) -> None:
        loc = f"{v['file']}:{v['line']}"
        detail = (
            "has no documentation comment"
            if v["reason"] == "no-doc-block"
            else "has a doc comment but no @deprecated command"
        )
        msg = f"[[deprecated]] at {loc} {detail}"
        if github:
            print(f"::error file={v['file']},line={v['line']}::{msg}")
        else:
            print(msg)

    print(f"== [[deprecated]] without @deprecated doc: {len(violations)} total ==\n")
    for v in violations:
        emit(v)

    if args.baseline:
        print(f"\n{len(new)} not in baseline.")
        return 1 if new else 0

    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
