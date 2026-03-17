#!/usr/bin/env python3
"""Emit GitHub Actions annotations from a Doxygen warnings log."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

LOCATION_RE = re.compile(
    r"^(?P<path>.+?):(?P<line>\d+):(?:(?P<col>\d+):)?\s*(?P<message>.+)$"
)
SEVERITY_PREFIX_RE = re.compile(r"^(?P<severity>warning|error):\s*", re.IGNORECASE)


def escape_annotation(value: str) -> str:
    return value.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")


def escape_annotation_property(value: str) -> str:
    return (
        value.replace("%", "%25")
        .replace("\r", "%0D")
        .replace("\n", "%0A")
        .replace(":", "%3A")
        .replace(",", "%2C")
    )


def annotation_severity(message: str) -> str:
    match = SEVERITY_PREFIX_RE.match(message)
    if match is not None and match.group("severity").lower() == "error":
        return "error"
    return "warning"


def annotation_message(message: str) -> str:
    return escape_annotation(SEVERITY_PREFIX_RE.sub("", message, count=1))


def emit_annotation(raw_entry: str) -> bool:
    entry = raw_entry.strip()
    if not entry:
        return False

    lines = entry.splitlines()
    line = lines[0].strip()
    continuation = [continuation_line.strip() for continuation_line in lines[1:]]

    match = LOCATION_RE.match(line)
    if match is None:
        severity = annotation_severity(line)
        message = annotation_message("\n".join([line, *continuation]))
        print(f"::{severity} title=doxygen::{message}")
        return True

    severity = annotation_severity(match.group("message"))
    path = Path(match.group("path"))
    if path.is_absolute():
        try:
            path = path.resolve().relative_to(Path.cwd().resolve())
        except ValueError:
            path = path.resolve()
    line_no = match.group("line")
    col_no = match.group("col")
    message_lines = [match.group("message"), *continuation]
    message = annotation_message("\n".join(message_lines))

    location = f"file={escape_annotation_property(str(path))},line={line_no}"
    if col_no is not None:
        location += f",col={col_no}"

    print(f"::{severity} {location},title=doxygen::{message}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("logfile", type=Path)
    args = parser.parse_args()

    if not args.logfile.exists():
        print(f"No Doxygen warnings log found at {args.logfile}.", file=sys.stderr)
        return 0

    count = 0
    current_entry: list[str] = []
    with args.logfile.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue

            if LOCATION_RE.match(stripped):
                if current_entry:
                    count += int(emit_annotation("\n".join(current_entry)))
                current_entry = [stripped]
                continue

            if current_entry:
                current_entry.append(stripped)
            else:
                current_entry = [stripped]

    if current_entry:
        count += int(emit_annotation("\n".join(current_entry)))

    print(f"Emitted {count} Doxygen annotation(s).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
