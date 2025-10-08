#!/usr/bin/env python3
import argparse
import json
import pathlib
import sys


def format_one(path: pathlib.Path, indent: int, sort_keys: bool) -> int:
    """
    Returns:
      0 if unchanged,
      1 if rewritten,
      2 if invalid JSON
    """
    try:
        raw = path.read_text(encoding="utf-8")
    except Exception as e:
        print(f"[format-json] Could not read {path}: {e}", file=sys.stderr)
        return 2

    try:
        data = json.loads(raw)
    except json.JSONDecodeError as e:
        print(f"[format-json] Invalid JSON in {path}: {e}", file=sys.stderr)
        return 2

    # Use separators to avoid trailing spaces after commas/colons.
    formatted = (
        json.dumps(
            data,
            indent=indent,
            ensure_ascii=False,
            sort_keys=sort_keys,
        )
        + "\n"
    )  # ensure single trailing newline

    if formatted != raw:
        try:
            path.write_text(formatted, encoding="utf-8")
        except Exception as e:
            print(f"[format-json] Could not write {path}: {e}", file=sys.stderr)
            return 2
        print(f"[format-json] Rewrote {path}")
        return 1

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Format JSON files in-place.")
    parser.add_argument(
        "--indent", type=int, default=2, help="Indent width (default: 2)"
    )
    parser.add_argument(
        "--sort-keys",
        action="store_true",
        help="Sort object keys for stable ordering (off by default)",
    )
    parser.add_argument("files", nargs="+")

    args = parser.parse_args()

    status_max = 0
    for f in args.files:
        p = pathlib.Path(f)
        # Skip non-files (pre-commit can pass deleted paths)
        if not p.is_file():
            continue
        rc = format_one(p, indent=args.indent, sort_keys=args.sort_keys)
        status_max = max(status_max, rc)

    # Return 0 if all good (even if files were rewritten).
    # Return 2 if any invalid JSON encountered to block the commit.
    return 0 if status_max < 2 else 2


if __name__ == "__main__":
    sys.exit(main())
