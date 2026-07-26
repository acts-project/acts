#!/usr/bin/env python3
"""Diff two public-API-surface snapshots and optionally gate on additions.

Consumes two JSON files produced by CI/public_api_surface.py (each carries a
``symbols`` list of stable per-symbol keys) and reports which public symbols a
PR adds or removes.

Gate: adding public API grows the committed surface, so it requires explicit
maintainer sign-off. If the head snapshot adds symbols and ``--allow-additions``
is not set, this exits non-zero. A maintainer records approval by applying the
designated label to the PR (the workflow maps the label to ``--allow-additions``).

Removals are breaking changes; they are reported but not gated here (use the
existing "Breaking change" label / SemVer-major process).

    CI/public_api_diff.py --base base.json --head head.json \
        --label-name "Public API" --allow-additions false
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def load_symbols(path: str) -> set[str]:
    data = json.loads(Path(path).read_text())
    return set(data.get("symbols", []))


def fmt_list(title: str, items: list[str], cap: int = 60) -> list[str]:
    lines = [f"<details><summary>{title} ({len(items)})</summary>", ""]
    for s in items[:cap]:
        lines.append(f"- `{s}`")
    if len(items) > cap:
        lines.append(f"- … and {len(items) - cap} more")
    lines += ["", "</details>"]
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="base snapshot JSON")
    ap.add_argument("--head", required=True, help="head snapshot JSON")
    ap.add_argument("--allow-additions", default="false",
                    help="'true' to permit additions (maintainer label applied)")
    ap.add_argument("--label-name", default="Public API",
                    help="label name shown in the failure message")
    ap.add_argument("--markdown", help="write Markdown report here ('-' for stdout)")
    args = ap.parse_args()

    allow = str(args.allow_additions).strip().lower() in ("true", "1", "yes")

    base = load_symbols(args.base)
    head = load_symbols(args.head)
    added = sorted(head - base)
    removed = sorted(base - head)

    lines = ["## Public API surface diff", ""]
    if not added and not removed:
        lines.append("No change to the public API surface. ✅")
    else:
        lines.append(f"**+{len(added)}** added, **−{len(removed)}** removed "
                     f"(net {len(added) - len(removed):+d}).")
        lines.append("")
        if added:
            gate = ("✅ approved via label" if allow
                    else f"⛔ needs the **{args.label_name}** label")
            lines.append(f"Additions grow the committed public API — {gate}.")
            lines += fmt_list("Added public symbols", added)
        if removed:
            lines.append("Removals are **breaking changes** "
                         "(SemVer-major / `Breaking change` label).")
            lines += fmt_list("Removed public symbols", removed)
    md = "\n".join(lines) + "\n"

    if args.markdown == "-" or args.markdown is None:
        print(md)
    elif args.markdown:
        Path(args.markdown).write_text(md)

    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as fh:
            fh.write(md)

    if added and not allow:
        print(f"::error::This PR adds {len(added)} public API symbol(s). "
              f"A maintainer must apply the '{args.label_name}' label to accept "
              f"the enlarged public surface.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
