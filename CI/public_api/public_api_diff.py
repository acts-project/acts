#!/usr/bin/env python3
"""Diff two public-API-surface snapshots and classify the change.

Consumes two JSON files from CI/public_api/public_api_surface.py. Each carries:
  * ``symbols``   -- name-level keys for non-function entities
                     (types, concepts, aliases, variables, enums)
  * ``callables`` -- source-callable function/method signatures -> return type,
                     with defaulted arguments already expanded into their
                     shorter callable forms.

Classification (source-level API only; ABI is out of scope):
  * ADDED    -- names/signatures present in head but not base
                (new type, new overload, added *defaulted* argument, ...)
  * BREAKING -- present in base but not head, or a changed return type:
                removals, renames (old name gone), removing an argument,
                adding a *non-defaulted* argument, or retyping a parameter
                (all drop the old callable form), and return-type changes.

Emits Markdown (for the PR/job summary) and a machine-readable JSON that a
labeling job can act on. Optionally fails the job (``--fail-on``).

    CI/public_api/public_api_diff.py --base base.json --head head.json \
        --json classification.json --markdown - --fail-on none
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import unicodedata
from pathlib import Path

NONFUNC_PREFIXES = ("type ", "concept ", "aliases ", "variables ", "enums ")


def load(path: str) -> dict:
    return json.loads(Path(path).read_text())


# The JSON consumed here is produced by an unprivileged job a PR fully
# controls (see api-surface-report.yml's header comment): a PR could skip
# public_api_surface.py's own extraction entirely and hand-craft this JSON,
# so a crafted C++20 Unicode identifier could carry bidi-override /
# zero-width characters into the Markdown this script renders. This is the
# one script guaranteed to run from a trusted (base-branch) checkout, so
# strip them here too rather than relying on the producer.
def _sanitize(text: str) -> str:
    return "".join(c for c in text if unicodedata.category(c) not in ("Cf", "Cc"))


def _sanitize_snapshot(d: dict) -> dict:
    return {
        "symbols": [_sanitize(s) for s in d.get("symbols", [])],
        "callables": {
            _sanitize(k): _sanitize(v) for k, v in d.get("callables", {}).items()
        },
        "fields": {_sanitize(k): _sanitize(v) for k, v in d.get("fields", {}).items()},
    }


def classify(base: dict, head: dict) -> dict:
    base = _sanitize_snapshot(base)
    head = _sanitize_snapshot(head)
    # non-function, name-level entities
    b_names = {s for s in base.get("symbols", []) if s.startswith(NONFUNC_PREFIXES)}
    h_names = {s for s in head.get("symbols", []) if s.startswith(NONFUNC_PREFIXES)}
    added_names = sorted(h_names - b_names)
    removed_names = sorted(b_names - h_names)

    # function/method callable signatures
    b_call = base.get("callables", {})
    h_call = head.get("callables", {})
    b_keys, h_keys = set(b_call), set(h_call)
    added_forms = sorted(h_keys - b_keys)
    removed_forms = sorted(b_keys - h_keys)
    ret_changed = sorted(
        f"{k}: {b_call[k]} -> {h_call[k]}"
        for k in (b_keys & h_keys)
        if b_call[k] != h_call[k]
    )

    # public data members (fields): removal / retype is breaking, add is not
    b_fld = base.get("fields", {})
    h_fld = head.get("fields", {})
    added_fields = sorted(set(h_fld) - set(b_fld))
    removed_fields = sorted(set(b_fld) - set(h_fld))
    field_retyped = sorted(
        f"{k}: {b_fld[k]} -> {h_fld[k]}"
        for k in (set(b_fld) & set(h_fld))
        if b_fld[k] != h_fld[k]
    )

    added = added_names + added_forms + added_fields
    breaking = (
        removed_names + removed_forms + ret_changed + removed_fields + field_retyped
    )
    return {
        "added": added,
        "breaking": breaking,
        "added_names": added_names,
        "added_signatures": added_forms,
        "added_fields": added_fields,
        "removed_names": removed_names,
        "removed_signatures": removed_forms,
        "return_type_changes": ret_changed,
        "removed_fields": removed_fields,
        "field_type_changes": field_retyped,
        "added_count": len(added),
        "breaking_count": len(breaking),
        "has_additions": bool(added),
        "has_breaking": bool(breaking),
    }


def _details(title: str, items: list[str], cap: int = 60) -> list[str]:
    out = [f"<details><summary>{title} ({len(items)})</summary>", ""]
    out += [f"- `{s}`" for s in items[:cap]]
    if len(items) > cap:
        out.append(f"- … and {len(items) - cap} more")
    return out + ["", "</details>"]


def render_markdown(c: dict) -> str:
    lines = ["## Public API surface diff", ""]
    if not c["has_additions"] and not c["has_breaking"]:
        return "\n".join(lines + ["No change to the public API surface. ✅", ""]) + "\n"

    lines.append(
        f"**+{c['added_count']} added**, " f"**{c['breaking_count']} breaking**."
    )
    lines.append("")
    if c["has_breaking"]:
        lines.append("### ⚠️ Breaking API changes (source-level)")
        if c["removed_names"]:
            lines += _details(
                "Removed types / aliases / enums / variables", c["removed_names"]
            )
        if c["removed_signatures"]:
            lines += _details(
                "Removed or changed call signatures", c["removed_signatures"]
            )
        if c["return_type_changes"]:
            lines += _details("Return-type changes", c["return_type_changes"])
        if c.get("removed_fields"):
            lines += _details("Removed public data members", c["removed_fields"])
        if c.get("field_type_changes"):
            lines += _details("Retyped public data members", c["field_type_changes"])
        lines.append("")
    if c["has_additions"]:
        lines.append("### ➕ Added public API")
        if c["added_names"]:
            lines += _details(
                "New types / aliases / enums / variables / concepts", c["added_names"]
            )
        if c["added_signatures"]:
            lines += _details(
                "New call signatures (incl. defaulted-arg overloads)",
                c["added_signatures"],
            )
        if c.get("added_fields"):
            lines += _details("New public data members", c["added_fields"])
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--base", required=True)
    ap.add_argument("--head", required=True)
    ap.add_argument("--json", help="write classification JSON here")
    ap.add_argument("--markdown", help="write Markdown here ('-' for stdout)")
    ap.add_argument(
        "--summary",
        action="store_true",
        help="also append the Markdown to $GITHUB_STEP_SUMMARY",
    )
    ap.add_argument(
        "--fail-on",
        choices=["none", "additions", "breaking", "any"],
        default="none",
        help="exit non-zero when this category is present",
    )
    ap.add_argument(
        "--allow-additions",
        default="false",
        help="'true' suppresses failure on additions (e.g. maintainer label present)",
    )
    args = ap.parse_args()

    c = classify(load(args.base), load(args.head))
    md = render_markdown(c)

    if args.markdown == "-" or args.markdown is None:
        print(md)
    elif args.markdown:
        Path(args.markdown).write_text(md)
    if args.json:
        Path(args.json).write_text(json.dumps(c, indent=2) + "\n")

    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if args.summary and summary:
        with open(summary, "a") as fh:
            fh.write(md)

    allow_add = str(args.allow_additions).strip().lower() in ("true", "1", "yes")
    fail = False
    if args.fail_on in ("breaking", "any") and c["has_breaking"]:
        print(
            f"::error::This PR makes {c['breaking_count']} breaking public API "
            f"change(s).",
            file=sys.stderr,
        )
        fail = True
    if args.fail_on in ("additions", "any") and c["has_additions"] and not allow_add:
        print(
            f"::error::This PR adds {c['added_count']} public API symbol(s); "
            f"a maintainer must accept the enlarged surface.",
            file=sys.stderr,
        )
        fail = True
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
