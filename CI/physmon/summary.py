#!/usr/bin/env python3
"""Render the physmon results JSON into the summary table, as markdown and HTML.

Input: physmon_results.json — a list of comparison rows, each with:
    {
      "group": "GBTS seeding",            # rows with the same group and subgroup
      "subgroup": "Single muon" | null,   # are listed under one heading
      "title": "CKF finding performance",
      "report": "html/.../foo.html",
      "info_only": false,             # true: shown, never fails the run
      "steps": [
        {"kind": "upstream"|"producer"|"plot"|"histcmp",
         "name": "trackfitting_kf",   # display label
         "exit_code": 0,
         "log": "logs/trackfitting_kf.log" | null,
         "cascade": false}            # true if earlier step in the row failed first
      ]
    }

Exit code is non-zero iff any non-info row's overall status is failing —
this is what surfaces physmon failures to snakemake / CI.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional, Tuple

# Herald serves the physmon artifact of a CI run; the summary links point there.
HERALD_URL = "https://acts-herald.app.cern.ch/view/{repo}/runs/{run_id}/artifacts/{artifact_name}/{path}"
IS_CI = "GITHUB_ACTIONS" in os.environ

# Step kinds, in display order.
STEP_KINDS = ["upstream", "producer", "plot", "histcmp"]
RESULT_HEADERS = ["", "Comparison", "Upstream", "Producer", "Plot", "Histcmp", "Report"]
CASCADE_NOTE = (
    "An asterisk (e.g. 🔴*) marks a cascade failure: "
    "an earlier step in the same row failed first."
)


# A table cell: its text and, optionally, the path it links to.
Cell = Tuple[str, Optional[str]]


def _markdown(cell: Cell) -> str:
    """The cell as markdown; on CI its link goes through Herald."""
    text, path = cell
    if not path:
        return text
    if IS_CI:
        path = HERALD_URL.format(
            repo=os.environ["GITHUB_REPOSITORY"],
            run_id=os.environ["GITHUB_RUN_ID"],
            artifact_name="physmon",
            path=path,
        )
    return f"[{text}]({path})"


def _html(cell: Cell) -> str:
    """The cell as HTML; summary.html is next to the files it links to."""
    text, path = cell
    return f'<a href="{path}">{text}</a>' if path else text


def _failed(step: dict) -> bool:
    """A step failed on any non-zero exit code, also the -1 of an unreadable witness."""
    return step["exit_code"] != 0


def _row_failed(row: dict) -> bool:
    """An info-only row never fails; any other row fails with one failed step."""
    return not row.get("info_only") and any(_failed(s) for s in row["steps"])


def _icon(rows: list[dict]) -> str:
    """Red if one of the rows failed, blue if all of them are info-only."""
    if any(_row_failed(row) for row in rows):
        return "🔴"
    return "🔵" if all(row.get("info_only") for row in rows) else "✅"


def _step_cell(step: dict | None) -> Cell:
    """The step's icon, linked to its log; a dash if the row has no such step."""
    if step is None:
        return ("—", None)
    icon = "🔴" if _failed(step) else "✅"
    return (f"{icon}*" if step.get("cascade") else icon, step.get("log"))


def _cells(row: dict) -> list[Cell]:
    steps = {s["kind"]: s for s in reversed(row["steps"])}
    return [
        (_icon([row]), None),
        (row["title"], None),
        *(_step_cell(steps.get(kind)) for kind in STEP_KINDS),
        ("report", row["report"]),
    ]


def _grouped(rows: list[dict]) -> list[tuple[str, str | None, list[dict]]]:
    """The rows by group and subgroup, each in the order of its first row."""
    groups: dict[str, dict[str | None, list[dict]]] = {}
    for row in rows:
        subgroups = groups.setdefault(row["group"], {})
        subgroups.setdefault(row.get("subgroup"), []).append(row)
    return [
        (group, subgroup, members)
        for group, subgroups in groups.items()
        for subgroup, members in subgroups.items()
    ]


def _memory_lines(memory: list[dict]) -> list[str]:
    """A warning that names the jobs which used more memory than they reserve."""
    over = [m for m in memory if m["max_pss_mb"] > m["reserved_mb"]]
    if not over:
        return []
    return [
        "",
        "## Memory",
        "",
        "These jobs used more memory than they reserve. Raise `mem_mb` in "
        "CI/physmon/config/physmon.yaml for a producer, or HISTCMP_MEM_MB / "
        "PLOT_MEM_MB in CI/physmon/Snakefile, so that the jobs running together "
        "stay within the memory budget:",
        "",
        *(
            f"- {m['job']}: {m['max_pss_mb']} MB used, {m['reserved_mb']} MB reserved"
            for m in over
        ),
    ]


def render_markdown(rows: list[dict], memory: list[dict] | None = None) -> str:
    def line(cells: list[str]) -> str:
        return "| " + " | ".join(cells) + " |"

    empty = [""] * (len(RESULT_HEADERS) - 2)
    lines = ["# physmon summary", "", line(RESULT_HEADERS)]
    lines.append(line([":---:", ":---"] + [":---:"] * len(empty)))
    for group, subgroup, members in _grouped(rows):
        name = f"**{group}**" + (f" — {subgroup}" if subgroup else "")
        lines.append(line([_icon(members), name] + empty))
        lines += [line([_markdown(cell) for cell in _cells(row)]) for row in members]
    lines += ["", f"*{CASCADE_NOTE}*".replace("🔴*", "🔴\\*")]
    lines += _memory_lines(memory or [])
    return "\n".join(lines) + "\n"


def render_html(rows: list[dict]) -> str:
    def line(cells: list[str], tag: str = "td") -> str:
        return "<tr>" + "".join(f"<{tag}>{cell}</{tag}>" for cell in cells) + "</tr>"

    empty = [""] * (len(RESULT_HEADERS) - 2)
    lines = [
        "<!DOCTYPE html>",
        '<html><head><title>physmon summary</title><meta charset="UTF-8"><style>',
        "body{font-family:sans-serif;margin:2em} table{border-collapse:collapse}",
        "th,td{border:1px solid #ccc;padding:.3em .6em;text-align:center}",
        "td:nth-child(2){text-align:left}",
        "</style></head><body><h1>physmon summary</h1><table>",
        line(RESULT_HEADERS, "th"),
    ]
    for group, subgroup, members in _grouped(rows):
        name = f"<b>{group}</b>" + (f" — {subgroup}" if subgroup else "")
        lines.append(line([_icon(members), name] + empty))
        lines += [line([_html(cell) for cell in _cells(row)]) for row in members]
    lines += ["</table>", f"<p><em>{CASCADE_NOTE}</em></p>", "</body></html>"]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "results", type=Path, help="physmon_results.json from build_results"
    )
    parser.add_argument(
        "--memory", type=Path, help="memory.json from build_results (optional)"
    )
    parser.add_argument("--md", type=Path, required=True)
    parser.add_argument("--html", type=Path, required=True)
    args = parser.parse_args()

    rows = json.loads(args.results.read_text())
    memory = json.loads(args.memory.read_text()) if args.memory else None

    args.md.write_text(render_markdown(rows, memory))
    args.html.write_text(render_html(rows))

    failed = [row for row in rows if _row_failed(row)]
    for row in failed:
        name = ", ".join(filter(None, [row["group"], row.get("subgroup")]))
        print(f"physmon: FAIL {name} | {row['title']}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
