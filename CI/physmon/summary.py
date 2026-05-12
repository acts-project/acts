#!/usr/bin/env python3
"""Render the physmon results JSON into markdown and HTML summary tables.

Input: physmon_results.json — a list of comparison rows, each with:
    {
      "title": "...",
      "report": "html/.../foo.html",
      "info_only": false,             # true for "Comparison ..." (gx2f_vs_kf)
      "steps": [
        {"kind": "upstream"|"generator"|"plot"|"histcmp",
         "name": "trackfitting_kf",   # display label
         "rc": 0,
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

# Same Herald artifact-URL pattern the legacy summary used for CI links.
HERALD_URL = "https://acts-herald.app.cern.ch/view/{repo}/runs/{run_id}/artifacts/{artifact_name}/{path}"
IS_CI = "GITHUB_ACTIONS" in os.environ

# Step kinds, in display order.
STEP_KINDS = ["upstream", "generator", "plot", "histcmp"]
STEP_HEADERS = {
    "upstream":  "Upstream",
    "generator": "Generator",
    "plot":      "Plot",
    "histcmp":   "Histcmp",
}


def _herald(path: str) -> str:
    if not IS_CI:
        return path
    return HERALD_URL.format(
        repo=os.environ["GITHUB_REPOSITORY"],
        run_id=os.environ["GITHUB_RUN_ID"],
        artifact_name="physmon",
        path=path,
    )


def _status_for(row: dict) -> tuple[str, bool]:
    """Return (icon, failed) for the row's overall status."""
    if row.get("info_only"):
        return "🔵", False
    rc = 0
    for s in row["steps"]:
        rc |= max(s["rc"], 0)
    if rc == 0:
        return "✅", False
    return "🔴", True


def _cell_md(step: dict | None) -> str:
    """Markdown cell for one step (or em-dash if step doesn't exist for this row)."""
    if step is None:
        return "—"
    icon = "✅" if step["rc"] == 0 else "🔴"
    if step.get("cascade"):
        icon = f"{icon}*"
    log = step.get("log")
    if log:
        return f"[{icon}]({_herald(log)})"
    return icon


def _cell_html(step: dict | None) -> str:
    if step is None:
        return "—"
    icon = "✅" if step["rc"] == 0 else "🔴"
    if step.get("cascade"):
        icon = f"{icon}*"
    log = step.get("log")
    if log:
        return f'<a href="{_herald(log)}" title="{step.get("name","")}">{icon}</a>'
    return icon


def _step_by_kind(row: dict) -> dict[str, dict]:
    by_kind: dict[str, dict] = {}
    for s in row["steps"]:
        by_kind.setdefault(s["kind"], s)
    return by_kind


def _md_escape(s: str) -> str:
    """Escape '|' for use inside markdown table cells. Catalog titles use '|' as a separator."""
    return s.replace("|", "\\|")


def render_markdown(rows: list[dict]) -> str:
    lines = ["# physmon summary", ""]
    headers = ["Status", "Comparison"] + [STEP_HEADERS[k] for k in STEP_KINDS] + ["Report"]
    aligns  = [":---:", ":---"] + [":---:"] * len(STEP_KINDS) + [":---:"]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(aligns) + " |")
    for row in rows:
        icon, _ = _status_for(row)
        by_kind = _step_by_kind(row)
        cells = [icon, _md_escape(row["title"])]
        for k in STEP_KINDS:
            cells.append(_cell_md(by_kind.get(k)))
        report_url = _herald(row["report"])
        cells.append(f"[report]({report_url})")
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    lines.append("*An asterisk (e.g. 🔴\\*) marks a cascade failure: an earlier step in the same row failed first.*")
    return "\n".join(lines) + "\n"


def render_html(rows: list[dict]) -> str:
    out = [
        "<!DOCTYPE html>",
        "<html><head><title>physmon summary</title><meta charset=\"UTF-8\">",
        "<style>",
        "body{font-family:sans-serif;margin:2em}",
        "table{border-collapse:collapse}",
        "th,td{border:1px solid #ccc;padding:.3em .6em;text-align:center}",
        "td:nth-child(2){text-align:left}",
        "</style>",
        "</head><body>",
        "<h1>physmon summary</h1>",
        "<table>",
        "<tr><th>Status</th><th>Comparison</th>"
        + "".join(f"<th>{STEP_HEADERS[k]}</th>" for k in STEP_KINDS)
        + "<th>Report</th></tr>",
    ]
    for row in rows:
        icon, _ = _status_for(row)
        by_kind = _step_by_kind(row)
        cells = [icon, row["title"]] + [_cell_html(by_kind.get(k)) for k in STEP_KINDS]
        cells.append(f'<a href="{_herald(row["report"])}">report</a>')
        out.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
    out.append("</table>")
    out.append("<p><em>An asterisk (e.g. 🔴*) marks a cascade failure: an earlier step in the same row failed first.</em></p>")
    out.append("</body></html>")
    return "\n".join(out) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path, help="physmon_results.json from build_results")
    parser.add_argument("--md", type=Path, required=True)
    parser.add_argument("--html", type=Path, required=True)
    args = parser.parse_args()

    rows = json.loads(args.results.read_text())

    args.md.write_text(render_markdown(rows))
    args.html.write_text(render_html(rows))

    any_failed = False
    for row in rows:
        _, failed = _status_for(row)
        if failed:
            any_failed = True
            print(f"physmon: FAIL {row['title']}", file=sys.stderr)
    return 1 if any_failed else 0


if __name__ == "__main__":
    sys.exit(main())
