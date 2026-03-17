#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "typer",
#     "rich",
#     "pydantic",
#     "pyyaml",
# ]
# ///
"""
Run cppcheck on source files and emit GitHub Actions annotations.

By default, only files changed relative to a base git ref are analysed
(PR mode).  Pass --all to analyse every file in compile_commands.json.

Sub-commands:
  analyze   Collect targets, run cppcheck, and save XML results.
  annotate  Emit GitHub Actions annotations from a results XML file.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from multiprocessing import cpu_count
from pathlib import Path, PurePosixPath
from typing import Annotated

import typer
import yaml
from pydantic import BaseModel, Field
from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from rich.text import Text

app = typer.Typer()
console = Console(stderr=True, width=None if sys.stderr.isatty() else 120)

SOURCE_SUFFIXES = {".cpp", ".cxx", ".cc", ".c"}


# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------


class FilterConfig(BaseModel):
    exclude_path_regexes: list[str] = Field(default_factory=list)
    exclude_check_regexes: list[str] = Field(default_factory=list)
    exclude_message_regexes: list[str] = Field(default_factory=list)
    severity: str = "warning"


class ParsedDiagnostic(BaseModel):
    check: str
    message: str
    path: Path
    line: int = 0
    col: int = 0
    abs_path: str = ""
    rel_path: str = ""


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------


def get_source_root() -> Path:
    # This script lives at <repo>/CI/cppcheck/run_cppcheck_pr.py
    return Path(__file__).resolve().parent.parent.parent


def load_filter_config(filter_config: Path | None) -> FilterConfig:
    if filter_config is not None and filter_config.exists():
        data = yaml.safe_load(filter_config.read_text()) or {}
        return FilterConfig.model_validate(data)
    return FilterConfig()


def load_compdb(build_dir: Path) -> list[dict]:
    compdb_path = build_dir / "compile_commands.json"
    return json.loads(compdb_path.read_text())


def filter_compdb(
    entries: list[dict], target_files: set[Path], build_dir: Path
) -> list[dict]:
    deps_dir = (build_dir / "_deps").resolve()
    result = []
    for entry in entries:
        filepath = Path(entry.get("file", ""))
        if not filepath.is_absolute():
            filepath = Path(entry.get("directory", "")) / filepath
        filepath = filepath.resolve()
        if filepath in target_files and not filepath.is_relative_to(deps_dir):
            result.append(entry)
    return result


# ---------------------------------------------------------------------------
#  Target collection
# ---------------------------------------------------------------------------


def _is_git_repo(source_root: Path) -> bool:
    return (source_root / ".git").exists()


def _get_changed_files_git(base_ref: str) -> list[str]:
    result = subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=ACM", f"{base_ref}...HEAD"],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def _get_changed_files_jj(base_ref: str) -> list[str]:
    result = subprocess.run(
        ["jj", "diff", "--from", base_ref, "--name-only"],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def get_changed_files(base_ref: str, source_root: Path) -> list[str]:
    if _is_git_repo(source_root):
        console.print("Detected git repository")
        return _get_changed_files_git(base_ref)
    console.print("No .git found, using jj")
    return _get_changed_files_jj(base_ref)


def is_path_excluded(path: str, exclude_path_regexes: list[str]) -> bool:
    return any(re.search(p, path) for p in exclude_path_regexes)


def resolve_source_targets(
    repo_paths: list[str],
    source_root: Path,
    exclude_path_regexes: list[str] | None = None,
) -> set[Path]:
    """Resolve repo-relative paths to absolute source file paths."""
    if exclude_path_regexes is None:
        exclude_path_regexes = []

    targets: set[Path] = set()
    for f in repo_paths:
        ext = PurePosixPath(f).suffix.lower()
        if ext not in SOURCE_SUFFIXES:
            continue
        abs_path = (source_root / f).resolve()
        if not abs_path.exists():
            console.print(f"  SKIP {f} (file no longer exists)")
            continue
        if is_path_excluded(str(abs_path), exclude_path_regexes):
            console.print(f"  SKIP {f} (excluded by filter)")
            continue
        targets.add(abs_path)
    return targets


# ---------------------------------------------------------------------------
#  cppcheck execution
# ---------------------------------------------------------------------------


def run_cppcheck(
    compdb_path: Path,
    output_xml: Path,
    jobs: int,
    cppcheck_bin: str,
    verbose: bool = False,
) -> None:
    cmd = [
        cppcheck_bin,
        f"--project={compdb_path}",
        "--enable=warning,performance,portability,style",
        "--std=c++20",
        "--suppress=missingIncludeSystem",
        "--suppress=unmatchedSuppression",
        "--suppress=checkersReport",
        "--xml",
        "--xml-version=2",
        f"-j{jobs}",
        "--error-exitcode=0",
    ]

    if verbose:
        console.print(f"Running: {' '.join(str(c) for c in cmd)}")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    # cppcheck writes XML to stderr when --xml is given
    xml_content = result.stderr

    output_xml.parent.mkdir(parents=True, exist_ok=True)
    output_xml.write_text(xml_content)

    if verbose and result.stdout:
        console.print(result.stdout)

    console.print(f"Wrote cppcheck results to {output_xml}")


def write_empty_results(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<results version="2"><cppcheck version=""/><errors/></results>\n'
    )
    console.print(f"Wrote empty results to {output}")


# ---------------------------------------------------------------------------
#  Parsing cppcheck XML output
# ---------------------------------------------------------------------------


def parse_cppcheck_xml(xml_path: Path) -> list[ParsedDiagnostic]:
    try:
        tree = ET.parse(xml_path)
    except ET.ParseError as e:
        console.print(f"[red]Failed to parse cppcheck XML: {e}[/]")
        return []

    root = tree.getroot()
    errors_elem = root.find("errors")
    if errors_elem is None:
        return []

    results: list[ParsedDiagnostic] = []
    for error in errors_elem.findall("error"):
        error_id = error.get("id", "unknown")
        msg = error.get("msg", "")

        location = error.find("location")
        if location is None:
            continue

        file_path = location.get("file", "")
        if not file_path:
            continue

        try:
            line = int(location.get("line", "0"))
        except ValueError:
            line = 0

        try:
            col = int(location.get("column", "0"))
        except ValueError:
            col = 0

        results.append(
            ParsedDiagnostic(
                check=f"cppcheck-{error_id}",
                message=msg,
                path=Path(file_path),
                line=line,
                col=col,
            )
        )
    return results


# ---------------------------------------------------------------------------
#  Annotation
# ---------------------------------------------------------------------------


def normalize_path(filepath: Path, source_root: Path) -> str:
    return str(filepath.relative_to(source_root, walk_up=True))


def is_excluded(diag: ParsedDiagnostic, config: FilterConfig) -> str | None:
    for pattern in config.exclude_path_regexes:
        if re.search(pattern, diag.abs_path):
            return f"path matches exclude pattern '{pattern}'"
    for pattern in config.exclude_check_regexes:
        if re.search(pattern, diag.check):
            return f"check matches exclude pattern '{pattern}'"
    for pattern in config.exclude_message_regexes:
        if re.search(pattern, diag.message):
            return f"message matches exclude pattern '{pattern}'"
    return None


def emit_annotations(
    diagnostics: list[ParsedDiagnostic],
    source_root: Path,
    config: FilterConfig,
    severity: str,
    verbose: bool = False,
    github_annotate: bool = True,
) -> int:
    severity = config.severity or severity

    for diag in diagnostics:
        diag.abs_path = str(diag.path)
        diag.rel_path = normalize_path(diag.path, source_root)

    # Deduplicate
    seen: dict[tuple, ParsedDiagnostic] = {}
    for diag in diagnostics:
        key = (diag.rel_path, diag.line, diag.col, diag.check)
        if key not in seen:
            seen[key] = diag
    unique = list(seen.values())

    remaining: list[ParsedDiagnostic] = []
    excluded_count = 0
    for d in unique:
        reason = is_excluded(d, config)
        if reason is not None:
            excluded_count += 1
            if verbose:
                console.print(
                    f"  [dim]EXCLUDE[/] {d.rel_path}:{d.line}:{d.col}"
                    f" [{d.check}] {d.message} -- {reason}"
                )
        else:
            remaining.append(d)

    console.print(
        f"After filtering: {len(remaining)} remaining, {excluded_count} excluded."
    )

    by_file: dict[str, list[ParsedDiagnostic]] = {}
    for diag in remaining:
        by_file.setdefault(diag.rel_path, []).append(diag)

    file_cache: dict[str, str] = {}

    for rel_path, file_diags in sorted(by_file.items()):
        for diag in sorted(file_diags, key=lambda d: d.line):
            if github_annotate:
                body = diag.message
                body = body.replace("%", "%25")
                body = body.replace("\r", "%0D")
                body = body.replace("\n", "%0A")
                print(
                    f"::{severity} file={diag.rel_path},line={diag.line},col={diag.col},title={diag.check}::{body}"
                )

            if diag.line > 0 and diag.abs_path:
                if diag.abs_path not in file_cache:
                    try:
                        file_cache[diag.abs_path] = Path(diag.abs_path).read_text(
                            errors="replace"
                        )
                    except OSError:
                        file_cache[diag.abs_path] = ""

                source = file_cache[diag.abs_path]
                if source:
                    ctx = 3
                    start = max(1, diag.line - ctx)
                    end = diag.line + ctx
                    syntax = Syntax(
                        source,
                        lexer="cpp",
                        line_numbers=True,
                        line_range=(start, end),
                        highlight_lines={diag.line},
                    )
                    title = Text(diag.check, style="bold red")
                    panel = Panel(
                        syntax,
                        title=title,
                        title_align="left",
                        border_style="red",
                    )
                    console.print(
                        Text.assemble((f"{rel_path}:{diag.line}:{diag.col}", "dim"))
                    )
                    console.print(
                        Panel(
                            Text(diag.message, style="yellow"),
                            border_style="dim",
                            title="message",
                            title_align="left",
                        )
                    )
                    console.print(panel)
                    continue

            console.print(
                f"[bold red]{diag.check}[/] {diag.rel_path}:{diag.line}:{diag.col}"
            )
            console.print(f"  [yellow]{diag.message}[/]")

    if remaining:
        console.print(
            f"\n[bold red]{len(remaining)} cppcheck diagnostic(s) "
            f"across {len(by_file)} file(s).[/]"
        )
        return 1
    console.print("[bold green]No cppcheck diagnostics (after filtering).[/]")
    return 0


# ---------------------------------------------------------------------------
#  Commands
# ---------------------------------------------------------------------------


@app.command()
def analyze(
    build_dir: Annotated[Path, typer.Argument(help="Build directory")],
    output_xml: Annotated[
        Path, typer.Argument(help="Output path for cppcheck XML results")
    ],
    base_ref: Annotated[
        str | None,
        typer.Option(
            "--base-ref",
            "-b",
            help="Git ref to diff against (required unless --all is given)",
        ),
    ] = None,
    all: Annotated[
        bool,
        typer.Option("--all", help="Analyse all files in compile_commands.json"),
    ] = False,
    source_root: Annotated[
        Path | None,
        typer.Option(help="Source root (default: derived from script location)"),
    ] = None,
    filter_config: Annotated[
        Path | None, typer.Option(help="Path to filter.yml")
    ] = Path(__file__)
    .resolve()
    .parent
    / "filter.yml",
    jobs: Annotated[
        int, typer.Option("-j", "--jobs", help="Parallel jobs")
    ] = cpu_count(),
    cppcheck: Annotated[str | None, typer.Option(help="cppcheck binary")] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Collect targets but skip cppcheck execution"),
    ] = False,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Print cppcheck command before running"),
    ] = False,
) -> None:
    """Collect targets, run cppcheck, and save XML results."""
    if cppcheck is None:
        cppcheck = shutil.which("cppcheck")
        if cppcheck is None:
            raise typer.BadParameter("cppcheck not found on PATH")

    if source_root is None:
        source_root = get_source_root()
    source_root = source_root.resolve()
    build_dir = build_dir.resolve()

    config = load_filter_config(filter_config)
    compdb_entries = load_compdb(build_dir)
    console.print(f"Loaded {len(compdb_entries)} entries from compile_commands.json")

    if all:
        if dry_run:
            console.print(
                f"[bold]Dry run:[/] would analyse all {len(compdb_entries)} compdb entries"
            )
            raise typer.Exit(0)
        run_cppcheck(
            build_dir / "compile_commands.json", output_xml, jobs, cppcheck, verbose
        )
    else:
        if base_ref is None:
            raise typer.BadParameter("--base-ref is required unless --all is given")

        changed = get_changed_files(base_ref, source_root)
        console.print(f"Found {len(changed)} changed files")

        target_files = resolve_source_targets(
            changed, source_root, config.exclude_path_regexes
        )
        console.print(f"Resolved {len(target_files)} source file targets")

        if not target_files:
            console.print("[bold green]No source files to analyse.[/]")
            write_empty_results(output_xml)
            raise typer.Exit(0)

        filtered = filter_compdb(compdb_entries, target_files, build_dir)
        console.print(f"Filtered to {len(filtered)} compdb entries")

        if not filtered:
            console.print("[bold green]No compdb entries match changed files.[/]")
            write_empty_results(output_xml)
            raise typer.Exit(0)

        if dry_run:
            console.print(f"[bold]Dry run:[/] would analyse {len(filtered)} file(s)")
            raise typer.Exit(0)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as tmp_f:
            json.dump(filtered, tmp_f)
            tmp_compdb = Path(tmp_f.name)

        try:
            run_cppcheck(tmp_compdb, output_xml, jobs, cppcheck, verbose)
        finally:
            tmp_compdb.unlink(missing_ok=True)

    raise typer.Exit(0)


@app.command()
def annotate(
    results: Annotated[Path, typer.Argument(help="cppcheck XML results file")],
    source_root: Annotated[
        Path | None,
        typer.Option(help="Source root (default: derived from script location)"),
    ] = None,
    filter_config: Annotated[
        Path | None, typer.Option(help="Path to filter.yml")
    ] = Path(__file__)
    .resolve()
    .parent
    / "filter.yml",
    exclude_path: Annotated[
        list[str] | None,
        typer.Option(help="Additional path regex to exclude (repeatable)"),
    ] = None,
    severity: Annotated[
        str, typer.Option(help="Annotation severity (error/warning)")
    ] = "warning",
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose",
            "-v",
            help="Print each excluded diagnostic with the reason",
        ),
    ] = False,
    github_annotate: Annotated[
        bool,
        typer.Option(
            "--github-annotate/--no-github-annotate",
            help="Emit ::error/::warning workflow commands (default: auto-detect from GITHUB_ACTIONS env var)",
        ),
    ] = os.environ.get("GITHUB_ACTIONS", "").lower()
    == "true",
) -> None:
    """Emit GitHub Actions annotations from a cppcheck XML results file."""
    if source_root is None:
        source_root = get_source_root()
    source_root = source_root.resolve()

    config = load_filter_config(filter_config)
    if exclude_path:
        config.exclude_path_regexes.extend(exclude_path)

    diagnostics = parse_cppcheck_xml(results)
    console.print(f"Found {len(diagnostics)} total diagnostic(s).")

    code = emit_annotations(
        diagnostics,
        source_root,
        config,
        severity,
        verbose=verbose,
        github_annotate=github_annotate,
    )
    raise typer.Exit(code)


if __name__ == "__main__":
    app()
