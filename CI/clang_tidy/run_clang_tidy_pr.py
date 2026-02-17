#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer",
#     "rich",
#     "pydantic",
#     "pyyaml",
# ]
# ///
"""
Run clang-tidy on source files and emit GitHub Actions annotations.

By default, only files changed relative to a base git ref are analysed
(PR mode).  Pass --all to analyse every file in compile_commands.json.

Sub-commands:
  analyze   Collect targets, run clang-tidy, and optionally persist fixes.
  annotate  Emit GitHub Actions annotations from a fixes YAML.
"""

import asyncio
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import cpu_count
from pathlib import Path, PurePosixPath
from typing import Annotated

import typer
import yaml
from pydantic import BaseModel, Field
from rich.console import Console, Group
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TaskID,
)
from rich.syntax import Syntax
from rich.text import Text

app = typer.Typer()
console = Console(stderr=True, width=None if sys.stderr.isatty() else 120)

SOURCE_SUFFIXES = {".cpp", ".cxx", ".cc", ".c"}
HEADER_SUFFIXES = {".hpp", ".hxx", ".hh", ".h", ".ipp"}


# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------


class FilterConfig(BaseModel):
    exclude_path_regexes: list[str] = Field(default_factory=list)
    exclude_check_regexes: list[str] = Field(default_factory=list)
    exclude_message_regexes: list[str] = Field(default_factory=list)
    severity: str = "error"


class Replacement(BaseModel):
    model_config = {"populate_by_name": True}

    file_path: str = Field(default="", alias="FilePath")
    offset: int = Field(default=0, alias="Offset")
    length: int = Field(default=0, alias="Length")
    replacement_text: str = Field(default="", alias="ReplacementText")


class DiagMessage(BaseModel):
    model_config = {"populate_by_name": True}

    message: str = Field(default="", alias="Message")
    file_path: str = Field(default="", alias="FilePath")
    file_offset: int | None = Field(default=None, alias="FileOffset")
    replacements: list[Replacement] = Field(default_factory=list, alias="Replacements")


class Diagnostic(BaseModel):
    model_config = {"populate_by_name": True}

    name: str = Field(default="unknown", alias="DiagnosticName")
    message: DiagMessage = Field(default_factory=DiagMessage, alias="DiagnosticMessage")
    analyzed_files: list[str] = Field(default_factory=list, alias="AnalyzedFiles")


class FixesFile(BaseModel):
    model_config = {"populate_by_name": True}

    diagnostics: list[Diagnostic] = Field(default_factory=list, alias="Diagnostics")
    analyzed_file: str = Field(default="", alias="AnalyzedFile")


class ParsedDiagnostic(BaseModel):
    check: str
    message: str
    path: Path
    line: int = 0
    col: int = 0
    abs_path: str = ""
    rel_path: str = ""
    suggestion: str = ""
    analyzed_files: list[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------


def get_source_root() -> Path:
    # This script lives at <repo>/CI/clang_tidy/run_clang_tidy_pr.py
    return Path(__file__).resolve().parent.parent.parent


def load_filter_config(filter_config: Path | None) -> FilterConfig:
    if filter_config is not None and filter_config.exists():
        data = yaml.safe_load(filter_config.read_text()) or {}
        return FilterConfig.model_validate(data)
    return FilterConfig()


def load_compdb(build_dir: Path) -> set[Path]:
    compdb_path = build_dir / "compile_commands.json"
    entries = json.loads(compdb_path.read_text())

    deps_dir = (build_dir / "_deps").resolve()
    files: set[Path] = set()
    for entry in entries:
        filepath = Path(entry.get("file", ""))
        if not filepath.is_absolute():
            filepath = Path(entry.get("directory", "")) / filepath
        filepath = filepath.resolve()
        if not filepath.is_relative_to(deps_dir):
            files.add(filepath)
    return files


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


def find_header_tu(
    header_abs_path: Path, source_root: Path, compdb_files: set[Path]
) -> Path | None:
    """Find the generated TU for a header by suffix-matching against the compdb.

    acts_compile_headers() generates TUs like:
      <build_dir>/<rel_path>/Foo.hpp.cpp
    We compute the header's path relative to the repo root and look for a
    compdb entry ending with that relative path + ``.cpp``.
    """
    try:
        rel = header_abs_path.relative_to(source_root)
    except ValueError:
        rel = Path(os.path.relpath(header_abs_path, source_root))
    suffix = "/" + str(rel).replace(os.sep, "/") + ".cpp"
    for f in compdb_files:
        if str(f).replace(os.sep, "/").endswith(suffix):
            return f
    return None


def is_path_excluded(path: str, exclude_path_regexes: list[str]) -> bool:
    return any(re.search(p, path) for p in exclude_path_regexes)


def resolve_targets(
    repo_paths: list[str],
    source_root: Path,
    compdb_files: set[Path],
    exclude_path_regexes: list[str] | None = None,
) -> list[Path]:
    """Classify repo-relative paths into sources and headers, resolve headers
    to their generated TUs, and return the list of clang-tidy targets that
    exist in compile_commands.json."""
    if exclude_path_regexes is None:
        exclude_path_regexes = []

    sources: list[str] = []
    headers: list[str] = []
    for f in repo_paths:
        ext = PurePosixPath(f).suffix.lower()
        if ext in SOURCE_SUFFIXES:
            sources.append(f)
        elif ext in HEADER_SUFFIXES:
            headers.append(f)

    console.print(f"  {len(sources)} source files, {len(headers)} headers")

    targets: set[Path] = set()

    for src in sources:
        abs_path = (source_root / src).resolve()
        if not abs_path.exists():
            console.print(f"  SKIP source {src} (file no longer exists)")
            continue
        if is_path_excluded(str(abs_path), exclude_path_regexes):
            console.print(f"  SKIP source {src} (excluded by filter)")
            continue
        if abs_path in compdb_files:
            targets.add(abs_path)
        else:
            console.print(f"  SKIP source {src} (not in compile_commands.json)")

    for hdr in headers:
        abs_path = (source_root / hdr).resolve()
        if not abs_path.exists():
            console.print(f"  SKIP header {hdr} (file no longer exists)")
            continue
        if is_path_excluded(str(abs_path), exclude_path_regexes):
            console.print(f"  SKIP header {hdr} (excluded by filter)")
            continue
        tu = find_header_tu(abs_path, source_root, compdb_files)
        if tu is not None:
            targets.add(tu)
        else:
            console.print(f"  SKIP header {hdr} (no TU in compile_commands.json)")

    console.print(f"Total targets: {len(targets)}")
    return sorted(targets)


def collect_changed_targets(
    base_ref: str,
    source_root: Path,
    compdb_files: set[Path],
    exclude_path_regexes: list[str] | None = None,
) -> list[Path]:
    changed = get_changed_files(base_ref, source_root)
    console.print(f"Found {len(changed)} changed files")

    return resolve_targets(
        changed, source_root, compdb_files, exclude_path_regexes
    )


def collect_targets_from_fixes(
    fixes_file: Path,
    source_root: Path,
    compdb_files: set[Path],
    exclude_path_regexes: list[str] | None = None,
) -> list[Path]:
    """Extract file paths from a previous fixes YAML and resolve them as targets."""
    fixes = FixesFile.model_validate(yaml.safe_load(fixes_file.read_text()) or {})
    paths: set[str] = set()
    for diag in fixes.diagnostics:
        fp = diag.message.file_path
        if fp:
            paths.add(fp)

    console.print(f"Found {len(paths)} unique file(s) in {fixes_file}")

    # Normalize to repo-relative paths for resolve_targets
    repo_paths: list[str] = []
    for p in paths:
        try:
            rel = os.path.relpath(p, source_root)
        except ValueError:
            rel = p
        repo_paths.append(rel)

    return resolve_targets(
        repo_paths, source_root, compdb_files, exclude_path_regexes
    )


def collect_all_targets(
    compdb_files: set[Path],
    exclude_path_regexes: list[str] | None = None,
) -> list[Path]:
    if exclude_path_regexes is None:
        exclude_path_regexes = []
    targets = sorted(
        f
        for f in compdb_files
        if f.suffix in SOURCE_SUFFIXES
        and not is_path_excluded(str(f), exclude_path_regexes)
    )
    console.print(f"Total targets: {len(targets)}")
    return targets


def _macos_sysroot() -> str | None:
    """Return the macOS SDK path, or None if not on macOS / xcrun fails."""
    if sys.platform != "darwin":
        return None
    try:
        result = subprocess.run(
            ["xcrun", "--show-sdk-path"],
            capture_output=True,
            text=True,
            check=True,
        )
        sdk = result.stdout.strip()
        if sdk:
            console.print(f"Using macOS sysroot: {sdk}")
            return sdk
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return None


# ---------------------------------------------------------------------------
#  clang-tidy execution
# ---------------------------------------------------------------------------


async def run_clang_tidy_on_targets(
    targets: list[Path],
    build_dir: Path,
    fixes_dir: Path,
    jobs: int,
    clang_tidy: str,
    verbose: bool = False,
    trace_includes: bool = False,
) -> None:
    fixes_dir.mkdir(parents=True, exist_ok=True)
    sem = asyncio.Semaphore(jobs)

    sysroot = _macos_sysroot()
    extra_args: list[str] = []
    if sysroot:
        extra_args += [f"--extra-arg=-isysroot", f"--extra-arg={sysroot}"]
    if trace_includes:
        extra_args.append("--extra-arg=-H")

    progress = Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    )

    async def analyse(file: Path, idx: int, task_id: TaskID) -> None:
        yaml_path = fixes_dir / f"{idx}.yaml"
        cmd = [
            clang_tidy,
            "-p",
            str(build_dir),
            str(file),
            "--quiet",
            f"--export-fixes={yaml_path}",
            "-header-filter=.*",  # export fixes for headers, even when analyzing a TU
            *extra_args,
        ]
        async with sem:
            if verbose or trace_includes:
                proc = await asyncio.create_subprocess_exec(
                    *cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.STDOUT,
                )
                stdout, _ = await proc.communicate()
                output = stdout.decode(errors="replace").strip()
                if output:
                    console.print(f"\n[bold]{file}[/]")
                    console.print(output)
            else:
                proc = await asyncio.create_subprocess_exec(
                    *cmd,
                    stdout=asyncio.subprocess.DEVNULL,
                    stderr=asyncio.subprocess.DEVNULL,
                )
                await proc.wait()

        # Inject the analyzed file into the YAML so we can trace diagnostics
        # back to the TU that produced them.
        if yaml_path.exists():
            data = yaml.safe_load(yaml_path.read_text()) or {}
            data["AnalyzedFile"] = str(file)
            yaml_path.write_text(
                yaml.dump(data, default_flow_style=False, sort_keys=False)
            )

        console.log(file)
        progress.advance(task_id)

    with progress:
        task_id = progress.add_task("clang-tidy", total=len(targets))
        async with asyncio.TaskGroup() as tg:
            for idx, file in enumerate(targets):
                tg.create_task(analyse(file, idx, task_id))


# ---------------------------------------------------------------------------
#  Fixes merging
# ---------------------------------------------------------------------------


def write_empty_fixes(output: Path) -> None:
    """Write an empty fixes YAML so downstream steps always find the file."""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.dump({"Diagnostics": []}, default_flow_style=False))
    console.print(f"Wrote empty fixes to {output}")


def merge_fixes_yaml(fixes_dir: Path, output: Path) -> None:
    """Merge and deduplicate all per-TU export-fixes YAML files into a single file.

    When the same diagnostic appears from multiple TUs, the ``AnalyzedFiles``
    list on the merged diagnostic collects all originating files.
    """
    yaml_files = sorted(fixes_dir.glob("*.yaml"))

    # Collect diagnostics, tagging each with its originating analyzed file.
    all_diagnostics: list[Diagnostic] = []
    for yf in yaml_files:
        fixes = FixesFile.model_validate(yaml.safe_load(yf.read_text()) or {})
        analyzed_file = fixes.analyzed_file
        for d in fixes.diagnostics:
            if analyzed_file:
                d.analyzed_files = [analyzed_file]
            all_diagnostics.append(d)

    # Deduplicate by (file, offset, check), merging analyzed_files lists.
    seen: dict[tuple, Diagnostic] = {}
    for d in all_diagnostics:
        key = (d.message.file_path, d.message.file_offset, d.name)
        if key in seen:
            existing = seen[key]
            for af in d.analyzed_files:
                if af not in existing.analyzed_files:
                    existing.analyzed_files.append(af)
        else:
            seen[key] = d

    unique = list(seen.values())
    deduped = len(all_diagnostics) - len(unique)
    if deduped:
        console.print(f"Deduplicated {deduped} diagnostic(s).")

    merged = {"Diagnostics": [d.model_dump(by_alias=True) for d in unique]}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.dump(merged, default_flow_style=False, sort_keys=False))
    console.print(f"Wrote {len(unique)} diagnostic(s) to {output}")


# ---------------------------------------------------------------------------
#  Annotation
# ---------------------------------------------------------------------------


def offset_to_line_col(text: str, offset: int) -> tuple[int, int]:
    """Convert a byte offset into (1-based line, 1-based column)."""
    line = text[:offset].count("\n") + 1
    last_nl = text.rfind("\n", 0, offset)
    col = offset - last_nl  # works even when last_nl == -1
    return line, col


def build_suggestion(
    source: str,
    diag_file: str,
    replacements: list[Replacement],
) -> str:
    """Build a clang-tidy-style suggestion from replacements.

    Produces output like::

        42 |   ContextType(const T& value) : m_data{value} {}
           |   ^
           |   explicit

    Only considers replacements that target the same file as the diagnostic.
    Returns an empty string when no applicable replacements exist.
    """
    applicable = [r for r in replacements if r.file_path == diag_file]
    if not applicable or not source:
        return ""

    lines = source.splitlines()
    parts: list[str] = []

    for repl in sorted(applicable, key=lambda r: r.offset):
        line_no, col = offset_to_line_col(source, repl.offset)
        if line_no < 1 or line_no > len(lines):
            continue

        source_line = lines[line_no - 1]
        line_num_width = len(str(line_no))
        gutter = " " * line_num_width

        parts.append(f"{line_no} | {source_line}")

        # col is 1-based; we need (col - 1) spaces to reach the caret
        indent = " " * (col - 1)
        if repl.length > 1:
            # Show ^~~~ spanning the replaced range
            span = "^" + "~" * (repl.length - 1)
        else:
            span = "^"
        parts.append(f"{gutter} | {indent}{span}")

        # Show the replacement text (may be multi-line)
        fix_lines = repl.replacement_text.splitlines() or [""]
        for fl in fix_lines:
            parts.append(f"{gutter} | {indent}{fl}")

    return "\n".join(parts)


def parse_fixes_yaml(path: Path) -> list[ParsedDiagnostic]:
    fixes = FixesFile.model_validate(yaml.safe_load(path.read_text()) or {})

    results: list[ParsedDiagnostic] = []
    file_cache: dict[Path, str] = {}

    for diag in fixes.diagnostics:
        msg = diag.message
        if not msg.file_path:
            continue

        filepath = Path(msg.file_path)
        line = 0
        col = 0
        if msg.file_offset is not None and msg.file_offset >= 0:
            if filepath not in file_cache:
                try:
                    file_cache[filepath] = filepath.read_text(errors="replace")
                except OSError:
                    file_cache[filepath] = ""
            source = file_cache[filepath]
            line, col = offset_to_line_col(source, msg.file_offset)
        else:
            source = ""

        suggestion = build_suggestion(source, msg.file_path, msg.replacements)

        results.append(
            ParsedDiagnostic(
                check=diag.name,
                message=msg.message,
                path=filepath,
                line=line,
                col=col,
                suggestion=suggestion,
                analyzed_files=list(diag.analyzed_files),
            )
        )
    return results


def normalize_path(filepath: Path, source_root: Path) -> str:
    """Make an absolute path repo-relative.  Falls back to the original for
    paths outside source_root (e.g. system headers).

    Uses os.path.relpath because pathlib's relative_to does not support
    paths outside the root."""
    try:
        return os.path.relpath(filepath, source_root).replace(os.sep, "/")
    except ValueError:
        # On Windows, relpath raises ValueError for paths on different drives.
        return str(filepath).replace(os.sep, "/")


def is_excluded(diag: ParsedDiagnostic, config: FilterConfig) -> str | None:
    """Check whether a diagnostic should be excluded.

    Returns ``None`` if the diagnostic passes all filters, or a human-readable
    reason string describing why it was excluded.
    """
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
    """Normalize paths, deduplicate, filter, and emit GH Actions annotations.
    Returns 1 if any diagnostics remain, 0 otherwise."""
    severity = config.severity or severity

    for diag in diagnostics:
        diag.abs_path = str(diag.path)
        diag.rel_path = normalize_path(diag.path, source_root)

    # Deduplicate, merging analyzed_files lists
    seen_map: dict[tuple, ParsedDiagnostic] = {}
    for diag in diagnostics:
        key = (diag.rel_path, diag.line, diag.col, diag.check)
        if key in seen_map:
            existing = seen_map[key]
            for af in diag.analyzed_files:
                if af not in existing.analyzed_files:
                    existing.analyzed_files.append(af)
            if verbose:
                console.print(
                    f"  [dim]DEDUP[/] {diag.rel_path}:{diag.line}:{diag.col}"
                    f" [{diag.check}] {diag.message}"
                )
        else:
            seen_map[key] = diag
    unique = list(seen_map.values())

    if verbose and len(diagnostics) != len(unique):
        console.print(f"Deduplicated {len(diagnostics) - len(unique)} diagnostic(s).")
    diagnostics = unique

    remaining: list[ParsedDiagnostic] = []
    excluded_count = 0
    for d in diagnostics:
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

    # Group by file for rich output
    by_file: dict[str, list[ParsedDiagnostic]] = {}
    for diag in remaining:
        by_file.setdefault(diag.rel_path, []).append(diag)

    file_cache: dict[str, str] = {}

    for rel_path, file_diags in sorted(by_file.items()):
        for diag in sorted(file_diags, key=lambda d: d.line):
            # GH Actions annotation on stdout
            if github_annotate:
                body = diag.message
                if diag.suggestion:
                    body += "\n" + diag.suggestion
                body = body.replace("%", "%25")
                body = body.replace("\r", "%0D")
                body = body.replace("\n", "%0A")
                print(
                    f"::{severity} file={diag.rel_path},line={diag.line},col={diag.col},title={diag.check}::{body}",
                )

            # Rich panel on stderr for the CI log
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
                    renderables: list = [
                        Text.assemble(
                            (f"{rel_path}:{diag.line}:{diag.col}", "dim"),
                        ),
                        syntax,
                    ]
                    renderables.append(
                        Panel(
                            Text(diag.message, style="yellow"),
                            border_style="dim",
                            title="message",
                            title_align="left",
                        )
                    )
                    if diag.suggestion:
                        renderables.append(
                            Panel(
                                Syntax(
                                    diag.suggestion,
                                    lexer="text",
                                    line_numbers=False,
                                    theme="ansi_dark",
                                ),
                                border_style="green",
                                title="suggestion",
                                title_align="left",
                            )
                        )
                    if diag.analyzed_files:
                        af_lines = "\n".join(diag.analyzed_files)
                        renderables.append(
                            Panel(
                                Text(af_lines, style="cyan"),
                                border_style="dim",
                                title=f"analyzed from ({len(diag.analyzed_files)})",
                                title_align="left",
                            )
                        )
                    panel = Panel(
                        Group(*renderables),
                        title=title,
                        title_align="left",
                        border_style="red",
                    )
                    console.print(panel)
                    continue

            # Fallback: no source available
            console.print(
                f"[bold red]{diag.check}[/] {diag.rel_path}:{diag.line}:{diag.col}"
            )
            console.print(f"  [yellow]{diag.message}[/]")
            if diag.analyzed_files:
                console.print(
                    f"  [dim]analyzed from:[/] [cyan]{', '.join(diag.analyzed_files)}[/]"
                )

    if remaining:
        console.print(
            f"\n[bold red]{len(remaining)} clang-tidy diagnostic(s) "
            f"across {len(by_file)} file(s).[/]"
        )
        return 1
    console.print("[bold green]No clang-tidy diagnostics (after filtering).[/]")
    return 0


def annotate_from_fixes_dir(
    fixes_dir: Path,
    source_root: Path,
    config: FilterConfig,
    severity: str,
    verbose: bool = False,
    github_annotate: bool = True,
) -> int:
    yaml_files = sorted(fixes_dir.glob("*.yaml"))
    if not yaml_files:
        console.print("No export-fixes YAML files found.")
        return 0

    console.print(f"Parsing {len(yaml_files)} YAML file(s)...")

    all_diagnostics: list[ParsedDiagnostic] = []
    for yf in yaml_files:
        all_diagnostics.extend(parse_fixes_yaml(yf))

    console.print(f"Found {len(all_diagnostics)} total diagnostic(s).")
    return emit_annotations(
        all_diagnostics,
        source_root,
        config,
        severity,
        verbose=verbose,
        github_annotate=github_annotate,
    )


def annotate_from_fixes_file(
    fixes_file: Path,
    source_root: Path,
    config: FilterConfig,
    severity: str,
    verbose: bool = False,
    github_annotate: bool = True,
) -> int:
    console.print(f"Parsing {fixes_file}...")
    all_diagnostics = parse_fixes_yaml(fixes_file)
    console.print(f"Found {len(all_diagnostics)} total diagnostic(s).")
    return emit_annotations(
        all_diagnostics,
        source_root,
        config,
        severity,
        verbose=verbose,
        github_annotate=github_annotate,
    )


# ---------------------------------------------------------------------------
#  Commands
# ---------------------------------------------------------------------------


@app.command()
def analyze(
    build_dir: Annotated[Path, typer.Argument(help="Build directory")],
    output_fixes: Annotated[
        Path, typer.Argument(help="Output path for merged fixes YAML")
    ],
    files: Annotated[
        list[Path] | None,
        typer.Argument(
            help="Explicit file paths to analyse (bypasses --base-ref / --all)"
        ),
    ] = None,
    base_ref: Annotated[
        str | None,
        typer.Option(
            "--base-ref",
            "-b",
            help="Git ref to diff against (required unless --all or files are given)",
        ),
    ] = None,
    all: Annotated[
        bool,
        typer.Option("--all", help="Analyse all files instead of only changed ones"),
    ] = False,
    source_root: Annotated[
        Path | None,
        typer.Option(help="Source root (default: git toplevel)"),
    ] = None,
    fixes_dir: Annotated[
        Path | None,
        typer.Option(help="Directory for export-fixes YAML (default: temp dir)"),
    ] = None,
    from_fixes: Annotated[
        Path | None,
        typer.Option(help="Re-analyse files from a previous fixes YAML"),
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
    clang_tidy: Annotated[str | None, typer.Option(help="clang-tidy binary")] = None,
    list_targets: Annotated[
        bool,
        typer.Option("--list-targets", help="Print resolved targets and exit"),
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Collect targets but skip clang-tidy execution"),
    ] = False,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Print clang-tidy output for each file"),
    ] = False,
    trace_includes: Annotated[
        bool,
        typer.Option(
            "--trace-includes", help="Pass -H to the compiler to dump the include tree"
        ),
    ] = False,
) -> None:
    """Collect targets, run clang-tidy, and optionally persist merged fixes."""
    if not dry_run and not list_targets:
        if clang_tidy is None:
            clang_tidy = shutil.which("clang-tidy")
            if clang_tidy is None:
                raise typer.BadParameter("clang-tidy not found on PATH")
    assert clang_tidy is not None

    if source_root is None:
        source_root = get_source_root()
    source_root = source_root.resolve()
    build_dir = build_dir.resolve()

    # Load filter config early so path excludes apply during target selection
    config = load_filter_config(filter_config)

    compdb_files = load_compdb(build_dir)
    console.print(f"Loaded {len(compdb_files)} entries from compile_commands.json")

    # 1) Collect targets
    if from_fixes is not None:
        targets = collect_targets_from_fixes(
            from_fixes,
            source_root,
            compdb_files,
            config.exclude_path_regexes,
        )
    elif files:
        file_list = [str(f) for f in files]
        console.print(f"{len(file_list)} file(s) specified as arguments")
        targets = resolve_targets(
            file_list,
            source_root,
            compdb_files,
            config.exclude_path_regexes,
        )
    elif all:
        targets = collect_all_targets(compdb_files, config.exclude_path_regexes)
    else:
        if base_ref is None:
            raise typer.BadParameter(
                "--base-ref is required unless --all, --from-fixes, or files are given"
            )
        targets = collect_changed_targets(
            base_ref, source_root, compdb_files, config.exclude_path_regexes
        )

    if list_targets:
        for t in targets:
            typer.echo(t)
        raise typer.Exit(0)

    if not targets:
        console.print("[bold green]No targets to analyse.[/]")
        write_empty_fixes(output_fixes)
        raise typer.Exit(0)

    if dry_run:
        console.print(f"[bold]Dry run:[/] would analyse {len(targets)} file(s)")
        raise typer.Exit(0)

    # 2) Run clang-tidy
    with tempfile.TemporaryDirectory(prefix="clang-tidy-fixes-") as tmp:
        if fixes_dir is None:
            fixes_dir = Path(tmp)

        asyncio.run(
            run_clang_tidy_on_targets(
                targets,
                build_dir,
                fixes_dir,
                jobs,
                clang_tidy,
                verbose=verbose,
                trace_includes=trace_includes,
            )
        )

        # 3) Persist merged fixes
        merge_fixes_yaml(fixes_dir, output_fixes)

    raise typer.Exit(0)


@app.command()
def annotate(
    fixes: Annotated[Path, typer.Argument(help="Merged fixes YAML file")],
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
    ] = "error",
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose",
            "-v",
            help="Print each excluded and deduplicated diagnostic with the reason",
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
    """Emit GitHub Actions annotations from a previously saved fixes YAML."""
    if source_root is None:
        source_root = get_source_root()
    source_root = source_root.resolve()

    config = load_filter_config(filter_config)
    if exclude_path:
        config.exclude_path_regexes.extend(exclude_path)

    code = annotate_from_fixes_file(
        fixes,
        source_root,
        config,
        severity,
        verbose=verbose,
        github_annotate=github_annotate,
    )
    raise typer.Exit(code)


if __name__ == "__main__":
    app()
