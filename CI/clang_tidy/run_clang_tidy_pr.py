#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer",
#     "rich",
#     "pyyaml",
# ]
# ///
"""
Run clang-tidy on source files and emit GitHub Actions annotations.

By default, only files changed relative to a base git ref are analysed
(PR mode).  Pass --all to analyse every file in compile_commands.json.
"""

import asyncio
import json
import os
import re
import shutil
import subprocess
import tempfile
from multiprocessing import cpu_count
from pathlib import Path, PurePosixPath
from typing import Annotated

import typer
import yaml
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeRemainingColumn

console = Console(stderr=True)

SOURCE_SUFFIXES = {".cpp", ".cxx", ".cc", ".c"}
HEADER_SUFFIXES = {".hpp", ".hxx", ".hh", ".h", ".ipp"}


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------


def get_source_root() -> Path:
    # This script lives at <repo>/CI/clang_tidy/run_clang_tidy_pr.py
    return Path(__file__).resolve().parent.parent.parent


def load_compdb(build_dir: Path) -> set[Path]:
    compdb_path = build_dir / "compile_commands.json"
    entries = json.loads(compdb_path.read_text())

    files: set[Path] = set()
    for entry in entries:
        filepath = Path(entry.get("file", ""))
        if not filepath.is_absolute():
            filepath = Path(entry.get("directory", "")) / filepath
        files.add(filepath.resolve())
    return files


# ---------------------------------------------------------------------------
#  Target collection
# ---------------------------------------------------------------------------


def get_changed_files(base_ref: str) -> list[str]:
    result = subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=ACM", f"{base_ref}...HEAD"],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def header_to_tu(header_repo_path: str, build_dir: Path) -> Path:
    """
    Map a repo-relative header path to the generated TU path.

    acts_compile_headers() generates:
      ${CMAKE_CURRENT_BINARY_DIR}/<rel_dir>/<name>.cpp
    which, for a header at <module>/include/<rest>/Foo.hpp (or .ipp),
    becomes:
      <build_dir>/<module>/include/<rest>/Foo.hpp.cpp
    """
    return (build_dir / (header_repo_path + ".cpp")).resolve()


def resolve_targets(
    repo_paths: list[str],
    build_dir: Path,
    source_root: Path,
    compdb_files: set[Path],
) -> list[Path]:
    """Classify repo-relative paths into sources and headers, resolve headers
    to their generated TUs, and return the list of clang-tidy targets that
    exist in compile_commands.json."""
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
        if abs_path in compdb_files:
            targets.add(abs_path)
        else:
            console.print(f"  SKIP source {src} (not in compile_commands.json)")

    for hdr in headers:
        tu = header_to_tu(hdr, build_dir)
        if tu in compdb_files:
            targets.add(tu)
        else:
            console.print(f"  SKIP header {hdr} (TU {tu} not in compile_commands.json)")

    console.print(f"Total targets: {len(targets)}")
    return sorted(targets)


def collect_changed_targets(
    base_ref: str, build_dir: Path, source_root: Path
) -> list[Path]:
    compdb_files = load_compdb(build_dir)
    console.print(f"Loaded {len(compdb_files)} entries from compile_commands.json")

    changed = get_changed_files(base_ref)
    console.print(f"Found {len(changed)} changed files")

    return resolve_targets(changed, build_dir, source_root, compdb_files)


def collect_all_targets(build_dir: Path) -> list[Path]:
    compdb_files = load_compdb(build_dir)
    console.print(f"Loaded {len(compdb_files)} entries from compile_commands.json")
    targets = sorted(f for f in compdb_files if f.suffix in SOURCE_SUFFIXES)
    console.print(f"Total targets: {len(targets)}")
    return targets


# ---------------------------------------------------------------------------
#  clang-tidy execution
# ---------------------------------------------------------------------------


async def run_clang_tidy_on_targets(
    targets: list[Path],
    build_dir: Path,
    fixes_dir: Path,
    jobs: int,
    clang_tidy: str,
) -> None:
    fixes_dir.mkdir(parents=True, exist_ok=True)
    sem = asyncio.Semaphore(jobs)

    progress = Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    )

    async def analyse(file: Path, idx: int, task_id: int) -> None:
        yaml_path = fixes_dir / f"{idx}.yaml"
        cmd = [
            clang_tidy,
            "-p",
            str(build_dir),
            str(file),
            "--quiet",
            f"--export-fixes={yaml_path}",
        ]
        async with sem:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.DEVNULL,
                stderr=asyncio.subprocess.DEVNULL,
            )
            await proc.wait()

        console.log(file)
        progress.advance(task_id)

    with progress:
        task_id = progress.add_task("clang-tidy", total=len(targets))
        async with asyncio.TaskGroup() as tg:
            for idx, file in enumerate(targets):
                tg.create_task(analyse(file, idx, task_id))


# ---------------------------------------------------------------------------
#  Annotation
# ---------------------------------------------------------------------------


def offset_to_line_col(text: str, offset: int) -> tuple[int, int]:
    """Convert a byte offset into (1-based line, 1-based column)."""
    line = text[:offset].count("\n") + 1
    last_nl = text.rfind("\n", 0, offset)
    col = offset - last_nl  # works even when last_nl == -1
    return line, col


def parse_fixes_yaml(path: Path) -> list[dict]:
    data = yaml.safe_load(path.read_text())
    if data is None:
        return []

    diagnostics = data.get("Diagnostics", [])
    results: list[dict] = []
    file_cache: dict[Path, str] = {}

    for diag in diagnostics:
        check = diag.get("DiagnosticName", "unknown")
        msg_block = diag.get("DiagnosticMessage", diag)
        message = msg_block.get("Message", "")
        raw_path = msg_block.get("FilePath", "")
        file_offset = msg_block.get("FileOffset", None)

        if not raw_path:
            continue

        filepath = Path(raw_path)
        line = 0
        col = 0
        if file_offset is not None and file_offset >= 0:
            if filepath not in file_cache:
                try:
                    file_cache[filepath] = filepath.read_text(errors="replace")
                except OSError:
                    file_cache[filepath] = ""
            line, col = offset_to_line_col(file_cache[filepath], file_offset)

        results.append(
            {
                "check": check,
                "message": message,
                "path": filepath,
                "line": line,
                "col": col,
            }
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


def is_excluded(diag: dict, config: dict) -> bool:
    path = diag["path"]
    check = diag["check"]
    message = diag["message"]

    for pattern in config.get("exclude_path_regexes", []):
        if re.search(pattern, path):
            return True
    for pattern in config.get("exclude_check_regexes", []):
        if re.search(pattern, check):
            return True
    for pattern in config.get("exclude_message_regexes", []):
        if re.search(pattern, message):
            return True
    return False


def annotate_fixes(
    fixes_dir: Path,
    source_root: Path,
    filter_config: Path | None,
    severity: str,
) -> int:
    config: dict = {}
    if filter_config is not None and filter_config.exists():
        config = yaml.safe_load(filter_config.read_text()) or {}

    severity = config.get("severity", severity)

    yaml_files = sorted(fixes_dir.glob("*.yaml"))
    if not yaml_files:
        console.print("No export-fixes YAML files found.")
        return 0

    console.print(f"Parsing {len(yaml_files)} YAML file(s)...")

    all_diagnostics: list[dict] = []
    for yf in yaml_files:
        all_diagnostics.extend(parse_fixes_yaml(yf))

    console.print(f"Found {len(all_diagnostics)} total diagnostic(s).")

    for diag in all_diagnostics:
        diag["path"] = normalize_path(diag["path"], source_root)

    # Deduplicate
    seen: set[tuple] = set()
    unique: list[dict] = []
    for diag in all_diagnostics:
        key = (diag["path"], diag["line"], diag["col"], diag["check"])
        if key not in seen:
            seen.add(key)
            unique.append(diag)
    all_diagnostics = unique

    remaining = [d for d in all_diagnostics if not is_excluded(d, config)]
    excluded_count = len(all_diagnostics) - len(remaining)

    console.print(
        f"After filtering: {len(remaining)} remaining, {excluded_count} excluded."
    )

    # Emit GH Actions annotations on stdout
    out = Console()  # stdout
    for diag in remaining:
        msg = diag["message"].replace("%", "%25")
        msg = msg.replace("\r", "%0D")
        msg = msg.replace("\n", "%0A")
        out.print(
            f"::{severity} file={diag['path']},line={diag['line']},col={diag['col']}"
            f"::[{diag['check']}] {msg}",
            highlight=False,
        )

    if remaining:
        console.print(f"[bold red]{len(remaining)} clang-tidy diagnostic(s) found.[/]")
        return 1
    console.print("[bold green]No clang-tidy diagnostics (after filtering).[/]")
    return 0


# ---------------------------------------------------------------------------
#  Main entry point
# ---------------------------------------------------------------------------


def main(
    build_dir: Annotated[Path, typer.Argument(help="Build directory")],
    base_ref: Annotated[
        str | None,
        typer.Argument(help="Git ref to diff against (required unless --all)"),
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
    files: Annotated[
        Path | None,
        typer.Option(help="File listing targets (one per line), bypasses git diff / --all"),
    ] = None,
    filter_config: Annotated[
        Path | None, typer.Option(help="Path to filter.yml")
    ] = None,
    jobs: Annotated[
        int, typer.Option("-j", "--jobs", help="Parallel jobs")
    ] = cpu_count(),
    clang_tidy: Annotated[str | None, typer.Option(help="clang-tidy binary")] = None,
    severity: Annotated[
        str, typer.Option(help="Annotation severity (error/warning)")
    ] = "error",
    list_targets: Annotated[
        bool,
        typer.Option("--list-targets", help="Print resolved targets and exit"),
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Collect targets but skip clang-tidy execution"),
    ] = False,
) -> None:
    if not dry_run and not list_targets:
        if clang_tidy is None:
            clang_tidy = shutil.which("clang-tidy")
            if clang_tidy is None:
                raise typer.BadParameter("clang-tidy not found on PATH")

    if source_root is None:
        source_root = get_source_root()
    source_root = source_root.resolve()
    build_dir = build_dir.resolve()

    # 1) Collect targets
    if files is not None:
        file_list = [l.strip() for l in files.read_text().splitlines() if l.strip()]
        compdb_files = load_compdb(build_dir)
        console.print(f"Loaded {len(compdb_files)} entries from compile_commands.json")
        console.print(f"Read {len(file_list)} file(s) from {files}")
        targets = resolve_targets(file_list, build_dir, source_root, compdb_files)
    elif all:
        targets = collect_all_targets(build_dir)
    else:
        if base_ref is None:
            raise typer.BadParameter(
                "base_ref is required unless --all or --files is given"
            )
        targets = collect_changed_targets(base_ref, build_dir, source_root)

    if list_targets:
        for t in targets:
            typer.echo(t)
        raise typer.Exit(0)

    if not targets:
        console.print("[bold green]No targets to analyse.[/]")
        raise typer.Exit(0)

    if dry_run:
        console.print(f"[bold]Dry run:[/] would analyse {len(targets)} file(s)")
        raise typer.Exit(0)

    # 2) Run clang-tidy
    with tempfile.TemporaryDirectory(prefix="clang-tidy-fixes-") as tmp:
        if fixes_dir is None:
            fixes_dir = Path(tmp)

        asyncio.run(
            run_clang_tidy_on_targets(targets, build_dir, fixes_dir, jobs, clang_tidy)
        )

        # 3) Annotate
        code = annotate_fixes(fixes_dir, source_root, filter_config, severity)

    raise typer.Exit(code)


if __name__ == "__main__":
    typer.run(main)
