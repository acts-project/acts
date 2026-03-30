#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer",
#     "rich",
# ]
# ///
"""
Utility to bump and make consistent all references to base Docker images and dependency versions.

This script can:
1. Bump Docker image tags (e.g., ubuntu2404:83 -> ubuntu2404:84)
2. Update spack-container versions in .devcontainer/Dockerfile
3. Make all references consistent across the codebase
"""

import re
import sys
import difflib
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.syntax import Syntax
from rich import box

app = typer.Typer(
    help="Bump and make consistent Docker image tags and dependency versions",
    add_completion=False,
)
console = Console()


class VersionBumper:
    """Handles version bumping for Docker images and dependencies."""

    # Patterns for different version formats
    PATTERNS = {
        # Matches: ubuntu2404:83, ubuntu2404_gnn:83, etc.
        "docker_tag": re.compile(
            r"((?:registry\.cern\.ch/)?ghcr\.io/acts-project/[a-z0-9_]+):(\d+)"
        ),
        # Matches: spack-container:18.0.0_linux-ubuntu24.04_gcc-13.3.0
        "spack_container": re.compile(
            r"(ghcr\.io/acts-project/spack-container):(\d+\.\d+\.\d+)(_[a-z0-9._-]+)"
        ),
        # Matches DEPENDENCY_TAG in two formats:
        # 1. GitLab CI: DEPENDENCY_TAG: v18.0.0
        # 2. GitHub Actions: DEPENDENCY_TAG: ... default: 'v18.0.0'
        "dependency_tag": re.compile(
            r"(DEPENDENCY_TAG:(?:\s*['\"]?|(?:[^\n]|\n(?!\s{0,2}\w))*?default:\s*['\"]))(v\d+\.\d+\.\d+)(['\"]?)",
            re.MULTILINE,
        ),
    }

    # File patterns to search
    FILE_PATTERNS = [
        ".gitlab-ci.yml",
        ".github/workflows/*.yml",
        ".github/workflows/*.yaml",
        ".github/actions/**/action.yml",
        ".devcontainer/Dockerfile",
        "CI/**/*.sh",
        "docs/**/*.md",
    ]

    def __init__(self, repo_root: Path):
        self.repo_root = repo_root

    def find_files(self) -> list[Path]:
        """Find all files that might contain version references."""
        files = set()
        for pattern in self.FILE_PATTERNS:
            files.update(self.repo_root.glob(pattern))

        # Filter out build directories and dependencies
        return [
            f
            for f in sorted(files)
            if f.is_file() and "build" not in f.parts and "_deps" not in f.parts
        ]

    def find_versions(self, content: str, pattern_name: str) -> set[str]:
        """Find all versions matching the given pattern."""
        pattern = self.PATTERNS[pattern_name]
        matches = pattern.findall(content)

        if pattern_name == "docker_tag":
            # Return unique tags (just the number part)
            return {match[1] for match in matches}
        elif pattern_name == "spack_container":
            # Return full version strings
            return {match[1] for match in matches}
        elif pattern_name == "dependency_tag":
            # Return full version strings (e.g., v18.0.0)
            return {match[1] for match in matches}

        return set()

    def scan_versions(self) -> dict:
        """Scan the repository for all current versions."""
        versions = {
            "docker_tags": set(),
            "spack_container_versions": set(),
            "dependency_tags": set(),
            "files_with_docker_tags": [],
            "files_with_spack_container": [],
            "files_with_dependency_tags": [],
        }

        files = self.find_files()
        console.print(f"[dim]Scanning {len(files)} files...[/dim]")

        for file_path in files:
            try:
                content = file_path.read_text()

                # Check for Docker tags
                docker_tags = self.find_versions(content, "docker_tag")
                if docker_tags:
                    versions["docker_tags"].update(docker_tags)
                    versions["files_with_docker_tags"].append(file_path)

                # Check for spack-container versions
                spack_versions = self.find_versions(content, "spack_container")
                if spack_versions:
                    versions["spack_container_versions"].update(spack_versions)
                    versions["files_with_spack_container"].append(file_path)

                # Check for dependency tags
                dependency_tags = self.find_versions(content, "dependency_tag")
                if dependency_tags:
                    versions["dependency_tags"].update(dependency_tags)
                    versions["files_with_dependency_tags"].append(file_path)

            except Exception as e:
                console.print(
                    f"[yellow]Warning: Error reading {file_path}: {e}[/yellow]"
                )

        return versions

    def bump_docker_tag(
        self,
        file_path: Path,
        old_tags: set[str],
        new_tag: str,
        dry_run: bool = False,
        show_diff: bool = False,
    ) -> tuple[int, str, str]:
        """Bump Docker image tags in a file. Returns (replacements, old_content, new_content)."""
        content = file_path.read_text()
        pattern = self.PATTERNS["docker_tag"]

        replacements = 0

        def replace_tag(match):
            nonlocal replacements
            old_tag = match.group(2)
            if old_tag in old_tags and old_tag != new_tag:
                replacements += 1
                return f"{match.group(1)}:{new_tag}"
            return match.group(0)

        new_content = pattern.sub(replace_tag, content)

        if replacements > 0:
            if not dry_run:
                file_path.write_text(new_content)

            prefix = "[dim][DRY RUN][/dim] " if dry_run else ""
            rel_path = file_path.relative_to(self.repo_root)
            console.print(
                f"{prefix}[green]✓[/green] Updated {replacements} occurrence(s) in [cyan]{rel_path}[/cyan]"
            )

            if show_diff:
                self._show_diff(file_path, content, new_content)

        return replacements, content, new_content

    def bump_spack_container(
        self,
        file_path: Path,
        old_versions: set[str],
        new_version: str,
        dry_run: bool = False,
        show_diff: bool = False,
    ) -> tuple[int, str, str]:
        """Bump spack-container version in a file. Returns (replacements, old_content, new_content)."""
        content = file_path.read_text()
        pattern = self.PATTERNS["spack_container"]

        replacements = 0

        def replace_version(match):
            nonlocal replacements
            old_version = match.group(2)
            if old_version in old_versions and old_version != new_version:
                replacements += 1
                return f"{match.group(1)}:{new_version}{match.group(3)}"
            return match.group(0)

        new_content = pattern.sub(replace_version, content)

        if replacements > 0:
            if not dry_run:
                file_path.write_text(new_content)

            prefix = "[dim][DRY RUN][/dim] " if dry_run else ""
            rel_path = file_path.relative_to(self.repo_root)
            console.print(
                f"{prefix}[green]✓[/green] Updated {replacements} occurrence(s) in [cyan]{rel_path}[/cyan]"
            )

            if show_diff:
                self._show_diff(file_path, content, new_content)

        return replacements, content, new_content

    def bump_dependency_tag(
        self,
        file_path: Path,
        old_versions: set[str],
        new_version: str,
        dry_run: bool = False,
        show_diff: bool = False,
    ) -> tuple[int, str, str]:
        """Bump dependency tag in a file. Returns (replacements, old_content, new_content)."""
        content = file_path.read_text()
        pattern = self.PATTERNS["dependency_tag"]

        replacements = 0

        def replace_version(match):
            nonlocal replacements
            old_version = match.group(2)
            if old_version in old_versions and old_version != new_version:
                replacements += 1
                return f"{match.group(1)}{new_version}{match.group(3)}"
            return match.group(0)

        new_content = pattern.sub(replace_version, content)

        if replacements > 0:
            if not dry_run:
                file_path.write_text(new_content)

            prefix = "[dim][DRY RUN][/dim] " if dry_run else ""
            rel_path = file_path.relative_to(self.repo_root)
            console.print(
                f"{prefix}[green]✓[/green] Updated {replacements} occurrence(s) in [cyan]{rel_path}[/cyan]"
            )

            if show_diff:
                self._show_diff(file_path, content, new_content)

        return replacements, content, new_content

    def _show_diff(self, file_path: Path, old_content: str, new_content: str):
        """Display a unified diff of the changes."""
        rel_path = file_path.relative_to(self.repo_root)
        diff = difflib.unified_diff(
            old_content.splitlines(),
            new_content.splitlines(),
            fromfile=str(rel_path),
            tofile=str(rel_path),
            lineterm="",
        )

        diff_lines = list(diff)
        if diff_lines:
            console.print()
            diff_text = "\n".join(diff_lines)
            syntax = Syntax(diff_text, "diff", theme="monokai", line_numbers=False)
            console.print(syntax)

    def bump_all_docker_tags(
        self,
        old_tags: set[str],
        new_tag: str,
        dry_run: bool = False,
        show_diff: bool = False,
    ) -> int:
        """Bump all Docker tags across the repository."""
        files = self.find_files()
        total_replacements = 0

        for file_path in files:
            try:
                replacements, _, _ = self.bump_docker_tag(
                    file_path, old_tags, new_tag, dry_run, show_diff
                )
                total_replacements += replacements
            except Exception as e:
                console.print(f"[red]Error processing {file_path}: {e}[/red]")

        return total_replacements

    def bump_all_spack_containers(
        self,
        old_versions: set[str],
        new_version: str,
        dry_run: bool = False,
        show_diff: bool = False,
    ) -> int:
        """Bump all spack-container versions across the repository."""
        files = self.find_files()
        total_replacements = 0

        for file_path in files:
            try:
                replacements, _, _ = self.bump_spack_container(
                    file_path, old_versions, new_version, dry_run, show_diff
                )
                total_replacements += replacements
            except Exception as e:
                console.print(f"[red]Error processing {file_path}: {e}[/red]")

        return total_replacements

    def bump_all_dependency_tags(
        self,
        old_versions: set[str],
        new_version: str,
        dry_run: bool = False,
        show_diff: bool = False,
    ) -> int:
        """Bump all dependency tags across the repository."""
        files = self.find_files()
        total_replacements = 0

        for file_path in files:
            try:
                replacements, _, _ = self.bump_dependency_tag(
                    file_path, old_versions, new_version, dry_run, show_diff
                )
                total_replacements += replacements
            except Exception as e:
                console.print(f"[red]Error processing {file_path}: {e}[/red]")

        return total_replacements


@app.command()
def scan(
    repo_root: Annotated[
        Path, typer.Option(help="Root directory of the repository")
    ] = Path.cwd(),
):
    """
    Scan repository for current versions.

    This command scans the codebase and displays all Docker image tags
    and spack-container versions currently in use.
    """
    bumper = VersionBumper(repo_root)
    versions = bumper.scan_versions()

    console.print()

    # Docker tags
    if versions["docker_tags"]:
        console.print("[bold]Docker image tags:[/bold]")
        for tag in sorted(versions["docker_tags"]):
            console.print(f"  [cyan]{tag}[/cyan]")
        console.print(
            f"[dim]  Found in {len(versions['files_with_docker_tags'])} file(s)[/dim]"
        )
    else:
        console.print("[bold]Docker image tags:[/bold] [yellow]none found[/yellow]")

    console.print()

    # Spack-container versions
    if versions["spack_container_versions"]:
        console.print("[bold]Spack-container versions:[/bold]")
        for version in sorted(versions["spack_container_versions"]):
            console.print(f"  [cyan]{version}[/cyan]")
        console.print(
            f"[dim]  Found in {len(versions['files_with_spack_container'])} file(s)[/dim]"
        )
    else:
        console.print(
            "[bold]Spack-container versions:[/bold] [yellow]none found[/yellow]"
        )

    console.print()

    # Dependency tags
    if versions["dependency_tags"]:
        console.print("[bold]Dependency tags:[/bold]")
        for version in sorted(versions["dependency_tags"]):
            console.print(f"  [cyan]{version}[/cyan]")
        console.print(
            f"[dim]  Found in {len(versions['files_with_dependency_tags'])} file(s)[/dim]"
        )
    else:
        console.print("[bold]Dependency tags:[/bold] [yellow]none found[/yellow]")


@app.command()
def bump_docker_tag(
    new_tag: Annotated[str, typer.Argument(help="New Docker tag to use (e.g., 84)")],
    repo_root: Annotated[
        Path, typer.Option(help="Root directory of the repository")
    ] = Path.cwd(),
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Preview changes without modifying files")
    ] = False,
    show_diff: Annotated[
        bool, typer.Option("--diff", help="Show diff of changes")
    ] = True,
):
    """
    Bump Docker image tags across the repository.

    This command finds all Docker images like 'ubuntu2404:83' and updates
    all tag numbers to the new value.

    Example:
        bump_versions.py bump-docker-tag 84
        bump_versions.py bump-docker-tag 84 --dry-run --diff
    """
    bumper = VersionBumper(repo_root)

    # Scan for existing tags
    versions = bumper.scan_versions()
    if not versions["docker_tags"]:
        console.print("[red]Error: No Docker tags found in repository[/red]")
        raise typer.Exit(1)

    # Replace ALL found tags
    old_tags = versions["docker_tags"]
    console.print(f"[dim]Found Docker tags: {', '.join(sorted(old_tags))}[/dim]")

    # Check if new tag is same as all old tags
    if len(old_tags) == 1 and new_tag in old_tags:
        console.print(f"[yellow]Tag is already {new_tag}, no changes needed[/yellow]")
        raise typer.Exit(0)

    console.print(f"[dim]Will replace ALL tags with: {new_tag}[/dim]")

    console.print()
    if dry_run:
        console.print(
            Panel(
                f"[bold]DRY RUN:[/bold] Bumping Docker tags [cyan]{', '.join(sorted(old_tags))}[/cyan] → [cyan]{new_tag}[/cyan]",
                border_style="yellow",
            )
        )
    else:
        console.print(
            Panel(
                f"Bumping Docker tags [cyan]{', '.join(sorted(old_tags))}[/cyan] → [cyan]{new_tag}[/cyan]",
                border_style="green",
            )
        )

    console.print()
    total = bumper.bump_all_docker_tags(old_tags, new_tag, dry_run, show_diff)

    console.print()
    if total > 0:
        if dry_run:
            console.print(
                f"[green]✓[/green] Would update [bold]{total}[/bold] occurrence(s)"
            )
            console.print("[dim]Run without --dry-run to apply changes[/dim]")
        else:
            console.print(
                f"[green]✓[/green] Updated [bold]{total}[/bold] occurrence(s)"
            )
    else:
        console.print(f"[yellow]No occurrences found[/yellow]")


@app.command()
def bump_spack(
    new_version: Annotated[
        str, typer.Argument(help="New spack-container version (e.g., 19.0.0)")
    ],
    repo_root: Annotated[
        Path, typer.Option(help="Root directory of the repository")
    ] = Path.cwd(),
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Preview changes without modifying files")
    ] = False,
    show_diff: Annotated[
        bool, typer.Option("--diff/--no-diff", help="Show diff of changes")
    ] = True,
):
    """
    Bump spack-container version and DEPENDENCY_TAG across the repository.

    This command updates the spack-container image version used in
    .devcontainer/Dockerfile and other configuration files, as well as
    the DEPENDENCY_TAG in GitHub Actions and GitLab CI. It replaces
    all found versions with the new version.

    Example:
        bump_versions.py bump-spack 19.0.0
        bump_versions.py bump-spack 19.0.0 --dry-run --diff
    """
    bumper = VersionBumper(repo_root)

    # Scan for existing versions
    versions = bumper.scan_versions()
    if not versions["spack_container_versions"]:
        console.print(
            "[red]Error: No spack-container versions found in repository[/red]"
        )
        raise typer.Exit(1)

    # Replace ALL found spack-container versions
    old_spack_versions = versions["spack_container_versions"]
    console.print(
        f"[dim]Found spack-container versions: {', '.join(sorted(old_spack_versions))}[/dim]"
    )

    # Also check for dependency tags
    old_dependency_versions = versions.get("dependency_tags", set())
    if old_dependency_versions:
        console.print(
            f"[dim]Found dependency tags: {', '.join(sorted(old_dependency_versions))}[/dim]"
        )

    # Check if new version is same as all old versions (both spack and dependency)
    new_dependency_version = f"v{new_version}"
    spack_is_current = (
        len(old_spack_versions) == 1 and new_version in old_spack_versions
    )
    dependency_is_current = len(old_dependency_versions) <= 1 and (
        not old_dependency_versions or new_dependency_version in old_dependency_versions
    )

    if spack_is_current and dependency_is_current:
        console.print(
            f"[yellow]All versions are already {new_version} (DEPENDENCY_TAG: {new_dependency_version}), no changes needed[/yellow]"
        )
        raise typer.Exit(0)

    console.print(f"[dim]Will replace ALL versions with: {new_version}[/dim]")

    console.print()
    if dry_run:
        console.print(
            Panel(
                f"[bold]DRY RUN:[/bold] Bumping spack-container and DEPENDENCY_TAG [cyan]{', '.join(sorted(old_spack_versions))}[/cyan] → [cyan]{new_version}[/cyan]",
                border_style="yellow",
            )
        )
    else:
        console.print(
            Panel(
                f"Bumping spack-container and DEPENDENCY_TAG [cyan]{', '.join(sorted(old_spack_versions))}[/cyan] → [cyan]{new_version}[/cyan]",
                border_style="green",
            )
        )

    console.print()

    # Bump spack-container versions
    console.print("[bold]Updating spack-container versions:[/bold]")
    total_spack = bumper.bump_all_spack_containers(
        old_spack_versions, new_version, dry_run, show_diff
    )

    # Bump dependency tags with "v" prefix
    if old_dependency_versions:
        console.print()
        console.print("[bold]Updating DEPENDENCY_TAG:[/bold]")
        new_dependency_version = f"v{new_version}"
        total_dependency = bumper.bump_all_dependency_tags(
            old_dependency_versions, new_dependency_version, dry_run, show_diff
        )
    else:
        total_dependency = 0

    console.print()
    total = total_spack + total_dependency
    if total > 0:
        if dry_run:
            console.print(
                f"[green]✓[/green] Would update [bold]{total}[/bold] occurrence(s) ([bold]{total_spack}[/bold] spack-container, [bold]{total_dependency}[/bold] DEPENDENCY_TAG)"
            )
            console.print("[dim]Run without --dry-run to apply changes[/dim]")
        else:
            console.print(
                f"[green]✓[/green] Updated [bold]{total}[/bold] occurrence(s) ([bold]{total_spack}[/bold] spack-container, [bold]{total_dependency}[/bold] DEPENDENCY_TAG)"
            )
    else:
        console.print(f"[yellow]No occurrences found[/yellow]")


if __name__ == "__main__":
    app()
