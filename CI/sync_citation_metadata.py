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
Synchronize citation metadata from CITATION.cff to .zenodo.json and AUTHORS.

This script maintains CITATION.cff as the single source of truth for citation
metadata and generates .zenodo.json and AUTHORS files from it.

Usage:
    sync_citation_metadata.py              # Update both files
    sync_citation_metadata.py --check      # Check if files are in sync
"""

import json
import sys
from pathlib import Path
from typing import Annotated

import typer
import yaml
from rich.console import Console
from rich.panel import Panel

app = typer.Typer(
    help="Synchronize citation metadata from CITATION.cff",
    add_completion=False,
)
console = Console()


def format_author_name(author: dict) -> str:
    """
    Format author name from CFF format.

    Handles:
    - Standard: {given-names} {family-names}
    - With particle: {given-names} {name-particle} {family-names}
    - Missing given-names: {family-names}
    """
    parts = []
    if "given-names" in author:
        parts.append(author["given-names"])
    if "name-particle" in author:
        parts.append(author["name-particle"])
    parts.append(author["family-names"])
    return " ".join(parts)


def extract_orcid_id(orcid: str) -> str:
    """
    Extract bare ORCID ID from URL or return as-is.

    Transforms:
    - https://orcid.org/0000-0002-2298-3605 -> 0000-0002-2298-3605
    - 0000-0002-2298-3605 -> 0000-0002-2298-3605
    """
    if not orcid.startswith("https://orcid.org/"):
        raise ValueError(f"Invalid ORCID format: {orcid}")
    return orcid.replace("https://orcid.org/", "")


def cff_author_to_zenodo_creator(author: dict) -> dict:
    """
    Convert CFF author to Zenodo creator format.

    CFF format:
        given-names: FirstName
        family-names: LastName
        affiliation: Institution
        orcid: https://orcid.org/XXXX-XXXX-XXXX-XXXX

    Zenodo format:
        name: FirstName LastName
        affiliation: Institution
        orcid: XXXX-XXXX-XXXX-XXXX
    """
    creator = {
        "affiliation": author.get("affiliation", ""),
        "name": format_author_name(author),
    }

    if "orcid" in author:
        creator["orcid"] = extract_orcid_id(author["orcid"])

    return creator


def cff_authors_to_authors_list(authors: list) -> str:
    """
    Generate AUTHORS file content from CFF authors.

    Format:
        The following people have contributed to the project (in alphabetical order):

        - FirstName LastName, Affiliation
        - ...

        Not associated with scientific/academic organisations:

        - FirstName LastName
        - ...

        See also the contributors list on github:
        https://github.com/acts-project/acts/graphs/contributors
    """
    # Separate authors with and without affiliations
    affiliated = [a for a in authors if a.get("affiliation")]
    unaffiliated = [a for a in authors if not a.get("affiliation")]

    # Sort each group by family name
    affiliated.sort(key=lambda a: a["family-names"])
    unaffiliated.sort(key=lambda a: a["family-names"])

    lines = [
        "The following people have contributed to the project (in alphabetical order):\n"
    ]

    # Add affiliated authors
    for author in affiliated:
        name = format_author_name(author)
        affiliation = author["affiliation"]
        lines.append(f"- {name}, {affiliation}")

    # Add unaffiliated section if there are any
    if unaffiliated:
        lines.append("\nNot associated with scientific/academic organisations:\n")
        for author in unaffiliated:
            name = format_author_name(author)
            lines.append(f"- {name}")

    # Add footer
    lines.append("\nSee also the contributors list on github:")
    lines.append("https://github.com/acts-project/acts/graphs/contributors")

    return "\n".join(lines) + "\n"


def generate_zenodo_json(cff_data: dict, existing_zenodo: dict) -> dict:
    """
    Generate .zenodo.json content from CITATION.cff.

    Preserves static fields from existing .zenodo.json and updates:
    - creators (from authors)
    - version
    - title (formatted as "acts-project/acts: v{version}")
    """
    # Start with existing data to preserve static fields
    zenodo_data = existing_zenodo.copy()

    # Update from CITATION.cff
    zenodo_data["creators"] = [
        cff_author_to_zenodo_creator(author) for author in cff_data["authors"]
    ]
    zenodo_data["version"] = cff_data["version"]
    zenodo_data["title"] = f"acts-project/acts: {cff_data['version']}"

    return zenodo_data


@app.command()
def generate(
    citation_file: Annotated[Path, typer.Option(help="Path to CITATION.cff")] = Path(
        "CITATION.cff"
    ),
    zenodo_file: Annotated[Path, typer.Option(help="Path to .zenodo.json")] = Path(
        ".zenodo.json"
    ),
    authors_file: Annotated[Path, typer.Option(help="Path to AUTHORS")] = Path(
        "AUTHORS"
    ),
    check: Annotated[
        bool,
        typer.Option("--check", help="Check if files are in sync"),
    ] = False,
):
    """
    Generate .zenodo.json and AUTHORS from CITATION.cff.

    In default mode, updates the files.
    In --check mode, verifies files are in sync and exits with code 1 if not.
    """
    # Read CITATION.cff
    if not citation_file.exists():
        console.print(f"[red]Error: {citation_file} not found[/red]")
        raise typer.Exit(1)

    try:
        with open(citation_file, "r", encoding="utf-8") as f:
            cff_data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        console.print(f"[red]Error parsing {citation_file}: {e}[/red]")
        raise typer.Exit(1)

    # Validate required fields
    if "authors" not in cff_data:
        console.print(f"[red]Error: 'authors' field missing in {citation_file}[/red]")
        raise typer.Exit(1)
    if "version" not in cff_data:
        console.print(f"[red]Error: 'version' field missing in {citation_file}[/red]")
        raise typer.Exit(1)

    # Read existing .zenodo.json
    existing_zenodo = {}
    if zenodo_file.exists():
        try:
            with open(zenodo_file, "r", encoding="utf-8") as f:
                existing_zenodo = json.load(f)
        except json.JSONDecodeError as e:
            console.print(f"[red]Error parsing {zenodo_file}: {e}[/red]")
            raise typer.Exit(1)
    else:
        console.print(
            f"[yellow]Warning: {zenodo_file} not found, will create new file[/yellow]"
        )

    # Generate new content
    new_zenodo_data = generate_zenodo_json(cff_data, existing_zenodo)
    new_zenodo_content = (
        json.dumps(new_zenodo_data, indent=2, ensure_ascii=False) + "\n"
    )

    new_authors_content = cff_authors_to_authors_list(cff_data["authors"])

    # Check mode: compare with existing files
    if check:
        files_out_of_sync = []

        # Check .zenodo.json
        if zenodo_file.exists():
            existing_zenodo_content = zenodo_file.read_text(encoding="utf-8")
            if existing_zenodo_content != new_zenodo_content:
                files_out_of_sync.append(str(zenodo_file))
        else:
            files_out_of_sync.append(str(zenodo_file))

        # Check AUTHORS
        if authors_file.exists():
            existing_authors_content = authors_file.read_text(encoding="utf-8")
            if existing_authors_content != new_authors_content:
                files_out_of_sync.append(str(authors_file))
        else:
            files_out_of_sync.append(str(authors_file))

        if files_out_of_sync:
            console.print()
            console.print(
                Panel(
                    "[red]Citation metadata files are out of sync![/red]\n\n"
                    + "\n".join(f"  - {f}" for f in files_out_of_sync)
                    + "\n\n"
                    + f"Run: [cyan]python {sys.argv[0]} generate[/cyan] to update them.",
                    title="Pre-commit Check Failed",
                    border_style="red",
                )
            )
            raise typer.Exit(1)
        else:
            console.print("[green]✓[/green] Citation metadata files are in sync")
            raise typer.Exit(0)

    # Write mode: update files
    console.print()
    console.print(
        Panel(
            f"Generating citation metadata from [cyan]{citation_file}[/cyan]",
            border_style="green",
        )
    )
    console.print()

    # Write .zenodo.json
    try:
        zenodo_file.write_text(new_zenodo_content, encoding="utf-8")
        console.print(f"[green]✓[/green] Updated [cyan]{zenodo_file}[/cyan]")
    except Exception as e:
        console.print(f"[red]Error writing {zenodo_file}: {e}[/red]")
        raise typer.Exit(1)

    # Write AUTHORS
    try:
        authors_file.write_text(new_authors_content, encoding="utf-8")
        console.print(f"[green]✓[/green] Updated [cyan]{authors_file}[/cyan]")
    except Exception as e:
        console.print(f"[red]Error writing {authors_file}: {e}[/red]")
        raise typer.Exit(1)

    console.print()
    console.print("[green]✓[/green] Citation metadata synchronized successfully")


if __name__ == "__main__":
    app()
