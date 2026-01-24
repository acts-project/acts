#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "typer",
#   "httpx",
#   "rich",
# ]
# ///

import asyncio
import concurrent.futures
import hashlib
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Annotated

import httpx
import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TaskID,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

console = Console()
app = typer.Typer()


def find_geant4_config() -> Path:
    """Find geant4-config in PATH."""
    result = shutil.which("geant4-config")
    if not result:
        console.print("[red]Error: geant4-config not found in PATH[/red]")
        raise typer.Exit(1)
    return Path(result)


def parse_datasets(config_path: Path) -> list[dict]:
    """Parse dataset information from geant4-config script."""
    with open(config_path) as f:
        content = f.read()

    # Find the dataset_list line
    for line in content.splitlines():
        if line.strip().startswith("dataset_list="):
            break
    else:
        console.print("[red]Error: Could not find dataset_list in geant4-config[/red]")
        raise typer.Exit(1)

    # Extract the awk script dataset string
    # Format: NAME|ENVVAR|PATH|FILENAME|MD5;...
    # The string is within quotes before the ", array," part
    start = line.find('"') + 1
    # Find the end quote before ", array,"
    end = line.find('", array,')
    if end == -1:
        end = line.rfind('"')
    dataset_string = line[start:end]

    datasets = []
    for entry in dataset_string.split(";"):
        if not entry.strip():
            continue
        parts = entry.split("|")
        if len(parts) >= 5:
            datasets.append(
                {
                    "name": parts[0],
                    "envvar": parts[1],
                    "path": parts[2],
                    "filename": parts[3],
                    "md5": parts[4],
                }
            )

    return datasets


def get_dataset_url(config_path: Path) -> str:
    """Get the base URL for datasets from geant4-config."""
    with open(config_path) as f:
        content = f.read()

    # Find the dataset_url line
    for line in content.splitlines():
        if line.strip().startswith("dataset_url="):
            # Extract URL from line like: dataset_url="https://cern.ch/geant4-data/datasets"
            start = line.find('"') + 1
            end = line.rfind('"')
            return line[start:end]

    # Fallback to default URL
    return "https://cern.ch/geant4-data/datasets"


def verify_md5(filepath: Path, expected_md5: str) -> bool:
    """Verify MD5 checksum of a file."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest() == expected_md5


def extract_and_install(
    tarball_path: Path, temp_dir: Path, dest_dir: Path
) -> tuple[bool, str]:
    """Extract tarball and install to final location (runs in process pool)."""
    try:
        # Extract
        with tarfile.open(tarball_path, "r:gz") as tar:
            tar.extractall(temp_dir, filter="data")

        # Install to final location
        dataset_dir_name = dest_dir.name
        src_dir = temp_dir / dataset_dir_name

        dest_dir.parent.mkdir(parents=True, exist_ok=True)

        if dest_dir.exists():
            shutil.rmtree(dest_dir)

        shutil.move(str(src_dir), str(dest_dir))

        # Clean up tarball
        tarball_path.unlink()

        return True, f"Successfully installed {dataset_dir_name}"
    except Exception as e:
        return False, f"Failed to extract/install: {e}"


async def download_dataset(
    client: httpx.AsyncClient,
    dataset: dict,
    base_url: str,
    temp_dir: Path,
    progress: Progress,
    executor: concurrent.futures.ProcessPoolExecutor,
) -> tuple[bool, str]:
    """Download and verify a single dataset."""
    filename = dataset["filename"]
    url = f"{base_url}/{filename}"
    dest_path = temp_dir / filename

    task_id = progress.add_task(f"[cyan]{filename}", total=None)

    try:
        # Download
        async with client.stream("GET", url, follow_redirects=True) as response:
            response.raise_for_status()
            total = int(response.headers.get("content-length", 0))
            progress.update(task_id, total=total)

            with open(dest_path, "wb") as f:
                downloaded = 0
                async for chunk in response.aiter_bytes(chunk_size=8192):
                    f.write(chunk)
                    downloaded += len(chunk)
                    progress.update(task_id, completed=downloaded)

        # Verify MD5 in process pool (non-blocking)
        progress.update(task_id, description=f"[yellow]{filename} (verifying)")
        loop = asyncio.get_event_loop()
        md5_valid = await loop.run_in_executor(
            executor,
            verify_md5,
            dest_path,
            dataset["md5"],
        )

        if not md5_valid:
            progress.update(task_id, description=f"[red]{filename} (MD5 mismatch)")
            return False, f"MD5 mismatch for {filename}"

        # Extract and install in process pool (non-blocking)
        progress.update(task_id, description=f"[yellow]{filename} (extracting)")
        dest_dir = Path(dataset["path"])
        success, msg = await loop.run_in_executor(
            executor,
            extract_and_install,
            dest_path,
            temp_dir,
            dest_dir,
        )

        if success:
            progress.update(task_id, description=f"[green]{filename} (installed)")
            return True, f"Successfully installed {dataset['name']}"
        else:
            progress.update(task_id, description=f"[red]{filename} (failed)")
            return False, msg

    except Exception as e:
        progress.update(task_id, description=f"[red]{filename} (failed)")
        return False, f"Failed to download {filename}: {e}"


async def download_all_datasets(
    datasets: list[dict],
    base_url: str,
    temp_dir: Path,
    max_concurrent: int,
    dry_run: bool = False,
    force: bool = False,
) -> int:
    """Download all datasets with limited concurrency.

    Returns:
        Number of failures (0 if all successful)
    """
    # Filter out already installed datasets
    if force:
        datasets_to_install = datasets
    else:
        datasets_to_install = [ds for ds in datasets if not Path(ds["path"]).exists()]

    if not datasets_to_install:
        console.print("[green]All datasets already installed[/green]")
        return 0

    if dry_run:
        console.print(
            f"[yellow]DRY RUN: Would download {len(datasets_to_install)} datasets:[/yellow]"
        )
        for ds in datasets_to_install:
            console.print(f"  [cyan]•[/cyan] {ds['name']} ({ds['filename']})")
            console.print(f"    URL: {base_url}/{ds['filename']}")
            console.print(f"    Destination: {ds['path']}")
            console.print(f"    MD5: {ds['md5']}")
        return 0

    console.print(f"[cyan]Downloading {len(datasets_to_install)} datasets...[/cyan]")

    progress = Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        DownloadColumn(),
        TransferSpeedColumn(),
        TimeRemainingColumn(),
        console=console,
    )

    # Use process pool for extraction
    with concurrent.futures.ProcessPoolExecutor() as executor:
        async with httpx.AsyncClient(timeout=1800.0) as client:
            with progress:
                # Use semaphore to limit concurrent downloads
                semaphore = asyncio.Semaphore(max_concurrent)

                async def bounded_download(dataset):
                    async with semaphore:
                        return await download_dataset(
                            client, dataset, base_url, temp_dir, progress, executor
                        )

                results = await asyncio.gather(
                    *[bounded_download(ds) for ds in datasets_to_install],
                    return_exceptions=True,
                )

    # Print summary
    console.print()
    successes = sum(1 for r in results if not isinstance(r, Exception) and r[0])
    failures = len(results) - successes

    if failures == 0:
        console.print(f"[green]✓ Successfully installed {successes} datasets[/green]")
    else:
        console.print(
            f"[yellow]⚠ Installed {successes} datasets, {failures} failed[/yellow]"
        )
        for result in results:
            if isinstance(result, Exception) or not result[0]:
                msg = str(result) if isinstance(result, Exception) else result[1]
                console.print(f"[red]  • {msg}[/red]")

    return failures


@app.command()
def main(
    max_concurrent: Annotated[
        int, typer.Option("--jobs", "-j", help="Maximum concurrent downloads")
    ] = 4,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            help="Show what would be downloaded without actually downloading",
        ),
    ] = False,
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            help="Redownload and reinstall datasets even if they already exist",
        ),
    ] = False,
    config: Annotated[
        Path | None, typer.Option("--config", help="Path to geant4-config script")
    ] = None,
) -> None:
    """Download Geant4 datasets in parallel."""
    # Find geant4-config
    if config:
        config_path = config
        if not config_path.exists():
            console.print(f"[red]Error: {config_path} does not exist[/red]")
            raise typer.Exit(1)
    else:
        config_path = find_geant4_config()
    console.print(f"[cyan]Found geant4-config at: {config_path}[/cyan]")

    # Parse datasets
    datasets = parse_datasets(config_path)
    console.print(f"[cyan]Found {len(datasets)} datasets[/cyan]")

    # Get base URL
    base_url = get_dataset_url(config_path)
    console.print(f"[cyan]Base URL: {base_url}[/cyan]")
    console.print()

    # Create temp directory
    with tempfile.TemporaryDirectory(prefix="geant4-downloads-") as temp_dir:
        # Download datasets
        failures = asyncio.run(
            download_all_datasets(
                datasets, base_url, Path(temp_dir), max_concurrent, dry_run, force
            )
        )

    # Exit with error code if any downloads failed
    if failures > 0:
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
