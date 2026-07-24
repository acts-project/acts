#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "typer",
#   "rich",
#   "httpx",
# ]
# ///
"""
Manage the physmon reference files.

The reference histograms are not committed to the repository. Instead,
`CI/physmon/reference.sha256` records the SHA256 of each reference file, and the
contents live as content-addressed blobs in an OCI registry. This script
resolves that manifest.

Run it with `uv run --no-project CI/physmon/reference.py <command>`.
"""

import concurrent.futures
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import threading
from pathlib import Path
from typing import Annotated

import httpx
import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn

# soft_wrap keeps long paths on one line, which matters when this runs in a CI
# log with a narrow assumed width
console = Console(stderr=True, soft_wrap=True)
app = typer.Typer(
    help=__doc__,
    no_args_is_help=True,
    add_completion=False,
    pretty_exceptions_show_locals=False,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = REPO_ROOT / "CI" / "physmon" / "reference.sha256"
DEFAULT_REFDIR = REPO_ROOT / "CI" / "physmon" / "reference"
DOCS = "docs/pages/contributing/physmon.md"
DEFAULT_REGISTRY = os.environ.get(
    "ACTS_PHYSMON_REGISTRY", "ghcr.io/acts-project/physmon-references"
)

# Blobs are stored raw, so the layer digest is the SHA256 of the file itself.
# Being explicit keeps oras from tarring the input.
MEDIA_TYPE = "application/octet-stream"

CHUNK = 1024 * 1024

RegistryOption = Annotated[
    str, typer.Option("--registry", help="OCI repository holding the blobs")
]
ManifestOption = Annotated[
    Path, typer.Option("--manifest", help="Manifest listing the reference hashes")
]
RefdirOption = Annotated[
    Path, typer.Option("--reference-dir", help="Where the reference files go")
]
JobsOption = Annotated[int, typer.Option("--jobs", "-j", help="Concurrent downloads")]


# --------------------------------------------------------------------------- #
# manifest
# --------------------------------------------------------------------------- #

# The manifest is `sha256sum` format: "<hex>  <relative path>". It is plain text
# so that reference updates show up as a reviewable diff, and so that
# `sha256sum -c` works on it directly from the reference directory.


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def read_manifest(path: Path) -> dict[str, str]:
    if not path.exists():
        console.print(f"[red]No manifest at {path}[/red]")
        raise typer.Exit(1)

    entries: dict[str, str] = {}
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            digest, relpath = line.split(None, 1)
        except ValueError:
            console.print(f"[red]{path}:{lineno}: malformed line: {line!r}[/red]")
            raise typer.Exit(1)
        digest = digest.lstrip("*")
        if len(digest) != 64 or not all(c in "0123456789abcdef" for c in digest):
            console.print(f"[red]{path}:{lineno}: not a sha256: {digest!r}[/red]")
            raise typer.Exit(1)
        entries[relpath.strip()] = digest
    return entries


def render_manifest(entries: dict[str, str]) -> str:
    return "".join(f"{entries[k]}  {k}\n" for k in sorted(entries))


def write_manifest(path: Path, entries: dict[str, str]) -> None:
    path.write_text(render_manifest(entries))


def manifest_from_dir(directory: Path) -> dict[str, str]:
    entries = {
        path.relative_to(directory).as_posix(): sha256_file(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }
    if not entries:
        console.print(f"[red]No files found under {directory}[/red]")
        raise typer.Exit(1)
    return entries


def manifest_tag(entries: dict[str, str]) -> str:
    """Tag derived from the manifest content, so re-publishing is idempotent."""
    return "m-" + hashlib.sha256(render_manifest(entries).encode()).hexdigest()[:16]


# --------------------------------------------------------------------------- #
# registry
# --------------------------------------------------------------------------- #


def split_registry(registry: str) -> tuple[str, str]:
    host, _, name = registry.partition("/")
    if not name:
        console.print(f"[red]Malformed registry reference: {registry!r}[/red]")
        raise typer.Exit(1)
    return host, name


# Registry failures land on people who have never opened this file, so they get
# an explanation of what to do rather than an httpx traceback.


_reported: set[str] = set()
_report_lock = threading.Lock()


def abort(headline: str, *detail: str) -> typer.Exit:
    """Print an explanation and return the exception to raise.

    Downloads run concurrently, so a registry-wide problem hits every worker at
    once. Print each explanation only the first time it comes up.
    """
    with _report_lock:
        if headline in _reported:
            return typer.Exit(1)
        _reported.add(headline)
        console.print(f"[red]{headline}[/red]")
        for line in detail:
            console.print(line)
    return typer.Exit(1)


def env_token() -> tuple[str | None, str | None]:
    """The credential to use, and the variable it came from."""
    for var in ("GITHUB_TOKEN", "ACTS_PHYSMON_TOKEN"):
        value = os.environ.get(var)
        if value:
            return value, var
    return None, None


def auth_failure(host: str, name: str, status: int) -> typer.Exit:
    _, var = env_token()
    if var is None:
        cause = f"{host} refused an anonymous pull; the package should be public."
    else:
        cause = f"The token in {var} was rejected."

    return abort(
        f"Cannot pull from {host}/{name} (HTTP {status})",
        "",
        cause,
        "To pull with credentials, export a token with the 'read:packages' scope as",
        "GITHUB_TOKEN or ACTS_PHYSMON_TOKEN. 'gh auth token' does not carry that scope",
        "by default; add it with 'gh auth refresh --scopes read:packages'.",
        f"See '{DOCS}'.",
    )


def missing_blob(host: str, name: str, digest: str) -> typer.Exit:
    return abort(
        f"{host}/{name} has no blob sha256:{digest}",
        "",
        "The manifest names a reference file that was never published. References have",
        "to be uploaded before the manifest naming them is committed; see 'How do I",
        f"update the reference files?' in '{DOCS}'.",
    )


def network_failure(host: str, exc: httpx.RequestError) -> typer.Exit:
    return abort(
        f"Could not reach {host}: {exc}",
        "",
        "The reference histograms are not committed, so physmon cannot run without them.",
        "Set ACTS_PHYSMON_NO_FETCH=1 to use a reference directory you populated yourself.",
    )


def registry_token(client: httpx.Client, host: str, name: str) -> str:
    """Get a pull token.

    Public packages allow anonymous pulls, but anonymous requests are rate
    limited aggressively enough that CI trips over it, so use GITHUB_TOKEN when
    one is available.
    """
    token, _ = env_token()
    auth = ("x-access-token", token) if token else None
    try:
        response = client.get(
            f"https://{host}/token",
            params={"service": host, "scope": f"repository:{name}:pull"},
            auth=auth,
            timeout=30,
        )
        response.raise_for_status()
    except httpx.HTTPStatusError as exc:
        status = exc.response.status_code
        if status in (401, 403):
            raise auth_failure(host, name, status) from None
        raise abort(f"{host} returned HTTP {status} for a pull token") from None
    except httpx.RequestError as exc:
        raise network_failure(host, exc) from None

    payload = response.json()
    # The spec calls the field "token"; "access_token" is the OAuth2 spelling
    # some registries answer with instead
    bearer = payload.get("token") or payload.get("access_token")
    if not bearer:
        raise abort(f"{host} did not return a token: {payload!r}")
    return bearer


def cache_dir() -> Path:
    override = os.environ.get("ACTS_PHYSMON_CACHE")
    if override:
        return Path(override)
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "acts" / "physmon-references"


def download_blob(
    client: httpx.Client, host: str, name: str, digest: str, dest: Path
) -> None:
    url = f"https://{host}/v2/{name}/blobs/sha256:{digest}"

    # Two reference files can legitimately have identical content, in which case
    # several workers download the same digest at once. Stage each attempt under
    # its own name so they cannot overwrite each other.
    handle, staged = tempfile.mkstemp(dir=dest.parent, prefix=f"{digest}.")
    tmp = Path(staged)
    try:
        try:
            # fdopen owns the descriptor, so take it before anything can raise
            with os.fdopen(handle, "wb") as f:
                with client.stream("GET", url, timeout=300) as response:
                    response.raise_for_status()
                    for chunk in response.iter_bytes(CHUNK):
                        f.write(chunk)
        except httpx.HTTPStatusError as exc:
            status = exc.response.status_code
            # A token that the token endpoint handed out happily can still be
            # refused here, when it does not carry the scope the blob needs
            if status in (401, 403):
                raise auth_failure(host, name, status) from None
            if status == 404:
                raise missing_blob(host, name, digest) from None
            raise abort(f"Fetching sha256:{digest} failed with HTTP {status}") from None
        except httpx.RequestError as exc:
            raise network_failure(host, exc) from None

        got = sha256_file(tmp)
        if got != digest:
            console.print(f"[red]Digest mismatch: expected {digest}, got {got}[/red]")
            raise typer.Exit(1)
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)


def materialize(
    entries: dict[str, str], output_dir: Path, registry: str, jobs: int
) -> int:
    """Populate output_dir with the files named by the manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir()
    cache.mkdir(parents=True, exist_ok=True)

    def up_to_date(relpath: str, digest: str) -> bool:
        target = output_dir / relpath
        return target.exists() and sha256_file(target) == digest

    outstanding = {k: v for k, v in entries.items() if not up_to_date(k, v)}
    if not outstanding:
        console.print(f"{len(entries)} reference file(s) already present")
        return 0

    # Only talk to the registry if something is actually missing from the cache
    needs_network = any(not (cache / d).exists() for d in outstanding.values())
    host, name = split_registry(registry)

    with httpx.Client(follow_redirects=True) as client:
        if needs_network:
            # httpx drops the Authorization header when a redirect crosses
            # origins, which is what blob reads do when they hand off to a CDN
            token = registry_token(client, host, name)
            client.headers["Authorization"] = f"Bearer {token}"

        def one(item: tuple[str, str]) -> None:
            relpath, digest = item
            blob = cache / digest
            if not (blob.exists() and sha256_file(blob) == digest):
                download_blob(client, host, name, digest, blob)
            target = output_dir / relpath
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(blob, target)

        with Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TimeRemainingColumn(),
            console=console,
            transient=True,
        ) as progress:
            task = progress.add_task("Fetching references", total=len(outstanding))
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
                for _ in pool.map(one, outstanding.items()):
                    progress.advance(task)

    console.print(
        f"{len(outstanding)} reference file(s) fetched, "
        f"{len(entries) - len(outstanding)} already present"
    )
    return len(outstanding)


def require(tool: str) -> None:
    if shutil.which(tool) is None:
        console.print(
            f"[red]'{tool}' is required for this command but is not on PATH[/red]"
        )
        raise typer.Exit(1)


def oras_push(registry: str, directory: Path, entries: dict[str, str]) -> str:
    """Push the full reference set as one tagged artifact.

    Pushing every file (rather than only the new ones) keeps every blob
    referenced by a manifest, so the registry never considers them garbage. The
    registry skips uploading blobs it already has, so unchanged files cost
    nothing beyond an existence check.
    """
    require("oras")
    ref = f"{registry}:{manifest_tag(entries)}"

    args = ["oras", "push", "--format", "json", ref]
    args += [f"{relpath}:{MEDIA_TYPE}" for relpath in sorted(entries)]

    console.print(f"Pushing {len(entries)} file(s) to [bold]{ref}[/bold]")
    result = subprocess.run(
        args, cwd=directory, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        console.print(result.stdout)
        console.print(f"[red]{result.stderr}[/red]")
        console.print(f"[red]oras push failed (exit {result.returncode})[/red]")
        raise typer.Exit(1)

    # Guard against oras ever packing the files instead of storing them raw: the
    # layer digests must be the file hashes, or fetch-by-hash cannot work.
    try:
        layers = {
            layer["digest"].removeprefix("sha256:")
            for layer in json.loads(result.stdout).get("layers", [])
        }
    except (json.JSONDecodeError, KeyError, TypeError):
        console.print(
            "[yellow]Could not parse oras output, skipping digest check[/yellow]"
        )
        return ref

    missing = set(entries.values()) - layers
    if layers and missing:
        console.print(
            "[red]oras did not store the files as raw blobs; "
            f"{len(missing)} expected digest(s) absent from the manifest[/red]"
        )
        raise typer.Exit(1)
    return ref


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


@app.command()
def pull(
    manifest: ManifestOption = DEFAULT_MANIFEST,
    reference_dir: RefdirOption = DEFAULT_REFDIR,
    registry: RegistryOption = DEFAULT_REGISTRY,
    jobs: JobsOption = 8,
) -> None:
    """Populate the reference directory from the manifest."""
    materialize(read_manifest(manifest), reference_dir, registry, jobs)


@app.command()
def verify(
    manifest: ManifestOption = DEFAULT_MANIFEST,
    reference_dir: RefdirOption = DEFAULT_REFDIR,
) -> None:
    """Check an existing reference directory against the manifest."""
    entries = read_manifest(manifest)

    problems = []
    for relpath, digest in sorted(entries.items()):
        path = reference_dir / relpath
        if not path.exists():
            problems.append(f"missing: {relpath}")
        elif sha256_file(path) != digest:
            problems.append(f"modified: {relpath}")

    if reference_dir.exists():
        present = {
            p.relative_to(reference_dir).as_posix()
            for p in reference_dir.rglob("*")
            if p.is_file()
        }
        problems += [f"untracked: {p}" for p in sorted(present - set(entries))]

    for line in problems:
        console.print(f"[red]{line}[/red]")
    if problems:
        console.print(f"[red]{len(problems)} problem(s) against {manifest}[/red]")
        raise typer.Exit(1)
    console.print(f"[green]{len(entries)} reference file(s) match {manifest}[/green]")


def read_results(path: Path) -> list[dict[str, str]]:
    """Read histcmp_results.csv, tolerating the pre-manifest 3-column format."""
    rows = []
    with path.open() as f:
        for row in csv.reader(f):
            if not row:
                continue
            row = row + [""] * (5 - len(row))
            title, html_path, ec, ref_path, monitored_path = row[:5]
            rows.append(
                {
                    "title": title,
                    "html_path": html_path,
                    "ec": ec,
                    "ref_path": ref_path,
                    "monitored_path": monitored_path,
                }
            )
    return rows


@app.command()
def candidate(
    results: Annotated[
        Path, typer.Option("--results", help="histcmp_results.csv from the run")
    ],
    data_dir: Annotated[
        Path, typer.Option("--data-dir", help="Directory holding the run's outputs")
    ],
    output: Annotated[
        Path, typer.Option("--output", help="Where to write the manifest")
    ],
    manifest: ManifestOption = DEFAULT_MANIFEST,
    replace_all: Annotated[
        bool,
        typer.Option("--all", help="Replace every entry, not just failing comparisons"),
    ] = False,
) -> None:
    """Build the manifest that would apply if a run's failures were accepted.

    The changed set comes from the histcmp exit codes, not from comparing
    hashes: ROOT files embed a creation timestamp, so every physmon run produces
    a different hash for every file even when the physics is unchanged.
    """
    entries = read_manifest(manifest)

    updated, skipped = [], []
    for row in read_results(results):
        if not row["ref_path"]:
            continue  # comparison between two monitored outputs, no reference
        if row["ec"] == "0" and not replace_all:
            continue
        source = data_dir / (row["monitored_path"] or row["ref_path"])
        if not source.exists():
            skipped.append(f"{row['ref_path']} (no output at {source})")
            continue
        entries[row["ref_path"]] = sha256_file(source)
        updated.append(row["ref_path"])

    write_manifest(output, entries)

    for relpath in sorted(updated):
        console.print(f"candidate update: {relpath}")
    for line in sorted(skipped):
        console.print(f"[yellow]candidate skipped: {line}[/yellow]")
    console.print(f"Wrote {output} ({len(updated)} of {len(entries)} entries updated)")


@app.command("import")
def import_references(
    reference_dir: RefdirOption = DEFAULT_REFDIR,
    manifest: ManifestOption = DEFAULT_MANIFEST,
    registry: RegistryOption = DEFAULT_REGISTRY,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Hash and write the manifest, upload nothing"),
    ] = False,
) -> None:
    """Seed the registry from an existing reference directory.

    This is the one-off migration step: it uploads the reference files that used
    to be committed and writes the manifest that replaces them.
    """
    entries = manifest_from_dir(reference_dir)
    console.print(f"{len(entries)} file(s) under {reference_dir}")

    if dry_run:
        console.print("[yellow]Dry run: not contacting the registry[/yellow]")
    else:
        console.print(
            f"Published [bold]{oras_push(registry, reference_dir, entries)}[/bold]"
        )

    write_manifest(manifest, entries)
    console.print(f"[green]Wrote {manifest}[/green]")


@app.command()
def update(
    run_id: Annotated[
        str, typer.Argument(help="Builds run holding the physmon artifact")
    ],
    repo: Annotated[
        str, typer.Option("--repo", help="Repository the run belongs to")
    ] = os.environ.get("GITHUB_REPOSITORY", "acts-project/acts"),
    output: Annotated[
        Path | None,
        typer.Option("--output", help="Write the manifest here instead of stdout"),
    ] = None,
    registry: RegistryOption = DEFAULT_REGISTRY,
    jobs: JobsOption = 8,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Report what would happen, upload nothing")
    ] = False,
) -> None:
    """Publish the references a physmon run implies, and print the new manifest."""
    require("gh")

    with tempfile.TemporaryDirectory() as tmp:
        artifact = Path(tmp) / "artifact"
        artifact.mkdir()
        console.print(f"Downloading physmon artifact from run {run_id}")
        result = subprocess.run(
            [
                "gh",
                "run",
                "download",
                run_id,
                "--repo",
                repo,
                "--name",
                "physmon",
                "--dir",
                str(artifact),
            ],
            check=False,
        )
        if result.returncode != 0:
            console.print(
                f"[red]Could not download the 'physmon' artifact from run {run_id}. "
                "GitHub deletes artifacts after a retention period, so the run may "
                "need to be repeated.[/red]"
            )
            raise typer.Exit(1)

        candidate_manifest = artifact / "reference-candidate.sha256"
        if not candidate_manifest.exists():
            console.print(
                f"[red]{candidate_manifest.name} is not in the artifact. The run "
                "predates the manifest-based reference system, or physmon did not "
                "finish.[/red]"
            )
            raise typer.Exit(1)
        entries = read_manifest(candidate_manifest)

        # Assemble the complete set: new content comes from the run's outputs,
        # everything else is already in the registry.
        staging = Path(tmp) / "staging"
        staging.mkdir(parents=True, exist_ok=True)
        data = artifact / "data"
        adopted = 0
        for relpath, digest in entries.items():
            source = data / relpath
            if source.exists() and sha256_file(source) == digest:
                target = staging / relpath
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
                adopted += 1
        console.print(
            f"{adopted} file(s) taken from the run, {len(entries) - adopted} reused"
        )

        if dry_run:
            console.print("[yellow]Dry run: not contacting the registry[/yellow]")
        else:
            # The unchanged files are already in the registry, but the push
            # re-references every blob so none of them can be garbage collected
            materialize(entries, staging, registry, jobs)
            console.print(
                f"Published [bold]{oras_push(registry, staging, entries)}[/bold]"
            )

    if output:
        write_manifest(output, entries)
        console.print(f"[green]Wrote {output}[/green]")
    else:
        print(render_manifest(entries), end="")


if __name__ == "__main__":
    app()
