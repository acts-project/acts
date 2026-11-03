#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["boto3==1.42.0"]
# ///
"""Restore, publish, and prune bounded ccache snapshots; keep compilation lookups local.

Overlapping publishers may leave an older valid snapshot as latest.
Readers are always anonymous. Remote cache failures never fail the build.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import tarfile
import tempfile
import uuid


def report(message):
    print(message, flush=True)
    if path := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(path, "a") as output:
            output.write(f"- {message}\n")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_manifest(client, bucket, prefix):
    try:
        response = client.get_object(Bucket=bucket, Key=prefix + "latest.json")
    except Exception as error:
        if getattr(error, "response", {}).get("Error", {}).get("Code") in (
            "NoSuchKey",
            "404",
        ):
            return None
        raise
    with response["Body"] as stream:
        manifest = json.loads(stream.read(16385))
    if (
        manifest["format"] != 1
        or not re.fullmatch(
            re.escape(prefix) + r"snapshots/[0-9]+-[0-9]+-[0-9a-f]{32}\.tar",
            manifest["key"],
        )
        or not re.fullmatch(r"[0-9a-f]{64}", manifest["sha256"])
        or any(
            type(manifest[field]) is not int or manifest[field] < 1
            for field in ("run_number", "run_attempt")
        )
    ):
        raise ValueError("Invalid ccache snapshot manifest")
    return manifest


def extract(archive, destination):
    with tarfile.open(archive, "r:") as source:
        for member in source:
            path = PurePosixPath(member.name)
            if (
                path.is_absolute()
                or ".." in path.parts
                or not (member.isfile() or member.isdir())
            ):
                raise ValueError("Unsafe ccache archive entry")
            source.extract(member, path=destination, filter="data")


def restore(client, bucket, prefix, cache):
    # The job uses a fresh, dedicated cache directory. Never combine partial/stale
    # restores with a new snapshot, or discard a caller's existing local cache.
    cache.mkdir(parents=True, exist_ok=True)
    if any(cache.iterdir()):
        raise ValueError("Restore requires an empty ccache directory")
    manifest = read_manifest(client, bucket, prefix)
    if manifest is None:
        report("ccache snapshot: no published snapshot; starting empty")
        return
    with tempfile.TemporaryDirectory(dir=cache.parent) as directory:
        staging = Path(directory)
        archive = staging / "snapshot.tar"
        client.download_file(bucket, manifest["key"], str(archive))
        if digest(archive) != manifest["sha256"]:
            raise ValueError("ccache snapshot checksum mismatch")
        unpacked = staging / "unpacked"
        unpacked.mkdir()
        extract(archive, unpacked)
        cache.rmdir()
        unpacked.rename(cache)
        report(
            f"ccache restored main run {manifest['run_number']}, "
            f"attempt {manifest['run_attempt']}"
        )


def publish(client, bucket, prefix, cache, run_number, run_id, attempt):
    # Best-effort freshness check, not a lock: another publisher may finish
    # after this read. An older snapshot is still safe for ccache to use.
    order = (run_number, attempt)
    current = read_manifest(client, bucket, prefix)
    if current and (current["run_number"], current["run_attempt"]) >= order:
        report("ccache publication skipped: an equal or newer main run exists")
        return
    with tempfile.TemporaryDirectory(dir=cache.parent) as directory:
        archive = Path(directory) / "snapshot.tar"
        # ccache already compresses entries. Avoid a second compression pass.
        with tarfile.open(archive, "w:") as output:
            for entry in sorted(cache.iterdir()):
                output.add(entry, arcname=entry.name)
        manifest = {
            "format": 1,
            "key": f"{prefix}snapshots/{run_id}-{attempt}-{uuid.uuid4().hex}.tar",
            "sha256": digest(archive),
            "run_number": run_number,
            "run_attempt": attempt,
            "commit": os.environ.get("GITHUB_SHA", ""),
        }
        client.upload_file(str(archive), bucket, manifest["key"])
        # A failed/cancelled upload cannot replace the currently usable snapshot.
        client.put_object(
            Bucket=bucket,
            Key=prefix + "latest.json",
            Body=json.dumps(manifest).encode(),
            ContentType="application/json",
            CacheControl="no-cache",
        )
        report(
            f"ccache snapshot published ({archive.stat().st_size / 1024**2:.1f} MiB)"
        )


def list_snapshots(client, bucket, prefix):
    # Finish listing before reading pointers. Later uploads are not candidates.
    variants = {}
    pattern = re.compile(
        re.escape(prefix) + r".+/snapshots/[0-9]+-[0-9]+-[0-9a-f]{32}\.tar"
    )
    for page in client.get_paginator("list_objects_v2").paginate(
        Bucket=bucket, Prefix=prefix
    ):
        for item in page.get("Contents", []):
            key = item["Key"]
            if pattern.fullmatch(key):
                variant = key.rsplit("snapshots/", 1)[0]
                variants.setdefault(variant, {})[key] = item["Size"]

    return variants


def prune_variant(client, bucket, variant, objects, dry_run):
    manifest = read_manifest(client, bucket, variant)
    if manifest is None:
        raise ValueError("Missing latest manifest")
    protected = {manifest["key"]}
    keys = sorted(objects)
    for offset in range(0, len(keys), 1000):
        # Protect both observed pointers if publication overlaps cleanup.
        current = read_manifest(client, bucket, variant)
        if current is None:
            raise ValueError("Missing latest manifest")
        protected.add(current["key"])
        candidates = [
            key for key in keys[offset : offset + 1000] if key not in protected
        ]
        if not candidates:
            continue
        if dry_run:
            yield candidates
            continue
        result = client.delete_objects(
            Bucket=bucket,
            Delete={"Objects": [{"Key": key} for key in candidates]},
        )
        deleted = [item["Key"] for item in result.get("Deleted", [])]
        yield deleted
        if result.get("Errors") or set(deleted) != set(candidates):
            raise RuntimeError("Incomplete deletion")


def prune(client, bucket, prefix, dry_run=False):
    variants = list_snapshots(client, bucket, prefix)
    total = sum(sum(objects.values()) for objects in variants.values())
    removed_count = removed_bytes = failures = 0
    for variant, objects in sorted(variants.items()):
        try:
            for deleted in prune_variant(client, bucket, variant, objects, dry_run):
                removed_count += len(deleted)
                removed_bytes += sum(objects[key] for key in deleted)
        except Exception as error:
            failures += 1
            report(f"ccache pruning failed for {variant} ({type(error).__name__})")
    verb = "would delete" if dry_run else "deleted"
    report(
        f"ccache pruning: {verb} {removed_count} archives "
        f"({removed_bytes / 1024**2:.1f} MiB); "
        f"retained {(total - removed_bytes) / 1024**2:.1f} MiB of listed archives"
    )
    if failures:
        raise RuntimeError(f"ccache pruning encountered {failures} failures")


def make_client(writable):
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    options = dict(
        connect_timeout=10,
        read_timeout=60,
        retries={"mode": "standard", "max_attempts": 2},
        s3={"addressing_style": "path"},
        # CERN's S3 gateway does not need optional streaming checksum trailers.
        request_checksum_calculation="when_required",
        response_checksum_validation="when_required",
    )
    if not writable:
        options["signature_version"] = UNSIGNED
    return boto3.client(
        "s3",
        endpoint_url=os.environ["CACHE_ENDPOINT"],
        region_name=os.environ["CACHE_REGION"],
        config=Config(**options),
    )


def prune_main(dry_run):
    if not (
        os.environ.get("GITHUB_REPOSITORY") == "acts-project/acts"
        and os.environ.get("GITHUB_REF") == "refs/heads/main"
        and os.environ.get("GITHUB_EVENT_NAME") in ("schedule", "workflow_dispatch")
        and os.environ.get("CACHE_PREFIX")
        == "acts-sccache/ccache-snapshots/acts-project/acts/"
    ):
        report("ccache pruning refused: requires an authorized main cleanup run")
        return 1
    try:
        prune(
            make_client(True),
            os.environ["CACHE_BUCKET"],
            os.environ["CACHE_PREFIX"],
            dry_run,
        )
    except Exception as error:
        report(f"ccache pruning failed ({type(error).__name__})")
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=["restore", "publish", "prune"])
    parser.add_argument("--dry-run", action="store_true", help="Preview pruning only")
    args = parser.parse_args()
    if args.dry_run and args.operation != "prune":
        parser.error("--dry-run requires prune")
    if args.operation == "prune":
        return prune_main(args.dry_run)
    writable = args.operation == "publish"
    if writable and not (
        os.environ.get("GITHUB_REF") == "refs/heads/main"
        and (
            os.environ.get("GITHUB_EVENT_NAME") == "push"
            or (
                os.environ.get("CACHE_ALLOW_MAIN_NON_PUSH") == "true"
                and os.environ.get("GITHUB_EVENT_NAME")
                in ("schedule", "workflow_dispatch")
            )
        )
    ):
        report("ccache publication skipped: only authorized main runs may publish")
        return 0
    try:
        client = make_client(writable)
        bucket = os.environ["CACHE_BUCKET"]
        prefix = os.environ["CACHE_PREFIX"].rstrip("/") + "/"
        cache = Path(os.environ["CCACHE_DIR"])
        if writable:
            publish(
                client,
                bucket,
                prefix,
                cache,
                int(os.environ["GITHUB_RUN_NUMBER"]),
                int(os.environ["GITHUB_RUN_ID"]),
                int(os.environ["GITHUB_RUN_ATTEMPT"]),
            )
        else:
            restore(client, bucket, prefix, cache)
    except Exception as error:
        # Upload credentials are never included in diagnostics.
        print(
            f"::warning::ccache snapshot {args.operation} failed ({type(error).__name__}); "
            "continuing without remote cache updates",
            flush=True,
        )
        report(f"ccache snapshot {args.operation}: failed ({type(error).__name__})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
