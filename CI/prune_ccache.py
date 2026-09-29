#!/usr/bin/env python3
"""Retain the newest main compiler cache per variant and archive format.

Dry-run by default. Only SHA-suffixed ccache keys are eligible; an uploaded,
non-empty replacement must exist before an older archive can be deleted.
"""

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import subprocess

MAIN_REF = "refs/heads/main"
CACHE_KEY = re.compile(r"^(ccache-.+-r\d+(?:-.+)?)-[0-9a-f]{40}$")


def superseded_caches(caches):
    """Return (replacement, obsolete archives) pairs, preserving cache versions."""
    groups = defaultdict(list)
    for cache in caches:
        match = CACHE_KEY.fullmatch(cache["key"])
        if (
            cache["ref"] == MAIN_REF
            and match
            and cache.get("version")
            and cache["size_in_bytes"] > 0
        ):
            groups[(match[1], cache["version"])].append(cache)
    result = []
    for group in groups.values():
        ordered = sorted(group, key=lambda c: (c["created_at"], c["id"]), reverse=True)
        if len(ordered) > 1:
            result.append((ordered[0], ordered[1:]))
    return result


def list_caches(repo, key="ccache-"):
    # Fetch all pages before any deletion, so changing page boundaries cannot
    # cause entries to be skipped. -X GET prevents -f switching to POST.
    output = subprocess.check_output(
        [
            "gh",
            "api",
            "--method",
            "GET",
            "--paginate",
            "--slurp",
            f"repos/{repo}/actions/caches",
            "-f",
            f"ref={MAIN_REF}",
            "-f",
            f"key={key}",
            "-f",
            "per_page=100",
        ],
        text=True,
    )
    return [cache for page in json.loads(output) for cache in page["actions_caches"]]


def prune(repo, caches, apply=False):
    total = 0
    for replacement, obsolete in superseded_caches(caches):
        if apply:
            # GitHub may have evicted our replacement since the initial listing.
            # Recheck it before pruning this variant. A new upload is never in
            # the deletion list, which contains only IDs from the initial read.
            available = list_caches(repo, replacement["key"])
            if not any(
                c["id"] == replacement["id"] and c["size_in_bytes"] > 0
                for c in available
            ):
                print(f"Keep older caches: replacement {replacement['id']} disappeared")
                continue
        print(f"Keep {replacement['id']}: {replacement['key']}")
        for cache in obsolete:
            print(
                f"{'Delete' if apply else 'Would delete'} {cache['id']}: {cache['key']}"
            )
            if apply:
                subprocess.run(
                    [
                        "gh",
                        "api",
                        "--method",
                        "DELETE",
                        f"repos/{repo}/actions/caches/{cache['id']}",
                    ],
                    check=True,
                )
            total += cache["size_in_bytes"]
    print(f"{'Removed' if apply else 'Reclaimable'} archive bytes: {total:,}")
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"))
    parser.add_argument(
        "--apply", action="store_true", help="Delete superseded archives"
    )
    parser.add_argument(
        "--snapshot", type=Path, help="Dry-run a saved cache API response"
    )
    args = parser.parse_args()
    if args.apply and args.snapshot:
        parser.error("--apply requires a fresh API listing, not --snapshot")
    if not args.snapshot and not args.repo:
        parser.error("--repo or GITHUB_REPOSITORY is required")
    caches = (
        json.loads(args.snapshot.read_text())["actions_caches"]
        if args.snapshot
        else list_caches(args.repo)
    )
    prune(args.repo, caches, args.apply)


if __name__ == "__main__":
    main()
