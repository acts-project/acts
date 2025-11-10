#!/usr/bin/env python3

# This file is part of the ACTS project.
#
# Copyright (C) 2025 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Download and validate GNN model files for CI testing.

This script downloads model tar files, verifies the integrity of the downloads
using hardcoded SHA256 hashes, and extracts them to a standard location.
"""

import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

# Model definitions with URLs, SHA256 hashes, and extraction directories
MODELS = {
    "metric-learning-onnx": {
        "url": "https://acts.web.cern.ch/ci/gnn/onnx_models_v01.tar",
        "sha256": "335d829439d9a5ae99ffc8f0bbf6a62119a140475e779e2ed00de21fdebb3cb4",
        "extract_to": "metric_learning",
    },
    "metric-learning-torch": {
        "url": "https://acts.web.cern.ch/ci/gnn/torchscript_models_v01.tar",
        "sha256": "1185060ce697bbc96c9dc32b85e5f0eb4db1f64a645c0fc4d2cb2731cb2ef3dc",
        "extract_to": "metric_learning",
    },
    "odd-module-map": {
        "url": "https://acts.web.cern.ch/ci/gnn/odd_module_map_v01.tar",
        "sha256": "59f0457f0043bac8594e9f5a3310a709244de980a7b0c206d7d0d95f15455d73",
        "extract_to": "odd_module_map",
    },
}


def download_file(url: str, dest_path: Path) -> None:
    """Download a file from URL to destination path."""
    print(f"Downloading {url} to {dest_path}")
    urllib.request.urlretrieve(url, dest_path)
    print(f"  Downloaded {dest_path.stat().st_size} bytes")


def compute_file_hash(file_path: Path) -> str:
    """Compute SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        # Read in chunks to handle large files efficiently
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def verify_hash(file_path: Path, expected_hash: str) -> bool:
    """Verify file hash matches expected hash."""
    print(f"Verifying hash of {file_path}")
    actual_hash = compute_file_hash(file_path)
    print(f"  Computed SHA256: {actual_hash}")

    if actual_hash == expected_hash:
        print(f"  ✓ Hash verification passed")
        return True
    else:
        print(f"  ✗ Hash mismatch!", file=sys.stderr)
        print(f"    Expected: {expected_hash}", file=sys.stderr)
        print(f"    Got:      {actual_hash}", file=sys.stderr)
        return False


def extract_tar(tar_path: Path, extract_dir: Path) -> None:
    """Extract tar file to specified directory."""
    print(f"Extracting {tar_path} to {extract_dir}")
    extract_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(tar_path) as tar:
        # Get list of members for reporting
        members = tar.getmembers()
        print(f"  Extracting {len(members)} files...")
        # Use data filter for Python 3.12+ compatibility
        tar.extractall(extract_dir, filter="data")

    print(f"  ✓ Extraction complete")


def main() -> int:
    """Main entry point: download, verify, and extract all models."""
    # Determine output directory (ci_models/ relative to repository root)
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent.parent  # CI/dependencies -> CI -> root
    output_dir = repo_root / "ci_models"

    print(f"CI Model Downloader")
    print(f"Output directory: {output_dir}")
    print(f"=" * 70)

    # Create temporary download directory
    download_dir = output_dir / ".downloads"
    download_dir.mkdir(parents=True, exist_ok=True)

    all_success = True

    for model_name, model_info in MODELS.items():
        print(f"\nProcessing {model_name}...")
        print("-" * 70)

        url = model_info["url"]
        expected_hash = model_info["sha256"]
        extract_to = model_info["extract_to"]

        # Determine file paths
        tar_filename = Path(url).name
        tar_path = download_dir / tar_filename
        final_extract_dir = output_dir / extract_to

        try:
            # Download tar file
            download_file(url, tar_path)

            # Verify hash
            print(f"Expected SHA256: {expected_hash}")
            if not verify_hash(tar_path, expected_hash):
                print(f"✗ Failed to verify {model_name}", file=sys.stderr)
                all_success = False
                continue

            # Extract
            extract_tar(tar_path, final_extract_dir)

            print(f"✓ Successfully processed {model_name}")

        except Exception as e:
            print(f"✗ Error processing {model_name}: {e}", file=sys.stderr)
            all_success = False
            continue

    print("\n" + "=" * 70)
    if all_success:
        print("✓ All models downloaded, verified, and extracted successfully")
        return 0
    else:
        print("✗ Some models failed to process", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
