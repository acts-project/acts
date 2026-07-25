#!/usr/bin/env bash
# Diagnostic wrapper around delocate for the macOS PyPI wheels.
#
# Enabling the Arrow/Parquet plugin for wheels (#5659) makes the wheel link
# openssl, after which delocate fails with:
#   DelocationError: Already planning to copy library with same basename as:
#   libssl.3.dylib
# even though the spack environment contains a single openssl (verified from
# the lockfile). delocate is therefore resolving libssl.3.dylib to two distinct
# paths (e.g. the spack view symlink vs the store realpath). This wrapper dumps
# the mach-o linkage and rpaths of every bundled binary so we can see those two
# paths, then runs delocate as usual (still expected to fail until the linkage
# is normalized).
#
# cibuildwheel invokes this as CIBW_REPAIR_WHEEL_COMMAND_MACOS with:
#   $1 = wheel, $2 = dest_dir, $3 = delocate_archs
set -euo pipefail

wheel="$1"
dest_dir="$2"
archs="$3"

inspect_dir="$(mktemp -d)"
unzip -q "$wheel" -d "$inspect_dir"

echo "::group::mach-o linkage of bundled binaries"
find "$inspect_dir" \( -name '*.so' -o -name '*.dylib' \) -print | sort | while read -r f; do
    echo "--- ${f#"$inspect_dir"/}"
    otool -L "$f" || true
    echo "  rpaths:"
    otool -l "$f" | awk '/LC_RPATH/{r=1} r&&/ path /{print "    " $2; r=0}' || true
done
echo "::endgroup::"

echo "Running delocate-wheel (expected to fail until the linkage is fixed)"
delocate-wheel --require-archs "$archs" -w "$dest_dir" -v "$wheel"
