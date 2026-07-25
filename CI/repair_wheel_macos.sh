#!/usr/bin/env bash
# Diagnostic + candidate-fix wrapper around delocate for the macOS PyPI wheels.
#
# Enabling Arrow/Parquet for wheels (#5659) makes the macOS wheel link openssl.
# The wheel's sole openssl consumer is libActsPluginArrow.dylib, which links it
# by absolute spack-store path while the same single openssl is also reachable
# via the spack view symlink (deps/view/lib) in the rpaths. delocate then fails:
#   DelocationError: Already planning to copy library with same basename as:
#   libssl.3.dylib
# i.e. it plans to copy the same physical libssl.3.dylib twice under two paths.
#
# Candidate fix: a newer delocate realpath-canonicalizes dependency paths before
# planning copies, collapsing the store/view symlink pair. We upgrade delocate
# and retry. Regardless of outcome we preserve the un-repaired wheel and a full
# dependency listing so the linkage can be analyzed offline if the fix misses.
#
# cibuildwheel invokes this as CIBW_REPAIR_WHEEL_COMMAND_MACOS with:
#   $1 = wheel, $2 = dest_dir, $3 = delocate_archs
set -euo pipefail

wheel="$1"
dest_dir="$2"
archs="$3"

# Preserve fallback artifacts (the un-repaired wheel + dependency tree) for
# offline analysis. DIAG_OUTDIR is provided by the diagnostics workflow.
outdir="${DIAG_OUTDIR:-$PWD/diag-out}"
mkdir -p "$outdir"
cp "$wheel" "$outdir/"

echo "::group::delocate version (before upgrade)"
delocate-wheel --version || true
echo "::endgroup::"

echo "::group::delocate-listdeps --all --depending"
delocate-listdeps --all --depending "$wheel" 2>&1 | tee "$outdir/listdeps.txt" || true
echo "::endgroup::"

echo "::group::mach-o linkage of bundled binaries"
inspect_dir="$(mktemp -d)"
unzip -q "$wheel" -d "$inspect_dir"
find "$inspect_dir" \( -name '*.so' -o -name '*.dylib' \) -print | sort | while read -r f; do
    echo "--- ${f#"$inspect_dir"/}"
    otool -L "$f" || true
done
echo "::endgroup::"

echo "::group::Upgrade delocate and retry"
# Candidate fix: newer delocate canonicalizes symlinked dependency paths.
python3 -m pip install --quiet --upgrade delocate || true
delocate-wheel --version || true
echo "::endgroup::"

delocate-wheel --require-archs "$archs" -w "$dest_dir" -v "$wheel"
