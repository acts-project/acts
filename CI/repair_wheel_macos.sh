#!/usr/bin/env bash
# Repair wrapper for the macOS PyPI wheels, invoked as
# CIBW_REPAIR_WHEEL_COMMAND_MACOS with: $1 = wheel, $2 = dest_dir, $3 = archs.
#
# Root cause (diagnosed with delocate-listdeps --all --depending):
#   * libActsPluginArrow.dylib (static Arrow) links spack's openssl:
#       .../spack/opt/spack/.../openssl-<hash>/lib/libssl.3.dylib
#   * spack's libthrift (an Arrow dependency) links the python.org framework's
#     openssl instead:
#       /Library/Frameworks/Python.framework/Versions/<X>/lib/libssl.3.dylib
# delocate then sees two different libssl.3.dylib (and libcrypto.3.dylib) with
# the same basename and aborts:
#   DelocationError: Already planning to copy library with same basename as:
#   libssl.3.dylib
#
# Workaround: before delocation, repoint every spack library that references the
# framework openssl at the spack openssl, so the whole graph resolves to a single
# copy that delocate can vendor. The proper fix belongs in ci-dependencies (build
# thrift/arrow against spack's openssl, or without openssl at all — ACTS Parquet
# uses no TLS); this keeps the wheels building until then.
#
# This modifies the spack store on the (ephemeral) CI runner. delocate re-signs
# the copies it vendors into the wheel, so the invalidated store-lib signatures
# do not matter, and the wheel's tests run against the self-contained wheel.
set -euo pipefail

wheel="$1"
dest_dir="$2"
archs="$3"

fw_re='Python\.framework/.*/lib(ssl|crypto)\.3\.dylib'

inspect_dir="$(mktemp -d)"
unzip -q "$wheel" -d "$inspect_dir"

# Locate spack's openssl via the plugin that links it by absolute path.
arrow_dylib="$(find "$inspect_dir" -name 'libActsPluginArrow.dylib' | head -1)"
spack_ssl="$(otool -L "$arrow_dylib" 2>/dev/null \
    | awk '/\/openssl-[^ ]*\/lib\/libssl\.3\.dylib /{print $1; exit}')"
if [ -z "${spack_ssl:-}" ]; then
    echo "ERROR: could not locate spack libssl from ${arrow_dylib}; linkage:" >&2
    otool -L "$arrow_dylib" >&2 || true
    exit 1
fi
spack_ssl_dir="$(dirname "$spack_ssl")"
spack_root="${spack_ssl_dir%/opt/spack/*}"
echo "spack openssl dir: ${spack_ssl_dir}"

echo "::group::Repoint framework openssl references to spack openssl"
find "$spack_root" -name '*.dylib' -type f | while read -r lib; do
    otool -L "$lib" 2>/dev/null | awk 'NR>1{print $1}' | grep -E "$fw_re" | while read -r fw; do
        name="$(basename "$fw")"
        echo "  ${lib}: ${fw} -> ${spack_ssl_dir}/${name}"
        install_name_tool -change "$fw" "${spack_ssl_dir}/${name}" "$lib" || true
    done
done
echo "::endgroup::"

if ! delocate-wheel --require-archs "$archs" -w "$dest_dir" -v "$wheel"; then
    echo "::group::delocate still failing — post-rewrite dependency tree"
    delocate-listdeps --all --depending "$wheel" || true
    echo "::endgroup::"
    exit 1
fi
