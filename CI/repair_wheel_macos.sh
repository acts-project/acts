#!/usr/bin/env bash
# Repair wrapper for the macOS PyPI wheels, invoked as
# CIBW_REPAIR_WHEEL_COMMAND_MACOS with: $1 = wheel, $2 = dest_dir, $3 = archs.
#
# Root cause (diagnosed with delocate-listdeps --all --depending, and confirmed
# against spack's thrift package.py and thrift's own build/cmake/
# DefineOptions.cmake):
#   * libActsPluginArrow.dylib (static Arrow) links spack's openssl:
#       .../spack/opt/spack/.../openssl-<hash>/lib/libssl.3.dylib
#   * spack's libthrift (an Arrow dependency) concretizes `~openssl` (its spack
#     "openssl" variant is off, so spack never wires in spack's openssl for
#     it). But thrift builds with build_system=cmake, and unlike thrift's
#     AutotoolsBuilder, its CMakeBuilder.cmake_args() never passes
#     -DWITH_OPENSSL. Thrift's own DefineOptions.cmake then runs a bare
#     find_package(OpenSSL) regardless of the variant and auto-links whatever
#     it finds on the runner's default search path — the python.org
#     framework's copy:
#       /Library/Frameworks/Python.framework/Versions/<X>/lib/libssl.3.dylib
# delocate then sees two different libssl.3.dylib (and libcrypto.3.dylib) with
# the same basename and aborts:
#   DelocationError: Already planning to copy library with same basename as:
#   libssl.3.dylib
#
# Workaround: before delocation, repoint libthrift's framework openssl
# references at spack's openssl, so the graph resolves to a single copy that
# delocate can vendor. Both are openssl 3.x, and ACTS' Parquet usage does not
# touch TLS. Scoped to libthrift specifically — the only known offender —
# rather than scanning/mutating the whole spack store, to keep this local.
#
# How long this is needed: until thrift's spack package (spack/spack-packages,
# repos/spack_repo/builtin/packages/thrift/package.py) wires its "openssl"
# variant into CMakeBuilder.cmake_args() the way it already does for
# AutotoolsBuilder — a one-line upstream fix
# (self.define_from_variant("WITH_OPENSSL", "openssl")). Once that lands and
# ci-dependencies picks up the updated spack-packages, thrift stops linking the
# framework openssl and the rewrite below finds nothing to do on its own; the
# script can then be deleted.
#
# This modifies the spack store on the (ephemeral) CI runner. delocate re-signs
# the copies it vendors into the wheel, so the invalidated store-lib signatures
# do not matter, and the wheel's tests run against the self-contained wheel.
set -euo pipefail

wheel="$1"
dest_dir="$2"
archs="$3"

# Any framework-provided openssl, regardless of soversion.
fw_re='/Library/Frameworks/Python\.framework/.*/lib(ssl|crypto)\.[0-9.]*dylib'

delocate() {
    if ! delocate-wheel --require-archs "$archs" -w "$dest_dir" -v "$wheel"; then
        echo "::group::delocate failed — dependency tree"
        delocate-listdeps --all --depending "$wheel" || true
        echo "::endgroup::"
        exit 1
    fi
}

inspect_dir="$(mktemp -d)"
trap 'rm -rf "$inspect_dir"' EXIT
unzip -q "$wheel" -d "$inspect_dir"

# Derive the spack store root from whatever absolute spack path the wheel's own
# binaries link against. If the wheel does not reference spack at all there is
# nothing to normalize.
spack_root=""
while IFS= read -r bin; do
    ref="$(otool -L "$bin" 2>/dev/null | awk '/\/opt\/spack\//{print $1; exit}')" || true
    if [ -n "${ref:-}" ]; then
        spack_root="${ref%/opt/spack/*}"
        break
    fi
done < <(find "$inspect_dir" \( -name '*.dylib' -o -name '*.so' \) -type f)

if [ -z "$spack_root" ]; then
    echo "wheel does not link the spack store; delocating as-is"
    delocate
    exit 0
fi
echo "spack root: ${spack_root}"

# Locate spack's own openssl, which is the copy we consolidate onto.
spack_ssl="$(find "${spack_root}/opt/spack" -path '*/openssl-*/lib/libssl.*.dylib' \
    -type f 2>/dev/null | head -1)"
if [ -z "${spack_ssl:-}" ]; then
    echo "no spack openssl found; delocating as-is"
    delocate
    exit 0
fi
spack_ssl_dir="$(dirname "$spack_ssl")"
echo "spack openssl dir: ${spack_ssl_dir}"

# libthrift is the only known offender (see root cause above); scope the
# rewrite to it instead of scanning/mutating the whole spack store.
libthrift="$(find "${spack_root}/opt/spack" -name 'libthrift*.dylib' \
    -type f 2>/dev/null | head -1)"
if [ -z "${libthrift:-}" ]; then
    echo "no libthrift in spack store; delocating as-is"
    delocate
    exit 0
fi

echo "::group::Repoint libthrift's framework openssl references to spack openssl"
# `|| true`: grep exits 1 when there is no framework reference (the expected
# steady state once the upstream thrift package is fixed), which would
# otherwise trip `set -e` under pipefail.
fw_refs="$(otool -L "$libthrift" 2>/dev/null | awk 'NR>1{print $1}' \
    | grep -E "$fw_re" || true)"
rewrote=0
if [ -n "$fw_refs" ]; then
    while IFS= read -r fw; do
        target="${spack_ssl_dir}/$(basename "$fw")"
        if [ ! -f "$target" ]; then
            echo "ERROR: ${libthrift} references ${fw}, but ${target} does not exist" >&2
            exit 1
        fi
        echo "  ${libthrift}: ${fw} -> ${target}"
        if ! install_name_tool -change "$fw" "$target" "$libthrift"; then
            echo "ERROR: install_name_tool failed to rewrite ${libthrift}" >&2
            exit 1
        fi
        rewrote=$((rewrote + 1))
    done <<< "$fw_refs"
fi
echo "rewrote ${rewrote} framework openssl reference(s) in libthrift"
echo "::endgroup::"

delocate
