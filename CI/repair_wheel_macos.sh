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
# copy that delocate can vendor. Both are openssl 3.x, and ACTS' Parquet usage
# does not touch TLS. The proper fix belongs in ci-dependencies (build
# thrift/arrow against spack's openssl, or without openssl at all); this keeps
# the wheels building until then, and becomes an automatic no-op once that lands.
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

echo "::group::Repoint framework openssl references to spack openssl"
rewrote=0
while IFS= read -r lib; do
    # `|| true`: grep exits 1 for the (vast majority of) libraries with no
    # framework reference, which would otherwise trip `set -e` under pipefail.
    fw_refs="$(otool -L "$lib" 2>/dev/null | awk 'NR>1{print $1}' \
        | grep -E "$fw_re" || true)"
    [ -n "$fw_refs" ] || continue
    while IFS= read -r fw; do
        target="${spack_ssl_dir}/$(basename "$fw")"
        if [ ! -f "$target" ]; then
            echo "ERROR: ${lib} references ${fw}, but ${target} does not exist" >&2
            exit 1
        fi
        echo "  ${lib}: ${fw} -> ${target}"
        if ! install_name_tool -change "$fw" "$target" "$lib"; then
            echo "ERROR: install_name_tool failed to rewrite ${lib}" >&2
            exit 1
        fi
        rewrote=$((rewrote + 1))
    done <<< "$fw_refs"
done < <(find "${spack_root}/opt/spack" -name '*.dylib' -type f)
echo "rewrote ${rewrote} framework openssl reference(s)"
echo "::endgroup::"

delocate
