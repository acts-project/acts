#!/bin/bash

# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Build the ACTS source archive that ships with pre-generated code.
#
# The archive holds the tracked content of a revision plus a `prebuilt-codegen/`
# directory, which cmake/ActsCodegen.cmake finds on its own. A build from the
# unpacked archive therefore never runs the code generators, and never downloads
# uv, a Python interpreter or the generator dependencies.
#
# Note that this is not the same as the source archive GitHub generates for a
# tag: that one cannot carry the generated code. Only the archive attached to
# the release as an asset does.
#
# Generate the code first, with CI/pregenerate_codegen.py.

set -e
set -u

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SOURCE_ROOT=$( dirname "$SCRIPT_DIR" )

output=""
prebuilt=""
revision="HEAD"
prefix=""

usage() {
    cat <<EOF
usage: $( basename "$0" ) --output <archive.tar.gz> --prebuilt <dir> [options]

  --output <path>     archive to write
  --prebuilt <dir>    pre-generated code, as produced by CI/pregenerate_codegen.py
  --revision <rev>    revision to take the tracked content from (default: HEAD)
  --prefix <name>     top-level directory inside the archive
                      (default: the output file name without .tar.gz)
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --output) output="$2"; shift 2 ;;
        --prebuilt) prebuilt="$2"; shift 2 ;;
        --revision) revision="$2"; shift 2 ;;
        --prefix) prefix="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "error: unknown argument '$1'" >&2; usage >&2; exit 1 ;;
    esac
done

if [ -z "$output" ] || [ -z "$prebuilt" ]; then
    echo "error: --output and --prebuilt are required" >&2
    usage >&2
    exit 1
fi

if [ ! -d "$prebuilt" ]; then
    echo "error: '$prebuilt' is not a directory" >&2
    exit 1
fi

count=$( find "$prebuilt" -type f | wc -l | tr -d ' ' )
if [ "$count" -eq 0 ]; then
    echo "error: '$prebuilt' holds no pre-generated files" >&2
    exit 1
fi

if [ -z "$prefix" ]; then
    prefix=$( basename "$output" )
    prefix=${prefix%.tar.gz}
    prefix=${prefix%.tgz}
fi

tmp=$( mktemp -d )
trap 'rm -rf "$tmp"' EXIT

# Tracked content only: a plain copy of the source directory would drag in build
# directories and whatever else is lying around, and this also honours
# `export-ignore` attributes. git archive puts it all under $prefix/ for us.
git -C "$SOURCE_ROOT" archive \
    --format=tar \
    --prefix="$prefix/" \
    --output="$tmp/base.tar" \
    "$revision"

# Append the generated code, laid out where ActsCodegen.cmake looks for it.
# Appending rather than unpacking and repacking is why this stays uncompressed
# until the end: tar can only extend an uncompressed archive.
mkdir -p "$tmp/stage/$prefix"
cp -R "$prebuilt" "$tmp/stage/$prefix/prebuilt-codegen"
tar -rf "$tmp/base.tar" -C "$tmp/stage" "$prefix/prebuilt-codegen"

gzip -c "$tmp/base.tar" > "$output"

echo "Wrote $output from $revision with $count pre-generated files"
