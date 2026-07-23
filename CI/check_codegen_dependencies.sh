#!/bin/bash
set -e
set -u

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")

# Must match the floor used by .github/workflows/update-pip-requirements.yml,
# otherwise this hook and the weekly workflow would rewrite the same lockfile
# with different pins.
PYTHON_VERSION=$("${SCRIPT_DIR}/python_versions.py" --floor)

input=$1
dir=$(dirname "$input")
output="$dir/requirements.txt"

uv python install "$PYTHON_VERSION"
# --universal keeps environment markers, so the lockfile stays installable on
# every Python from $PYTHON_VERSION upwards rather than just that one version.
uv pip compile \
  --universal \
  --python-version "$PYTHON_VERSION" \
  "$input" \
  -o "$output"
