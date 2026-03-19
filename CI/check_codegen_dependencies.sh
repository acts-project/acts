#!/bin/bash
set -e
set -u

PYTHON_VERSION=3.14

input=$1
dir=$(dirname "$input")
output="$dir/requirements.txt"

uv python install "$PYTHON_VERSION"
uv pip compile \
  --python-version "$PYTHON_VERSION" \
  "$input" \
  -o "$output"
