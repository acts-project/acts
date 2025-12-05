#!/bin/bash

set -e
set -u
set -o pipefail

# Assert MODEL_STORAGE environment variable is set
if [[ -z "${MODEL_STORAGE+x}" ]]; then
  echo "Error: MODEL_STORAGE environment variable is not set"
  exit 1
fi

function download {
  URL=$1
  HASH=$2

  tarname=$(basename "$URL")
  curl -SL "$URL" -o "$tarname"
  echo "$HASH $tarname" | sha256sum -c
  tar -xf "$tarname"
}

mkdir -p "${MODEL_STORAGE}"
cd "${MODEL_STORAGE}"

download \
  https://acts.web.cern.ch/ci/gnn/onnx_models_v01.tar \
  335d829439d9a5ae99ffc8f0bbf6a62119a140475e779e2ed00de21fdebb3cb4

download \
  https://acts.web.cern.ch/ci/gnn/torchscript_models_v01.tar \
  1185060ce697bbc96c9dc32b85e5f0eb4db1f64a645c0fc4d2cb2731cb2ef3dc

download \
  https://acts.web.cern.ch/ci/gnn/odd_module_map_v01.tar \
  59f0457f0043bac8594e9f5a3310a709244de980a7b0c206d7d0d95f15455d73
