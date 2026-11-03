#!/bin/bash

set -e
set -u
set -o pipefail

# Assert COLLIDERML_DATA_DIR environment variable is set
if [[ -z "${COLLIDERML_DATA_DIR+x}" ]]; then
  echo "Error: COLLIDERML_DATA_DIR environment variable is not set"
  exit 1
fi

function download {
  URL=$1
  HASH=$2

  tarname=$(basename "$URL")
  curl -SL "$URL" -o "$tarname"
  echo "$HASH $tarname" | sha256sum -c
  tar --no-same-owner -xf "$tarname"
}

mkdir -p "${COLLIDERML_DATA_DIR}"
cd "${COLLIDERML_DATA_DIR}"

download \
  https://acts.web.cern.ch/ci/colliderml_ttbar_pu0_ci_sample_v1.tar \
  942901db5c860cc8e95d620394970b99704593264a35d68c8b4661947c2054c9
