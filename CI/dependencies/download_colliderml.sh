#!/bin/bash

set -e
set -u
set -o pipefail

# Assert COLLIDERML_DATA environment variable is set
if [[ -z "${COLLIDERML_DATA+x}" ]]; then
  echo "Error: COLLIDERML_DATA environment variable is not set"
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

mkdir -p "${COLLIDERML_DATA}"
cd "${COLLIDERML_DATA}"

download \
  https://acts.web.cern.ch/ci/colliderml_ttbar_pu0_ci_sample_v1.tar \
  942901db5c860cc8e95d620394970b99704593264a35d68c8b4661947c2054c9
