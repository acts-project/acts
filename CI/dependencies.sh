#!/bin/bash
set -e
set -u


url=${1:-${DEPENDENCY_URL:-}}
destination=${2:-${DEPENDENCY_DIR:-}}

if [ -z "${url}" ]; then
    echo "url is not set"
    exit 1
fi

if [ -z "${destination}" ]; then
    echo "destination is not set"
    exit 1
fi

echo "URL: $url"
echo "DESTINATION: $destination"

# check curl location
CURL=$(command -v curl)
if [ -z "$CURL" ]; then
    echo "curl is not available"
    exit 1
fi

UNZSTD=$(command -v unzstd)
if [ -z "$UNZSTD" ]; then
    echo "unzstd is not available"
    exit 1
fi

TAR=$(command -v tar)
if [ -z "$TAR" ]; then
    echo "tar is not available"
    exit 1
fi

mkdir -p "${destination}"



$CURL \
  --retry 5 \
  --connect-timeout 2 \
  --location $url \
  | unzstd \
  | tar \
    -x \
    --strip-components=1 \
    --directory "${destination}"
