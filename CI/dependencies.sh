#!/bin/bash
set -e
set -u

function run() {
    set -x
    "$@"
    { set +x;  } 2> /dev/null
}

function set_env {
  key="$1"
  value="$2"

  echo "=> ${key}=${value}"

  if [ -n "${GITHUB_ACTIONS:-}" ]; then
    echo "${key}=${value}" >> $GITHUB_ENV
  else
    export ${key}=${value}
  fi
}

url=${1:-${DEPENDENCY_URL:-}}

if [ -n "${GITHUB_ACTIONS:-}" ]; then
    destination="${GITHUB_WORKSPACE}/dependencies"
elif [ -n "${GITLAB_CI:-}" ];then
    destination="${CI_PROJECT_DIR}/dependencies"
else
    destination=${2}
fi

set_env DEPENDENCY_DIR "${destination}"

if [ -z "${url}" ]; then
    echo "url is not set"
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


echo "check geant4 dataset"

# Patch up geant4-config data install script
orig_share=$(${destination}/bin/geant4-config --datasets|head -n1|perl -pe 's|.*?(\/.*)\/share.*|\1|')
echo "Original share: $orig_share"
orig_share_escaped=$(echo $orig_share|perl -pe 's|/|\\/|g')
echo "Original share escaped: $orig_share_escaped"
destination_escaped=$(echo "$destination"|perl -pe 's|/|\\/|g')
echo "Destination escaped: $destination_escaped"
run perl -pi.bak -e "s/$orig_share_escaped/$destination_escaped/g" ${destination}/bin/geant4-config

echo "Check CI mode"
if [ -n "${GITHUB_ACTIONS:-}" ]; then
  echo "Running in GitHub Actions"
  venv="${GITHUB_WORKSPACE}/venv"
fi
echo "not github"

if [ -n "${GITLAB_CI:-}" ];then
  echo "Running in GitLab CI"
  venv="${CI_PROJECT_DIR}/venv"
fi
echo "after gitlab"

if [ -n "${CI:-}" ];then
  run "${destination}/bin/python3" -m venv "${venv}"
  run "${venv}/bin/python3" -m pip install pyyaml jinja2
  set_env PATH "${venv}/bin:${destination}/bin/:${PATH}"
  set_env CMAKE_PREFIX_PATH "${destination}"
  set_env LD_LIBRARY_PATH "${destination}/lib"
  set_env ROOT_INCLUDE_PATH "${destination}/include"
  # Geant4 puts CLHEP in a subdirectory
  set_env ROOT_INCLUDE_PATH "${destination}/include/Geant4"
  # Pythia8 looks for settings in this directory
  set_env PYTHIA8DATA "${destination}/share/Pythia8/xmldoc"
fi
