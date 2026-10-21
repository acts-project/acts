#!/bin/bash
set -e
set -u

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
_spack_version=${SPACK_VERSION:-$(sed -n 's/^SPACK_VERSION=//p' "${SCRIPT_DIR}/versions.env")}

_spack_folder=$1

if [ ! -d "${_spack_folder}" ]; then
    echo "Cloning spack"
    time git clone -c feature.manyFiles=true https://github.com/spack/spack.git "${_spack_folder}"
    pushd "${_spack_folder}" > /dev/null
    echo "Checking out spack version ${_spack_version}"
    time git checkout "${_spack_version}"
    git config user.name 'CI'
    git config user.email '<>'
    popd > /dev/null
else
  echo "Using cached spack at ${_spack_folder}"
fi
