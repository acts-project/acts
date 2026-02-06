#!/bin/bash
set -e
set -u

_spack_version=${SPACK_VERSION:-develop}

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
