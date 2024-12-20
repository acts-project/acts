#!/bin/bash
set -e
set -u

_spack_version=${SPACK_VERSION:-develop}

if [ ! -d "spack" ]; then
    echo "Cloning spack"
    git clone --branch ${_spack_version} -c feature.manyFiles=true https://github.com/spack/spack.git
    pushd spack > /dev/null
    git config user.name 'CI'
    git config user.email '<>'
    echo "Applying patch for spack improvements"
    curl https://patch-diff.githubusercontent.com/raw/spack/spack/pull/47370.patch | git am
    curl https://patch-diff.githubusercontent.com/raw/spack/spack/pull/48236.patch | git am

    rm -rf .git

    echo "Populating the repository index"
    bin/spack list > /dev/null
    popd > /dev/null
fi
