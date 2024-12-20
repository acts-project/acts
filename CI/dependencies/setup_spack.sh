#!/bin/bash
set -e
set -u

_spack_version=${SPACK_VERSION:-develop}

if [ ! -d "spack" ]; then
    echo "Cloning spack"
    git clone --branch ${_spack_version} -c feature.manyFiles=true https://github.com/acts-project/spack.git
    pushd spack > /dev/null
    git config user.name 'CI'
    git config user.email '<>'

    rm -rf .git

    echo "Populating the repository index"
    bin/spack list > /dev/null
    popd > /dev/null
fi
