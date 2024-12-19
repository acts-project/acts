#!/bin/bash
set -e
set -u

if ! command -v spack &> /dev/null; then
    if [ ! -d "spack" ]; then
        echo "Cloning spack"
        git clone -c feature.manyFiles=true https://github.com/spack/spack.git
        pushd spack > /dev/null
        git config user.name github-actions[bot]
        git config user.email 41898282+github-actions[bot]@users.noreply.github.com
        # Apply patch for spack improvements
        curl https://patch-diff.githubusercontent.com/raw/spack/spack/pull/47370.patch | git am
        popd > /dev/null
    else
        echo "Updating spack"
        pushd spack > /dev/null
        git pull --rebase
        popd > /dev/null
    fi

    source "$(pwd)/spack/share/spack/setup-env.sh"
fi
