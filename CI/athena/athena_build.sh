#!/bin/bash

set -e

env_file=$PWD/.env
if [[ -f "$env_file"  ]]; then
    echo "Loading env from: $env_file"
    cat $env_file
    source $env_file
fi


function group {
    if [[ -v GITHUB_ACTIONS  ]]; then
        echo ""
        echo "::group::${@}"
        echo ""
    else
        n=${#1}
        s=$(python3 -c "print(\"-\" * $n)")
        echo ""
        echo $s
        echo "${@}"
        echo $s
        echo ""
    fi
}

group "Setup ACTS"

echo "LCG_RELEASE: ${LCG_RELEASE}"
echo "LCG_PLATFORM: ${LCG_PLATFORM}"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
ATHENA_GIT_REPO=https://gitlab.cern.ch/atlas/athena.git
ATHENA_RELEASE=Athena,master,latest


source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh || true
asetup none,gcc11 || true
lsetup cmake || true
lsetup "views ${LCG_RELEASE} ${LCG_PLATFORM}" || true

CCACHE=$(command -v ccache)
# ln -s $(command -v ccache) $CCACHE

# ls -al

$CCACHE -z

group "Configure ACTS"

cmake -S $PWD -B acts-build \
	-DCMAKE_INSTALL_PREFIX=$PWD/acts-install \
	-DACTS_BUILD_PLUGIN_JSON:BOOL=ON \
	-DCMAKE_CXX_COMPILER_LAUNCHER=$CCACHE \
	-DACTS_BUILD_FATRAS:BOOL=ON

group "Build ACTS"

cmake --build acts-build -- -j $(nproc)
$CCACHE -s

group "Install ACTS"
cmake --install acts-build

group "Setup Athena"

asetup ${ATHENA_RELEASE} || true
nightly=$(basename $ATLAS_RELEASE_BASE)
branch=$(echo $ATLAS_RELEASE_BASE | perl -pe 's/.*\/sw\/(\w*?)_(\w*?)_.*/\1/g')
project=$(echo $ATLAS_RELEASE_BASE | perl -pe 's/.*\/sw\/(\w*?)_(\w*?)_.*/\2/g')
ATHENA_REF="nightly/${branch}/${nightly}"

group "Cloning ${ATHENA_GIT_REPO} @ ${ATHENA_REF}"
git clone ${ATHENA_GIT_REPO} -b ${ATHENA_REF}

$CCACHE -z

export CMAKE_PREFIX_PATH="$PWD/acts-install:$CMAKE_PREFIX_PATH"

group "Configure Athena"

cmake -S athena/Projects/WorkDir -B athena-build \
	-DATLAS_PACKAGE_FILTER_FILE=$PWD/CI/athena/package_filters.txt \
	-DCMAKE_CXX_COMPILER_LAUNCHER=$CCACHE

group "Build Athena"

cmake --build athena-build -- -j $(nproc)
$CCACHE -s

group "Run Athena based tests"

source athena-build/x*/setup.sh
export LD_LIBRARY_PATH=$PWD/acts-install/lib:$LD_LIBRARY_PATH
mkdir run
cd run

ec=0
function runTest {
    group "Running ${@}"
    ${@}
    ec=$(($ec | $?))
}

runTest python3 ../athena/Tracking/Acts/ActsGeometry/test/ActsITkTest.py
runTest ../athena/AtlasTest/CITest/test/ActsKfRefitting.sh
# runTest ../athena/AtlasTest/CITest/test/ActsPersistifyEDM.sh
# runTest ../athena/AtlasTest/CITest/test/ActsValidateClusters.sh
# runTest ../athena/AtlasTest/CITest/test/ActsValidateOrthogonalSeeds.sh
# runTest ../athena/AtlasTest/CITest/test/ActsValidateSeeds.sh
# runTest ../athena/AtlasTest/CITest/test/ActsValidateSpacePoints.sh
# runTest ../athena/AtlasTest/CITest/test/ActsValidateTracks.sh
runTest ../athena/AtlasTest/CITest/test/ActsWorkflow.sh

exit $ec
