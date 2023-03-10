#!/bin/bash

set -e


function group {
    if [[ -v GITHUB_ACTIONS  ]]; then
        echo "::group::${1}"
    else
        n=${#1}
        s=$(python3 -c "print(\"-\" * $n)")
        echo ""
        echo $s
        echo "${1}"
        echo $s
        echo ""
    fi
}

group "Setup"

LCG_RELEASE=$1
LCG_PLATFORM=$2

echo "LCG_RELEASE: ${LCG_RELEASE}"
echo "LCG_PLATFORM: ${LCG_PLATFORM}"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
ATHENA_GIT_REPO=https://gitlab.cern.ch/atlas/athena.git
ATHENA_RELEASE=Athena,master,latest


source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh || true
asetup none,gcc11 || true
lsetup cmake || true
lsetup "views ${LCG_RELEASE} ${LCG_PLATFORM}" || true

ln -s $(command -v ccache) $PWD/ccache

$PWD/ccache -z

group "Build ACTS"

cmake -S $PWD -B acts-build \
	-DCMAKE_INSTALL_PREFIX=$PWD/acts-install \
	-DACTS_BUILD_PLUGIN_JSON:BOOL=ON \
	-DCMAKE_CXX_COMPILER_LAUNCHER=$PWD/ccache \
	-DACTS_BUILD_FATRAS:BOOL=ON
cmake --build acts-build -- -j $(nproc)
$PWD/ccache -s
cmake --install acts-build

asetup ${ATHENA_RELEASE} || true

nightly=$(basename $ATLAS_RELEASE_BASE)
branch=$(echo $ATLAS_RELEASE_BASE | perl -pe 's/.*\/sw\/(\w*?)_(\w*?)_.*/\1/g')
project=$(echo $ATLAS_RELEASE_BASE | perl -pe 's/.*\/sw\/(\w*?)_(\w*?)_.*/\2/g')
ATHENA_REF="nightly/${branch}/${nightly}"

group "Cloning ${ATHENA_GIT_REPO} @ ${ATHENA_REF}"
git clone ${ATHENA_GIT_REPO} -b ${ATHENA_REF}

$PWD/ccache -z

export CMAKE_PREFIX_PATH="$PWD/acts-install:$CMAKE_PREFIX_PATH"

cat >package_filters.txt <<EOL
+ .*Acts.*
+ .*InnerDetector/InDetRecAlgs/InDetPriVxFinder.*
+ .*InnerDetector/InDetRecTools/InDetRecToolInterfaces.*
+ .*Tracking/TrkExtrapolation/TrkExAlgs.*

- .*
... 
EOL

group "Build Athena"

cmake -S athena/Projects/WorkDir -B athena-build \
	-DATLAS_PACKAGE_FILTER_FILE=$PWD/package_filters.txt \
	-DCMAKE_CXX_COMPILER_LAUNCHER=$PWD/ccache


cmake --build athena-build -- -j $(nproc)
$PWD/ccache -s

group "Run Athena based tests"

source athena-build/x*/setup.sh
export LD_LIBRARY_PATH=$PWD/acts-install/lib:$LD_LIBRARY_PATH
mkdir run
cd run

function runTest {
    group "Running ${@}"
    ${@}
}

runTest python3 ../athena/Tracking/Acts/ActsGeometry/test/ActsITkTest.py
runTest ../athena/AtlasTest/CITest/test/ActsKfRefitting.sh
runTest ../athena/AtlasTest/CITest/test/ActsPersistifyEDM.sh
runTest ../athena/AtlasTest/CITest/test/ActsValidateClusters.sh
runTest ../athena/AtlasTest/CITest/test/ActsValidateOrthogonalSeeds.sh
runTest ../athena/AtlasTest/CITest/test/ActsValidateSeeds.sh
runTest ../athena/AtlasTest/CITest/test/ActsValidateSpacePoints.sh
runTest ../athena/AtlasTest/CITest/test/ActsValidateTracks.sh
runTest ../athena/AtlasTest/CITest/test/ActsWorkflow.sh

