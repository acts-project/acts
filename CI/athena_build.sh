#!/bin/bash

set -e

# nightly=$(basename $ATLAS_RELEASE_BASE)
# branch=$(echo $ATLAS_RELEASE_BASE | perl -pe 's/.*\/sw\/(\w*?)_(\w*?)_.*/\1/g')
# project=$(echo $ATLAS_RELEASE_BASE | perl -pe 's/.*\/sw\/(\w*?)_(\w*?)_.*/\2/g')
# tag="nightly/${branch}/${nightly}"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
ACTS_GIT_REPO=https://github.com/acts-project/acts.git
ACTS_BRANCH=main
ATHENA_GIT_REPO=https://gitlab.cern.ch/atlas/athena.git
ATHENA_RELEASE=Athena,master,latest


source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh || true
asetup none,gcc11 || true
lsetup cmake || true
lsetup "views ${LCG_RELEASE} x86_64-centos7-gcc11-opt" || true

