#!/bin/sh -ex
#
# setup LLVM compiler via CVMFS repository

# determine os release
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env_info.sh

_llvm_version=${ACTS_LLVM_VERSION:-9.0.0}

source /cvmfs/sft.cern.ch/lcg/releases/clang/${_llvm_version}-*/x86_64-${ACTS_OS}/setup.sh
hash -r

echo "CXX: $CXX"
