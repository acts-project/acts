#!/bin/sh -ex
#
# setup LLVM compiler via CVMFS repository

# determine os release
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env_info.sh

source /cvmfs/sft.cern.ch/lcg/releases/clang/7.0.0-05e9c/x86_64-${ACTS_OS}/setup.sh
hash -r
