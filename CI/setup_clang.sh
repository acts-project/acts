#!/bin/sh -ex
#
# setup LLVM 4.0 compiler via CVMFS CLIC repository

# determine os release
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env_info.sh

source /cvmfs/clicdp.cern.ch/compilers/llvm/6.0.0/x86_64-${ACTS_OS}/setup.sh
hash -r
