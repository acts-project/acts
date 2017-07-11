#!/bin/sh -ex
#
# setup requirements on Scientific Linux CERN 6 via cvmfs

platform=x86_64-slc6-gcc62-opt
view=/cvmfs/sft.cern.ch/lcg/views/LCG_88/${platform}
dd4hep=/cvmfs/sft.cern.ch/lcg/releases/DD4hep/00-20-b3d88/${platform}
gcc=/cvmfs/sft.cern.ch/lcg/external/gcc/6.2/${platform}

source ${view}/setup.sh
# additional variables that are not set automatically
export GCC_TOOLCHAIN=${gcc}
export PYTHIA8_INCLUDE_DIR="${view}/include"
export PYTHIA8_LIBRARY_DIR="${view}/lib"
# dd4hep config is missing the config file
pushd ${dd4hep} >/dev/null
source bin/thisdd4hep.sh
popd >/dev/null
