#!/bin/sh -ex
#
# setup LCG release 88 on Scientific Linux CERN 6 via cvmfs

platform=x86_64-slc6
release=LCG_88
lcg=/cvmfs/sft.cern.ch/lcg/views/${release}/${platform}-gcc62-opt
dd4hep=/cvmfs/sft.cern.ch/lcg/releases/${release}/DD4hep/00-20/${platform}-gcc62-opt
gcc=/cvmfs/sft.cern.ch/lcg/releases/${release}/gcc/6.2/${platform}

source ${lcg}/setup.sh
# dd4hep config is missing the config file
pushd ${dd4hep} >/dev/null
source bin/thisdd4hep.sh
popd >/dev/null
# additional variables
export GCC_TOOLCHAIN=${gcc}
