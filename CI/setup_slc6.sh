#!/bin/sh -ex
#
# setup requirements on Scientific Linux CERN 6 via cvmfs

platform=x86_64-slc6-gcc62-opt
view=/cvmfs/sft.cern.ch/lcg/views/LCG_88/${platform}
dd4hep=/cvmfs/sft.cern.ch/lcg/releases/DD4hep/00-20-b3d88/${platform}

source ${view}/setup.sh
# additional variables that are not set automatically
export BOOST_ROOT="${view}"
export EIGEN_INCLUDE_DIR="${view}/include/eigen3"
export PYTHIA8_INCLUDE_DIR="${view}/include"
export PYTHIA8_LIBRARY_DIR="${view}/lib"
# dd4hep config is missing the config file
pushd ${dd4hep} >/dev/null
source bin/thisdd4hep.sh
popd >/dev/null
