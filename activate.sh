acts_source=~/cern/source/acts/acts/dev3
acts_build=~/cern/build/acts/acts/dev3

. ~/cern/install/geant4/11.2.1/bin/geant4.sh
. ~/cern/install/root/6.30.06/bin/thisroot.sh
. ~/cern/install/DD4hep/01-27-02/bin/thisdd4hep.sh

# copied and modified from $acts_build/this_acts.sh
export PATH="${acts_build}/bin:${PATH}"
export LD_LIBRARY_PATH="${acts_build}/lib:${acts_build}/thirdparty/OpenDataDetector/factory:${LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="${acts_build}/lib:${acts_build}/thirdparty/OpenDataDetector/factory:${DYLD_LIBRARY_PATH}"
export DD4HEP_LIBRARY_PATH="${acts_build}/lib:${acts_build}/thirdparty/OpenDataDetector/factory:${DD4HEP_LIBRARY_PATH}"

# python deps
. ~/cern/venv/bin/activate

# acts python
. "${acts_build}/python/setup.sh"
export PYTHONPATH="${acts_source}/Examples/Scripts/Python:${PYTHONPATH}"

# root include path
export ROOT_INCLUDE_PATH="${ROOT_INCLUDE_PATH}":~/cern/install/podio/01-00-01/include
export ROOT_INCLUDE_PATH="${ROOT_INCLUDE_PATH}":~/cern/install/edm4hep/00-10-03/include

# cmake
export CMAKE_PREFIX_PATH=~/cern/install/json/3.11.3:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/tbb/2021.11.0:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/eigen/3.4.0:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/root/6.30.06:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/geant4/11.2.1:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/pythia/8312:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/podio/01-00-01:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/edm4hep/00-10-03:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/dd4hep/01-27-02:"${CMAKE_PREFIX_PATH}"
export CMAKE_PREFIX_PATH=~/cern/install/hepmc3/3.2.7:"${CMAKE_PREFIX_PATH}"

# ATLAS
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
