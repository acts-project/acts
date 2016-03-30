export ATS=$PWD
export EIGEN=/afs/cern.ch/sw/lcg/releases/eigen/3.2.7-292e1/x86_64-slc6-gcc49-opt/logs
# CMTCONFIG to define platform + build-type (opt vs dbg)
export CMTCONFIG=x86_64-slc6-gcc49-opt
# CMTPROJECTPATH as additional search path for heptools
export CMTPROJECTPATH=/afs/cern.ch/sw/Gaudi/releases:/afs/cern.ch/sw/lcg/releases:/afs/cern.ch/sw/lcg/app/releases
# set up the compiler
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh

# should you need anything in addition to Eigen and Gaudi, you can set up the LCG stack with this:
# NOTE that there is a bug in the setup script that prevents it to work with zsh.
source /afs/cern.ch/sw/lcg/views/LCG_83/x86_64-slc6-gcc49-opt/setup.sh

# Either you use the cmake dir from LHCb which has all the package configs for heptools:
# export CMAKE_PREFIX_PATH=$ATS/cmake:/afs/cern.ch/lhcb/software/releases/LBSCRIPTS/LBSCRIPTS_v8r5p7/LbUtils/cmake:$CMAKE_PREFIX_PATH
# Or you explicitly set the prefix path for the Gaudi version you are using (has to be synchronised with the main CMakeLists.txt):
export CMAKE_PREFIX_PATH=$ATS/cmake:/afs/cern.ch/sw/Gaudi/releases/GAUDI/GAUDI_v27r0/cmake:$CMAKE_PREFIX_PATH

export CMTPROJECTPATH=/afs/cern.ch/exp/fcc/sw/0.5/
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh

# add DD4hep
export inithere=$PWD
cd  /afs/cern.ch/exp/fcc/sw/0.5/DD4hep/20152711/x86_64-slc6-gcc49-opt
source bin/thisdd4hep.sh
cd $inithere

# export CMTPROJECTPATH=/afs/cern.ch/exp/fcc/sw/0.5/

#needed for detector building in DD4hep
export LD_LIBRARY_PATH=$ATS/build.x86_64-slc6-gcc49-opt/lib:$LD_LIBRARY_PATH
