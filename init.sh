export ACTS=$PWD
export EIGEN=/afs/cern.ch/sw/lcg/releases/eigen/3.2.7-292e1/x86_64-slc6-gcc49-opt/logs
# CMTCONFIG to define platform + build-type (opt vs dbg)
export CMTCONFIG=x86_64-slc6-gcc49-opt
# CMTPROJECTPATH as additional search path for heptools
export CMTPROJECTPATH=/afs/cern.ch/sw/Gaudi/releases:/afs/cern.ch/sw/lcg/releases:/afs/cern.ch/sw/lcg/app/releases
# set up the compiler
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh
#Gaudi
export CMAKE_PREFIX_PATH=$ATS/cmake:/afs/cern.ch/sw/Gaudi/releases/GAUDI/GAUDI_v27r0/cmake:$CMAKE_PREFIX_PATH
export CMTPROJECTPATH=/afs/cern.ch/exp/fcc/sw/0.5/
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh

# later implement that used only when needed
# Root needed for geoDisplay
source /afs/cern.ch/sw/lcg/releases/LCG_83/ROOT/6.06.00/x86_64-slc6-gcc49-opt/bin/thisroot.sh

# DD4hep needed for DD4hep Plugin and example - decouple later
export inithere=$PWD
cd  /afs/cern.ch/exp/fcc/sw/0.5/DD4hep/20152711/x86_64-slc6-gcc49-opt
source bin/thisdd4hep.sh
cd $inithere
#needed for detector building in DD4hep
export LD_LIBRARY_PATH=$ACTS/build.x86_64-slc6-gcc49-opt/lib:$LD_LIBRARY_PATH



