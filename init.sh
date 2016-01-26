export ATS=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw
export EIGEN=/afs/cern.ch/sw/lcg/releases/eigen/3.2.7-292e1/x86_64-slc6-gcc49-opt/logs
export GAUDI=/afs/cern.ch/sw/Gaudi/releases/GAUDI/GAUDI_v26r4/cmake

source /afs/cern.ch/lhcb/software/releases/LBSCRIPTS/LBSCRIPTS_v8r4p3/InstallArea/scripts/LbLogin.sh --cmtconfig x86_64-slc6-gcc49-opt

export CMAKE_PREFIX_PATH=$GAUDI:/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/cmake

export CMTPROJECTPATH=/afs/cern.ch/exp/fcc/sw/0.5/
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh
