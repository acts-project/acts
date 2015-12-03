export GAUDI=$ATLAS/TrackingSW/ats-Gaudi/GAUDI/GAUDI_v25r2
export ATS=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw
export EIGEN=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/a-tracking-sw/Eigen3

if [[ "x$GAUDI" == "x" ]]; then
echo "Need to set the GAUDI environment variable to the path of the Gaudi software directory (contains GaudiKernel/)."
return 1
else
echo "Gaudi   :    $GAUDI"
fi

if [[ "x$EIGEN" == "x" ]]; then
echo "Need to set the EIGEN environment variable to the path of the eigen software directory. containes algebra library."
return 1
else
echo "eigen   :    $EIGEN"
fi

# set up CMake:
export PATH=/afs/cern.ch/sw/lcg/contrib/CMake/2.8.12.2/Linux-i386/bin:$PATH

export CMAKE_PREFIX_PATH=/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/GAUDI/GAUDI_v25r2/cmake:/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/GAUDI/GAUDI_v25r2/InstallArea/x86_64-slc6-gcc48-opt:/afs/cern.ch/work/j/jhrdinka/ATLAS/TrackingSW/ats-Gaudi/Eigen3/lib/cmake/eigen3:/afs/cern.ch/sw/lcg/releases:$PWD/..

export CMTCONFIG=x86_64-slc6-gcc48-opt

# set up the compilers
export PATH=/afs/cern.ch/lhcb/software/releases/LBSCRIPTS/LBSCRIPTS_v8r0/InstallArea/scripts:$PATH

source /afs/cern.ch/sw/lcg/contrib/gcc/4.8.1/x86_64-slc6/setup.sh

