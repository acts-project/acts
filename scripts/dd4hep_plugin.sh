source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6-gcc49-opt/setup.sh
export PATH=/afs/cern.ch/sw/lcg/contrib/CMake/3.5.2/Linux-x86_64/bin/:${PATH}

export DD4hep_DIR=/afs/cern.ch/exp/fcc/sw/0.7/DD4hep/20161003/x86_64-slc6-gcc49-opt
export ROOT_DIR=/afs/cern.ch/sw/lcg/releases/LCG_83/ROOT/6.06.00/x86_64-slc6-gcc49-opt/

source $ROOT_DIR/bin/thisroot.sh

# DD4hep needed for DD4hep Plugin and example - decouple later
export inithere=$PWD
cd  $DD4hep_DIR
source bin/thisdd4hep.sh
cd $inithere

#export LD_LIBRARY_PATH=/afs/cern.ch/work/j/jhrdinka/ACTS/ACTS/a-common-tracking-sw/build/installed/lib64:/afs/cern.ch/user/j/jhrdinka/DD4hep/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=`pwd`/build/installed/lib64:DD4hep_DIR/lib:$LD_LIBRARY_PATH
