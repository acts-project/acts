#!/bin/sh -ex
#
# setup LLVM 4.0 compiler via CVMFS CLIC repository

# determine os release
if [ "$(cat /etc/redhat-release | grep 'Scientific Linux CERN SLC release 6')" ]; then
  os=slc6
elif [ "$(cat /etc/centos-release | grep 'CentOS Linux release 7')" ]; then
  os=centos7
else
  echo "Unknown OS" 1>&2
  exit 1
fi

source /cvmfs/clicdp.cern.ch/compilers/llvm/4.0.1/x86_64-${os}/setup.sh
hash -r
