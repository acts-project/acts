# setup LCG release 95 via cvmfs

# determine os release
if [ -n "$BASH_SOURCE" ]; then
  this_script=$BASH_SOURCE
elif [ -n "$ZSH_VERSION" ]; then
  setopt function_argzero
  this_script=$0
else
  echo 1>&2 "Unsupported shell. Please use bash, or zsh."
  exit 2
fi

dir="$( cd "$( dirname "${this_script}" )" && pwd )"

source $dir/env_info.sh ${this_script}

release=LCG_95
platform=x86_64-${ACTS_OS}-${ACTS_COMPILER}
lcg=/cvmfs/sft.cern.ch/lcg/views/${release}/${platform}

source ${lcg}/setup.sh
# extra variables required to build acts
export DD4hep_DIR=${lcg}
#~ export HepMC3_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-${ACTS_OS}-${ACTS_COMPILER}
