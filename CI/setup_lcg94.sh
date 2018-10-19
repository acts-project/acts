# setup LCG release 94 via cvmfs

# determine os release
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env_info.sh

release=LCG_94
platform=x86_64-${ACTS_OS}-${ACTS_COMPILER}
lcg=/cvmfs/sft.cern.ch/lcg/views/${release}/${platform}

source ${lcg}/setup.sh
# extra variables required to build acts
export DD4hep_DIR=${lcg}
