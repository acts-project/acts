# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#
# This script is meant to configure the build/runtime environment of the
# Docker contaners that are used in the project's CI configuration.
#
# Usage: source .github/ci_setup.sh <platform name>
#

# The platform name.
PLATFORM_NAME=$1

# Make sure that the build and CTest would use all available cores.
export CMAKE_BUILD_PARALLEL_LEVEL=`nproc`
export CTEST_PARALLEL_LEVEL=${CMAKE_BUILD_PARALLEL_LEVEL}
export MAKEFLAGS="-j${CMAKE_BUILD_PARALLEL_LEVEL}"

# Set up the correct environment for the SYCL tests.
if [[ "${PLATFORM_NAME}" == *"SYCL"* ]]; then
   if [ -f "/opt/intel/oneapi/setvars.sh" ]; then
      OLD_CPATH=${CPATH}
      source /opt/intel/oneapi/setvars.sh --include-intel-llvm
      export CPATH=${OLD_CPATH}
   fi
fi
