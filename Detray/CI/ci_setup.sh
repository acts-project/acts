# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# This script is meant to configure the build/runtime environment of the
# Docker containers that are used in the project's CI configuration.
#
# Usage: source .github/ci_setup.sh <platform name>
#

# The platform name.
PLATFORM_NAME=$1

# Set up the correct environment for the SYCL tests.
# Do this also for the HIP based tests, since they are using the same image
if [[ "${PLATFORM_NAME}" = "SYCL" ]] || [[ "${PLATFORM_NAME}" = "HIP-AMD" ]]; then
   echo "Setting up oneapi env for ${PLATFORM_NAME}"
   if [[ -f "/opt/intel/oneapi/setvars.sh" ]]; then
      source /opt/intel/oneapi/setvars.sh --include-intel-llvm
   fi
fi

# Make sure that GNU Make and CTest would use all available cores.
export MAKEFLAGS="-j$(nproc)"
export CTEST_PARALLEL_LEVEL=$(nproc)
