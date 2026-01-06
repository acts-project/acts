# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building ModuleMapGraph as part of the Acts project")

# Declare where to get traccc from.
FetchContent_Declare(
    ModuleMapGraph
    SYSTEM
    ${ACTS_MODULEMAPGRAPH_SOURCE}
    OVERRIDE_FIND_PACKAGE
)

set(MMG_USE_MPI OFF)
set(MMG_USE_CUDA ON)
set(MMG_WITH_NUMPY_EXPORT OFF)
set(MMG_WITH_PYTORCH_EXPORT OFF)
set(CMAKE_CUDA_STANDARD 20)

# Now set up its build.
FetchContent_MakeAvailable(ModuleMapGraph)
