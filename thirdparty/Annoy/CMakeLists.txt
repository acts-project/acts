# This file is part of the Acts project.
#
# Copyright (C) 2024 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
include(FetchContent)

# Information message.
message(STATUS "Building Annoy as part of the ACTS project")

# Version of the Annoy library to fetch.
set(ANNOY_VERSION "v${_acts_annoy_version}")

# Declare the fetching of the Annoy library.
FetchContent_Declare(
    Annoy
    ${ACTS_ANNOY_SOURCE}
    PATCH_COMMAND
        "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR> ${GIT_EXECUTABLE} config
        user.name "CI" && "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR>
        ${GIT_EXECUTABLE} config user.email "<>" && "${CMAKE_COMMAND}" -E chdir
        <SOURCE_DIR> ${GIT_EXECUTABLE} am
        ${CMAKE_CURRENT_SOURCE_DIR}/0001-Modify-annoy.patch
)

# Make the fetched content available.
# Annoy contains a CMakeLists.txt file, so it will be added as a subdirectory.
FetchContent_MakeAvailable(Annoy)

target_compile_options(
    Annoy
    INTERFACE "-Wno-old-style-cast" "-Wno-zero-as-null-pointer-constant"
)
