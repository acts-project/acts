# This file is part of the Acts project.
#
# Copyright (C) 2024 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
include(FetchContent)

# Minimum CMake version.
cmake_minimum_required(VERSION 3.14) # FetchContent_MakeAvailable was introduced in 3.14

# Information message.
message(STATUS "Building Annoy as part of the ACTS project")

# Version of the Annoy library to fetch.
set(ANNOY_VERSION "v${_acts_annoy_version}")

# Git repository for the Annoy library.
set(ACTS_ANNOY_GIT_REPOSITORY "https://github.com/spotify/annoy.git"
   CACHE STRING "Git repository to take Annoy from")

# Declare the fetching of the Annoy library.
FetchContent_Declare(
   Annoy
   GIT_REPOSITORY "${ACTS_ANNOY_GIT_REPOSITORY}"
   GIT_TAG "${ANNOY_VERSION}"
   CONFIGURE_COMMAND ""
   BUILD_COMMAND ""
   INSTALL_COMMAND ""
   PATCH_COMMAND "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR> ${GIT_EXECUTABLE} config user.name "CI" &&
                  "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR> ${GIT_EXECUTABLE} config user.email "<>" &&
                  "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR> ${GIT_EXECUTABLE} am ${CMAKE_CURRENT_SOURCE_DIR}/0001-Modify-annoy.patch
)

# Make the fetched content available.
# Annoy contains a CMakeLists.txt file, so it will be added as a subdirectory.
FetchContent_MakeAvailable(Annoy)
