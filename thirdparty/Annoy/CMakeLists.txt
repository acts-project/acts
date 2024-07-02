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
   PATCH_COMMAND "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR> ${GIT_EXECUTABLE} config user.name "Jeremy Couthures" &&
                  "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR> ${GIT_EXECUTABLE} config user.email "jeremy.couthures@lapp.in2p3.fr" &&
                  "${CMAKE_COMMAND}" -E chdir <SOURCE_DIR> ${GIT_EXECUTABLE} am ${CMAKE_CURRENT_SOURCE_DIR}/0001-Modify-annoy.patch
)

# Make the fetched content available.
# Sets the variables <content-name>_SOURCE_DIR and <content-name>_BINARY_DIR and <content-name>_POPULATED
# in the same way as FetchContent_GetProperties() would.
# Annoy contains a CMakeLists.txt file, so it will be added as a subdirectory.
FetchContent_MakeAvailable(Annoy)

# Print source dir and binary dir
message(STATUS "Annoy source dir: ${annoy_SOURCE_DIR}")
message(STATUS "Annoy binary dir: ${annoy_BINARY_DIR}")

# Interface library for Annoy.
# add_library(Annoy::annoy INTERFACE IMPORTED GLOBAL)
# # Set the include directory to the src directory of the fetched content
# target_include_directories(Annoy::annoy INTERFACE ${annoy_SOURCE_DIR}/src)
# add_dependencies(Annoy::annoy annoy)