# This file is part of the Acts project.
#
# Copyright (C) 2022 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.11)
include(FetchContent)

message(STATUS "Building FRNN as part of the ACTS project")

set(ACTS_FRNN_GIT_REPOSITORY "https://github.com/lxxue/FRNN"
    CACHE STRING "Git repository to take FRNN from")
set(ACTS_FRNN_GIT_TAG "3e370d8d9073d4e130363faf87d2370598b5fbf2"
    CACHE STRING "Version of FRNN to build")

mark_as_advanced(ACTS_FRNN_GIT_REPOSITORY ACTS_FRNN_GIT_TAG)

FetchContent_Declare(frnncontent
    GIT_REPOSITORY "${ACTS_FRNN_GIT_REPOSITORY}"
    GIT_TAG "${ACTS_FRNN_GIT_TAG}" )

# FRNN does not provide a CMakeLists.txt, so we use a custom one. Because of this,
# we have to implement the populate step manually
FetchContent_GetProperties(frnncontent)
if(NOT frnncontent_POPULATED)
  FetchContent_Populate(frnncontent)
  configure_file(CMakeLists.txt.in "${frnncontent_SOURCE_DIR}/CMakeLists.txt" COPYONLY)
  add_subdirectory(${frnncontent_SOURCE_DIR} ${frnncontent_BINARY_DIR})
endif()


