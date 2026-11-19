# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.25)
project(metadata_registry_test LANGUAGES NONE)
include("${DETRAY_SOURCE}/cmake/detray-metadata.cmake")
set(CMAKE_INSTALL_INCLUDEDIR include)
add_library(detectors INTERFACE)
detray_register_metadata(detectors ID telescope HEADER telescope.hpp
    TYPE "example::telescope<@ALGEBRA@, example::rectangle>"
)
add_library(external INTERFACE)
# The provider's location and C++ spelling need not resemble its stable ID.
detray_register_metadata(external ID itk HEADER atlas/tracker.hpp
    TYPE "atlas::tracker<@ALGEBRA@>"
)
if(CASE STREQUAL "duplicate_local")
    detray_register_metadata(external ID itk HEADER another.hpp TYPE "other<@ALGEBRA@>")
elseif(CASE STREQUAL "duplicate_provider")
    detray_register_metadata(detectors ID itk HEADER another.hpp TYPE "other<@ALGEBRA@>")
elseif(CASE STREQUAL "invalid_id")
    detray_register_metadata(external ID "bad-id" HEADER bad.hpp TYPE "other<@ALGEBRA@>")
elseif(CASE STREQUAL "missing_algebra")
    detray_register_metadata(external ID bad HEADER bad.hpp TYPE "other<float>")
elseif(CASE STREQUAL "missing_provider")
    detray_generate_metadata_registry(detectors nonexistent)
endif()
detray_generate_metadata_registry(detectors external)
if(CASE STREQUAL "unavailable")
    detray_validate_metadata_ids(detectors unknown)
endif()
detray_validate_metadata_ids(detectors itk telescope)
get_target_property(ids detectors DETRAY_AVAILABLE_METADATA_IDS)
if(NOT ids STREQUAL "telescope;itk")
    message(FATAL_ERROR "Unexpected IDs: ${ids}")
endif()
install(TARGETS detectors external EXPORT providers)
install(EXPORT providers DESTINATION lib/cmake)

install(
    FILES
        "${DETRAY_SOURCE}/cmake/detray-metadata.cmake"
        "${DETRAY_SOURCE}/cmake/detray-metadata-registry.hpp.in"
    DESTINATION lib/cmake
)
