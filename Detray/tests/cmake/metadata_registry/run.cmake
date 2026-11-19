# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

file(REMOVE_RECURSE "${TEST_BINARY_DIR}")
execute_process(
    COMMAND
        "${CMAKE_COMMAND}" -S "${CMAKE_CURRENT_LIST_DIR}" -B
        "${TEST_BINARY_DIR}" "-DDETRAY_SOURCE=${DETRAY_SOURCE}" "-DCASE=${CASE}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
)
if(CASE STREQUAL "valid")
    if(NOT result EQUAL 0)
        message(
            FATAL_ERROR
            "Registry configuration failed: ${output}\n${error}"
        )
    endif()
    execute_process(
        COMMAND
            "${CMAKE_COMMAND}" --install "${TEST_BINARY_DIR}" --prefix
            "${TEST_BINARY_DIR}/install"
        RESULT_VARIABLE result
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Registry install failed")
    endif()
    # Consume the exported targets in a separate project. Registration properties
    # must survive installation, as must the provider dependency.
    file(MAKE_DIRECTORY "${TEST_BINARY_DIR}/consumer")
    file(
        WRITE
        "${TEST_BINARY_DIR}/consumer/CMakeLists.txt"
        "
cmake_minimum_required(VERSION 3.25)
project(consumer LANGUAGES NONE)
include(\"${TEST_BINARY_DIR}/install/lib/cmake/providers.cmake\")
include(\"${TEST_BINARY_DIR}/install/lib/cmake/detray-metadata.cmake\")
detray_validate_metadata_ids(detectors itk telescope)
get_target_property(ids detectors DETRAY_AVAILABLE_METADATA_IDS)
get_target_property(type external DETRAY_METADATA_itk_TYPE)
get_target_property(deps detectors INTERFACE_LINK_LIBRARIES)
if(NOT ids STREQUAL \"telescope;itk\" OR
   NOT type STREQUAL \"atlas::tracker<@ALGEBRA@>\" OR NOT deps STREQUAL \"external\")
    message(FATAL_ERROR \"Lost exported registry properties\")
endif()
"
    )
    execute_process(
        COMMAND
            "${CMAKE_COMMAND}" -S "${TEST_BINARY_DIR}/consumer" -B
            "${TEST_BINARY_DIR}/consumer-build"
        RESULT_VARIABLE result
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Installed registry consumption failed")
    endif()
else()
    if(result EQUAL 0 OR NOT error MATCHES "${EXPECTED_ERROR}")
        message(
            FATAL_ERROR
            "Expected '${EXPECTED_ERROR}', got: ${output}\n${error}"
        )
    endif()
endif()
