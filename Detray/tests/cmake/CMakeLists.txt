# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

foreach(
    case
    IN
    ITEMS
        valid
        duplicate_local
        duplicate_provider
        invalid_id
        missing_algebra
        missing_provider
        unavailable
)
    set(expected "")
    if(case MATCHES "^duplicate")
        set(expected "Duplicate metadata ID")
    elseif(case STREQUAL "invalid_id")
        set(expected "requires ID, HEADER and TYPE")
    elseif(case STREQUAL "missing_algebra")
        set(expected "TYPE must contain")
    elseif(case STREQUAL "unavailable")
        set(expected "is unavailable")
    elseif(case STREQUAL "missing_provider")
        set(expected "Unknown metadata provider")
    endif()
    add_test(
        NAME detray_metadata_registry_${case}
        COMMAND
            "${CMAKE_COMMAND}"
            "-DDETRAY_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}/../.."
            "-DTEST_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/${case}"
            "-DCASE=${case}" "-DEXPECTED_ERROR=${expected}" -P
            "${CMAKE_CURRENT_SOURCE_DIR}/metadata_registry/run.cmake"
    )
endforeach()
