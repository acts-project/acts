# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

include_guard(GLOBAL)

# Register a metadata family on a provider target. TYPE is a C++ type expression
# with @ALGEBRA@ as its single substitution point (extra parameters are allowed).
function(detray_register_metadata target)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "" "ID;HEADER;TYPE" "")
    if(
        ARG_UNPARSED_ARGUMENTS
        OR ARG_KEYWORDS_MISSING_VALUES
        OR NOT ARG_ID MATCHES "^[a-z][a-z0-9_]*$"
        OR NOT ARG_HEADER
        OR NOT ARG_TYPE
    )
        message(
            FATAL_ERROR
            "detray_register_metadata requires ID, HEADER and TYPE"
        )
    endif()
    if(NOT ARG_TYPE MATCHES "@ALGEBRA@")
        message(FATAL_ERROR "Metadata TYPE must contain @ALGEBRA@")
    endif()
    get_target_property(_ids ${target} DETRAY_METADATA_IDS)
    if(ARG_ID IN_LIST _ids)
        message(FATAL_ERROR "Duplicate metadata ID '${ARG_ID}'")
    endif()
    set_property(
        TARGET ${target}
        APPEND
        PROPERTY DETRAY_METADATA_IDS "${ARG_ID}"
    )
    set_property(
        TARGET ${target}
        PROPERTY DETRAY_METADATA_${ARG_ID}_HEADER "${ARG_HEADER}"
    )
    set_property(
        TARGET ${target}
        PROPERTY DETRAY_METADATA_${ARG_ID}_TYPE "${ARG_TYPE}"
    )
    set_property(
        TARGET ${target}
        APPEND
        PROPERTY
            EXPORT_PROPERTIES
                DETRAY_METADATA_IDS
                DETRAY_METADATA_${ARG_ID}_HEADER
                DETRAY_METADATA_${ARG_ID}_TYPE
    )
endfunction()

# Assemble all providers once, after registrations are complete.
function(detray_generate_metadata_registry target)
    set(_ids "")
    set(DETRAY_REGISTRY_INCLUDES "")
    set(DETRAY_REGISTRY_ALIASES "")
    set(DETRAY_REGISTRY_LIST "")
    foreach(_provider IN ITEMS ${target} ${ARGN})
        if(NOT TARGET ${_provider})
            message(FATAL_ERROR "Unknown metadata provider '${_provider}'")
        endif()
        get_target_property(_provider_ids ${_provider} DETRAY_METADATA_IDS)
        if(NOT _provider_ids)
            message(
                FATAL_ERROR
                "Metadata provider '${_provider}' has no registrations"
            )
        endif()
        foreach(_id IN LISTS _provider_ids)
            if(_id IN_LIST _ids)
                message(
                    FATAL_ERROR
                    "Duplicate metadata ID '${_id}' across providers"
                )
            endif()
            list(APPEND _ids "${_id}")
            get_target_property(
                _header
                ${_provider}
                DETRAY_METADATA_${_id}_HEADER
            )
            get_target_property(_type ${_provider} DETRAY_METADATA_${_id}_TYPE)
            string(REPLACE "@ALGEBRA@" "algebra_t" _type "${_type}")
            string(TOUPPER "${_id}" _macro)
            string(
                APPEND
                DETRAY_REGISTRY_INCLUDES
                "#include <${_header}>\n#define DETRAY_HAS_${_macro}_METADATA\n"
            )
            string(
                APPEND
                DETRAY_REGISTRY_ALIASES
                "template <typename algebra_t> using ${_id}_metadata = ${_type};\n"
            )
            string(APPEND DETRAY_REGISTRY_LIST " X(${_id})")
        endforeach()
        if(NOT _provider STREQUAL target)
            target_link_libraries(${target} INTERFACE ${_provider})
        endif()
    endforeach()
    set_property(
        TARGET ${target}
        PROPERTY DETRAY_AVAILABLE_METADATA_IDS "${_ids}"
    )
    set_property(
        TARGET ${target}
        APPEND
        PROPERTY EXPORT_PROPERTIES DETRAY_AVAILABLE_METADATA_IDS
    )
    set(_include_dir
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include"
    )
    configure_file(
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/detray-metadata-registry.hpp.in"
        "${_include_dir}/detray/detectors/metadata_registry.hpp"
        @ONLY
    )
    target_include_directories(
        ${target}
        INTERFACE $<BUILD_INTERFACE:${_include_dir}>
    )
    install(
        FILES "${_include_dir}/detray/detectors/metadata_registry.hpp"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/detray/detectors"
    )
endfunction()

# Consumers select expensive instantiations independently of availability.
function(detray_validate_metadata_ids target)
    get_target_property(_available ${target} DETRAY_AVAILABLE_METADATA_IDS)
    foreach(_id IN LISTS ARGN)
        if(NOT _id IN_LIST _available)
            message(
                FATAL_ERROR
                "Metadata ID '${_id}' is unavailable. Registered metadata IDs: ${_available}"
            )
        endif()
    endforeach()
endfunction()
