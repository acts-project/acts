# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building 'detray::detectors' component")

# Set the common C++ flags.
include(detray-compiler-options-cpp)

# Set up a detector library, which contains the different detector types that
# are currently available by default.
file(
    GLOB _detray_detector_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/detectors/*.hpp"
)

detray_add_library( detray_detectors detectors ${_detray_detectors_headers})
target_link_libraries(
    detray_detectors
    INTERFACE vecmem::core detray::core detray::io
)

unset(_detray_detector_headers)

# Compile all headers one-by-one.
detray_compile_headers(
    detectors
    detray::detectors
    "include/detray/detectors/*.hpp"
    "include/detray/detectors/**/*.hpp"
)
if(TARGET detray_detectors_HEADERS)
    # The detector headers need an algebra plugin to compile.
    target_link_libraries(
        detray_detectors_HEADERS
        PRIVATE detray::core_array detray::io_array
    )
endif()

# Location of the generated metadata headers. This directory always exists, so
# that clients can always include the aggregate "generated_metadata.hpp" header.
set(DETRAY_GENERATED_METADATA_DIR
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include"
)
set(OUTPUT_DIR "${DETRAY_GENERATED_METADATA_DIR}/detray/detectors")
file(MAKE_DIRECTORY ${OUTPUT_DIR})
message(STATUS "Detector metadata header location: ${OUTPUT_DIR}")

# Set up how clients should find the generated headers.
target_include_directories(
    detray_detectors
    INTERFACE
        $<BUILD_INTERFACE:${DETRAY_GENERATED_METADATA_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(
    DIRECTORY "${DETRAY_GENERATED_METADATA_DIR}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    OPTIONAL
)

# Content of the aggregate header: one include and one feature macro per
# generated metadata header.
set(DETRAY_GENERATED_METADATA_CONTENT "")
# X-macro list of all generated metadata names, e.g. " X(itk)"
set(DETRAY_GENERATED_METADATA_LIST "")

# Create the metadata headers
if(NOT "${DETRAY_GENERATE_METADATA}" STREQUAL "")
    message(STATUS "Configuring detray detector types...")

    find_package(Python COMPONENTS Interpreter REQUIRED)

    # Generate the custom metadata headers
    foreach(GENERATOR ${DETRAY_GENERATE_METADATA})
        set(GENERATED_HEADER "${OUTPUT_DIR}/${GENERATOR}.hpp")
        set(GENERATOR_SCRIPT
            "${PROJECT_SOURCE_DIR}/detectors/python/${GENERATOR}.py"
        )

        add_custom_command(
            OUTPUT "${GENERATED_HEADER}"
            COMMAND
                ${CMAKE_CURRENT_SOURCE_DIR}/generate_metadata.sh -m
                ${GENERATOR_SCRIPT} -o ${OUTPUT_DIR} -p ${CMAKE_BINARY_DIR}
            DEPENDS
                "${GENERATOR_SCRIPT}"
                "${PROJECT_SOURCE_DIR}/detectors/generate_metadata.sh"
                "${CMAKE_BINARY_DIR}/python/setup.sh"
        )

        add_custom_target(
            ${GENERATOR}_generation
            ALL
            DEPENDS ${GENERATED_HEADER}
        )

        add_dependencies(detray_detectors ${GENERATOR}_generation)

        # Register the header in the aggregate header, e.g.
        # "<name>_metadata" -> DETRAY_HAS_<name>_METADATA and X(name)
        string(MAKE_C_IDENTIFIER "${GENERATOR}" _detray_md_macro)
        string(TOUPPER "${_detray_md_macro}" _detray_md_macro)
        string(
            APPEND DETRAY_GENERATED_METADATA_CONTENT
            "#define DETRAY_HAS_${_detray_md_macro}\n"
        )
        if(
            EXISTS
                "${CMAKE_CURRENT_SOURCE_DIR}/include/detray/detectors/${GENERATOR}.hpp"
        )
            # A static header with the same name exists: Do not include the
            # generated header, because it defines the same metadata type again
            message(
                STATUS
                "Static header for '${GENERATOR}' exists: "
                "Not added to 'generated_metadata.hpp'"
            )
        else()
            string(REGEX REPLACE "_metadata$" "" _detray_md_name "${GENERATOR}")
            string(
                APPEND DETRAY_GENERATED_METADATA_CONTENT
                "#include \"detray/detectors/${GENERATOR}.hpp\"\n"
            )
            string(
                APPEND DETRAY_GENERATED_METADATA_LIST
                " X(${_detray_md_name})"
            )
            unset(_detray_md_name)
        endif()
        unset(_detray_md_macro)

        message(STATUS "Creating metadata: ${GENERATOR}")
    endforeach()
endif()

# Always create the aggregate header (it is empty if no metadata is generated)
string(
    APPEND DETRAY_GENERATED_METADATA_CONTENT
    "#define DETRAY_FOR_EACH_GENERATED_METADATA(X)${DETRAY_GENERATED_METADATA_LIST}\n"
)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/detray/detectors/generated_metadata.hpp.in"
    "${OUTPUT_DIR}/generated_metadata.hpp"
    @ONLY
)
