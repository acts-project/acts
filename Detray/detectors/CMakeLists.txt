# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

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

# Create the metadata headers
if(NOT "${DETRAY_GENERATE_METADATA}" STREQUAL "")
    message(STATUS "Configuring detray detector types...")

    find_package(Python COMPONENTS Interpreter REQUIRED)

    set(OUTPUT_DIR
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/detray/detectors"
    )
    file(MAKE_DIRECTORY ${OUTPUT_DIR})

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

        add_custom_target(${GENERATOR}_generation DEPENDS ${GENERATED_HEADER})

        # Set up how clients should find its headers.
        target_include_directories(
            detray_detectors
            INTERFACE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include>
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

        install(
            DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
            OPTIONAL
        )

        add_dependencies(detray_detectors ${GENERATOR}_generation)

        message(STATUS "Creating metadata: ${GENERATOR}")
    endforeach()
endif()
