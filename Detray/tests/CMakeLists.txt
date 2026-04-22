# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the common C++ flags.
include(detray-compiler-options-cpp)

include_directories(
    SYSTEM
    $<TARGET_PROPERTY:covfie::core,INTERFACE_INCLUDE_DIRECTORIES>
)

# Detray unit and integration test framework (depends on google test)
if(DETRAY_BUILD_TESTING)
    message(STATUS "Building 'detray::test_framework' component")

    file(
        GLOB _detray_test_headers
        RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
        "include/detray/test/framework/*.hpp"
    )

    add_library(detray_test_framework INTERFACE "${_detray_test_headers}")
    add_library(detray::test_framework ALIAS detray_test_framework)

    target_include_directories(
        detray_test_framework
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )

    target_link_libraries(
        detray_test_framework
        INTERFACE GTest::gtest vecmem::core detray::core detray::detectors
    )

    unset(_detray_test_headers)
endif()

# Common test dependencies such as test detectors, B-fields and track generators
# (depends on covfie)
if(DETRAY_BUILD_TEST_UTILS OR DETRAY_BUILD_BENCHMARKS)
    message(STATUS "Building 'detray::test_common' component")

    file(
        GLOB _detray_test_common_headers
        RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
        "include/detray/test/common/*.hpp"
        "include/detray/test/common/factories/*.hpp"
        "include/detray/test/common/event_generator/*.hpp"
    )

    detray_add_library(detray_test_common test_common "${_detray_test_common_headers}")

    target_link_libraries(
        detray_test_common
        INTERFACE vecmem::core covfie::core detray::core detray::detectors
    )

    unset(_detray_test_common_headers)
endif()

# Lightweight test utils library without external dependencies
if(DETRAY_BUILD_TEST_UTILS OR DETRAY_BUILD_VALIDATION_TOOLS)
    message(STATUS "Building 'detray::test_utils' component")

    file(
        GLOB _detray_test_utils_headers
        RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
        "include/detray/test/utils/*.hpp"
    )

    detray_add_library(detray_test_utils test_utils ${_detray_test_utils_headers})

    target_link_libraries(detray_test_utils INTERFACE vecmem::core detray::core)

    unset(_detray_test_utils_headers)
endif()

# Build extended validation tools (actsvg and csv dependencies)
if(DETRAY_BUILD_VALIDATION_TOOLS)
    message(STATUS "Building 'detray::validation_utils' component")

    file(
        GLOB _detray_validation_headers
        RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
        "include/detray/test/validation/*.hpp"
    )

    detray_add_library(detray_validation_utils validation_utils ${_detray_validation_headers})

    target_link_libraries(
        detray_validation_utils
        INTERFACE
            vecmem::core
            detray::core
            detray::test_utils
            detray::csv_io
            detray::svgtools
    )

    unset(_detray_validation_headers)
endif()

# Include all of the code-holding sub-directories.

# Build the benchmarks
if(DETRAY_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# Build test executables
if(DETRAY_BUILD_TESTING)
    # Build host and device specific test code
    add_subdirectory(include/detray/test/cpu)

    if(DETRAY_BUILD_CUDA OR DETRAY_BUILD_HIP OR DETRAY_BUILD_SYCL)
        add_subdirectory(include/detray/test/device)
    endif()

    # Build the unittests
    if(DETRAY_BUILD_UNITTESTS)
        add_subdirectory(unit_tests)
    endif()

    # Build the integration tests
    if(DETRAY_BUILD_INTEGRATIONTESTS)
        add_subdirectory(integration_tests)
    endif()
endif()

# Build the command line tools (also required for the host jacobian
# validation integration test)
if(
    (DETRAY_BUILD_INTEGRATIONTESTS AND DETRAY_BUILD_HOST)
    OR DETRAY_BUILD_CLI_TOOLS
)
    add_subdirectory(tools)
endif()
