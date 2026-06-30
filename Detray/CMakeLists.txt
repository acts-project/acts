# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set up the project.
cmake_minimum_required(VERSION 3.21)
project(detray VERSION 0.111.0 LANGUAGES CXX)

# Set up the used C++ standard(s).
set(CMAKE_CXX_STANDARD 20 CACHE STRING "The (host) C++ standard to use")
set(CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Disable (host) C++ extensions")
set(CMAKE_CUDA_STANDARD 20 CACHE STRING "The (CUDA) C++ standard to use")
set(CMAKE_CUDA_EXTENSIONS FALSE CACHE BOOL "Disable (CUDA) C++ extensions")
set(CMAKE_SYCL_STANDARD 20 CACHE STRING "The (SYCL) C++ standard to use")
set(CMAKE_HIP_STANDARD 20 CACHE STRING "The (HIP) C++ standard to use")
set(CMAKE_HIP_EXTENSIONS FALSE CACHE BOOL "Disable (HIP) C++ extensions")

if(PROJECT_IS_TOP_LEVEL)
    if(${CMAKE_CXX_STANDARD} LESS 20)
        message(
            SEND_ERROR
            "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}, but detray requires C++>=20"
        )
    endif()
endif()

# Explicitly set the output directory for the binaries. Such that if this
# project is included by another project, the main project's configuration would
# win out.
include(GNUInstallDirs)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}"
    CACHE PATH
    "Directory for the built binaries"
)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}"
    CACHE PATH
    "Directory for the built libraries"
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}"
    CACHE PATH
    "Directory for the built static libraries"
)

set(DETRAY_PYTHON_INSTALL_DIR "python/detray")

# Include the Detray CMake code.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(detray-functions)

#
# Define build options
#
include(CMakeDependentOption)

# General
option(DETRAY_USE_SYSTEM_LIBS "Use system libraries be default" FALSE)
option(
    DETRAY_FAIL_ON_WARNINGS
    "Make the build fail on compiler/linker warnings"
    FALSE
)
set(DETRAY_SET_LOGGING "INFO" CACHE STRING "Set the log level")
set_property(
    CACHE DETRAY_SET_LOGGING
    PROPERTY STRINGS "NONE" "WARNING" "INFO" "VERBOSE" "DEBUG"
)
option(DETRAY_BUILD_CLI_TOOLS "Build the command line tools of Detray" OFF)
# Which detector metadata headers to generate during the build
set(DETRAY_GENERATE_METADATA
    ""
    CACHE STRING
    "Semicolon-separated list of metadata generator scripts to run in this build"
)

# Device compilation

# Check if CUDA is available.
include(CheckLanguage)
check_language(CUDA)
option(
    DETRAY_BUILD_CUDA
    "Build the CUDA sources included in detray"
    ${CMAKE_CUDA_COMPILER}
)

option(DETRAY_BUILD_SYCL "Build the SYCL sources included in detray" OFF)

# Check if HIP is available
option(DETRAY_BUILD_HIP "Build the HIP sources included in detray" OFF)

cmake_dependent_option(
    DETRAY_BUILD_HOST
    "Build the host sources included in detray"
    ON
    "DETRAY_BUILD_CUDA OR DETRAY_BUILD_SYCL OR DETRAY_BUILD_HIP"
    ON
)

# Plugins / Algebra backends
option(DETRAY_SVG_DISPLAY "Build ActSVG display module" OFF)
option(DETRAY_ARRAY_PLUGIN "Build Eigen math plugin" ON)
option(DETRAY_EIGEN_PLUGIN "Build Eigen math plugin" OFF)
option(DETRAY_FASTOR_PLUGIN "Build Fastor math plugin" OFF)
option(DETRAY_SMATRIX_PLUGIN "Build ROOT/SMatrix math plugin" OFF)
option(DETRAY_VC_AOS_PLUGIN "Build Vc based AoS math plugin" OFF)
option(DETRAY_VC_SOA_PLUGIN "Build Vc based SoA math plugin" OFF)

# Python related options
option(
    DETRAY_BUILD_PYTHON_BINDINGS
    "Build python bindings for all enabled components"
    OFF
)

# Test and benchmark suite
include(CTest)

option(DETRAY_BUILD_ALL_TESTS "Build unit and integrations tests of Detray" OFF)
option(DETRAY_BUILD_UNITTESTS "Build the unit tests of Detray" OFF)
option(
    DETRAY_BUILD_INTEGRATIONTESTS
    "Build the integration tests of Detray"
    OFF
)
option(DETRAY_BUILD_TEST_UTILS "Build the test utility library of Detray" OFF)
option(DETRAY_BUILD_VALIDATION_TOOLS "Build the validation tools of Detray" OFF)
option(DETRAY_BUILD_BENCHMARKS "Build the benchmark tests" OFF)
option(DETRAY_BUILD_TUTORIALS "Build the tutorial executables of Detray" OFF)

#
# Resolve build options and option inter-dependencies
#

# Set the internal log level from user input
set(DETRAY_LOG_LVL -1 CACHE INTERNAL "Disable logging")
if(DETRAY_SET_LOGGING STREQUAL "WARN")
    set(DETRAY_LOG_LVL 0 CACHE INTERNAL "Print warnings and errors" FORCE)
elseif(DETRAY_SET_LOGGING STREQUAL "INFO")
    set(DETRAY_LOG_LVL 1 CACHE INTERNAL "Print general information" FORCE)
elseif(DETRAY_SET_LOGGING STREQUAL "VERBOSE")
    set(DETRAY_LOG_LVL 2 CACHE INTERNAL "Print detailed information" FORCE)
elseif(DETRAY_SET_LOGGING STREQUAL "DEBUG")
    set(DETRAY_LOG_LVL 3 CACHE INTERNAL "Print expert information" FORCE)
endif()

# Need test utils in CLI tools for the example detector generation
if(DETRAY_BUILD_CLI_TOOLS)
    set(DETRAY_BUILD_TEST_UTILS ON)
endif()

# Check CUDA and SYCL and HIP C++ standards
if(${DETRAY_BUILD_CUDA} AND ${CMAKE_CUDA_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_CUDA_STANDARD=${CMAKE_CUDA_STANDARD}, but detray requires C++>=20"
    )
endif()

if(${DETRAY_BUILD_SYCL} AND ${CMAKE_SYCL_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_SYCL_STANDARD=${CMAKE_SYCL_STANDARD}, but detray requires C++>=20"
    )
endif()

if(${DETRAY_BUILD_HIP} AND ${CMAKE_HIP_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_HIP_STANDARD=${CMAKE_HIP_STANDARD}, but detray requires C++>=20"
    )
endif()

# Disable the [[no_unique_address]] annotation if we are on CUDA>=13.0.
# TODO: Update this statement once a fix in CUDA is available.
set(DETRAY_INTERNAL_USE_NO_UNIQUE_ADDRESS_ANNOTATION ON)

if(DETRAY_BUILD_CUDA)
    # Figure out the version of CUDA being used.
    find_package(CUDAToolkit REQUIRED)

    if(CUDAToolkit_VERSION_MAJOR GREATER_EQUAL 13)
        message(
            STATUS
            "Disabling [[no_unique_address]] annotation due to CUDA incompatibility"
        )
        set(DETRAY_INTERNAL_USE_NO_UNIQUE_ADDRESS_ANNOTATION OFF)
    endif()
endif()

# Convenience option to build tests
if(DETRAY_BUILD_ALL_TESTS)
    set(DETRAY_BUILD_UNITTESTS ON)
    set(DETRAY_BUILD_INTEGRATIONTESTS ON)
endif()

# Alias flag to enable tests in detray (pulls in google test and triggers the build of executables that depend on it)
if(BUILD_TESTING AND (DETRAY_BUILD_UNITTESTS OR DETRAY_BUILD_INTEGRATIONTESTS))
    set(DETRAY_BUILD_TESTING ON)
endif()

# Validation utilities are needed for integration tests and tutorials
if(DETRAY_BUILD_TESTING OR DETRAY_BUILD_TUTORIALS)
    set(DETRAY_BUILD_VALIDATION_TOOLS ON)
endif()

# Test utilities are needed for validation utilities
if(DETRAY_BUILD_TESTING OR DETRAY_BUILD_VALIDATION_TOOLS)
    set(DETRAY_BUILD_TEST_UTILS ON)
endif()

# Svg display is needed for the validation utilities
if(DETRAY_BUILD_VALIDATION_TOOLS)
    set(DETRAY_SVG_DISPLAY ON)
endif()

# Allow the compilation of test with code sanitizers
cmake_dependent_option(
    DETRAY_ENABLE_SANITIZER
    "Compile tests with sanitizers"
    OFF
    "BUILD_TESTING AND DETRAY_BUILD_TESTING"
    OFF
)

#
# Configure dependencies
#

# Set up VecMem.
option(DETRAY_SETUP_VECMEM "Set up the VecMem target(s) explicitly" TRUE)
option(
    DETRAY_USE_SYSTEM_VECMEM
    "Pick up an existing installation of VecMem from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_VECMEM)
    set(DETRAY_VECMEM_VERSION 1.23.0)
    set(DETRAY_VECMEM_DIGEST
        c4b87f88856b5273b89373d880bdc0fd
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_VECMEM)
        find_package(vecmem ${DETRAY_VECMEM_VERSION} REQUIRED)
    else()
        add_subdirectory(extern/vecmem)
    endif()
    # Make the "VecMem language code" available for the whole project.
    list(APPEND CMAKE_MODULE_PATH "${VECMEM_LANGUAGE_DIR}")
endif()

# Set up JSON for I/O
option(
    DETRAY_SETUP_NLOHMANN
    "Set up the nlohmann::json target(s) explicitly"
    TRUE
)
option(
    DETRAY_USE_SYSTEM_NLOHMANN
    "Pick up an existing installation of nlohman::json from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_NLOHMANN)
    set(DETRAY_NLOHMANN_VERSION 3.11.3)
    set(DETRAY_NLOHMANN_DIGEST
        2074caa675f8097d9b03c0f4976ffc3410170937
        CACHE STRING
        "SHA1 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_NLOHMANN)
        find_package(nlohmann_json ${DETRAY_NLOHMANN_VERSION} REQUIRED)
    else()
        add_subdirectory(extern/nlohmann_json)
    endif()
endif()

# Set up Eigen3.
option(
    DETRAY_SETUP_EIGEN3
    "Set up the Eigen3 target(s) explicitly"
    ${DETRAY_EIGEN_PLUGIN}
)
option(
    DETRAY_USE_SYSTEM_EIGEN3
    "Pick up an existing installation of Eigen3 from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_EIGEN3)
    set(DETRAY_EIGEN_VERSION 3.4.0)
    set(DETRAY_EIGEN_DIGEST
        132dde48fe2b563211675626d29f1707
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_EIGEN3)
        find_package(Eigen3 REQUIRED)
    else()
        add_subdirectory(extern/eigen3)
    endif()
endif()

# Set up Fastor.
option(
    DETRAY_SETUP_FASTOR
    "Set up the Fastor target(s) explicitly"
    ${DETRAY_FASTOR_PLUGIN}
)
option(
    DETRAY_USE_SYSTEM_FASTOR
    "Pick up an existing installation of Fastor from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_FASTOR)
    set(DETRAY_FASTOR_VERSION 0.6.4)
    set(DETRAY_FASTOR_DIGEST
        0644a0bf6337e0ffc9e0173be757aa45
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_FASTOR)
        find_package(Fastor ${DETRAY_FASTOR_VERSION} REQUIRED)
    else()
        add_subdirectory(extern/fastor)
    endif()
endif()

# Set up Vc.
if(DETRAY_VC_SOA_PLUGIN)
    set(DETRAY_VC_AOS_PLUGIN ON)
endif()

if(DETRAY_VC_AOS_PLUGIN)
    set(DETRAY_VC_PLUGIN ON)
endif()

option(DETRAY_SETUP_VC "Set up the Vc target(s) explicitly" ${DETRAY_VC_PLUGIN})
option(
    DETRAY_USE_SYSTEM_VC
    "Pick up an existing installation of Vc from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_VC)
    set(DETRAY_VC_VERSION 1.4.5)
    set(DETRAY_VC_DIGEST
        03831cbf0921a10322d8baf08001cbf5
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_VC)
        find_package(Vc ${DETRAY_VC_VERSION} REQUIRED)
    else()
        add_subdirectory(extern/vc)
    endif()
endif()

# Set up ACTSVG for displaying
option(
    DETRAY_SETUP_ACTSVG
    "Set up the actsvg target(s) explicitly"
    ${DETRAY_SVG_DISPLAY}
)
option(
    DETRAY_USE_SYSTEM_ACTSVG
    "Pick up an existing installation of actsvg from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_ACTSVG)
    set(DETRAY_ACTSVG_VERSION 0.4.57)
    set(DETRAY_ACTSVG_DIGEST
        5136e50ff2d78f0b1ff83a56be6a1954
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_ACTSVG)
        find_package(
            actsvg
            ${DETRAY_ACTSVG_VERSION}
            REQUIRED
            COMPONENTS core meta
        )
    else()
        add_subdirectory(extern/actsvg)
    endif()
endif()

# Set up GoogleTest.
option(
    DETRAY_SETUP_GOOGLETEST
    "Set up the GoogleTest target(s) explicitly"
    ${DETRAY_BUILD_TESTING}
)
option(
    DETRAY_USE_SYSTEM_GOOGLETEST
    "Pick up an existing installation of GoogleTest from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_GOOGLETEST)
    set(DETRAY_GOOGLETEST_VERSION 1.17.0)
    set(DETRAY_GOOGLETEST_DIGEST
        b6f100bc2a5853a48046aa168ececf84
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_GOOGLETEST)
        find_package(GTest REQUIRED)
    else()
        add_subdirectory(extern/googletest)
    endif()
endif()

# Set up Google Benchmark.
option(
    DETRAY_SETUP_BENCHMARK
    "Set up the Google Benchmark target(s) explicitly"
    ${DETRAY_BUILD_BENCHMARKS}
)
option(
    DETRAY_USE_SYSTEM_BENCHMARK
    "Pick up an existing installation of Google Benchmark from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_BENCHMARK)
    set(DETRAY_BENCHMARK_VERSION 1.9.5)
    set(DETRAY_BENCHMARK_DIGEST
        12c6c0c228fc07106c62634222bd2541
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_BENCHMARK)
        find_package(benchmark ${DETRAY_BENCHMARK_VERSION} REQUIRED)
    else()
        add_subdirectory(extern/benchmark)
    endif()
endif()

# Set up covfie.
if(DETRAY_BUILD_TESTING OR DETRAY_BUILD_BENCHMARKS OR DETRAY_BUILD_CLI_TOOLS)
    set(DETRAY_BUILD_COVFIE ON)
endif()
option(
    DETRAY_SETUP_COVFIE
    "Set up the covfie target(s) explicitly"
    ${DETRAY_BUILD_COVFIE}
)
option(
    DETRAY_USE_SYSTEM_COVFIE
    "Pick up an existing installation of covfie from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_COVFIE)
    set(DETRAY_COVFIE_VERSION 0.15.4)
    set(DETRAY_COVFIE_DIGEST
        80de9b43644c59d8239fb8cc184c6276
        CACHE STRING
        "MD5 digest of the downloaded zip"
    )

    if(DETRAY_USE_SYSTEM_COVFIE)
        find_package(covfie ${DETRAY_COVFIE_VERSION} REQUIRED)
    else()
        add_subdirectory(extern/covfie)
    endif()
endif()

# Set up pybind11.
option(
    DETRAY_SETUP_PYBIND11
    "Build python bindings for all enabled components"
    ${DETRAY_BUILD_PYTHON_BINDINGS}
)
option(
    DETRAY_USE_SYSTEM_PYBIND11
    "Pick up an existing installation of pybind11 from the build environment"
    ${DETRAY_USE_SYSTEM_LIBS}
)
if(DETRAY_SETUP_PYBIND11)
    set(DETRAY_PYBIND11_VERSION 3.0.1)

    find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
    if(DETRAY_USE_SYSTEM_PYBIND11)
        find_package(pybind11 CONFIG REQUIRED)
    else()
        add_subdirectory(extern/pybind11)
    endif()
endif()

#
# Set up all of the libraries of the project.
#

add_subdirectory(core)
add_subdirectory(detectors)
add_subdirectory(io)
add_subdirectory(plugins)
add_subdirectory(python)

# Test utils and validation tools can also be required standalone
# (e.g. in ACTS detray plugin)
if(
    DETRAY_BUILD_TESTING
    OR DETRAY_BUILD_TEST_UTILS
    OR DETRAY_BUILD_VALIDATION_TOOLS
    OR DETRAY_BUILD_BENCHMARKS
    OR DETRAY_BUILD_CLI_TOOLS
)
    add_subdirectory(tests)
endif()

# Set up the tutorial(s).
if(DETRAY_BUILD_TUTORIALS)
    add_subdirectory(tutorials)
endif()

# Set up the packaging of the project.
include(detray-packaging)
