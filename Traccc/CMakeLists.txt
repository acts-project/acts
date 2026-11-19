# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set up the project.
cmake_minimum_required(VERSION 3.25)
project(traccc VERSION 1.6.0 LANGUAGES CXX)

# Standard CMake include(s).
include(GNUInstallDirs)

# Nothing reads TRACCC_CUSTOM_SCALARTYPE, so a value for it would be ignored.
if(DEFINED TRACCC_CUSTOM_SCALARTYPE)
    message(
        FATAL_ERROR
        "TRACCC_CUSTOM_SCALARTYPE has no effect. The scalar type of traccc "
        "is set with DETRAY_CUSTOM_SCALARTYPE."
    )
endif()

if(PROJECT_IS_TOP_LEVEL)
    message(
        DEPRECATION
        "Building traccc as a standalone project is deprecated, and support "
        "for it will be removed in a future version. Build traccc as part of "
        "ACTS instead, with ACTS_BUILD_TRACCC=ON and the ACTS_TRACCC_* options."
    )

    # Set up the used C++ standard(s).
    set(CMAKE_CXX_STANDARD 23 CACHE STRING "The (host) C++ standard to use")
    set(CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Disable (host) C++ extensions")
    set(CMAKE_CUDA_STANDARD 20 CACHE STRING "The (CUDA) C++ standard to use")
    set(CMAKE_CUDA_EXTENSIONS FALSE CACHE BOOL "Disable (CUDA) C++ extensions")
    set(CMAKE_SYCL_STANDARD 20 CACHE STRING "The (SYCL) C++ standard to use")
    set(CMAKE_HIP_STANDARD 20 CACHE STRING "The (HIP) C++ standard to use")

    # Set the CUDA architecture to build code for.
    set(CMAKE_CUDA_ARCHITECTURES
        "75"
        CACHE STRING
        "CUDA architectures to build device code for"
    )

    # Set the HIP architecture to build code for.
    set(CMAKE_HIP_ARCHITECTURES
        "gfx1101;gfx90a"
        CACHE STRING
        "HIP architectures to build device code for"
    )

    # Flag controlling whether warnings should make the build fail.
    option(
        TRACCC_FAIL_ON_WARNINGS
        "Make the build fail on compiler/linker warnings"
        FALSE
    )

    # Explicitly set the output directory for the binaries. Such that if this
    # project is included by another project, the main project's configuration would
    # win out.
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

    # Temporary setting for the traccc device log level, until it can be removed.
    set(TRACCC_DEVICE_LOG_LVL
        "NONE"
        CACHE STRING
        "Log level for traccc and detray device code"
    )

    # Flags controlling which parts of traccc to build.
    option(TRACCC_BUILD_CUDA "Build the traccc::cuda library" FALSE)
    option(
        TRACCC_BUILD_CUDA_UTILS
        "Build the traccc::cuda_utils library"
        ${TRACCC_BUILD_CUDA}
    )
    option(TRACCC_BUILD_HIP "Build the traccc::hip library" FALSE)
    option(
        TRACCC_BUILD_HIP_UTILS
        "Build the traccc::hip_utils library"
        ${TRACCC_BUILD_HIP}
    )
    option(TRACCC_BUILD_SYCL "Build the traccc::sycl library" FALSE)
    option(
        TRACCC_BUILD_SYCL_UTILS
        "Build the traccc::sycl_utils library"
        ${TRACCC_BUILD_SYCL}
    )
    option(
        TRACCC_BUILD_ALPAKA
        "Build the Alpaka sources included in traccc"
        FALSE
    )
    option(
        TRACCC_BUILD_IO
        "Build the IO module (needed by examples, performance, testing)"
        TRUE
    )
    option(TRACCC_BUILD_PERFORMANCE "Build the performance module" TRUE)
    option(TRACCC_BUILD_SIMULATION "Build the simulation module" TRUE)
    option(TRACCC_BUILD_TESTING "Build the (unit) tests of traccc" TRUE)
    option(TRACCC_BUILD_EXAMPLES "Build the examples of traccc" TRUE)

    # Flags controlling what traccc should use.
    option(TRACCC_USE_SYSTEM_LIBS "Use system libraries be default" FALSE)
    option(TRACCC_USE_SPACK_LIBS "Use Spack libraries by default" FALSE)
    option(TRACCC_USE_ROOT "Use ROOT in the build (if needed)" TRUE)
else()
    set(CMAKE_CXX_STANDARD "${ACTS_CXX_STANDARD}")

    # traccc is built as part of ACTS, and ACTS owns its configuration. So the
    # traccc settings are normal variables, set from the ACTS options. A cache
    # entry with the same name has no effect, so it is an error. Remove the
    # entry, so that it only reports a stale entry from an older build once.
    # The optional third argument names the ACTS setting(s) that control the
    # traccc setting, for the error message.
    macro(traccc_set_from_acts name value)
        if(DEFINED CACHE{${name}})
            if(${ARGC} GREATER 2)
                set(_traccc_set_from_acts_hint "ACTS sets it from ${ARGV2}.")
            else()
                set(_traccc_set_from_acts_hint "ACTS sets it itself.")
            endif()
            message(
                SEND_ERROR
                "${name} cannot be set when traccc is built as part of ACTS. "
                "${_traccc_set_from_acts_hint}"
            )
            unset(_traccc_set_from_acts_hint)
            unset(${name} CACHE)
        endif()
        set(${name} "${value}")
    endmacro()

    # The CUDA dependencies (vecmem::cuda, CUDA architectures) follow
    # ACTS_ENABLE_CUDA, so traccc CUDA code cannot be built without it.
    if(_acts_traccc_build_cuda_utils AND NOT ACTS_ENABLE_CUDA)
        message(
            FATAL_ERROR
            "The traccc CUDA code (ACTS_TRACCC_BUILD_CUDA or "
            "alpaka_ACC_GPU_CUDA_ENABLE) requires ACTS_ENABLE_CUDA=ON."
        )
    endif()

    # Settings with an ACTS_TRACCC_* option of their own.
    traccc_set_from_acts(
        TRACCC_BUILD_CUDA
        "${ACTS_TRACCC_BUILD_CUDA}"
        ACTS_TRACCC_BUILD_CUDA
    )
    traccc_set_from_acts(
        TRACCC_BUILD_CUDA_UTILS
        "${_acts_traccc_build_cuda_utils}"
        "ACTS_TRACCC_BUILD_CUDA and alpaka_ACC_GPU_CUDA_ENABLE"
    )
    traccc_set_from_acts(
        TRACCC_BUILD_HIP
        "${ACTS_TRACCC_BUILD_HIP}"
        ACTS_TRACCC_BUILD_HIP
    )
    traccc_set_from_acts(
        TRACCC_BUILD_HIP_UTILS
        "${_acts_traccc_build_hip_utils}"
        "ACTS_TRACCC_BUILD_HIP and alpaka_ACC_GPU_HIP_ENABLE"
    )
    traccc_set_from_acts(
        TRACCC_BUILD_SYCL
        "${ACTS_TRACCC_BUILD_SYCL}"
        ACTS_TRACCC_BUILD_SYCL
    )
    traccc_set_from_acts(
        TRACCC_BUILD_SYCL_UTILS
        "${_acts_traccc_build_sycl_utils}"
        "ACTS_TRACCC_BUILD_SYCL and alpaka_ACC_SYCL_ENABLE"
    )
    traccc_set_from_acts(
        TRACCC_BUILD_ALPAKA
        "${ACTS_TRACCC_BUILD_ALPAKA}"
        ACTS_TRACCC_BUILD_ALPAKA
    )
    traccc_set_from_acts(
        TRACCC_BUILD_IO
        "${ACTS_TRACCC_BUILD_IO}"
        ACTS_TRACCC_BUILD_IO
    )
    traccc_set_from_acts(
        TRACCC_BUILD_PERFORMANCE
        "${ACTS_TRACCC_BUILD_PERFORMANCE}"
        ACTS_TRACCC_BUILD_PERFORMANCE
    )
    traccc_set_from_acts(
        TRACCC_BUILD_TESTING
        "${ACTS_TRACCC_BUILD_TESTS}"
        ACTS_TRACCC_BUILD_TESTS
    )
    traccc_set_from_acts(
        TRACCC_BUILD_EXAMPLES
        "${ACTS_TRACCC_BUILD_EXAMPLES}"
        ACTS_TRACCC_BUILD_EXAMPLES
    )
    traccc_set_from_acts(
        TRACCC_ENABLE_NVTX_PROFILING
        "${ACTS_TRACCC_ENABLE_NVTX_PROFILING}"
        ACTS_TRACCC_ENABLE_NVTX_PROFILING
    )
    traccc_set_from_acts(
        TRACCC_SUPPORTED_DETECTORS
        "${ACTS_TRACCC_SUPPORTED_DETECTORS}"
        ACTS_TRACCC_SUPPORTED_DETECTORS
    )

    traccc_set_from_acts(
        TRACCC_BUILD_SIMULATION
        OFF
        "ACTS_TRACCC_BUILD_TESTS and ACTS_TRACCC_BUILD_EXAMPLES"
    )
    set_option_if(
        TRACCC_BUILD_SIMULATION
        TRACCC_BUILD_TESTING
        OR
        TRACCC_BUILD_EXAMPLES
    )

    # Settings that follow the general ACTS options.
    traccc_set_from_acts(
        TRACCC_USE_ROOT
        "${ACTS_BUILD_PLUGIN_ROOT}"
        ACTS_BUILD_PLUGIN_ROOT
    )
    traccc_set_from_acts(TRACCC_FAIL_ON_WARNINGS OFF ACTS_FAIL_ON_WARNINGS)
    set_option_if(TRACCC_FAIL_ON_WARNINGS ACTS_FAIL_ON_WARNINGS)

    # Dependencies. ACTS sets up all of them, except for the ones below.
    traccc_set_from_acts(TRACCC_SETUP_VECMEM FALSE)
    traccc_set_from_acts(
        TRACCC_USE_SYSTEM_VECMEM
        FALSE
        ACTS_USE_SYSTEM_VECMEM
    )
    # Only the examples and oneDPL need TBB. Build it if it is not available.
    traccc_set_from_acts(TRACCC_SETUP_TBB FALSE)
    if(TRACCC_BUILD_EXAMPLES OR TRACCC_BUILD_SYCL_UTILS)
        find_package(TBB ${_acts_tbb_version} CONFIG)
        if(NOT TBB_FOUND)
            traccc_set_from_acts(TRACCC_SETUP_TBB TRUE)
        endif()
    endif()
    traccc_set_from_acts(TRACCC_USE_SYSTEM_TBB FALSE)
    traccc_set_from_acts(TRACCC_SETUP_THRUST TRUE)
    # rocThrust is needed for HIP on AMD, and is taken from the ROCm installation.
    if(
        (TRACCC_BUILD_HIP AND NOT "${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia")
        OR (TRACCC_BUILD_ALPAKA AND alpaka_ACC_GPU_HIP_ENABLE)
    )
        traccc_set_from_acts(TRACCC_SETUP_ROCTHRUST TRUE)
    else()
        traccc_set_from_acts(TRACCC_SETUP_ROCTHRUST FALSE)
    endif()
    traccc_set_from_acts(TRACCC_USE_SYSTEM_ROCTHRUST TRUE)
    traccc_set_from_acts(TRACCC_SETUP_DPL "${TRACCC_BUILD_SYCL_UTILS}")
    traccc_set_from_acts(TRACCC_USE_SYSTEM_DPL FALSE)
    traccc_set_from_acts(
        TRACCC_SETUP_ALPAKA
        "${TRACCC_BUILD_ALPAKA}"
        ACTS_TRACCC_BUILD_ALPAKA
    )
    traccc_set_from_acts(
        TRACCC_USE_SYSTEM_ALPAKA
        "${ACTS_USE_SYSTEM_ALPAKA}"
        ACTS_USE_SYSTEM_ALPAKA
    )
    traccc_set_from_acts(TRACCC_SETUP_COVFIE FALSE)
    traccc_set_from_acts(
        TRACCC_USE_SYSTEM_COVFIE
        FALSE
        ACTS_USE_SYSTEM_COVFIE
    )
    traccc_set_from_acts(TRACCC_SETUP_ACTS FALSE)
    traccc_set_from_acts(
        TRACCC_SETUP_GOOGLETEST
        "${TRACCC_BUILD_TESTING}"
        ACTS_TRACCC_BUILD_TESTS
    )
    traccc_set_from_acts(TRACCC_USE_SYSTEM_GOOGLETEST FALSE)
endif()

# Include the traccc CMake code.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(traccc-functions)

# Honour ACTS_FORCE_ASSERTIONS in standalone builds inside the ACTS tree. When
# traccc is built via add_subdirectory from ACTS this is a no-op.
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/ActsForceAssertions.cmake)

# Use an absolute path to avoid changing package lookup via CMAKE_MODULE_PATH.
include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/ActsWarningFlags.cmake")

# Check CUDA and SYCL C++ standards
if(${TRACCC_BUILD_CUDA} AND ${CMAKE_CUDA_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_CUDA_STANDARD=${CMAKE_CUDA_STANDARD}, but traccc requires C++>=20"
    )
endif()

if(${TRACCC_BUILD_SYCL} AND ${CMAKE_SYCL_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_SYCL_STANDARD=${CMAKE_SYCL_STANDARD}, but traccc requires C++>=20"
    )
endif()

if(${TRACCC_BUILD_HIP} AND ${CMAKE_HIP_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_HIP_STANDARD=${CMAKE_HIP_STANDARD}, but traccc requires C++>=20"
    )
endif()

# If we are building with CUDA, we need to find the toolkit early because it
# contains CCCL and Thrust.
if(TRACCC_BUILD_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# The detector types supported by this build of traccc.
if(PROJECT_IS_TOP_LEVEL)
    set(TRACCC_SUPPORTED_DETECTORS
        "default_detector;odd_detector;telescope_detector"
        CACHE STRING
        "Semicolon-separated list of detector types to support in this build"
    )
endif()

# Set up build profiling for the project.
if(CTEST_USE_LAUNCHERS)
    # Find the bash and time executables.
    find_program(BASH_EXECUTABLE bash REQUIRED)
    find_program(TIME_EXECUTABLE time REQUIRED)

    # Decide what flag to use with the time executable to make it print verbose
    # information.
    if("${CMAKE_HOST_SYSTEM_NAME}" MATCHES "Darwin")
        set(TIME_VERBOSE_FLAG "-l")
    elseif("${CMAKE_HOST_SYSTEM_NAME}" MATCHES "Linux")
        set(TIME_VERBOSE_FLAG "-v")
    else()
        message(
            WARNING
            "Build profiling is only supported on Linux and macOS."
            "This build will likely fail."
        )
        set(TIME_VERBOSE_FLAG "")
    endif()

    # Configure the script that would intercept the build commands and save
    # profile logs for them.
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/traccc-ctest.sh.in"
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/traccc-ctest.sh"
        @ONLY
    )
    set(CMAKE_CTEST_COMMAND
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/traccc-ctest.sh"
    )

    # Clean up.
    unset(TIME_VERBOSE_FLAG)

    # Remove the performance log during a cleaning step.
    set_property(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        APPEND
        PROPERTY
            ADDITIONAL_MAKE_CLEAN_FILES
                "${CMAKE_CURRENT_BINARY_DIR}/traccc_build_performance.log"
    )

    # Let the user know what happened.
    message(
        STATUS
        "Saving traccc build performance logs using: ${CMAKE_CTEST_COMMAND}"
    )
endif()

# Set up VecMem.
option(TRACCC_SETUP_VECMEM "Set up the VecMem target(s) explicitly" TRUE)
option(
    TRACCC_USE_SYSTEM_VECMEM
    "Pick up an existing installation of VecMem from the build environment"
    ${TRACCC_USE_SYSTEM_LIBS}
)
if(TRACCC_SETUP_VECMEM)
    if(TRACCC_USE_SYSTEM_VECMEM)
        find_package(vecmem 1.25.0 REQUIRED)
    else()
        add_subdirectory(extern/vecmem)
    endif()
    # Make the "VecMem language code" available for the whole project.
    list(PREPEND CMAKE_MODULE_PATH "${VECMEM_LANGUAGE_DIR}")
endif()

# Set up TBB.
option(TRACCC_SETUP_TBB "Set up the TBB target(s) explicitly" TRUE)
if(TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
    set(TRACCC_USE_SYSTEM_TBB_DEFAULT ON)
else()
    set(TRACCC_USE_SYSTEM_TBB_DEFAULT OFF)
endif()
option(
    TRACCC_USE_SYSTEM_TBB
    "Pick up an existing installation of TBB from the build environment"
    ${TRACCC_USE_SYSTEM_TBB_DEFAULT}
)
unset(TRACCC_USE_SYSTEM_TBB_DEFAULT)
if(TRACCC_SETUP_TBB)
    if(TRACCC_USE_SYSTEM_TBB)
        find_package(TBB REQUIRED)
    else()
        add_subdirectory(extern/tbb)
    endif()
endif()

# Set up CCCL.
option(TRACCC_SETUP_THRUST "Set up the Thrust target(s) explicitly" TRUE)
if(TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS OR CUDAToolkit_FOUND)
    set(TRACCC_USE_SYSTEM_THRUST_DEFAULT ON)
else()
    set(TRACCC_USE_SYSTEM_THRUST_DEFAULT OFF)
endif()
option(
    TRACCC_USE_SYSTEM_THRUST
    "Pick up an existing installation of Thrust from the build environment"
    ${TRACCC_USE_SYSTEM_THRUST_DEFAULT}
)
unset(TRACCC_USE_SYSTEM_THRUST_DEFAULT)
if(TRACCC_SETUP_THRUST)
    if(TRACCC_USE_SYSTEM_THRUST)
        if(CUDAToolkit_FOUND)
            find_package(
                Thrust
                REQUIRED
                HINTS "${CUDAToolkit_LIBRARY_DIR}/cmake"
            )
        else()
            find_package(Thrust REQUIRED)
        endif()
    else()
        if(CUDAToolkit_FOUND)
            set(TRACCC_CCCL_DIAGNOSTIC WARNING)
            if("$ENV{CI}")
                set(TRACCC_CCCL_DIAGNOSTIC FATAL_ERROR)
            endif()
            message(
                ${TRACCC_CCCL_DIAGNOSTIC}
                "CUDA toolkit is available, but a bundled CCCL is to be used. This may lead to cross-version include path pollution and errors."
            )
            unset(TRACCC_CCCL_DIAGNOSTIC)
        endif()
        add_subdirectory(extern/cccl)
    endif()
endif()

# Set up rocThrust.
option(TRACCC_SETUP_ROCTHRUST "Set up the rocThrust target(s) explicitly" FALSE)
option(
    TRACCC_USE_SYSTEM_ROCTHRUST
    "Pick up an existing installation of rocThrust from the build environment"
    ${TRACCC_USE_SYSTEM_LIBS}
)
if(TRACCC_SETUP_ROCTHRUST)
    if(TRACCC_USE_SYSTEM_ROCTHRUST)
        find_package(rocthrust REQUIRED)
        # Dress up the rocthrust target a little.
        target_compile_definitions(
            roc::rocthrust
            INTERFACE THRUST_IGNORE_CUB_VERSION_CHECK
        )
    else()
        add_subdirectory(extern/rocThrust)
        # Dress up the rocthrust target a little.
        target_compile_definitions(
            rocthrust
            INTERFACE THRUST_IGNORE_CUB_VERSION_CHECK
        )
    endif()
endif()

# Set up DPL if SYCL, or Alpaka with its SYCL back-end, is built.
if(TRACCC_BUILD_SYCL OR (TRACCC_BUILD_ALPAKA AND alpaka_ACC_SYCL_ENABLE))
    set(TRACCC_SETUP_DPL_DEFAULT ON)
else()
    set(TRACCC_SETUP_DPL_DEFAULT OFF)
endif()
option(
    TRACCC_SETUP_DPL
    "Set up the DPL target(s) explicitly"
    ${TRACCC_SETUP_DPL_DEFAULT}
)
unset(TRACCC_SETUP_DPL_DEFAULT)
if(TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
    set(TRACCC_USE_SYSTEM_DPL_DEFAULT ON)
else()
    set(TRACCC_USE_SYSTEM_DPL_DEFAULT OFF)
endif()
option(
    TRACCC_USE_SYSTEM_DPL
    "Pick up an existing installation of DPL from the build environment"
    ${TRACCC_USE_SYSTEM_DPL_DEFAULT}
)
unset(TRACCC_USE_SYSTEM_DPL_DEFAULT)
if(TRACCC_SETUP_DPL)
    if(TRACCC_USE_SYSTEM_DPL)
        # OneDPL determines whether SYCL is supported by asking the C++ compiler
        # rather than the SYCL compiler, as a dedicated SYCL compiler is a non-
        # standard traccc idea. To override the default behaviour (i.e. OneDPL
        # testing the C++ compiler and finding that it does _not_ support SYCL)
        # we simply override the support flags. This is fragile, as the flags
        # are internal to OneDPL, but it's the best we can do. Note thr flag was
        # renamed in https://github.com/uxlfoundation/oneDPL/pull/1949.
        set(_sycl_support ON)
        set(SYCL_SUPPORT ON)
        find_package(oneDPL REQUIRED)
    else()
        add_subdirectory(extern/dpl)
    endif()
    # OneDPL attaches the `-fsycl` flag to the C++ compiler, which causes
    # it to incorrectly trigger some preprocessor definitions, thereby
    # loading SYCL-specific header files which do not exist for e.g. a
    # generic clang installation. Therefore, we have to manually wipe the
    # compile flags that OneDPL sets.
    set_target_properties(
        oneDPL
        PROPERTIES
            INTERFACE_COMPILE_OPTIONS ""
            INTERFACE_LINK_LIBRARIES TBB::tbb
    )
endif()

# Set up Alpaka.
option(TRACCC_SETUP_ALPAKA "Set up the Alpaka library" ${TRACCC_BUILD_ALPAKA})
if(TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
    set(TRACCC_USE_SYSTEM_ALPAKA_DEFAULT ON)
else()
    set(TRACCC_USE_SYSTEM_ALPAKA_DEFAULT OFF)
endif()
option(
    TRACCC_USE_SYSTEM_ALPAKA
    "Pick up an existing installation of Alpaka from the build environment"
    ${TRACCC_USE_SYSTEM_ALPAKA_DEFAULT}
)
unset(TRACCC_USE_SYSTEM_ALPAKA_DEFAULT)
if(TRACCC_SETUP_ALPAKA)
    # Default options for the Alpaka build.
    set(alpaka_ACC_CPU_B_SEQ_T_THREADS_ENABLE
        TRUE
        CACHE BOOL
        "Enable the serial backend of Alpaka"
    )
    if(TRACCC_USE_SYSTEM_ALPAKA)
        find_package(alpaka REQUIRED)
    else()
        add_subdirectory(extern/alpaka)
    endif()
endif()
# Alpaka drops members under NDEBUG that its (now forced) asserts still use.
if(ACTS_FORCE_ASSERTIONS AND TARGET alpaka::alpaka)
    get_target_property(_alpaka_target alpaka::alpaka ALIASED_TARGET)
    if(NOT _alpaka_target)
        set(_alpaka_target alpaka::alpaka)
    endif()
    set_property(
        TARGET ${_alpaka_target}
        APPEND
        PROPERTY INTERFACE_COMPILE_OPTIONS -UNDEBUG
    )
    unset(_alpaka_target)
endif()

# Set up covfie.
option(TRACCC_SETUP_COVFIE "Set up the covfie target(s) explicitly" TRUE)
option(
    TRACCC_USE_SYSTEM_COVFIE
    "Pick up an existing installation of covfie from the build environment"
    ${TRACCC_USE_SYSTEM_LIBS}
)
if(TRACCC_SETUP_COVFIE)
    if(TRACCC_USE_SYSTEM_COVFIE)
        find_package(covfie REQUIRED)
    else()
        add_subdirectory(extern/covfie)
    endif()
endif()

# Set up Acts. Acts (and, transitively, Detray) is colocated in-tree and is the
# single source of truth, so it is always built from that source rather than
# picked up from an external installation.
option(TRACCC_SETUP_ACTS "Set up the Acts target(s) explicitly" TRUE)
if(TRACCC_SETUP_ACTS)
    add_subdirectory(extern/acts)
endif()

# Set up GoogleTest.
include(CTest)
if(BUILD_TESTING AND TRACCC_BUILD_TESTING)
    set(TRACCC_DEFAULT_SETUP_GOOGLETEST TRUE)
endif()
option(
    TRACCC_SETUP_GOOGLETEST
    "Set up the GoogleTest target(s) explicitly"
    ${TRACCC_DEFAULT_SETUP_GOOGLETEST}
)
if(TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
    set(TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT ON)
else()
    set(TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT OFF)
endif()
option(
    TRACCC_USE_SYSTEM_GOOGLETEST
    "Pick up an existing installation of GoogleTest from the build environment"
    ${TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT}
)
unset(TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT)
if(TRACCC_SETUP_GOOGLETEST)
    if(TRACCC_USE_SYSTEM_GOOGLETEST)
        find_package(GTest REQUIRED)
    else()
        add_subdirectory(extern/googletest)
    endif()
endif()

option(
    TRACCC_ENABLE_NVTX_PROFILING
    "Use instrument functions to enable fine grained profiling"
    FALSE
)

# Build the traccc code.
add_subdirectory(core)
add_subdirectory(device/common)
if(TRACCC_BUILD_CUDA_UTILS)
    add_subdirectory(device/cuda_utils)
endif()
if(TRACCC_BUILD_CUDA)
    add_subdirectory(device/cuda)
endif()
if(TRACCC_BUILD_HIP_UTILS)
    add_subdirectory(device/hip_utils)
endif()
if(TRACCC_BUILD_HIP)
    add_subdirectory(device/hip)
endif()
if(TRACCC_BUILD_SYCL_UTILS)
    add_subdirectory(device/sycl_utils)
endif()
if(TRACCC_BUILD_SYCL)
    add_subdirectory(device/sycl)
endif()
if(TRACCC_BUILD_ALPAKA)
    add_subdirectory(device/alpaka)
endif()
if(TRACCC_BUILD_IO)
    add_subdirectory(io)
else()
    message(
        STATUS
        "traccc::io not built, traccc::performance and traccc::simulation are forcefully switched off."
    )
endif()

if(TRACCC_BUILD_PERFORMANCE)
    if(NOT TRACCC_BUILD_IO)
        message(
            FATAL_ERROR
            "traccc::io is disabled, but it is required to build the performance tests."
        )
    endif()
    add_subdirectory(performance)
endif()
if(TRACCC_BUILD_SIMULATION)
    if(NOT TRACCC_BUILD_IO)
        message(
            FATAL_ERROR
            "traccc::io is disabled, but it is required to build the simulation."
        )
    endif()
    add_subdirectory(simulation)
endif()

if(TRACCC_BUILD_EXAMPLES)
    # Find Boost.
    find_package(Boost CONFIG REQUIRED COMPONENTS program_options filesystem)
    if(NOT TRACCC_BUILD_IO)
        message(
            FATAL_ERROR
            "traccc::io is disabled, but it is required to build the examples."
        )
    endif()
    if(NOT TRACCC_BUILD_PERFORMANCE)
        message(
            FATAL_ERROR
            "traccc::performance is disabled, but it is required to build the examples."
        )
    endif()
    add_subdirectory(examples)
endif()

# Set up the test(s).
if(BUILD_TESTING AND TRACCC_BUILD_TESTING)
    if(NOT TRACCC_BUILD_IO)
        message(
            FATAL_ERROR
            "traccc::io is disabled, but it is required to build the tests."
        )
    endif()
    if(NOT TRACCC_BUILD_SIMULATION)
        message(
            FATAL_ERROR
            "traccc::simulation is disabled, but it is required to build the tests."
        )
    endif()
    if(NOT TRACCC_BUILD_PERFORMANCE)
        message(
            FATAL_ERROR
            "traccc::performance is disabled, but it is required to build the tests."
        )
    endif()
    add_subdirectory(tests)
endif()

# Set up the packaging of the project.
include(traccc-packaging)
