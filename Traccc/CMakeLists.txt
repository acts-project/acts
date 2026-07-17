# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the project.
cmake_minimum_required( VERSION 3.25 )
project( traccc VERSION 1.6.0 LANGUAGES CXX )

# Set up the used C++ standard(s).
set( CMAKE_CXX_STANDARD 20 CACHE STRING "The (host) C++ standard to use" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Disable (host) C++ extensions" )
set( CMAKE_CUDA_STANDARD 20 CACHE STRING "The (CUDA) C++ standard to use" )
set( CMAKE_CUDA_EXTENSIONS FALSE CACHE BOOL "Disable (CUDA) C++ extensions" )
set( CMAKE_SYCL_STANDARD 20 CACHE STRING "The (SYCL) C++ standard to use" )
set( CMAKE_HIP_STANDARD 20 CACHE STRING "The (HIP) C++ standard to use" )

if(${CMAKE_CXX_STANDARD} LESS 20)
   message(SEND_ERROR "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}, but traccc requires C++>=20")
endif()

# Set the CUDA architecture to build code for.
set( CMAKE_CUDA_ARCHITECTURES "75" CACHE STRING
   "CUDA architectures to build device code for" )

# Set the HIP architecture to build code for.
set( CMAKE_HIP_ARCHITECTURES "gfx1101;gfx90a" CACHE STRING
   "HIP architectures to build device code for" )

# Flag controlling whether warnings should make the build fail.
option( TRACCC_FAIL_ON_WARNINGS
   "Make the build fail on compiler/linker warnings" FALSE )

# Standard CMake include(s).
include( GNUInstallDirs )

# Explicitly set the output directory for the binaries. Such that if this
# project is included by another project, the main project's configuration would
# win out.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" CACHE PATH
   "Directory for the built binaries" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built libraries" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built static libraries" )

# Include the traccc CMake code.
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )
include( traccc-functions )

# Temporary setting for the traccc::scalar type, until it can be removed.
set( TRACCC_CUSTOM_SCALARTYPE "float" CACHE STRING
   "Scalar type to use in the TRACCC code" )

# Temporary setting for the traccc device log level, until it can be removed.
set( TRACCC_DEVICE_LOG_LVL "NONE" CACHE STRING
   "Log level for traccc and detray device code" )

# Flags controlling which parts of traccc to build.
option( TRACCC_BUILD_CUDA
   "Build the traccc::cuda library" FALSE )
option( TRACCC_BUILD_CUDA_UTILS
   "Build the traccc::cuda_utils library" ${TRACCC_BUILD_CUDA} )
option( TRACCC_BUILD_HIP "Build the HIP sources included in traccc" FALSE)
option( TRACCC_BUILD_SYCL
   "Build the traccc::sycl library" FALSE )
option( TRACCC_BUILD_SYCL_UTILS
   "Build the traccc::sycl_utils library" ${TRACCC_BUILD_SYCL} )
option( TRACCC_BUILD_ALPAKA "Build the Alpaka sources included in traccc"
   FALSE )
option( TRACCC_BUILD_IO "Build the IO module (needed by examples, performance, testing)" TRUE )
option( TRACCC_BUILD_PERFORMANCE "Build the performance module" TRUE )
option( TRACCC_BUILD_SIMULATION "Build the simulation module" TRUE )
option( TRACCC_BUILD_TESTING "Build the (unit) tests of traccc" TRUE )
option( TRACCC_BUILD_EXAMPLES "Build the examples of traccc" TRUE )

# Flags controlling what traccc should use.
option( TRACCC_USE_SYSTEM_LIBS "Use system libraries be default" FALSE )
option( TRACCC_USE_SPACK_LIBS "Use Spack libraries by default" FALSE )
option( TRACCC_USE_ROOT "Use ROOT in the build (if needed)" TRUE )

# Check CUDA and SYCL C++ standards
if(${TRACCC_BUILD_CUDA} AND ${CMAKE_CUDA_STANDARD} LESS 20)
   message(SEND_ERROR "CMAKE_CUDA_STANDARD=${CMAKE_CUDA_STANDARD}, but traccc requires C++>=20")
endif()

if(${TRACCC_BUILD_SYCL} AND ${CMAKE_SYCL_STANDARD} LESS 20)
   message(SEND_ERROR "CMAKE_SYCL_STANDARD=${CMAKE_SYCL_STANDARD}, but traccc requires C++>=20")
endif()

if(${TRACCC_BUILD_HIP} AND ${CMAKE_HIP_STANDARD} LESS 20)
   message(SEND_ERROR "CMAKE_HIP_STANDARD=${CMAKE_HIP_STANDARD}, but traccc requires C++>=20")
endif()

# In order to generate specializations for our kernels, we need to have both
# a working Python interpreter, and a list of supported detectors.
find_package (Python COMPONENTS Interpreter REQUIRED)
# The detector types supported by this build of traccc.
set( TRACCC_SUPPORTED_DETECTORS
   "default_detector;odd_detector;telescope_detector"
   CACHE STRING
   "Semicolon-separated list of detector types to support in this build" )

# Set up build profiling for the project.
if( CTEST_USE_LAUNCHERS )

   # Find the bash and time executables.
   find_program( BASH_EXECUTABLE bash REQUIRED )
   find_program( TIME_EXECUTABLE time REQUIRED )

   # Decide what flag to use with the time executable to make it print verbose
   # information.
   if( "${CMAKE_HOST_SYSTEM_NAME}" MATCHES "Darwin" )
      set( TIME_VERBOSE_FLAG "-l" )
   elseif( "${CMAKE_HOST_SYSTEM_NAME}" MATCHES "Linux" )
      set( TIME_VERBOSE_FLAG "-v" )
   else()
      message( WARNING "Build profiling is only supported on Linux and macOS."
                       "This build will likely fail." )
      set( TIME_VERBOSE_FLAG "" )
   endif()

   # Configure the script that would intercept the build commands and save
   # profile logs for them.
   configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/traccc-ctest.sh.in"
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/traccc-ctest.sh"
      @ONLY )
   set( CMAKE_CTEST_COMMAND
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/traccc-ctest.sh" )

   # Clean up.
   unset( TIME_VERBOSE_FLAG )

   # Remove the performance log during a cleaning step.
   set_property( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" APPEND PROPERTY
      ADDITIONAL_MAKE_CLEAN_FILES
      "${CMAKE_CURRENT_BINARY_DIR}/traccc_build_performance.log" )

   # Let the user know what happened.
   message( STATUS
      "Saving traccc build performance logs using: ${CMAKE_CTEST_COMMAND}" )
endif()

# Set up VecMem.
option( TRACCC_SETUP_VECMEM
   "Set up the VecMem target(s) explicitly" TRUE )
option( TRACCC_USE_SYSTEM_VECMEM
   "Pick up an existing installation of VecMem from the build environment"
   ${TRACCC_USE_SYSTEM_LIBS} )
if( TRACCC_SETUP_VECMEM )
   if( TRACCC_USE_SYSTEM_VECMEM )
      find_package( vecmem 1.25.0 REQUIRED )
   else()
      add_subdirectory( extern/vecmem )
   endif()
   # Make the "VecMem language code" available for the whole project.
   list( PREPEND CMAKE_MODULE_PATH "${VECMEM_LANGUAGE_DIR}" )
endif()

# Set up TBB.
option( TRACCC_SETUP_TBB
   "Set up the TBB target(s) explicitly" TRUE )
if (TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
   set(TRACCC_USE_SYSTEM_TBB_DEFAULT ON)
else()
   set(TRACCC_USE_SYSTEM_TBB_DEFAULT OFF)
endif()
option( TRACCC_USE_SYSTEM_TBB
   "Pick up an existing installation of TBB from the build environment"
   ${TRACCC_USE_SYSTEM_TBB_DEFAULT} )
unset(TRACCC_USE_SYSTEM_TBB_DEFAULT)
if( TRACCC_SETUP_TBB )
   if( TRACCC_USE_SYSTEM_TBB )
      find_package( TBB REQUIRED )
   else()
      add_subdirectory( extern/tbb )
   endif()
endif()

# Set up CCCL.
option( TRACCC_SETUP_THRUST
   "Set up the Thrust target(s) explicitly" TRUE )
if (TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
   set(TRACCC_USE_SYSTEM_THRUST_DEFAULT ON)
else()
   set(TRACCC_USE_SYSTEM_THRUST_DEFAULT OFF)
endif()
option( TRACCC_USE_SYSTEM_THRUST
   "Pick up an existing installation of Thrust from the build environment"
   ${TRACCC_USE_SYSTEM_THRUST_DEFAULT} )
unset(TRACCC_USE_SYSTEM_THRUST_DEFAULT)
if( TRACCC_SETUP_THRUST )
   if( TRACCC_USE_SYSTEM_THRUST )
      find_package( Thrust REQUIRED )
   else()
      add_subdirectory( extern/cccl )
   endif()
endif()

# Set up rocThrust.
option( TRACCC_SETUP_ROCTHRUST
   "Set up the rocThrust target(s) explicitly" FALSE )
option( TRACCC_USE_SYSTEM_ROCTHRUST
   "Pick up an existing installation of rocThrust from the build environment"
   ${TRACCC_USE_SYSTEM_LIBS} )
if( TRACCC_SETUP_ROCTHRUST )
   if( TRACCC_USE_SYSTEM_ROCTHRUST )
      find_package( rocthrust REQUIRED )
      # Dress up the rocthrust target a little.
      target_compile_definitions( roc::rocthrust INTERFACE
         THRUST_IGNORE_CUB_VERSION_CHECK )
   else()
      add_subdirectory( extern/rocThrust )
      # Dress up the rocthrust target a little.
      target_compile_definitions( rocthrust INTERFACE
         THRUST_IGNORE_CUB_VERSION_CHECK )
   endif()
endif()

# Set up DPL if SYCL is built.
option( TRACCC_SETUP_DPL
   "Set up the DPL target(s) explicitly" ${TRACCC_BUILD_SYCL} )
if (TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
   set(TRACCC_USE_SYSTEM_DPL_DEFAULT ON)
else()
   set(TRACCC_USE_SYSTEM_DPL_DEFAULT OFF)
endif()
option( TRACCC_USE_SYSTEM_DPL
   "Pick up an existing installation of DPL from the build environment"
   ${TRACCC_USE_SYSTEM_DPL_DEFAULT} )
unset(TRACCC_USE_SYSTEM_DPL_DEFAULT)
if( TRACCC_SETUP_DPL )
   if( TRACCC_USE_SYSTEM_DPL )
      # OneDPL determines whether SYCL is supported by asking the C++ compiler
      # rather than the SYCL compiler, as a dedicated SYCL compiler is a non-
      # standard traccc idea. To override the default behaviour (i.e. OneDPL
      # testing the C++ compiler and finding that it does _not_ support SYCL)
      # we simply override the support flags. This is fragile, as the flags
      # are internal to OneDPL, but it's the best we can do. Note thr flag was
      # renamed in https://github.com/uxlfoundation/oneDPL/pull/1949.
      set(_sycl_support ON)
      set(SYCL_SUPPORT ON)
      find_package( oneDPL REQUIRED )
   else()
      add_subdirectory( extern/dpl )
   endif()
   # OneDPL attaches the `-fsycl` flag to the C++ compiler, which causes
   # it to incorrectly trigger some preprocessor definitions, thereby
   # loading SYCL-specific header files which do not exist for e.g. a
   # generic clang installation. Therefore, we have to manually wipe the
   # compile flags that OneDPL sets.
   set_target_properties(oneDPL PROPERTIES
      INTERFACE_COMPILE_OPTIONS ""
      INTERFACE_LINK_LIBRARIES "" )
endif()

# Set up Alpaka.
option( TRACCC_SETUP_ALPAKA
   "Set up the Alpaka library" ${TRACCC_BUILD_ALPAKA})
if (TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
   set(TRACCC_USE_SYSTEM_ALPAKA_DEFAULT ON)
else()
   set(TRACCC_USE_SYSTEM_ALPAKA_DEFAULT OFF)
endif()
option( TRACCC_USE_SYSTEM_ALPAKA
   "Pick up an existing installation of Alpaka from the build environment"
   ${TRACCC_USE_SYSTEM_ALPAKA_DEFAULT} )
unset(TRACCC_USE_SYSTEM_ALPAKA_DEFAULT)
if( TRACCC_SETUP_ALPAKA )
   # Default options for the Alpaka build.
   set( alpaka_ACC_CPU_B_SEQ_T_THREADS_ENABLE TRUE CACHE BOOL
      "Enable the serial backend of Alpaka" )
   if( TRACCC_USE_SYSTEM_ALPAKA )
      find_package( alpaka REQUIRED )
   else()
      add_subdirectory( extern/alpaka )
   endif()
endif()

# Set up covfie.
option( TRACCC_SETUP_COVFIE
   "Set up the covfie target(s) explicitly" TRUE )
option( TRACCC_USE_SYSTEM_COVFIE
   "Pick up an existing installation of covfie from the build environment"
   ${TRACCC_USE_SYSTEM_LIBS} )
if( TRACCC_SETUP_COVFIE )
   if( TRACCC_USE_SYSTEM_COVFIE )
      find_package( covfie REQUIRED )
   else()
      add_subdirectory( extern/covfie )
   endif()
endif()

# Set up Acts.
option( TRACCC_SETUP_ACTS
   "Set up the Acts target(s) explicitly" TRUE )
if (TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
   set(TRACCC_USE_SYSTEM_ACTS_DEFAULT ON)
else()
   set(TRACCC_USE_SYSTEM_ACTS_DEFAULT OFF)
endif()
option( TRACCC_USE_SYSTEM_ACTS
   "Pick up an existing installation of Acts from the build environment"
   ${TRACCC_USE_SYSTEM_ACTS_DEFAULT} )
unset(TRACCC_USE_SYSTEM_ACTS_DEFAULT)
if( TRACCC_SETUP_ACTS )
   if( TRACCC_USE_SYSTEM_ACTS )
      find_package( Acts REQUIRED COMPONENTS PluginJson )
      find_package( detray REQUIRED )
   else()
      add_subdirectory( extern/acts )
   endif()
endif()

# Set up GoogleTest.
include( CTest )
if (BUILD_TESTING AND TRACCC_BUILD_TESTING)
set (TRACCC_DEFAULT_SETUP_GOOGLETEST TRUE)
endif ()
option( TRACCC_SETUP_GOOGLETEST
   "Set up the GoogleTest target(s) explicitly" ${TRACCC_DEFAULT_SETUP_GOOGLETEST} )
if (TRACCC_USE_SYSTEM_LIBS OR TRACCC_USE_SPACK_LIBS)
   set(TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT ON)
else()
   set(TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT OFF)
endif()
option( TRACCC_USE_SYSTEM_GOOGLETEST
   "Pick up an existing installation of GoogleTest from the build environment"
   ${TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT} )
unset(TRACCC_USE_SYSTEM_GOOGLETEST_DEFAULT)
if( TRACCC_SETUP_GOOGLETEST )
   if( TRACCC_USE_SYSTEM_GOOGLETEST )
      find_package( GTest REQUIRED )
   else()
      add_subdirectory( extern/googletest )
   endif()
endif()


option( TRACCC_ENABLE_NVTX_PROFILING
        "Use instrument functions to enable fine grained profiling" FALSE )

# Build the traccc code.
add_subdirectory( core )
add_subdirectory( device/common )
if( TRACCC_BUILD_CUDA_UTILS )
   add_subdirectory( device/cuda_utils )
endif()
if( TRACCC_BUILD_CUDA )
   add_subdirectory( device/cuda )
endif()
if( TRACCC_BUILD_HIP )
   add_subdirectory( device/hip )
endif()
if( TRACCC_BUILD_SYCL_UTILS )
   add_subdirectory( device/sycl_utils )
endif()
if( TRACCC_BUILD_SYCL )
   add_subdirectory( device/sycl )
endif()
if( TRACCC_BUILD_ALPAKA )
   add_subdirectory( device/alpaka )
endif()
if ( TRACCC_BUILD_IO )
   add_subdirectory( io )
else()
   message(STATUS "traccc::io not built, traccc::performance and traccc::simulation are forcefully switched off.")
endif()

if ( TRACCC_BUILD_PERFORMANCE )
   if ( NOT TRACCC_BUILD_IO )
      message(FATAL_ERROR "traccc::io is disabled, but it is required to build the performance tests.")
   endif()
   add_subdirectory( performance )
endif()
if ( TRACCC_BUILD_SIMULATION )
   if ( NOT TRACCC_BUILD_IO )
      message(FATAL_ERROR "traccc::io is disabled, but it is required to build the simulation.")
   endif()
   add_subdirectory( simulation )
endif()

if ( TRACCC_BUILD_EXAMPLES )
   # Find Boost.
   find_package( Boost CONFIG REQUIRED COMPONENTS program_options filesystem )
   if ( NOT TRACCC_BUILD_IO )
      message(FATAL_ERROR "traccc::io is disabled, but it is required to build the examples.")
   endif()
   if ( NOT TRACCC_BUILD_PERFORMANCE )
      message(FATAL_ERROR "traccc::performance is disabled, but it is required to build the examples.")
   endif()
   add_subdirectory( examples )
endif()

# Set up the test(s).
if( BUILD_TESTING AND TRACCC_BUILD_TESTING )
   if ( NOT TRACCC_BUILD_IO )
      message(FATAL_ERROR "traccc::io is disabled, but it is required to build the tests.")
   endif()
   if ( NOT TRACCC_BUILD_SIMULATION )
      message(FATAL_ERROR "traccc::simulation is disabled, but it is required to build the tests.")
   endif()
   if ( NOT TRACCC_BUILD_PERFORMANCE )
      message(FATAL_ERROR "traccc::performance is disabled, but it is required to build the tests.")
   endif()
   add_subdirectory( tests )
endif()

# Set up the packaging of the project.
include( traccc-packaging )
