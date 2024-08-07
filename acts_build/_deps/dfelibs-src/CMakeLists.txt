# to allow per-target cxx standard requirement
cmake_minimum_required(VERSION 3.8)

# Check if being used directly or via add_subdirectory
set(dfelibs_MASTER_PROJECT FALSE)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(dfelibs_MASTER_PROJECT TRUE)
endif()

project(dfelibs VERSION 20200416 LANGUAGES CXX)

# CMake include(s).
include(GNUInstallDirs)

# options are on by default if build directly, i.e. not via add_subdirectory
option(dfelibs_BUILD_EXAMPLES "Build examples" ${dfelibs_MASTER_PROJECT})
option(dfelibs_BUILD_UNITTESTS "Build unit tests" ${dfelibs_MASTER_PROJECT})
option(dfelibs_ENABLE_INSTALL "Enable library installation" ${dfelibs_MASTER_PROJECT})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# enable most warnings and treat them as errors
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")
set(CMAKE_CXX_EXTENSIONS off)

# define header-only library
add_library(dfelibs INTERFACE)
target_include_directories(dfelibs INTERFACE
   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(dfelibs INTERFACE cxx_std_14)
add_library(dfelibs::dfelibs ALIAS dfelibs)

if(dfelibs_ENABLE_INSTALL)
  install(DIRECTORY dfe DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(TARGETS dfelibs EXPORT dfelibs-exports)

  include(CMakePackageConfigHelpers)
  set(CMAKE_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/dfelibs-${PROJECT_VERSION}")
  install(EXPORT dfelibs-exports
    NAMESPACE "dfelibs::"
    FILE "dfelibs-config-targets.cmake"
    DESTINATION "${CMAKE_INSTALL_CMAKEDIR}")
  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/dfelibs-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_CMAKEDIR
    NO_CHECK_REQUIRED_COMPONENTS_MACRO)
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config-version.cmake"
    COMPATIBILITY "AnyNewerVersion")
  install( FILES
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_CMAKEDIR}")
endif()

# optional components
if(dfelibs_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
if(dfelibs_BUILD_UNITTESTS)
  enable_testing()
  add_subdirectory(unittests)
endif()

# source code formatting targets (only if build directly)
if(dfelibs_MASTER_PROJECT)
  include(ClangFormatTargets)
  add_format_targets(dfe/*.hpp examples/*.cpp unittests/*.cpp unittests/*.hpp)
endif()
