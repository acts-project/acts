cmake_minimum_required (VERSION 3.5)

project (ACTS LANGUAGES CXX)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# print all available options
if (DEFINED PRINT_OPTIONS)
  # silence cmake
  set (_print_options ${PRINT_OPTIONS})
  # include options file
  include (cmake/PrintOptions.cmake)
  unset (PRINT_OPTIONS CACHE)
  return ()
endif (DEFINED PRINT_OPTIONS)

# include support for colored output
include (cmake/Colors.cmake)

# Use GNU-style hierarchy for installing build products
include (GNUInstallDirs)

# include ACTS specific compilation options
include (cmake/ACTSCompilerOptions.cmake)

# include ACTS specific macros
include (cmake/ACTSFunctions.cmake)

# include versions for dependencies of ACTS
include (cmake/ACTSDependencies.cmake)

# Configure CCache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message (STATUS "${Yellow}use ccache${ColorReset}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# silence warning about missing RPATH on Mac OSX
set (CMAKE_MACOSX_RPATH 1)

# set the current version of ACTS
set (ACTS_VERSION dev)

# convert the version into a plain number to be used in the version header file
string (REPLACE "." "" ACTS_VERSION_NUMBER ${ACTS_VERSION})

# default parameter definitions
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "ACTS/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE FILEPATH "default track parameter definition")

# finding boost on lxplus does not work with standard boost-cmake look up
set (Boost_NO_BOOST_CMAKE ON CACHE BOOL "disable boost-cmake find mechanism")

# initialize empty list of built components
set (_supported_components)

# delegate
add_subdirectory (Core)
add_subdirectory (doc)
add_subdirectory (Plugins)
add_subdirectory (Examples)

# add unit tests
include (CTest)
add_subdirectory (Tests)

# prepare integration tests
# add custom target that runs the invokation shell script
option (BUILD_INTEGRATION_TESTS "build integration tests" OFF)

if (BUILD_INTEGRATION_TESTS)

  add_custom_target(integration_tests
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
  )

  # macro to collect integration test binaries from subdirectories
  macro(add_integration_test)
    file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    add_custom_target("Integration${ARGV0}" 
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMAND "${_relPath}/${ARGV0}")
    add_dependencies(integration_tests "Integration${ARGV0}")
  endmacro()

  add_subdirectory (IntegrationTests)

else (BUILD_INTEGRATION_TESTS)
  message (STATUS "${Blue}disable build of Integration Tests${ColorReset}")
endif (BUILD_INTEGRATION_TESTS)


# add version header
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/ACTSVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ACTS)

# add setup scripts
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh DESTINATION ${CMAKE_INSTALL_BINDIR})

# list components which are going to be installed
message (STATUS "going to install the following components:")
foreach (_scomp ${_supported_components})
  message (STATUS "   - ${_scomp}")
endforeach ()

# add CMake infrastructure
include (cmake/ACTSCreateConfig.cmake)
