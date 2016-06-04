cmake_minimum_required (VERSION 2.8)
project (ACTS)

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

# set the current version of ACTS
set (ACTS_VERSION 0.01.00)

# convert the version into a plain number to be used in the version header file
string (REPLACE "." "" ACTS_VERSION_NUMBER ${ACTS_VERSION})

# default parameter definitions
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "ACTS/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE FILEPATH "default track parameter definition")

# default message macros
set (ACTS_MSG_MACROS_PLUGIN "ACTS/Utilities/detail/DefaultMsgMacros.hpp" CACHE FILEPATH "default message macros")

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
enable_testing ()
add_subdirectory (Tests)

# add version header
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/ACTSVersion.h ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ACTS)

# add CMake infrastructure
include (cmake/ACTSCreateConfig.cmake)

