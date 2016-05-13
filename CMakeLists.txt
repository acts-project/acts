cmake_minimum_required (VERSION 2.8)
project (ACTS)

# Use GNU-style hierarchy for installing build products
include (GNUInstallDirs)

# set the current version of ACTS
set (ACTS_VERSION 0.01.00)

# convert the version into a plain number to be used in the version header file
string (REPLACE "." "" ACTS_VERSION_NUMBER ${ACTS_VERSION})

# default installation direction
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "/opt/ACTS/${ACTS_VERSION}" CACHE PATH "installation directory" FORCE)
endif (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# default parameter definitions
<<<<<<< HEAD
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "ACTS/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE FILEPATH "default track parameter definition")

# default message macros
set (ACTS_MSG_MACROS_PLUGIN "ACTS/Utilities/detail/DefaultMsgMacros.hpp" CACHE FILEPATH "default message macros")

# delegate
add_subdirectory (Core)
add_subdirectory (doc)
add_subdirectory (Plugins)

# add version header
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/ACTSVersion.h ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ACTS)

# add CMake infrastructure
include (cmake/ACTSCreateConfig.cmake)

