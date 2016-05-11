cmake_minimum_required (VERSION 2.6)
project (ACTS)

# Use GNU-style hierarchy for installing build products
include(GNUInstallDirs)

# set the current version of ACTS
set (ACTS_VERSION 0.01.00)

# convert the version into a plain number to be used in the version header file
string (REPLACE "." "" ACTS_VERSION_NUMBER ${ACTS_VERSION})

# default installation direction
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "/opt/ACTS/${ACTS_VERSION}" CACHE PATH "installation directory" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# default parameter definitions
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "ACTS/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE STRING "default track parameter definition")

# default message macros
set (ACTS_MSG_MACROS_PLUGIN "ACTS/Utilities/detail/DefaultMsgMacros.hpp" CACHE STRING "default message macros")

# set location of cmake configuration files
set(ConfigPackageLocation lib/cmake/ACTS)

add_subdirectory (Core)
add_subdirectory (doc)
add_subdirectory (Plugins)

# add version header
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ACTSVersion.h
               ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ACTS)

# add CMake infrastructure
include(cmake/ACTSCreateConfig.cmake)

