cmake_minimum_required (VERSION 2.6)
project (ACTS)

# set the current version of ACTS
set (ACTS_VERSION 0.0.1)

# default installation direction
set (CMAKE_INSTALL_PREFIX "/opt/ACTS/${ACTS_VERSION}" CACHE STRING "installation directory")

# specify compile flags
add_definitions (-std=c++14)
add_definitions (-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"ACTS/Utilities/detail/DefaultParameterDefinitions.h\")

add_subdirectory (Core)
add_subdirectory (doc)

include(GenerateExportHeader)
generate_export_header (ACTSCore)
set_property (TARGET ACTSCore PROPERTY VERSION ${ACTS_VERSION})

set(ConfigPackageLocation lib/cmake/ACTS)
include(CMakePackageConfigHelpers)
write_basic_package_version_file ("${ConfigPackageLocation}/ACTSConfigVersion.cmake" VERSION ${ACTS_VERSION} COMPATIBILITY AnyNewerVersion)
configure_file(cmake/ACTSConfig.cmake "${CMAKE_INSTALL_PREFIX}/${ConfigPackageLocation}/ACTSConfig.cmake" COPYONLY)
install(EXPORT ACTSTargets FILE ACTSTargets.cmake DESTINATION ${ConfigPackageLocation})

